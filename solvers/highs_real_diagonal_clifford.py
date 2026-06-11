from __future__ import annotations

import os
import time

import highspy
import numpy as np
from numpy.typing import NDArray


def get_coefficient_matrix(num_qubits: int) -> NDArray[np.float64]:
    """Load the real diagonal Clifford group and build the LP coefficient matrix.

    Each element of the real diagonal Clifford group contributes two columns to
    the matrix — one positive and one negative — so that non-negative LP variables
    can represent both positive and negative coefficients in the decomposition.

    Args:
        num_qubits: Number of qubits (1–6).

    Returns:
        Array of shape (2^num_qubits, 2 * num_clifford_gates) containing the
        coefficient matrix for the LP.

    Raises:
        ValueError: If num_qubits is not in the range [1, 6].
    """
    dataset_paths = {
        1: 'datasets/real_diagonals/one_qubit_real_diagonal_clifford_groups.txt',
        2: 'datasets/real_diagonals/two_qubit_real_diagonal_clifford_groups.txt',
        3: 'datasets/real_diagonals/three_qubit_real_diagonal_clifford_groups.txt',
        4: 'datasets/real_diagonals/four_qubit_real_diagonal_clifford_groups.txt',
        5: 'datasets/real_diagonals/five_qubit_real_diagonal_clifford_groups.txt',
        6: 'datasets/real_diagonals/six_qubit_real_diagonal_clifford_groups.txt',
    }
    if num_qubits not in dataset_paths:
        raise ValueError("Unsupported number of qubits. Supported values are 1 to 6.")

    clifford_group = np.loadtxt(dataset_paths[num_qubits])

    num_clifford_gates = clifford_group.shape[0]
    matrix_dim = 2 ** num_qubits  # number of diagonal elements
    Final = np.zeros((matrix_dim, 2 * num_clifford_gates), dtype=float)

    print('Getting real part of coefficient matrix')
    for j in range(matrix_dim):
        for i in range(num_clifford_gates):
            Final[j][i] = clifford_group[i][j]
            Final[j][i + num_clifford_gates] = -clifford_group[i][j]

    return Final


def create_multi_controlled_z_target_matrix(num_qubits: int) -> NDArray[np.float64]:
    """Create the target diagonal vector for a C^(n-1)Z (multi-controlled Z) gate.

    The gate applies a -1 phase to the all-ones basis state |11...1> and leaves
    all other basis states unchanged. The first (n-1) qubits are controls and
    the last qubit is the target.

    Args:
        num_qubits: Total number of qubits (must be >= 2).

    Returns:
        Array of length 2^num_qubits with +1 entries everywhere except the last,
        which is -1.

    Raises:
        ValueError: If num_qubits < 2.
    """
    if num_qubits < 2:
        raise ValueError("Multi-controlled Z gate requires at least 2 qubits")

    cz_target = np.ones(2 ** num_qubits, dtype=float)
    cz_target[-1] = -1
    return cz_target


def ccz_tensor_cz() -> NDArray[np.float64]:
    """Return the diagonal of the tensor product CCZ ⊗ CZ.

    Returns:
        Array of length 32 (= 8 * 4) with ±1 entries.
    """
    ccz_diag = np.ones(8, dtype=float)
    ccz_diag[7] = -1

    cz_diag = np.ones(4, dtype=float)
    cz_diag[3] = -1

    return np.kron(ccz_diag, cz_diag)


def ccz_tensor_ccz() -> NDArray[np.float64]:
    """Return the diagonal of the tensor product CCZ ⊗ CCZ.

    Returns:
        Array of length 64 (= 8 * 8) with ±1 entries.
    """
    ccz_diag = np.ones(8, dtype=float)
    ccz_diag[7] = -1

    return np.kron(ccz_diag, ccz_diag)


def cz_tensor_cz() -> NDArray[np.float64]:
    """Return the diagonal of the tensor product CZ ⊗ CZ.

    Returns:
        Array of length 16 (= 4 * 4) with ±1 entries.
    """
    cz_diag = np.ones(4, dtype=float)
    cz_diag[3] = -1

    return np.kron(cz_diag, cz_diag)


def multi_cz_unitary(
    n_qubits: int,
    gates: list[tuple[int, ...]],
    *,
    msb_top: bool = True,
    dtype: type = complex,
) -> NDArray:
    """Build the diagonal unitary matrix for a circuit of multi-controlled-Z gates.

    Each gate in the circuit applies a -1 phase to all basis states where every
    qubit in the gate tuple is in state |1>.

    Args:
        n_qubits: Total number of qubits.
        gates: Each element is a tuple of qubit indices for one gate.
            A 2-tuple gives a CZ, a 3-tuple gives a CCZ, etc.
            Order within a tuple does not matter.
        msb_top: If True (default), qubit 0 is the most-significant bit in the
            standard basis ordering |q0 q1 ... q_{n-1}>. Set False to use
            qubit 0 as the least-significant bit.
        dtype: NumPy dtype for the output matrix.

    Returns:
        Real diagonal matrix of shape (2^n_qubits, 2^n_qubits).

    Raises:
        ValueError: If a gate acts on fewer than 2 qubits or a qubit index is
            out of range.
    """
    dim = 1 << n_qubits
    diag = np.ones(dim, dtype=dtype)

    def bitpos(q: int) -> int:
        return (n_qubits - 1 - q) if msb_top else q

    state_ids = np.arange(dim, dtype=np.uint64)
    for gate in gates:
        if len(gate) < 2:
            raise ValueError(f"Gate {gate} must act on at least 2 qubits.")
        mask = 0
        for q in gate:
            if not (0 <= q < n_qubits):
                raise ValueError(f"Qubit index {q} out of range for n={n_qubits}.")
            mask |= (1 << bitpos(q))
        hits = (state_ids & mask) == mask
        diag[hits] *= -1

    return np.diag(np.diag(diag)).real


def run_highs_solver(target: str = 'ccz', num_threads: int = 1) -> None:
    """Solve the Clifford decomposition LP for a given target gate using HiGHS.

    Builds the LP from the pre-computed real diagonal Clifford group datasets,
    solves it with the HiGHS simplex method, and writes results to disk.

    Supported target gates:
        'cz', 'ccz', 'cccz', 'ccccz', 'cccccz',
        'cz_tensor_cz', 'ccz_tensor_cz', 'ccz_tensor_ccz',
        'circuit_cz_4', 'circuit_cz_5', 'circuit_cz_6'

    Args:
        target: Name of the target gate to decompose.
        num_threads: Number of threads passed to the HiGHS parallel option.
            Results are written to ``highs_solution_<num_threads>/``.

    Raises:
        ValueError: If ``target`` is not a supported gate name.

    Side effects:
        Creates directory ``highs_solution_<num_threads>/`` if it does not exist.
        Writes ``highs_solution_<target>`` (HiGHS raw output) and
        ``solution_summary_<target>.txt`` into that directory.
    """
    qubit_map = {
        'cz': 2,
        'ccz': 3,
        'cccz': 4,
        'ccccz': 5,
        'cccccz': 6,
        'ccz_tensor_cz': 5,
        'ccz_tensor_ccz': 6,
        'cz_tensor_cz': 4,
        'circuit_cz_4': 4,
        'circuit_cz_5': 5,
        'circuit_cz_6': 6,
    }
    if target not in qubit_map:
        raise ValueError(
            f"Unsupported target gate '{target}'. "
            f"Supported values are: {sorted(qubit_map)}."
        )
    num_qubits = qubit_map[target]
    print(f'Setting up HiGHS solver for {target} on {num_qubits} qubits')

    h = highspy.Highs()
    coefficient_matrix = get_coefficient_matrix(num_qubits)
    print('Done')
    num_vars = coefficient_matrix.shape[1]
    num_rows = coefficient_matrix.shape[0]

    inf = highspy.kHighsInf
    cost = np.ones(num_vars, dtype=np.double)
    lower = np.zeros(num_vars, dtype=np.double)
    upper = np.full(num_vars, inf, dtype=np.double)
    h.addCols(num_vars, cost, lower, upper, 0, 0, 0, 0)

    # Build the target vector
    tensor_targets = {
        'ccz_tensor_cz': ccz_tensor_cz,
        'ccz_tensor_ccz': ccz_tensor_ccz,
        'cz_tensor_cz': cz_tensor_cz,
    }
    circuit_targets = {
        'circuit_cz_4': lambda: multi_cz_unitary(4, [(0, 1, 2), (1, 2, 3)]),
        'circuit_cz_5': lambda: multi_cz_unitary(5, [(0, 1, 2), (1, 2, 3), (2, 3, 4)]),
        'circuit_cz_6': lambda: multi_cz_unitary(6, [(0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5)]),
    }
    if target in tensor_targets:
        target_vector = tensor_targets[target]()
    elif target in circuit_targets:
        target_vector = circuit_targets[target]()
    else:
        target_vector = create_multi_controlled_z_target_matrix(num_qubits)

    row_lower = target_vector.copy()
    row_upper = target_vector.copy()
    num_nz = num_rows * num_vars
    start = np.array([i * num_vars for i in range(num_rows)], dtype=np.int32)
    index = np.array([j for _ in range(num_rows) for j in range(num_vars)], dtype=np.int32)
    value = coefficient_matrix.flatten()
    h.addRows(num_rows, row_lower, row_upper, num_nz, start, index, value)

    h.setOptionValue("solver", "simplex")
    h.setOptionValue("parallel", "on")
    h.setOptionValue("threads", num_threads)

    print('Start solving')
    start_time = time.time()
    h.run()
    end_time = time.time()
    print('Done solving')

    solution = h.getSolution()
    variable_values = solution.col_value
    non_zero_count = int(np.sum(np.abs(variable_values) > 1e-8))

    solution_dir = f'highs_solution_{num_threads}'
    if not os.path.exists(solution_dir):
        os.makedirs(solution_dir)
        print(f"Created directory: {solution_dir}")

    info = h.getInfo()
    h.writeSolution(f'{solution_dir}/highs_solution_{target}', 0)
    model_status = h.getModelStatus()
    elapsed = end_time - start_time

    print('Model status = ', h.modelStatusToString(model_status))
    print('Optimal objective = ', info.objective_function_value)
    print('Optimal objective squared = ', info.objective_function_value ** 2)
    print(f"HiGHS solver execution time: {elapsed:.4f} seconds")
    print('Iteration count = ', info.simplex_iteration_count)
    print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
    print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
    print('Basis validity = ', h.basisValidityToString(info.basis_validity))
    print("Number of non-zero variables = ", non_zero_count)

    non_zero_indices = np.where(np.abs(variable_values) > 1e-8)[0]
    selected_columns = coefficient_matrix[:, non_zero_indices]

    output_filename = f"{solution_dir}/solution_summary_{target}.txt"
    with open(output_filename, "w") as f:
        f.write(f"Model Status: {h.modelStatusToString(model_status)}\n")
        f.write(f"Optimal objective: {info.objective_function_value}\n")
        f.write(f"Optimal objective squared: {info.objective_function_value ** 2}\n")
        f.write(f"HiGHS solver execution time: {elapsed:.4f} seconds\n")
        f.write(f"Iteration count: {info.simplex_iteration_count}\n")
        f.write(f"Primal solution status: {h.solutionStatusToString(info.primal_solution_status)}\n")
        f.write(f"Dual solution status: {h.solutionStatusToString(info.dual_solution_status)}\n")
        f.write(f"Basis validity: {h.basisValidityToString(info.basis_validity)}\n")
        f.write(f"Number of non-zero variables: {non_zero_count}\n")

        for i in range(selected_columns.shape[0]):
            row_str = " ".join([f"{val:.6f}" for val in selected_columns[i, :]])
            f.write(f"Row {i}: {row_str}\n")


if __name__ == "__main__":
    run_highs_solver('cccccz', num_threads=10)
