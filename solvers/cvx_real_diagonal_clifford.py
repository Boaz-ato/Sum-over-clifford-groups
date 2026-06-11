from __future__ import annotations

import os
import time

import cvxpy as cp
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
    all other basis states unchanged.

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


def run_cvxpy_solver(target: str = 'ccz') -> tuple[cp.Problem, NDArray[np.float64]] | None:
    """Solve the Clifford decomposition LP for a given target gate using CVXPY.

    Builds and solves the LP using the CLARABEL solver via the CVXPY modelling
    framework, then writes a summary of the solution to disk.

    Supported target gates:
        'cz', 'ccz', 'cccz', 'ccccz', 'cccccz',
        'ccz_tensor_cz', 'ccz_tensor_ccz'

    Args:
        target: Name of the target gate to decompose.

    Returns:
        A ``(problem, variable_values)`` tuple if a solution was found, or
        ``None`` if no solution exists.

    Raises:
        ValueError: If ``target`` is not a supported gate name.

    Side effects:
        Creates ``cvxpy_final_solutions/`` if it does not exist.
        Writes ``cvxpy_solution_<target>.txt`` and
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
    }
    if target not in qubit_map:
        raise ValueError(
            f"Unsupported target gate '{target}'. "
            f"Supported values are: {sorted(qubit_map)}."
        )
    num_qubits = qubit_map[target]
    print(f'Setting up CVXPY solver for {target} on {num_qubits} qubits')

    coefficient_matrix = get_coefficient_matrix(num_qubits)
    print('Done')

    num_vars = coefficient_matrix.shape[1]

    # Build target vector
    tensor_targets = {
        'ccz_tensor_cz': ccz_tensor_cz,
        'ccz_tensor_ccz': ccz_tensor_ccz,
    }
    if target in tensor_targets:
        target_vector = tensor_targets[target]()
    else:
        target_vector = create_multi_controlled_z_target_matrix(num_qubits)

    # Define and solve the LP
    x = cp.Variable(num_vars, nonneg=True)
    constraints = [coefficient_matrix @ x == target_vector]
    objective = cp.Minimize(cp.sum(x))
    problem = cp.Problem(objective, constraints)

    print('Start solving')
    start_time = time.time()
    problem.solve(solver=cp.CLARABEL, verbose=True)
    end_time = time.time()
    print('Done solving')

    if problem.status not in ["optimal", "optimal_inaccurate"]:
        print(f"Warning: Problem status is {problem.status}")

    variable_values = x.value
    if variable_values is None:
        print("No solution found!")
        return None

    non_zero_count = int(np.sum(np.abs(variable_values) > 1e-8))
    elapsed = end_time - start_time

    os.makedirs('cvxpy_final_solutions', exist_ok=True)

    print('Model status =', problem.status)
    print('Optimal objective =', problem.value)
    print('Optimal objective squared =', problem.value ** 2)
    print(f"CVXPY solver execution time: {elapsed:.4f} seconds")
    print("Number of non-zero variables =", non_zero_count)

    np.savetxt(f'cvxpy_final_solutions/cvxpy_solution_{target}.txt', variable_values)

    non_zero_indices = np.where(np.abs(variable_values) > 1e-8)[0]
    selected_columns = coefficient_matrix[:, non_zero_indices]

    output_filename = f"cvxpy_final_solutions/solution_summary_{target}.txt"
    with open(output_filename, "w") as f:
        f.write(f"Model Status: {problem.status}\n")
        f.write(f"Optimal objective: {problem.value}\n")
        f.write(f"Optimal objective squared: {problem.value ** 2}\n")
        f.write(f"CVXPY solver execution time: {elapsed:.4f} seconds\n")
        f.write(f"Number of non-zero variables: {non_zero_count}\n")
        f.write(f"Solver used: {problem.solver_stats.solver_name}\n")
        if hasattr(problem.solver_stats, 'num_iters'):
            f.write(f"Number of iterations: {problem.solver_stats.num_iters}\n")
        if hasattr(problem.solver_stats, 'solve_time'):
            f.write(f"Solver time: {problem.solver_stats.solve_time}\n")

        f.write("\nNon-zero Clifford groups:\n")
        for i in range(selected_columns.shape[0]):
            row_str = " ".join([f"{val:.6f}" for val in selected_columns[i, :]])
            f.write(f"Row {i}: {row_str}\n")

    print(f"\nNon-zero variables (threshold 1e-8):")
    for idx in non_zero_indices[:10]:
        print(f"x[{idx}] = {variable_values[idx]:.6f}")
    if len(non_zero_indices) > 10:
        print(f"... and {len(non_zero_indices) - 10} more")

    return problem, variable_values


if __name__ == "__main__":
    run_cvxpy_solver('cccccz')
