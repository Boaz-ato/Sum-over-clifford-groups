import highspy
import numpy as np
import time


def get_coefficient_matrix(num_qubits, num_clifford_gates, clifford_file_path):
    """
    Generate coefficient matrix for quantum Clifford group optimization.

    Args:
        num_qubits: Number of qubits
        num_clifford_gates: Number of Clifford gates
        clifford_file_path: Path to the text file containing Clifford group data

    Returns:
        Final: Coefficient matrix for the optimization problem
    """
    clifford_group = np.loadtxt(clifford_file_path, delimiter='\t')
    clifford_group = clifford_group.reshape((2 * num_clifford_gates, 2 ** num_qubits, 2 ** num_qubits))

    matrix_dim = 2 ** num_qubits
    Final = np.zeros((2 * matrix_dim * matrix_dim, 2 * num_clifford_gates), dtype=float)

    print('Getting real part of coefficient matrix')
    for j in range(matrix_dim * matrix_dim):
        k = j // matrix_dim
        l = j % matrix_dim

        count = 0
        for i in range(num_clifford_gates):
            Final[j][count] = clifford_group[2 * i][k][l]
            Final[j][count + 1] = -clifford_group[2 * i][k][l]
            count += 2

    print('Getting imaginary part of coefficient matrix')
    for j in range(matrix_dim * matrix_dim, 2 * matrix_dim * matrix_dim):
        matrix_index = j - matrix_dim * matrix_dim
        k = matrix_index // matrix_dim
        l = matrix_index % matrix_dim

        count = 0
        for i in range(num_clifford_gates):
            Final[j][count] = clifford_group[2 * i + 1][k][l]
            Final[j][count + 1] = -clifford_group[2 * i + 1][k][l]
            count += 2

    return Final


def create_multi_controlled_z_target_matrix(num_qubits):
    """
    Create target matrix for multi-controlled Z gate decomposition.
    For n qubits: first (n-1) qubits are controls, last qubit is target.

    Examples:
    - 2 qubits: CZ gate (qubit 1 controls qubit 2)
    - 3 qubits: CCZ gate (qubits 1,2 control qubit 3)
    - 4 qubits: CCCZ gate (qubits 1,2,3 control qubit 4)

    Args:
        num_qubits: Number of qubits in the system

    Returns:
        target_vector: Flattened target matrix [real_part, imag_part]
    """
    if num_qubits < 2:
        raise ValueError("Multi-controlled Z gate requires at least 2 qubits")

    dim = 2 ** num_qubits
    mcz_matrix = np.eye(dim, dtype=complex)


    all_ones_state = (1 << num_qubits) - 1
    mcz_matrix[all_ones_state, all_ones_state] = -1

    real_part = mcz_matrix.real.flatten()
    imag_part = mcz_matrix.imag.flatten()

    target_vector = np.concatenate([real_part, imag_part])

    return target_vector




def run_highs_solver(num_qubits=2, num_clifford_gates=128, clifford_file_path='Diagonal_Clifford_group_2qubit.txt'):


    h = highspy.Highs()
    print('Getting coefficient matrix for {} qubits and {} Clifford gates from {}'.format(num_qubits, num_clifford_gates, clifford_file_path))
    coefficient_matrix = get_coefficient_matrix(num_qubits, num_clifford_gates, clifford_file_path)
    print('Done')
    num_vars = coefficient_matrix.shape[1]
    num_rows = coefficient_matrix.shape[0]

    inf = highspy.kHighsInf
    cost = np.array([1]*num_vars, dtype=np.double)
    lower = np.array([0]*num_vars, dtype=np.double)
    upper = np.array([inf]*num_vars, dtype=np.double)
    num_nz = 0
    start = 0
    index = 0
    value = 0
    h.addCols(num_vars, cost, lower, upper, num_nz, start, index, value)

    # Add the rows, with the constraint matrix row-wise
    target_vector = create_multi_controlled_z_target_matrix(num_qubits)
    lower = target_vector.copy()
    upper = target_vector.copy()

    num_nz = num_rows*num_vars
    start = np.array([i*num_vars for i in range(num_rows)], dtype=np.int32)
    index = np.array([j for _ in range(num_rows) for j in range(num_vars)], dtype=np.int32)
    value = coefficient_matrix.flatten()
    h.addRows(num_rows, lower, upper, num_nz, start, index, value)
    print('Start solving')
    start_time = time.time()
    h.run()
    end_time = time.time()
    print('Done solving')
    info = h.getInfo()
    model_status = h.getModelStatus()
    print('Model status = ', h.modelStatusToString(model_status))
    print()
    print('Optimal objective = ', info.objective_function_value)
    print('Optimal objective squared =' , info.objective_function_value ** 2)
    print(f"HiGHS solver execution time: {end_time - start_time:.4f} seconds")
    print('Iteration count = ', info.simplex_iteration_count)
    print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
    print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
    print('Basis validity = ', h.basisValidityToString(info.basis_validity))





if __name__ == "__main__":

    run_highs_solver(4, 16384, 'Diagonal_Clifford_group_4qubit.txt')