import highspy
import numpy as np
import time
from generate_diagonal_components_of_real_diagonal_Clifford_groups import (
 one_qubit_diagonal_group,
two_qubit_diagonal_group,
three_qubit_diagonal_group,
four_qubit_diagonal_group,
five_qubit_diagonal_group,
six_qubit_diagonal_group
)


def get_coefficient_matrix(num_qubits):
    """
    Generate coefficient matrix for quantum Clifford group optimization.

    Args:
        num_qubits: Number of qubits

    Returns:
        Final: Coefficient matrix for the optimization problem
    """

    if num_qubits == 1:
        clifford_group = one_qubit_diagonal_group()
    elif num_qubits == 2:
        clifford_group = two_qubit_diagonal_group()
    elif num_qubits == 3:
        clifford_group = three_qubit_diagonal_group()
    elif num_qubits == 4:
        clifford_group = four_qubit_diagonal_group()
    elif num_qubits == 5:
        clifford_group = five_qubit_diagonal_group()
    elif num_qubits == 6:
        clifford_group = six_qubit_diagonal_group()
    else:
        raise ValueError("Unsupported number of qubits. Supported values are 1 to 6.")

    num_clifford_gates = clifford_group.shape[0]
    matrix_dim = 2 ** num_qubits
    Final = np.zeros((matrix_dim, 2 * num_clifford_gates), dtype=float)

    print('Getting real part of coefficient matrix')
    for j in range(matrix_dim):
        for i in range(num_clifford_gates):
            Final[j][i] = clifford_group[i][j]
            Final[j][i + num_clifford_gates] = -clifford_group[i][j]

    return Final


def create_multi_controlled_z_target_matrix(num_qubits):
    """
    Create target matrix for multi-controlled Z gate decomposition.
    For n qubits: first (n-1) qubits are controls, last qubit is target.

    Args:
        num_qubits: Number of qubits in the system

    Returns:
        target_vector: Flattened target matrix [real_part, imag_part]
    """
    if num_qubits < 2:
        raise ValueError("Multi-controlled Z gate requires at least 2 qubits")

    cz_target = np.ones(2 ** num_qubits, dtype=float)
    cz_target[-1] = -1

    return cz_target




def run_highs_solver(num_qubits):


    h = highspy.Highs()
    coefficient_matrix = get_coefficient_matrix(num_qubits)
    print('Done')
    num_vars = coefficient_matrix.shape[1]
    num_rows = coefficient_matrix.shape[0]

    inf = highspy.kHighsInf
    # The constraint matrix is defined with the rows below, but parameters
    # for an empty (column-wise) matrix must be passed
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
    print('Optimal objective = ', info.objective_function_value)
    print('Optimal objective squared =' , info.objective_function_value ** 2)
    print(f"HiGHS solver execution time: {end_time - start_time:.4f} seconds")
    print('Iteration count = ', info.simplex_iteration_count)
    print('Primal solution status = ', h.solutionStatusToString(info.primal_solution_status))
    print('Dual solution status = ', h.solutionStatusToString(info.dual_solution_status))
    print('Basis validity = ', h.basisValidityToString(info.basis_validity))





if __name__ == "__main__":

    #run_highs_solver()
    #final = get_coefficient_matrix(2, 128, 'Diagonal_Clifford_group_2qubit.txt')

    #print(final)
    #run_highs_solver(5, 1048576, 'Diagonal_Clifford_group_5qubit.txt')
    #run_highs_solver(4, 16384, 'Diagonal_Clifford_group_4qubit.txt')
    #run_highs_solver(3, 512, 'Diagonal_Clifford_group_3qubit.txt')
    run_highs_solver(4)