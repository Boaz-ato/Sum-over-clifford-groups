import numpy as np
from itertools import combinations


def one_qubit_diagonal_group():
    nmatriz = 2**1
    num_clifford_gates = 4
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)

    diagonal_group = []
    for i in range(1,5):
        P_power = np.linalg.matrix_power(P, i)
        diagonal_group.append(P_power)

    f = open(f'Diagonal_Clifford_group_1qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))


    return np.array(diagonal_group)


def two_qubit_diagonal_group():
    # Define matrices
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)
    I1 = np.eye(2, dtype=complex)  # Identity matrix 2x2
    CZ = np.array([[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, -1]], dtype=complex)

    # Define tensor products
    P1 = np.kron(P, I1)
    P2 = np.kron(I1, P)
    P1P2 = np.kron(P, P)

    # Create the diagonal group
    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 3):
                    # Calculate powers and multiply matrices
                    Aux = np.linalg.matrix_power(P1, i) @ np.linalg.matrix_power(P2, j) @ \
                          np.linalg.matrix_power(P1P2, k) @ np.linalg.matrix_power(CZ, l)
                    diagonal_group.append(Aux)

    nmatriz = 2 ** 2
    num_clifford_gates = len(diagonal_group)

    f = open(f'Diagonal_Clifford_group_2qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        # Write real part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        # Write imaginary part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))

    f.close()
    return np.array(diagonal_group)


def three_qubit_diagonal_group():
    # Define basic matrices
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)
    I1 = np.eye(2, dtype=complex)  # Identity matrix 2x2

    # Define tensor products for 3 qubits
    P1 = np.kron(P, np.kron(I1, I1))
    P2 = np.kron(I1, np.kron(P, I1))
    P3 = np.kron(I1, np.kron(I1, P))
    I = np.kron(I1, np.kron(I1, I1))

    # Define controlled-Z gates
    CZ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]], dtype=complex)
    CZ12 = np.kron(CZ, I1)
    CZ23 = np.kron(I1, CZ)

    # Define CZ13 gate (specific construction from MATLAB)
    CZ13 = np.diag([1, 1, 1, 1, 1, -1, 1, -1]).astype(complex)

    # Create the diagonal group
    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 3):
                    for m in range(1, 3):
                        for n in range(1, 3):
                            # Calculate powers and multiply matrices
                            Aux = np.linalg.matrix_power(P1, i) @ np.linalg.matrix_power(P2, j) @ \
                                  np.linalg.matrix_power(P3, k) @ np.linalg.matrix_power(CZ12, l) @ \
                                  np.linalg.matrix_power(CZ23, n) @ np.linalg.matrix_power(CZ13, m)
                            diagonal_group.append(Aux)

    nmatriz = 2 ** 3
    num_clifford_gates = len(diagonal_group)

    f = open('Diagonal_Clifford_group_3qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        # Write real part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        # Write imaginary part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))

    f.close()
    return np.array(diagonal_group)


def four_qubit_diagonal_group():
    # Define basic matrices
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)
    I1 = np.eye(2, dtype=complex)  # Identity matrix 2x2

    # Define tensor products for 4 qubits
    P1 = np.kron(P, np.kron(I1, np.kron(I1, I1)))
    P2 = np.kron(I1, np.kron(P, np.kron(I1, I1)))
    P3 = np.kron(I1, np.kron(I1, np.kron(P, I1)))
    P4 = np.kron(I1, np.kron(I1, np.kron(I1, P)))

    # Define controlled-Z gates for all pairs
    CZ = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]], dtype=complex)

    # CZ gates between adjacent qubits
    CZ12 = np.kron(CZ, np.kron(I1, I1))
    CZ23 = np.kron(I1, np.kron(CZ, I1))
    CZ34 = np.kron(I1, np.kron(I1, CZ))

    # CZ gates between non-adjacent qubits (constructed as diagonal matrices)
    # CZ13: qubits 1 and 3
    CZ13 = np.diag([1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1]).astype(complex)
    # CZ14: qubits 1 and 4
    CZ14 = np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1]).astype(complex)
    # CZ24: qubits 2 and 4
    CZ24 = np.diag([1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1]).astype(complex)

    # Create the diagonal group
    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 3):
                        for n in range(1, 3):
                            for o in range(1, 3):
                                for p in range(1, 3):
                                    for q in range(1, 3):
                                        for r in range(1, 3):
                                            # Calculate powers and multiply matrices
                                            Aux = np.linalg.matrix_power(P1, i) @ np.linalg.matrix_power(P2, j) @ \
                                                  np.linalg.matrix_power(P3, k) @ np.linalg.matrix_power(P4, l) @ \
                                                  np.linalg.matrix_power(CZ12, m) @ np.linalg.matrix_power(CZ23, n) @ \
                                                  np.linalg.matrix_power(CZ34, o) @ np.linalg.matrix_power(CZ13, p) @ \
                                                  np.linalg.matrix_power(CZ14, q) @ np.linalg.matrix_power(CZ24, r)
                                            diagonal_group.append(Aux)

    nmatriz = 2 ** 4
    num_clifford_gates = len(diagonal_group)

    f = open('Diagonal_Clifford_group_4qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        # Write real part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        # Write imaginary part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))

    f.close()
    return np.array(diagonal_group)


def five_qubit_diagonal_group():
    # Define basic matrices
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)
    I1 = np.eye(2, dtype=complex)

    # Define tensor products for 5 qubits
    P1 = np.kron(P, np.kron(I1, np.kron(I1, np.kron(I1, I1))))
    P2 = np.kron(I1, np.kron(P, np.kron(I1, np.kron(I1, I1))))
    P3 = np.kron(I1, np.kron(I1, np.kron(P, np.kron(I1, I1))))
    P4 = np.kron(I1, np.kron(I1, np.kron(I1, np.kron(P, I1))))
    P5 = np.kron(I1, np.kron(I1, np.kron(I1, np.kron(I1, P))))

    # Helper function to create CZ gate for specific qubit pairs
    def create_cz_gate(qubit1, qubit2, n_qubits):
        """Create CZ gate between qubit1 and qubit2 for n_qubits system"""
        dim = 2 ** n_qubits
        cz_matrix = np.eye(dim, dtype=complex)

        for i in range(dim):
            # Check if both qubits are in state |1>
            if ((i >> (n_qubits - qubit1)) & 1) and ((i >> (n_qubits - qubit2)) & 1):
                cz_matrix[i, i] = -1

        return cz_matrix

    # Create all CZ gates for 5 qubits (10 pairs total)
    cz_gates = []
    for i in range(1, 6):
        for j in range(i + 1, 6):
            cz_gates.append(create_cz_gate(i, j, 5))

    # Create the diagonal group
    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 5):
                        # For each CZ gate, use powers 1 and 2
                        for cz_powers in np.ndindex(*([2] * 10)):  # 10 CZ gates, each with power 1 or 2
                            Aux = np.linalg.matrix_power(P1, i) @ np.linalg.matrix_power(P2, j) @ \
                                  np.linalg.matrix_power(P3, k) @ np.linalg.matrix_power(P4, l) @ \
                                  np.linalg.matrix_power(P5, m)

                            for idx, power in enumerate(cz_powers):
                                Aux = Aux @ np.linalg.matrix_power(cz_gates[idx], power + 1)

                            diagonal_group.append(Aux)

    nmatriz = 2 ** 5
    num_clifford_gates = len(diagonal_group)

    f = open('Diagonal_Clifford_group_5qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        # Write real part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        # Write imaginary part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))

    f.close()
    return np.array(diagonal_group)


def six_qubit_diagonal_group():
    # Define basic matrices
    P = np.array([[1, 0], [0, complex(0, 1)]], dtype=complex)
    I1 = np.eye(2, dtype=complex)

    # Define tensor products for 6 qubits
    P1 = np.kron(P, np.kron(I1, np.kron(I1, np.kron(I1, np.kron(I1, I1)))))
    P2 = np.kron(I1, np.kron(P, np.kron(I1, np.kron(I1, np.kron(I1, I1)))))
    P3 = np.kron(I1, np.kron(I1, np.kron(P, np.kron(I1, np.kron(I1, I1)))))
    P4 = np.kron(I1, np.kron(I1, np.kron(I1, np.kron(P, np.kron(I1, I1)))))
    P5 = np.kron(I1, np.kron(I1, np.kron(I1, np.kron(I1, np.kron(P, I1)))))
    P6 = np.kron(I1, np.kron(I1, np.kron(I1, np.kron(I1, np.kron(I1, P)))))

    # Helper function to create CZ gate for specific qubit pairs
    def create_cz_gate(qubit1, qubit2, n_qubits):
        """Create CZ gate between qubit1 and qubit2 for n_qubits system"""
        dim = 2 ** n_qubits
        cz_matrix = np.eye(dim, dtype=complex)

        for i in range(dim):
            # Check if both qubits are in state |1>
            if ((i >> (n_qubits - qubit1)) & 1) and ((i >> (n_qubits - qubit2)) & 1):
                cz_matrix[i, i] = -1

        return cz_matrix

    # Create all CZ gates for 6 qubits (15 pairs total)
    cz_gates = []
    for i in range(1, 7):
        for j in range(i + 1, 7):
            cz_gates.append(create_cz_gate(i, j, 6))

    # Create the diagonal group
    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 5):
                        for n in range(1, 5):
                            # For each CZ gate, use powers 1 and 2
                            for cz_powers in np.ndindex(*([2] * 15)):  # 15 CZ gates, each with power 1 or 2
                                Aux = np.linalg.matrix_power(P1, i) @ np.linalg.matrix_power(P2, j) @ \
                                      np.linalg.matrix_power(P3, k) @ np.linalg.matrix_power(P4, l) @ \
                                      np.linalg.matrix_power(P5, m) @ np.linalg.matrix_power(P6, n)

                                for idx, power in enumerate(cz_powers):
                                    Aux = Aux @ np.linalg.matrix_power(cz_gates[idx], power + 1)

                                diagonal_group.append(Aux)

    nmatriz = 2 ** 6
    num_clifford_gates = len(diagonal_group)

    f = open('Diagonal_Clifford_group_6qubit.txt', 'w')

    for x in range(0, num_clifford_gates):
        group = diagonal_group[x]

        # Write real part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].real))
                else:
                    f.write('%lf\t' % (group[i][j].real))

        # Write imaginary part
        for i in range(0, nmatriz):
            for j in range(0, nmatriz):
                if (j == (nmatriz - 1)):
                    f.write('%lf\n' % (group[i][j].imag))
                else:
                    f.write('%lf\t' % (group[i][j].imag))

    f.close()
    return np.array(diagonal_group)




if __name__ == '__main__':
    one_qubit_diagonal_group()
    two_qubit_diagonal_group()
    three_qubit_diagonal_group()
    print("Generating 4-qubit diagonal group...")
    four_qubit_diagonal_group()
    print("4-qubit diagonal group generated.")

    print("Generating 5-qubit diagonal group...")
    five_qubit_diagonal_group()
    print("5-qubit diagonal group generated.")

    print("Generating 6-qubit diagonal group...")
    #six_qubit_diagonal_group()
    #print("6-qubit diagonal group generated.")


