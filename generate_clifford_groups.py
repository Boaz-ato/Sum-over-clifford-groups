import qecc as q
import numpy as np
import math
import cmath




def generate_clifford_groups(n = 1, with_phases = False):
    """
    Generates the Clifford group for a given number of qubits.
    The Clifford group is a group of quantum gates that map Pauli operators to Pauli operators.
    """

    n = int(n)
    iterator = q.clifford_group(n, True)
    if with_phases:

        phases = [
            complex(np.cos(np.pi/4), np.sin(np.pi/4)),
            complex(np.cos(np.pi/2), np.sin(np.pi/2)),
            complex(np.cos(3*np.pi/4), np.sin(3*np.pi/4)),
            complex(np.cos(np.pi), np.sin(np.pi)),
            complex(np.cos(5*np.pi/4), np.sin(5*np.pi/4)),
            complex(np.cos(3*np.pi/2), np.sin(3*np.pi/2)),
            complex(np.cos(7*np.pi/4), np.sin(7*np.pi/4)),
            complex(np.cos(2*np.pi), np.sin(2*np.pi))
        ]

    else:
        phases = [1]

    if n == 1:

        num_clifford_gates = 24
    elif n == 2:
        num_clifford_gates = 11520
    else:
        raise ValueError("n must be 1 or 2")



    if with_phases:
        f = open(f'Clifford_group_{n}_with_phases.txt', 'w')
    else:
        f = open(f'Clifford_group_{n}.txt', 'w')
    nmatriz = 2 ** n
    for x in range(0, num_clifford_gates):
        Matrix = next(iterator).as_unitary()
        for phase in phases:
            Matrix = Matrix * phase
            for i in range(0, nmatriz):
                for j in range(0, nmatriz):
                    if (j == (nmatriz - 1)):
                        f.write('%lf\n' % (Matrix[i][j].real))
                    else:
                        f.write('%lf\t' % (Matrix[i][j].real))

            for i in range(0, nmatriz):
                for j in range(0, nmatriz):
                    if (j == (nmatriz - 1)):
                        f.write('%lf\n' % (Matrix[i][j].imag))
                    else:
                        f.write('%lf\t' % (Matrix[i][j].imag))


if __name__ == "__main__":

    generate_clifford_groups(2, with_phases=True)