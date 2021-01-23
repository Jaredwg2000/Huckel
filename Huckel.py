"""Calculates Huckel pi orbital energies."""

import numpy as np


def printEvals(Beta):
    """Print a list of eigenvalues.

    Takes the Huckel matrix as input, finds the eigenvalues,
    rounds them and prints them, sorted.
    """
    rounded = []

    evals, evecs = np.linalg.eig(Beta)
    for i in range(len(evals)):
        # Need real part: was getting some complex numbers for buckyball
        # with imaginary part on the order of 10^-16. This is due to
        # some kind of floating point error in the algorithm that finds
        # the eigenvalues.
        rounded.append(round(np.real(evals[i]), 3))

    rounded.sort()
    # Create a dictionary with 1 key for each eigenvalue and value equal
    # to the degeneracy of that eigenvalue
    eigenvalues = {}
    for i in rounded:
        if i in eigenvalues:
            eigenvalues[i] += 1
        else:
            eigenvalues[i] = 1
    print(eigenvalues)


def LinearPolyene():
    """Pi orbital energies for linear polyenes."""
    # Input the length of the chain
    print("Length?")
    n = int(input())

    Beta = np.zeros([n, n])

    # Link every position to those adjacent to it.
    for i in range(n):
        if(i != 0):
            Beta[i, i - 1] = -1
        if(i != (n - 1)):
            Beta[i, i + 1] = -1

    # Print the eigenvalues.
    printEvals(Beta)


def CyclicPolyene():
    """Pi orbital energies of cyclic polyenes."""
    # Input the number of atoms in the ring.
    print("Number of atoms?")
    n = int(input())

    Beta = np.zeros([n, n])

    # Connect every atom to those it is adjacent to, looping around from n to 0
    for i in range(n):
        Beta[i, (i + 1) % n] = -1
        Beta[i, (i - 1) % n] = -1

    # Print the eigenvalues.
    printEvals(Beta)


def Polyhedron():
    """Pi orbital energies for polyhedra.

    This one was slightly more awkward, as there is often no easy way to
    find the connections, so they are all hard coded, either by finding
    loops on a graph of the polyhedron, or by literally hard-coding the
    Huckel matrix.
    """
    print("Which type? (C)ube; (T)etrahedron; (D)odecahedron; (O)ctahedron; " +
          "(I)cosahedron; (B)uckminsterfullerene")
    polyType = input()

    if(polyType == "C"):
        # Huckel matrix for a cube.
        Beta = np.array([[0, -1, -1, 0, -1, 0, 0, 0],
                         [-1, 0, 0, -1, 0, -1, 0, 0],
                         [-1, 0, 0, -1, 0, 0, -1, 0],
                         [0, -1, -1, 0, 0, 0, 0, -1],
                         [-1, 0, 0, 0, 0, -1, -1, 0],
                         [0, -1, 0, 0, -1, 0, 0, -1],
                         [0, 0, -1, 0, -1, 0, 0, -1],
                         [0, 0, 0, -1, 0, -1, -1, 0]])

        # Print the eigenvalues.
        printEvals(Beta)

    if(polyType == "T"):
        # Huckel matrix for a tetrahedron.
        Beta = np.identity(4) - np.ones([4, 4])

        # Print the eigenvalues.
        printEvals(Beta)

    if(polyType == "D"):
        # Huckel matrix for a dodecahedron. Tried to simplify it by finding
        # Cycles in the graph. Not sure if it was actually simpler than
        # Just writing out the matrix.
        Beta = np.zeros([20, 20])
        for i in range(5):
            Beta[i, (i + 1) % 5] = -1
            Beta[i, (i - 1) % 5] = -1
            Beta[i, (5 + 2*i)] = -1
            Beta[(5 + 2*i), i] = -1
            Beta[6 + 2*i, 15 + i] = -1
            Beta[15 + i, 6 + 2*i] = -1

            if(i != 4):
                Beta[15 + i, 16 + i] = -1
            else:
                Beta[15 + i, 11 + i] = -1

            if(i != 0):
                Beta[15 + i, 14 + i] = -1
            else:
                Beta[15 + i, 19 + i] = -1

        for i in range(10):
            if(i != 9):
                Beta[5 + i, 6 + i] = -1
            if(i != 0):
                Beta[5 + i, 4 + i] = -1

        # Print the eigenvalues
        printEvals(Beta)

    if(polyType == "O"):
        # Huckel matrix for an ortahedron.
        Beta = np.zeros([6, 6])
        for i in range(3):
            Beta[i, (i + 1) % 3] = -1
            Beta[(i + 1) % 3, i] = -1

            Beta[i + 3, ((i + 1) % 3)+3] = -1
            Beta[((i + 1) % 3) + 3, i + 3] = -1

            if(i != 0):
                Beta[i, 4] = -1
                Beta[4, i] = -1
            if(i != 1):
                Beta[i, 5] = -1
                Beta[5, i] = -1
            if(i != 2):
                Beta[i, 3] = -1
                Beta[3, i] = -1

        # Print the eigenvalues
        printEvals(Beta)

    if(polyType == "I"):
        # Huckel matrix for an icosahedron.
        Beta = np.zeros([12, 12])

        # See note on buckyball.
        with open("Icosahedron.dat", "r") as f:
            for line in f:
                splitL = line.split()
                Beta[int(splitL[0]) - 1, int(splitL[1]) - 1] = -1
                Beta[int(splitL[1]) - 1, int(splitL[0]) - 1] = -1

        # Print the eigenvalues
        printEvals(Beta)

    if(polyType == "B"):
        # Huckel matrix for a buckyball.
        Beta = np.zeros([60, 60])

        # I procrastrinated on doing this for quite a while. Did the
        # Other tasks in the meantime, so will use a datafile for this
        # instead of writing out however many lines it would take
        # otherwise.
        with open("Buckyball.dat", "r") as f:
            for line in f:
                splitL = line.split()
                Beta[int(splitL[0]) - 1, int(splitL[1]) - 1] = -1
                Beta[int(splitL[1]) - 1, int(splitL[0]) - 1] = -1

        # Print the eigenvalues
        printEvals(Beta)


def OtherFromData():
    """Return Huckel energies for a specified molecule.

    Takes a datafile containing the connections in a molecule,
    with each connection on a separate line, and returns the Huckel
    pi energies. Each bond need only be represented once in the
    datafile. A datafile for butadiene would be:

    1 2
    2 3
    3 4

    """
    print("Number of atoms:")
    noAtoms = int(input())
    Beta = np.zeros([noAtoms, noAtoms])
    print("Path:")
    path = input()

    with open(path, "r") as f:
        for line in f:
            splitL = line.split()
            Beta[int(splitL[0]) - 1, int(splitL[1]) - 1] = -1
            Beta[int(splitL[1]) - 1, int(splitL[0]) - 1] = -1

    # Print the eigenvalues
    printEvals(Beta)


# Ask the user what they want to do.
print("Which type of molecule? (L)inear Polyene; (C)yclic Polyene; " +
      "(P)olyhedron; (O)ther from datafile.")
input1 = input()

if(input1 == "L"):
    LinearPolyene()
elif(input1 == "C"):
    CyclicPolyene()
elif(input1 == "P"):
    Polyhedron()
elif(input1 == "O"):
    OtherFromData()
else:
    print("Invalid option")
