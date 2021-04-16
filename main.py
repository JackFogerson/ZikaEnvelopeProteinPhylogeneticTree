from Bio import Phylo
from io import StringIO
import numpy

#compares sequences and adds to matrix
def sequenceCompare():
    #original sequence
    sequenceBrazil =       'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAG'

    #Basic data Set
    sequenceVietnam =      'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAGGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequencePuertoRico =   'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACTGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequenceGuatamala =    'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTGGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequenceGuadeloupe =   'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceUSA =          'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceMexico =       'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceMartinique =   'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceYucatan =      'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceHaiti =        'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceFiji =         'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceFrenchGuiana = 'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceAngola =       'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceCuba =         'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceHonduras =     'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACNNGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '

    # add sequences to array and get length
    sequence_list = [sequenceBrazil, sequenceVietnam, sequencePuertoRico, sequenceGuatamala,sequenceGuadeloupe,
                     sequenceUSA, sequenceMexico, sequenceMartinique, sequenceYucatan, sequenceHaiti, sequenceFiji,
                     sequenceFrenchGuiana, sequenceAngola, sequenceCuba, sequenceHonduras]
    n = len(sequence_list)

    # make a numpy array of size n x n
    upgmaArray = numpy.zeros((n, n))

    # fill array with pairwise distances
    upgmaArray = makeArray(sequence_list, upgmaArray)
    print("UPGMA table for sequences:")
    print(upgmaArray)

    # use the array to get the tree contents
    treeContent = computeTree(upgmaArray)

    # add data to file and draw to tree
    tree = open("tree.dnd", "w+")
    tree.write(treeContent)
    drawTree(treeContent)

    # legend
    print("a-Brazil")
    print("b-Vietnam")
    print("c-Puerto Rico")
    print("d-Guatemala")
    print("e-Guadeloupe")
    print("f-USA")
    print("g-Mexico")
    print("h-Martinique")
    print("i-Yucatan")
    print("j-Haiti")
    print("k-Fiji")
    print("l-French Guiana")
    print("m-Angola")
    print("n-Cuba")
    print("o-Honduras")

    # takes tree file and constructs a tree using biopython


def drawTree(treeFile):
    tree = Phylo.read(StringIO(treeFile), "newick")
    Phylo.draw_ascii(tree)


# makes a upgma array using pairwise distances between sequences
def makeArray(list, array):
    for i, element1 in enumerate(list):
        for j, element2 in enumerate(list):
            if j >= i:
                # since matrix is mirrored, no need to enumerate j past i
                break
            # calculate pairwise distance
            distance = pairwise(element1, element2)
            # add pairwise distances to array
            array[i, j] = distance
            array[j, i] = distance
    return array


# use upgma array to implement upgma algorithm
def computeTree(array):
    # create labels for use in Array and tree
    labels = []
    for i in range(ord("A"), ord("O") + 1):
        labels.append(chr(i))

    # get lower triangular matrix
    arrayData = [
        [],
        [array[1][0]],
        [array[2][0], array[2][1]],
        [array[3][0], array[3][1], array[3][2]],
        [array[4][0], array[4][1], array[4][2], array[4][3]],
        [array[5][0], array[5][1], array[5][2], array[5][3], array[5][4]],
        [array[6][0], array[6][1], array[6][2], array[6][3], array[6][4], array[6][5]],
        [array[7][0], array[7][1], array[7][2], array[7][3], array[7][4], array[7][5], array[7][6]],
        [array[8][0], array[8][1], array[8][2], array[8][3], array[8][4], array[8][5], array[8][6], array[8][7]],
        [array[9][0], array[9][1], array[9][2], array[9][3], array[9][4], array[9][5], array[9][6], array[9][7], array[9][8]],
        [array[10][0], array[10][1], array[10][2], array[10][3], array[10][4], array[10][5], array[10][6], array[10][7], array[10][8], array[10][9]],
        [array[11][0], array[11][1], array[11][2], array[11][3], array[11][4], array[11][5], array[11][6], array[11][7], array[11][8], array[11][9], array[11][10]],
        [array[12][0], array[12][1], array[12][2], array[12][3], array[12][4], array[12][5], array[12][6], array[12][7], array[12][8], array[12][9], array[12][10], array[12][11]],
        [array[13][0], array[13][1], array[13][2], array[13][3], array[13][4], array[13][5], array[13][6], array[13][7], array[13][8], array[13][9], array[13][10], array[13][11], array[13][12]],
        [array[14][0], array[14][1], array[14][2], array[14][3], array[14][4], array[14][5], array[14][6], array[14][7], array[14][8], array[14][9], array[14][10], array[14][11], array[14][12], array[14][13]],
    ]

    # distance array for adding branch lengths
    distances = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # while labels are avaliable to be joined
    while len(labels) > 1:
        print("\nTable being analyzed:")
        for i in range(len(arrayData)):
            print(arrayData[i])
        # find lowest value
        x, y = lowValue(arrayData)

        # merges array entries for x,y by averaging data
        distancesx = (arrayData[x][y] / 2) - distances[x]
        distancesy = (arrayData[x][y] / 2) - distances[y]
        distances[x] = (arrayData[x][y] / 2)
        distances[y] = (arrayData[x][y] / 2)
        mergeArray(arrayData, x, y, labels, distances, distancesx, distancesy)

    return labels[0]


# returns pairwise distance for compared sequences
def pairwise(seq1, seq2):
    return sum(x != y for x, y in zip(seq1, seq2))


# finds location of lowest value in array
def lowValue(array):
    # initialize lowest value at infinity
    low = float("inf")
    x, y = -1, -1

    # iterate table for lowest value
    for i in range(len(array)):
        for j in range(len(array[i])):
            if array[i][j] < low:
                low = array[i][j]
                x, y = i, j
    return x, y


# merges array entries for x,y by averaging data
def mergeArray(array, x, y, labels, distance, distancex, distancey):
    # makes sure array is only worked on left side, swap if not
    if y < x:
        x, y = y, x
        tempy = distancey
        distancey = distancex
        distancex = tempy

    distancex = str(distancex)
    distancey = str(distancey)
    # create new label for merged data
    labels[x] = "(" + labels[x] + ":" + distancex + "," + labels[y] + ":" + distancey + ")"

    # recontruct the x row
    row = []
    for i in range(0, x):
        row.append((array[x][i] + array[y][i]) / 2)
    array[x] = row

    # reconstruct x column
    for i in range(x + 1, y):
        array[i][x] = (array[i][x] + array[y][i]) / 2

    for i in range(y + 1, len(array)):
        array[i][x] = (array[i][x] + array[i][y]) / 2
        # since data merged, delete leftover column
        del array[i][y]

    # since data merged, delete leftover row
    del array[y]

    # delete old label
    del labels[y]

    # delete old distance
    del distance[y]


if __name__ == '__main__':
    print("Comparing sequences of Zika Strain Envelope Protein...\n")
    sequenceCompare()
