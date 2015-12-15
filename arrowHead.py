import argparse
import matplotlib.pyplot as plt
import numpy
from PIL import Image
from PIL import ImageOps

def fn(a):
    if a:
        return a
    return 1

def readMat(input_file):
    print "Reading Matrix..."
    np = numpy.loadtxt(input_file)
    print "Finished Reading Matrix..."
    return np

def normalize(np):
    print "Normalizing Matrix..."
    sumRow = np.sum(axis=0)
    sumCol = np.sum(axis=1)
    row,col = np.shape
    for r in range(row):
        for c in range(col):
            np[r][c] = np[r][c]/fn(sumRow[r]*sumCol[c])
    print "Finished Normalizing Matrix..."

def display(np):
    print "About to display image"
    # Write image to file
    im = Image.fromarray(numpy.uint8(plt.cm.gist_earth(np)*255))
    im.save("img.jpg")

    # Display image
    plt.matshow(np, figure=1, cmap=plt.cm.gray)
    plt.show()

#def readMat(input_file):
#    fd = open(input_file, 'r')
#    # Read input file, split by newline, split by whitespace and then
#    # convert each to float
#    M = [map(lambda x:float(x), y) for y in [row.split() for row in
#                                             fd.read()[:-1].split('\n')]]
#    return (M, len(M), len(M[0]))
#

#def normalize(M, numRow, numCol):
#    # Get the sum of each row
#    sumRow = [fn(sum(x)) for x in M]
#
#    # Get the sum of each column
#    sumCol = [fn(cval) for cval in reduce(lambda x,y: [a+b for a,b in zip(x,y)], M)]
#
#    # Divide each entry by sum of its column
#    for row in range(numRow):
#        for col in range(numCol):
#            M[row][col] = M[row][col]/(sumRow[row]*sumCol[col])
#
#    # Divide each entry by sum of its column
#    #M = map(lambda x: [a/fn(b) for a,b in zip(x,sumCol)], M)
#    # Divide each entry by sum of its row
#    #M = [map(lambda x: x/fn(sumRow[row]), M[row]) for row in range(len(sumRow))]
#    return M
#

#def display(M, numRow, numCol):
#    plt.matshow(M, figure=1, cmap=plt.cm.gray)
#    plt.show()
#
def computeArrowHead(np):
    print "Computing ArrowHead..."
    A = numpy.zeros(np.shape)
    row,col = np.shape
    for r in range(row):
        for c in range(col):
            d = c - r
            if  r-d >= col or r-d<0:
                continue
            A[r][r+d] = (np[r][r-d]-np[r][r+d])/fn(np[r][r-d]+np[r][r+d])

    print "Finished Computing ArrowHead..."
    return A

def sgn(x):
    if x < 0:
        return -1
    elif x > 0:
        return 1
    return 0

def computeUsgn(A):
    print "Computing Usgn..."
    Usgn = numpy.zeros(A.shape)

    row,col = A.shape

    # Compute U[a][b]
    for d in range(1,row):
        for r in range(row):
            if r+d < col:
                Usgn[r][r+d] = Usgn[r][r+d-1] + numpy.sign(A[r:r+d/2+1,r+d]).sum()
    print "Finished Computing Usgn..."
    return Usgn

def computeLsgn(A):
    print "Computing Lsgn..."
    Lsgn = numpy.zeros(A.shape)

    row,col = A.shape

    # compute L[a][b]
    for d in range(1, row):
        for r in range(row):
            if r+d >= col or 2*(r+d)-r>=col:
                break
            Lsgn[r][r+d] = Lsgn[r][r+d-1] + numpy.sign(A[r+d, r+d:2*(r+d)-r+1]).sum()
    print "Finished Computing Lsgn..."
    return Lsgn

def sum_sign(A):
    ss = 0
    for x in numpy.nditer(A):
        if x < 0:
            ss = ss - 1
        elif x > 0:
            ss = ss + 1
    return ss

def normalizeS(S):
    return  S/S.max()

def computeSsign(Usgn, Lsgn):
    print "Computing Ssign"
    Ssign = Lsgn - Usgn
    print "Finished computing Ssign"
    return Ssign

def computeU(A):
    print "Computing U"
    U = numpy.zeros(A.shape)

    row,col = A.shape

    # Compute U[a][b]
    for d in range(1,row):
        for r in range(row):
            if r+d < col:
                U[r][r+d] = U[r][r+d-1] + A[r:r+d/2+1,r+d].sum()
    print "Finished Computing U..."
    return U

def computeL(A):
    print "Computing L..."
    L = numpy.zeros(A.shape)

    row,col = A.shape

    # compute L[a][b]
    for d in range(1, row):
        for r in range(row):
            if r+d >= col or 2*(r+d)-r>=col:
                break
            L[r][r+d] = L[r][r+d-1] + A[r+d, r+d:2*(r+d)-r+1].sum()
    print "Finished Computing L..."
    return L

def computeSsum(U, L):
    print "Computing Ssum..."
    Ssum = L - U
    print "Finished Computing Ssum..."
    return Ssum

def computeSvar(A, U, L):
    print "Computing Svar..."
    Sx = U + L
    countU = numpy.zeros(A.shape)
    countL = numpy.zeros(A.shape)
    Sx2 = numpy.zeros(A.shape)
    Svar = numpy.zeros(A.shape)
    sq = A**2

    (row,col) = A.shape
    for d in range(1,row):
        for r in range(row):
            if r+d >= col or 2*(r+d)-r >= col:
                break
            countU[r][r+d] = countU[r][r+d-1] + A[r:r+d/2+1, r+d].size
            countL[r][r+d] = countL[r][r+d-1] + A[r+d, r+d:2*(r+d)-r+1].size
            Sx2[r][r+d] = Sx2[r][r+d-1] + sq[r:r+d/2+1, r+d].sum() + sq[r+d, r+d:2*(r+d)-r+1].sum()

    for d in range(1, row):
        for r in range(row):
            if r+d >= col or 2*(r+d) >= col:
                break
            Svar[r][r+d] = (Sx2[r][r+d] -Sx[r][r+d])/(countU[r][r+d]+countL[r][r+d])

    print "Finished Computing Svar..."
    return (Svar,countU,countL)


def getCornerScore(Ssign, Ssum, Svar):
    print "Computing Scorner..."
    Scorner =  Ssign + Ssum + Svar
    print "Finished Computing Scorner..."
    return Scorner

def getAllMat(A):
    Usgn = computeUsgn(A)
    Lsgn = computeLsgn(A)
    U = computeU(A)
    L = computeL(A)
    Ssign = computeSsign(Usgn, Lsgn)
    Ssum = computeSsum(U,L)
    (Svar, countU, countL) = computeSvar(A, U, L)
    Ssign = normalizeS(Ssign)
    Ssum = normalizeS(Ssum)
    Svar = normalizeS(Svar)
    Scorner = getCornerScore(Ssign, Ssum, Svar)
    return (Usgn, Lsgn, countU, countL, Svar, Scorner)

def main(args):
    np = readMat(args.input_data)
    if args.is_normal is 'n':
        normalize(np)
    A = computeArrowHead(np)
    #display(A)
    (Usgn, Lsgn, countU, countL, Svar, Scorner) = getAllMat(A)
    display(Scorner)


if __name__== "__main__":
    parser = argparse.ArgumentParser(description="Arrowhead")
    parser.add_argument("--input_data")
    parser.add_argument("--is_normal")
    parser.add_argument("--t1")
    parser.add_argument("--t2")
    args = parser.parse_args()
    main(args)
