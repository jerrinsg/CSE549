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

def display(np, filename):
    print "About to display image"
    # Write image to file
    im = Image.fromarray(numpy.uint8(plt.cm.gist_earth(np)*255))
    im.save(filename)

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
    for r in range(row):
        Usgn[r][r]=A[r][r]

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
    for r in range(row):
        Lsgn[r][r]=A[r][r]

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
    for r in range(row):
        U[r][r]=A[r][r]

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
    for r in range(row):
        L[r][r]=A[r][r]
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
            if r+d >= col :
                continue
            countU[r][r+d] = countU[r][r+d-1] + A[r:r+d/2+1, r+d].size
            if 2*(r+d)-r >= col:
                continue
            countL[r][r+d] = countL[r][r+d-1] + A[r+d, r+d:2*(r+d)-r+1].size
            Sx2[r][r+d] = Sx2[r][r+d-1] + sq[r:r+d/2+1, r+d].sum() + sq[r+d, r+d:2*(r+d)-r+1].sum()

    for d in range(1, row):
        for r in range(row):
            if r+d >= col or 2*(r+d) >= col:
                break
            Svar[r][r+d] = (Sx2[r][r+d])/(countU[r][r+d]+countL[r][r+d]) -(Sx[r][r+d]/(countU[r][r+d]+countL[r][r+d]))**2

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
    Svarn = normalizeS(Svar)
    Scorner = getCornerScore(Ssign, Ssum, Svarn)

    newFn = numpy.vectorize(fn)
    NewCountU = newFn(countU)
    NewCountL = newFn(countL)
    MeanSgnU = Usgn/NewCountU;
    MeanSgnL = Lsgn/NewCountL;


    return (Usgn, Lsgn, countU, countL, Svar, Svarn, Scorner, MeanSgnU, MeanSgnL)

def main(args):
    np = readMat(args.input_data)

    apply_threshold1 = args.apply_threshold1
    if apply_threshold1 is 'y':
        t1=float(args.t1)
        t2=float(args.t2)
        t3=float(args.t3)

    apply_threshold2 = args.apply_threshold2
    if apply_threshold2 is 'y':
        t4=float(args.t4)
        t5=float(args.t5)

    if args.is_normal is 'n':
        normalize(np)

    A = computeArrowHead(np)

    display(A,"A.jpg")
    (Usgn, Lsgn, countU, countL, Svar, Svarn, Scorner, MeanSgnU, MeanSgnL) = getAllMat(A)

    if apply_threshold1 is 'y':
        print "First stage filtering begin :"
        for x, y in numpy.ndindex(Scorner.shape):
            if Svar[x][y]<t1 and MeanSgnU[x][y]<-t2 and MeanSgnL[x][y]>t3:
                continue 
            else:
                Scorner[x][y]=0

        print "First stage filtering end."

    if apply_threshold2 is 'y':
        print "Second stage filtering begin :"
        for x, y in numpy.ndindex(Scorner.shape):
            if MeanSgnU[x][y]<-t4 and MeanSgnL[x][y]>t5:
                continue
            else:
                Scorner[x][y]=0
        print "Second stage filtering end."

    display(Scorner,"Scorner.jpg")
    numpy.savetxt('ScornerData',Scorner)


if __name__== "__main__":
    parser = argparse.ArgumentParser(description="Arrowhead")
    parser.add_argument("--input_data")
    parser.add_argument("--is_normal")
    parser.add_argument("--t1")
    parser.add_argument("--t2")
    parser.add_argument("--t3")
    parser.add_argument("--t4")
    parser.add_argument("--t5")
    parser.add_argument("--apply_threshold1")
    parser.add_argument("--apply_threshold2")
    args = parser.parse_args()
    main(args)
