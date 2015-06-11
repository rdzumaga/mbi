import argparse
import os
import re
import sys
import unittest


# Create an empty matrix
def create_matrix(m, n):
    return [[0]*n for _ in xrange(m)]

# Read the BLOSUM50 table
def readBlosum(fname):

    #dictionary holding pairs of nucleotides and their BLOSUM matrix value, e.g.: ('T', 'A'): -4
    d = {}
    lines = open(fname, "rt").readlines()
    alpha = lines[0].rstrip('\n\r').split()
    assert(len(alpha) == len(lines)-1)
    for r in lines[1:]:
        r = r.rstrip('\n\r').split()
        a1 = r[0]
        for a2, score in zip(alpha, r[1:]):
            d[(a1, a2)] = int(score)
    return d


     






seq02= "GAAAGAT" #horizontal
seq01 = "GATGAA"#vertical
#seqRef= "GGCTCAATCA"
#seq= "ACCTAAGG"
#seq = 'AGCACACA'
#seqRef = 'ACACACTA'
seq02 = "FTFTALILLAVAV"
seq01 = "FTALLLAAV"


seq02= "CGTGAATTCAT" #horizontal
seq01 = "GACTTAC"#vertical

#---------------functions--------------------

def createScoreMatrix(seq, seqRef, penalty, match, mismatch):
    '''Create a matrix and fill it with values representing possible alignments of two sequences

    Best alignment can be found by locating a path in the matrix (when represenet as a 2D graph)
    with highest cumulative score. 
    '''
    rows=len(seq)+1
    cols=len(seqRef)+1
    
    #initialize the matrix with 0
    scoreMatrix = [[0 for col in range(cols)] for row in range(rows)]
   
    maxScore = 0
    bestPos   = None    # i and j index for matrix cell with highest score
    
    # Fill the scoring matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calcScore(seq, seqRef, scoreMatrix, i, j, penalty, match, mismatch)
            if score > maxScore:
                maxScore = score
                bestPos   = (i, j)

            scoreMatrix[i][j] = score

    assert bestPos is not None, 'position with the highest score not found'
    return scoreMatrix, bestPos


def calcScore(seq, seqRef, matrix, i, j, penalty, match, mismatch):
    '''Calculate score for a given position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    similarity=match if seq[i-1]==seqRef[j-1] else mismatch
  
    diagScore = matrix[i - 1][j - 1] + similarity
    upScore   = matrix[i - 1][j] + penalty
    leftScore = matrix[i][j - 1] + penalty

    return max(0, diagScore, upScore, leftScore)


def calcStep(steps, seq, seqRef, matrix, i, j, penalty, match, mismatch):
    '''Calculate score for a given position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    #print "steps: ", steps
    nextStep=[]
    possibilities=[]
    bestIndex=0
    similarity=match if seq[i-1]==seqRef[j-1] else mismatch
  
    diag = [matrix[i - 1][j - 1] + similarity, i-1, j-1]
    possibilities.append(diag)
	
    up = [matrix[i - 1][j] + penalty, i-1, j]
    possibilities.append(up)
    if(up[0]>diag[0]):
	bestIndex=1
		
    left = [matrix[i][j - 1] + penalty, i, j-1]
    possibilities.append(left)
	
    zero =[0, -1, -1]
    possibilities.append(zero)
	
    if left[0]>up[0] and left[0]>diag[0]:
        bestIndex=2
		
    if possibilities[bestIndex][0]<0:
        
	bestIndex=3

    #print "possi[best]=",bestIndex, " ___________", diag, up, left	
    best=possibilities[bestIndex]
    steps.append([best[0], best[1], best[2],  possibilities[0][0], possibilities[1][0], possibilities[2][0]])

    return steps, possibilities[bestIndex][0]

def calcMatrixStepByStep(step, seq, seqRef, penalty, match, mismatch):
    rows=len(seq)+1
    cols=len(seqRef)+1
    
    #initialize the matrix with 0
    scoreMatrix = [[0 for col in range(cols)] for row in range(rows)]
    steps=[]
	
    maxScore = 0
    bestPos   = None    # i and j index for matrix cell with highest score
    
    counter=0
    # Fill the scoring matrix.
    for i in range(1, rows):
        for j in range(1, cols):
	    counter+=1
			
	    #check if reached the indicated step nr
            if counter<=step:
                steps, score = calcStep(steps,seq, seqRef, scoreMatrix, i, j, penalty, match, mismatch)
                if score > maxScore:
                    maxScore = score
                    bestPos   = (i, j)
            else:
                break
            
	    scoreMatrix[i][j] = score
	    

    assert bestPos is not None, 'position with the highest score not found'
    return scoreMatrix, steps
	
	
def traceback(scoreMatrix, startPos, seq, seqRef, penalty, match, mismatch):
    '''Find the optimal path through the matrix representing the alignment.

    Starting from the best position (bottom right of a path), trace the whole path
    back up (top-left corner), finding thus the best local alignment. Each step of the path (matrix cell)
    corresponds to either a gap in a sequence (or both sequences) or a match/mismatch in the following way:
        diagonal (i-1, j-1) - match/mismatch
        up       (i-1, j  ) - gap in sequence 1
        left     (i  , j-1) - gap in sequence 2

    A step that should be taken is the one that leads to the predecessor cell
    '''

    END, DIAG, UP, LEFT = range(4)
    alignedSeq = []
    alignedSeqRef = []
    i, j         = startPos
    step         = nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch)
    
    while step != END:
        if step == DIAG:
            alignedSeq.append(seq[i - 1])
            alignedSeqRef.append(seqRef[j - 1])
            i -= 1
            j -= 1
        elif step == UP:
            alignedSeq.append(seq[i - 1])
            alignedSeqRef.append('-')
            i -= 1
        else:
            alignedSeq.append('-')
            alignedSeqRef.append(seqRef[j - 1])
            j -= 1
       
        step = nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch)
       
    return ''.join(reversed(alignedSeq)), ''.join(reversed(alignedSeqRef))


def nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch):
    if(i==0 or j==0):
        return 0

    score=scoreMatrix[i][j]
    diag = scoreMatrix[i - 1][j - 1]
    up   = scoreMatrix[i - 1][j]
    left = scoreMatrix[i][j - 1]

    similarity=match if seq[i-1]==seqRef[j-1] else mismatch

    if(score==diag+similarity):
        return 1
    if (score==up+penalty):
        return 2
    if (score==left+penalty):
        return 3

    return 0

"""    if diag >= up and diag >= left:     # Tie - DIAG step "wins".
        return 1 if diag != 0 else 0    # 1 signals a DIAG step. 0 signals the end.
    elif up > diag and up >= left:      # Tie -  UP step"wins".
        return 2 if up != 0 else 0      # UP step or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT step or end.
    else:
        raise ValueError('invalid move during traceback')"""


def createAlignmentString(alignedSeq, alignedSeqRef):
    '''Construct a special string showing identities, gaps, and mismatches.

    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.

    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignmentString = []
    for base1, base2 in zip(alignedSeq, alignedSeqRef):
        if base1 == base2:
            alignmentString.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignmentString.append(' ')
            gaps += 1
        else:
            alignmentString.append(':')
            mismatches += 1

    return ''.join(alignmentString), idents, gaps, mismatches


def print_matrix(matrix):
    '''Print the scoring matrix.

    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    for row in matrix:
        for col in row:
            print('{0:>4}'.format(col)),
        print


class ScoreMatrixTest(unittest.TestCase):
    '''Compare the matrix produced by create_score_matrix() with a known matrix.'''
    def test_matrix(self):
        # From Wikipedia (en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
        #                -   A   C   A   C   A   C   T   A
        knownMatrix = [[0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
                        [0,  2,  1,  2,  1,  2,  1,  0,  2],  # A
                        [0,  1,  1,  1,  1,  1,  1,  0,  1],  # G
                        [0,  0,  3,  2,  3,  2,  3,  2,  1],  # C
                        [0,  2,  2,  5,  4,  5,  4,  3,  4],  # A
                        [0,  1,  4,  4,  7,  6,  7,  6,  5],  # C
                        [0,  2,  3,  6,  6,  9,  8,  7,  8],  # A
                        [0,  1,  4,  5,  8,  8, 11, 10,  9],  # C
                        [0,  2,  3,  6,  7, 10, 10, 10, 12]]  # A

        global seq, seqRef
        seq = 'AGCACACA'
        seqRef = 'ACACACTA'
        rows = len(seq) + 1
        cols = len(seqRef) + 1

        matrixToTestest, bestPos = createScoreMatrix(rows, cols)
        self.assertEqual(knownMatrix, matrixToTest)

		
def printMatrix(matrix,seq, seqRef):
    #print ref sequence's nukleotides
    rows=len(seq)
    cols=len(seqRef)
    print
    print"         ",
    for n in range(cols):
        print seqRef[n], "  ",
    print
    i=0
    for row in matrix:
        if(row>0):
            print seq[i-1],
        else:
            print "    ",
        i+=1
        for col in row:
            #if col!=0:
            print ('{0:>4}'.format(col)),
            #else:
               # print "    ",
        print
		
#-4,5, -3	
seq02= "CGTGAATTCAT" #horizontal
seq01 = "GACTTAC"#vertical
def SmithWaterman(step=0, seq="GACTTAC", seqRef="CGTGAATTCAT", penalty=-4, match=5, mismatch=-3):
    rows=len(seq)+1
    cols=len(seqRef)+1
    matrix=[]
    
    if step>0 and step<len(seq)*len(seqRef)-1:
        matrix, steps=calcMatrixStepByStep(step, seq, seqRef, penalty, match, mismatch)
        printMatrix(matrix, seq, seqRef)
        return steps
	
    scoreMatrix, bestPos = createScoreMatrix(seq, seqRef, penalty, match, mismatch)
    
    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seqAligned, seqRefAligned = traceback(scoreMatrix, bestPos, seq, seqRef, penalty, match, mismatch)
    assert len(seqAligned) == len(seqRefAligned), 'aligned strings are not the same size'
    return [row[1:] for row in matrix[1:]], seqAligned, seqRefAligned





""""
#------------------------------------
try:
    t=1
    #parse_cmd_line()
except ValueError as err:
    print('error:', err)


# The scoring matrix contains an extra row and column for the gap (-), hence
# the +1 here.
rows = len(seq) + 1
cols = len(seqRef) + 1

# Initialize the scoring matrix.
scoreMatrix, bestPos = createScoreMatrix(rows, cols)
print "----------------------------------------"
print_matrix(scoreMatrix)

# Traceback. Find the optimal path through the scoring matrix. This path
# corresponds to the optimal local sequence alignment.
seqAligned, seqRefAligned = traceback(scoreMatrix, bestPos)
assert len(seqAligned) == len(seqRefAligned), 'aligned strings are not the same size'

# Pretty print the results. The printing follows the format of BLAST results
# as closely as possible.
alignmentStr, idents, gaps, mismatches = createAlignmentString(seqAligned, seqRefAligned)
alength = len(seqAligned)


print "***************stats**************"
print "seq=", seq, "rows=len(seq), seq is vertical and string in"
print "seqRef=", seqRef,"cols=len(seqRef), seqRef is horizontal and string ref"
print "seqAligned=",seqAligned
print "seqRefAligned=",seqRefAligned
print "*********************************"

print
print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
      alength, idents / alength, gaps, alength, gaps / alength))
print
for i in range(0, alength, 60):
    seqSlice = seqAligned[i:i+60]
    print('Query      {0:<4}  {1}  {2:<4}'.format(i + 1, seqSlice, i + len(seqSlice)))
    print('                 {0}'.format(alignmentStr[i:i+60]))
    seqRefAlignedSlice = seqRefAligned[i:i+60]
    print('Reference  {0:<4}  {1}  {2:<4}'.format(i + 1, seqRefAlignedSlice, i + len(seqRefAlignedSlice)))
    print()

strG = str(raw_input("Genome sequence: "))
"""

