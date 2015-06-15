"""!@package FastaSW
Find the best match between query sequence and reference sequences from a database
"""

import argparse
import os
import re
import sys

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


def isWithinBounds(i,j, k, rowColPath)	:
	for row, col in rowColPath:
		if row==i:
			if j>col+k or j<col-k:
				return False
		if col==j:
			if i>row+k or i<row-k:
				return False

	return True

def createScoreMatrix(seq, seqRef, penalty, match, mismatch, path, k):
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
            score = calcScore(seq, seqRef, scoreMatrix, i, j, penalty, match, mismatch, k,path)
            if score > maxScore:
                maxScore = score
                bestPos   = (i, j)

            scoreMatrix[i][j] = score
    assert bestPos is not None, 'position with the highest score not found'
    return scoreMatrix, bestPos


def calcScore(seq, seqRef, matrix, i, j, penalty, match, mismatch, k, path):
    '''Calculate score for a given position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    similarity=match if seq[i-1]==seqRef[j-1] else mismatch

    #check bounds
    diagScore = matrix[i - 1][j - 1] + similarity if isWithinBounds(i-1,j-1, k, path) else -1
    upScore   = matrix[i - 1][j] + penalty if isWithinBounds(i-1,j, k, path) else -1
    leftScore = matrix[i][j - 1] + penalty if isWithinBounds(i,j-1, k, path) else -1

    return max(0, diagScore, upScore, leftScore)


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
    score=0
    END, DIAG, UP, LEFT = range(4)
    alignedSeq = []
    alignedSeqRef = []

    i, j         = startPos

    step,newscore   = nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch,score )

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

        step,score = nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch,newscore)

	alignedSeqRefStartIndex=j
    return ''.join(reversed(alignedSeq)), ''.join(reversed(alignedSeqRef)), score, alignedSeqRefStartIndex


def nextStep(scoreMatrix, i, j, seq, seqRef, penalty, match, mismatch,score ):
    if(i==0 or j==0):
        return 0, score

    val=scoreMatrix[i][j]
    diag = scoreMatrix[i - 1][j - 1]
    up   = scoreMatrix[i - 1][j]
    left = scoreMatrix[i][j - 1]


    similarity=match if seq[i-1]==seqRef[j-1] else mismatch

    if(val==diag+similarity):
        return 1, score+scoreMatrix[i][j]
    if (val==up+penalty):
        return 2, score+scoreMatrix[i][j]
    if (val==left+penalty):
        return 3, score+scoreMatrix[i][j]

    return 0, score+scoreMatrix[i][j]


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


def printMatrix(matrix,seq, seqRef):
    #print ref sequence's nukleotides
    f=open('SW_Debug.txt', 'w+')
    rows=len(seq)
    cols=len(seqRef)
    print >>f
    print >>f, "         ",
    for n in range(cols):
        print >>f, seqRef[n], "  ",
    print >>f
    i=0
    for row in matrix:
        if(row>0):
            print >>f,  seq[i-1],
        else:
            print >>f, "    ",
        i+=1
        for col in row:
            #if col!=0:
            print >>f, ('{0:>4}'.format(col)),
            #else:
               # print "    ",
        print >>f


def SmithWaterman(seq, seqRef, path, k, penalty=-5, match=1, mismatch=-1):
    """
    !@brief	A modified version of Smith-Waterman algorithm for use with FastA
    """

    rows=len(seq)+1
    cols=len(seqRef)+1
    matrix, bestPos = createScoreMatrix(seq, seqRef, penalty, match, mismatch, path, k)

    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seqAligned, seqRefAligned, score, alignedSeqRefStartIndex = traceback(matrix, bestPos, seq, seqRef, penalty, match, mismatch)
    #printMatrix(matrix, seq, seqRef)#to File
    assert len(seqAligned) == len(seqRefAligned), 'aligned strings are not the same size'
    return matrix, seqAligned, seqRefAligned,score, alignedSeqRefStartIndex



