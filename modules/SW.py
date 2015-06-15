"""!@package SW
Compute local alignement of two DNA sequences
"""

import os
import re
import sys

# Create an empty matrix
def createMatrix(rows, cols):
    """
    initialize matrix with 0
    """
    return [[0]*cols for _ in xrange(rows)]
	
# Read the BLOSUM table
def readBlosum(fname):
    """!@brief	    Read from file the Blosum table
    """
	
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

def createScoreMatrix(seq, seqRef, penalty, match, mismatch):
    '''!@brief	Create a matrix and fill it with values representing possible alignments of two sequences

    Best alignment can be found by locating a path in the matrix with highest cumulative score.

    @param  seq	        String representing query DNA sequency
    @param  seqRef	DNA reference string, against which the query sequence will be compared
    @param  penalty	Gap penalty, must be negative or 0
    @param  match	The score to add for a match between a pair of nucleotides, must be >0
    @oparam mismatch	The score to add for a mismatch between a pair of nucleotides
	
    @retval             A Tuple; Score matrix and the (row,col)position of a cell with best score
                        
    '''
    
    rows=len(seq)+1
    cols=len(seqRef)+1

    #initialize the matrix with 0
    scoreMatrix = createMatrix(rows, cols)

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
    '''!@brief Calculate score for a given position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    similarity=match if seq[i-1]==seqRef[j-1] else mismatch

    diagScore = matrix[i - 1][j - 1] + similarity
    upScore   = matrix[i - 1][j] + penalty
    leftScore = matrix[i][j - 1] + penalty

    return max(0, diagScore, upScore, leftScore)

def createMatrixStepByStep(step, seq, seqRef, penalty, match, mismatch):
    """!@brief    Calculate the score matrix until idicated step is reached

    Parameters:
    @param  step       Number a steps to take when calculating the algorithm
    @param  seq	        String representing query DNA sequency
    @param  seqRef	DNA reference string, against which the query sequence will be compared
    @param  matrix	Score matrix
    @param  penalty	Gap penalty, must be negative or 0
    @param  match	The score to add for a match between a pair of nucleotides, must be >0
    @oparam mismatch	The score to add for a mismatch between a pair of nucleotides
	
    @retval             A Tuple; partially calulated score matrix, a list of steps.
                        A step consists of [nextBestPredecessorIndex, nextBestPredecessorRow, nextBestPredecessorCol, scoreFromDiagonal, scoreFromUp, ScoreFromLeft].
                        nextBestStepIndex takes on the following values
                            0 - diagonal predecessor
                            1 - left predecessor
                            2 - up predecessor
                            3 - no predecessor    
    """
    rows=len(seq)+1
    cols=len(seqRef)+1

    #initialize the matrix with 0
    scoreMatrix = createMatrix(rows, cols)
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

def calcStep(steps, seq, seqRef, matrix, i, j, penalty, match, mismatch):
    '''!@briefCalculate score for a given position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    best index can take on following values:
	0 - diag
	1 - left
	2 - up
	3 - zero
	
    Parameters:
    @param  steps       Number a steps to take when calculating the algorithm
    @param  seq	        String representing query DNA sequency
    @param  seqRef	DNA reference string, against which the query sequence will be compared
    @param  matrix	Score matrix
    @param  i           row index for cell from which a step will be taken
    @param  j           col index for cell from which a step will be taken
    @param  penalty	Gap penalty, must be negative or 0
    @param  match	The score to add for a match between a pair of nucleotides, must be >0
    @oparam mismatch	The score to add for a mismatch between a pair of nucleotides
	
    @retval             A Tuple; the list of steps with appended new step, score calculated for the nest step.
                        The appended step consists of [nextBestPredecessorIndex, nextBestPredecessorRow, nextBestPredecessorCol, scoreFromDiagonal, scoreFromUp, ScoreFromLeft]
    '''
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

    best=possibilities[bestIndex]
    steps.append([bestIndex, best[1], best[2],  possibilities[0][0], possibilities[1][0], possibilities[2][0]])

    return steps, possibilities[bestIndex][0]
	
def traceback(scoreMatrix, startPos, seq, seqRef, penalty, match, mismatch):
    '''!@briefFind the optimal path through the matrix representing the alignment.

    Starting from the best position (bottom right of a path), trace the whole path
    back up (top-left corner), finding thus the best local alignment. Each step of the path (matrix cell)
    corresponds to either a gap in a sequence (or both sequences) or a match/mismatch in the following way:
        diagonal (i-1, j-1) - match/mismatch
        up       (i-1, j  ) - gap in sequence 1
        left     (i  , j-1) - gap in sequence 2

    A step that should be taken is the one that leads to the predecessor cell
    Parameters:
    @param  scoreMatrix	Matrix with each cell scored
    @param  seq	        String representing query DNA sequency
    @param  seqRef	DNA reference string, against which the query sequence will be compared
    @param  penalty	Gap penalty, must be negative or 0
    @param  match	The score to add for a match between a pair of nucleotides, must be >0
    @oparam mismatch	The score to add for a mismatch between a pair of nucleotides
	
    @retval  A Tuple of strings: (aligned query sequence,aligned reference sequence)
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
    """!@brief	Calculate next step in the tracedback path
	
	Return:
		1 - diagonal step
		2 - up step
		3 - left step
		0 - reached the end of path
    """	
	
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

def createAlignmentString(alignedSeq, alignedSeqRef):
    '''!@brief	Construct a special string showing identities, gaps, and mismatches.

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

def printMatrix(matrix,seq, seqRef):
    """!@brief	Helper function for printing the score matrix and letters of both sequences on the console
    """
    
    #print ref sequence's nukleotides
    rows=len(seq)
    cols=len(seqRef)
    print
    print"         ",
    for n in range(cols):
        print seqRef[n], "  ",
    print

    #print a letter from query sequence and the corresponding values
    i=0
    for row in matrix:
        if(row>0):
            print seq[i-1],
        else:
            print "    ",
        i+=1
        for col in row:
            print ('{0:>4}'.format(col)),
        print

def SmithWaterman(step=-1, seq="GACTTAC", seqRef="CGTGAATTCAT", penalty=-4, match=5, mismatch=-3):
    """!@brief	Method calculating local alignment of two sequences.
    Parameters:
    @param  step	number of steps to take when calculating the score matrix. Each step move from cell F[i][j] to cell F[i][j+1] or F[i+1][0] when reached the end of a row. If step<0, the whole algorithm is executed, along with string alignements
    @param  seq	        string representing query DNA sequency
    @param  seqRef	DNA reference string, against which the query sequence will be compared
    @param  penalty	gap penalty, must be negative or 0
    @param  match	The score to add for a match between a pair of nucleotides, must be >0
    @oparam mismatch	The score to add for a mismatch between a pair of nucleotides
	
    @retval Array of steps taken, if step parameter>0, otherwise, computed score matrix (without the 'null' row and col). When returning an array of steps, one step consists of [nextBestPredecessorIndex, nextBestPredecessorRow, nextBestPredecessorCol, scoreFromDiagonal, scoreFromUp, ScoreFromLeft]
    """
	
    rows=len(seq)+1
    cols=len(seqRef)+1
    matrix=[]

    if step>=0 and step<len(seq)*len(seqRef)+1:
        matrix, steps=createMatrixStepByStep(step, seq, seqRef, penalty, match, mismatch)
        return steps

    matrix, bestPos = createScoreMatrix(seq, seqRef, penalty, match, mismatch)

    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seqAligned, seqRefAligned = traceback(matrix, bestPos, seq, seqRef, penalty, match, mismatch)
    assert len(seqAligned) == len(seqRefAligned), 'aligned strings are not the same size'
    return [row[1:] for row in matrix[1:]], seqAligned, seqRefAligned
