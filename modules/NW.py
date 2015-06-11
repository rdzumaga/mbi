
#seqRef="GAATTC"
#seq="GATTA"
#penalty=-5

# Create an empty matrix
def createMatrix(rows, cols):
    return [[0]*cols for _ in xrange(rows)]

# Read BLOSUM table
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

def calcMatrix(seqVertical, seqHorizontal, blosum, penalty):
    rows = len(seqVertical)+1 
    cols = len(seqHorizontal)+1 
    
    F = createMatrix(rows, cols)
	
    for i in range(0, rows):
        F[i][0] = i * penalty
    for j in range(0, cols):
        F[0][j] = j * penalty
 
    for i in range(1, rows):
        for j in range(1, cols):
            match = F[i-1][j-1] + blosum[(seqVertical[i-1], seqHorizontal[j-1])]
            delete = F[i-1][j] + penalty
            insert = F[i][j-1] + penalty

            F[i][j] = max(match, delete, insert)
 
    return F

def calcMatrixStepByStep(seqVertical, seqHorizontal, blosum, penalty, step):
    rows = len(seqVertical)+1 
    cols = len(seqHorizontal)+1
    #rows=step/len(seqHorizontal)+1
    #cols=step%len(seqHorizontal)
    print "step=", step, "rows=", rows, "cols=", cols
    
    F = createMatrix(rows, cols)

    steps=[]
	
    for i in range(0, rows):
        F[i][0] = i * penalty
    for j in range(0, cols):
        F[0][j] = j * penalty
    
    counter=0
    for i in range(1, rows):
        for j in range(1, cols):
            counter+=1
            
            #check if reached the indicated step nr
            if counter<=step:
                #match or delete or insert =(value, row, col)
                possibilities=[]
                indexOfBest=0
                #match (diag)
                diag=[F[i-1][j-1] + blosum[(seqVertical[i-1], seqHorizontal[j-1])],i-1, j-1]
                possibilities.append(diag)
               
                #delete (up)
                up=[F[i-1][j] + penalty, i-j, j]
                possibilities.append(up)
                if up[0]>diag[0]:
                    indexOfBest=1

                #insert (left)
                left=[F[i][j-1] + penalty, i, j-1]
                possibilities.append(left)
                if left[0]>up[0] and left[0]>diag[0]:
                    indexOfBest=2

                best=possibilities[indexOfBest]
                steps.append([best[0], best[1], best[2],  possibilities[0][0], possibilities[1][0], possibilities[2][0]])
                F[i][j] = best[0]
            else:
                break
    
    return F, steps
    
def traceback(seq, seqRef, scoreMatrix, startPos, blosum,penalty):
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
    step         = nextStep(seq, seqRef, scoreMatrix, i, j,blosum, penalty)
    
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
       
        step = nextStep(seq, seqRef, scoreMatrix, i, j,blosum, penalty)
       
    return ''.join(reversed(alignedSeq)), ''.join(reversed(alignedSeqRef))


def nextStep(seq, seqRef, scoreMatrix, i, j, blosum, penalty):
    score=scoreMatrix[i][j]
    diag = scoreMatrix[i - 1][j - 1]
    up   = scoreMatrix[i - 1][j]
    left = scoreMatrix[i][j - 1]
    similarity=blosum[(seq[i-1], seqRef[j-1])]

    if(score==diag+similarity):
        return 1
    if (score==up+penalty):
        return 2
    if (score==left+penalty):
        return 3

    return 0

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

def needlemanWunsch(step=0, seq="GATTA", seqRef="GAATTC", penalty=-5):
    blosum=readBlosum("blosum.txt")
    matrix=[]
                                 
    if step>0 and step<len(seq)*len(seqRef)-1:
        matrix, steps=calcMatrixStepByStep(seq, seqRef, blosum, penalty, step)
        return matrix, steps, "", ""
        
    matrix =calcMatrix(seq, seqRef, blosum, penalty)
    rightBottomCell=(len(seq), len(seqRef))
    seqAligned, seqRefAligned = traceback(seq, seqRef, matrix, rightBottomCell, blosum, penalty)
    return matrix, [], seqAligned, seqRefAligned
    




