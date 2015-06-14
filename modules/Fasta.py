import graph
import copy
import FastaSW
seqRef1="GAATTC"
seq1="GATTA"

seqRefk="AACACTTTTCAAT"
seqk="ACTTATCA"


#seq="GCATCGGC"
#seq="ACATCTCCAGCG"
seqRefa="CCATCGCCATCGG"
seqa="CTCGCACATGG"

seqRef="ACCTGTTAAC"
seq="GCCTTGTAA"



#ktup lenght; should be 4,5 or 6 for DNA, preferably 6
k=2
# some value to filter out diagonals with lowest scores (below cutoff means low score)
cutoff=10



def getTuplesList(str):
	tuples=set()
	tuplesDict={}

	for i in range(0, len(str)-k+1):
		tuple=str[i:i+k]

		if tuple not in tuples:
			tuples.add(tuple)
			tuplesDict[tuple]=[]

		tuplesDict[tuple].append(i)
	return tuples,tuplesDict

def createDiagonalDict0():
	#TODO: verify that
	#diagonalNum=rows+cols-2*k+1

	diagonalIndexDict={}
	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=0
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=0
	return diagonalIndexDict

def createDiagonalDict(seq, seqRef):
	#TODO: verify that
	#diagonalNum=rows+cols-2*k+1

	rows=len(seq)
	cols=len(seqRef)
	diagonalIndexDict={}

	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=[graph.DiagonalRun(row, str(row))]
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=[graph.DiagonalRun(-col, str(-col))]
	return diagonalIndexDict

def createDiagonalDictFrom(bestDiagonals):
	newDict=dict(bestDiagonals)
	for diag in bestDiagonals:
		newDict[diag]=0
	return newDict

def createHotspotsDict():
	#TODO: verify that
	#diagonalNum=rows+cols-2*k+1

	diagonalIndexDict={}
	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=[]
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=[]
	return diagonalIndexDict
def initDiagonals(seq, seqRef):
	diags=[]
	rows=len(seq)
	cols=len(seqRef)

	for row in range(0, rows-k+1):
		diagonal=[graph.DiagonalRun(row, str(row)+".0")]
		diags.append(diagonal)
	for col in range(1, cols-k+1):
		diagonal=[graph.DiagonalRun(-col, str(row)+".0")]
		diags.append(diagonal)
	return diags

def calcDiagonalSums0(tuplesRef, tuplesRefDict):
	diagonalSums=createDiagonalDict()
	hotspotRows=createHotspotsDict()
	#scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]

	for row in range(0, rows-k+1):
		tuple=seq[row:row+k]
		if tuple in tuplesRef:
			for col in tuplesRefDict[tuple]:
				offset=row-col

				diagonalSums[offset]+=1
				hotspotRows[offset].append(row)

	#get rid of entries with sum==0
	answer=dict(diagonalSums)
	for diag in diagonalSums:
		if diagonalSums[diag]==0:
			del answer[diag]
			del hotspotRows[diag]

	print"-------------------"
	print "diagonalSums="
	print answer
	print "-------------------"
	print "hotspotRows"
	print hotspotRows
	print"-----------------------"
	print "tuplesRefDict (hotspotCols):"
	print tuplesRefDict
	return answer, hotspotRows

def calcDiagonalSums(seq, seqRef,tuplesRef, tuplesRefDict):
	#diagonalSums=createDiagonalDict0()
	diagonalDict=createDiagonalDict(seq, seqRef)

	#hotspotRows=createHotspotsDict()
	#scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]
	rows=len(seq)
	cols=len(seqRef)

	for row in range(0, rows-k+1):
		tuple=seq[row:row+k]
		if tuple in tuplesRef:
			for col in tuplesRefDict[tuple]:
				offset=row-col
				diagonalDict[offset][0].add(row, col)

	#get rid of entries with sum==0
	goodDiagonalDict=dict(diagonalDict)
	for diag in diagonalDict:
		if diagonalDict[diag][0].value==0:
			del goodDiagonalDict[diag]

	return goodDiagonalDict

def scoreDiagonals0(diagonalSumDict, hotspotRows):
	"""
	In order to evaluate each diagonal run, FASTA gives each hot spot a positive score,
	and the space between consecutive hot spots in a run is given a negative score that
	decreases with the increasing distance. The score of the diagonal run is the sum of the
	hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
	runs under this evaluating scheme.

	"""
	rescoredDiagonals=dict(diagonalSumDict)
	reward=20 #this should be some kind of positive value

	for diag in diagonalSumDict:
		print "diag=",diag

		firstRow=0 if diag<=0 else diag
		gapPenalty=-reward/2 #this should be some kind of negative value
		interspotPenaltySum=0
		sum=0

		#for each row, score the cell and update diagonalSums
		for row in range(firstRow,rows-k+1):
			col=row-diag

			# check if col within boundaries
			if(col<cols):
				hotspot=row in hotspotRows[diag]

				#if for this col there is a gap between two hotspots:
				#(sum>0 indicates at least 1 hotspot was found before)
				if(hotspot==False and sum>0):
					#sum penalties for each consecutive gap
					interspotPenaltySum+=gapPenalty

					#make sure penalty stays a negativ value (in case there were too many consecutive gaps)
					if(gapPenalty<0):
						#each consecutive gap has less of a penalty
						gapPenalty+=1

				elif hotspot:
					#sum the reward and add the penalty for the gaps between this hotspot and the previous one
					sum=sum+reward+interspotPenaltySum

					#since a hotspot was found, reset the interspotPenaltySum and revert to the initial gapPenalty value
					interspotPenaltySum=0
					gapPenalty=-reward/2
				#print "row=", row, "sum=", sum, "gapPenalty=", interspotPenaltySum

		rescoredDiagonals[diag]=sum
		sum=0

	if(len(rescoredDiagonals)>10):
		rescoredDiagonals=getTopDiagonals(getTopDiagonals)
	return rescoredDiagonals

def scoreDiagonals(diagonalDict, seq, seqRef):
	"""
	In order to evaluate each diagonal run, FASTA gives each hot spot a positive score,
	and the space between consecutive hot spots in a run is given a negative score that
	decreases with the increasing distance. The score of the diagonal run is the sum of the
	hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
	runs under this evaluating scheme.

	"""
	rows=len(seq)
	cols=len(seqRef)
	rescoredDiagonals=dict(diagonalDict)
	reward=20 #this should be some kind of positive value

	#group all consecutive hotspots into regions
	for diag in diagonalDict:
		regionNr=1
		regions=[]
		hspots=sorted(diagonalDict[diag][0].hotspots)
		region=graph.DiagonalRun(diag,str(diag)+"."+str(regionNr))
		rescoredDiagonals[diag]=[]
		for i in range(len(hspots)):
			if(i>0):
				if hspots[i][0]!=hspots[i-1][0]+1:
					rescoredDiagonals[diag].append(copy.copy(region))
					regionNr+=1
					region=graph.DiagonalRun(diag,str(diag)+"."+str(regionNr))
					regionNr+=1

			region.add(hspots[i][0], hspots[i][1])
		rescoredDiagonals[diag].append(copy.copy(region))


	#TODO check if regions would have a better score when combined with others, inspite of gaps

		""""firstRow=0 if diag<=0 else diag
		gapPenalty=-reward/2 #this should be some kind of negative value
		interspotPenaltySum=0
		sum=0

		#for each row, score the cell and update diagonalSums
		for row in range(firstRow,rows-k+1):
			regionNr=1
			col=row-diag
			newDiagRun=graph.DiagonalRun(diag, str(diag)+"."+str(regionNr))

			# check if col within boundaries
			if(col<cols):
				hotspot=row in diagonalDict[diag][0].hotspots

				#if for this col there is a gap between two hotspots:
				#(sum>0 indicates at least 1 hotspot was found before)
				if(hotspot==False and sum>0):
					recoredDiagonals[diag][0]=copy.copy(newDiagRun)
					newDiagRun.hotspot.append(row,col)
					#sum penalties for each consecutive gap
					interspotPenaltySum+=gapPenalty

					#make sure penalty stays a negativ value (in case there were too many consecutive gaps)
					if(gapPenalty<0):
						#each consecutive gap has less of a penalty
						gapPenalty+=1

				elif hotspot:
					#sum the reward and add the penalty for the gaps between this hotspot and the previous one
					newDiagRun.value=newDiagRun.value+reward+interspotPenaltySum
					newDiagRun.hotspots.append(row,col)

					sumWithGaps=sum+reward+interspotPenaltySum
					if(sumWithGaps>sum):
						sum=sumWithGaps
					sum=sum+reward+interspotPenaltySum

					#since a hotspot was found, reset the interspotPenaltySum and revert to the initial gapPenalty value
					interspotPenaltySum=0
					gapPenalty=-reward/2
				#print "row=", row, "sum=", sum, "gapPenalty=", interspotPenaltySum

			regionNr+=1
		rescoredDiagonals[diag][0].value=sum
		sum=0"""

	if(len(rescoredDiagonals)>10):
		rescoredDiagonals=getTopDiagonals(rescoredDiagonals)
	return rescoredDiagonals

def getTopDiagonals(diagonals):
	regions=listAllRegions(diagonals)

	#get updated diagonals dictionary
	bestDiagonals=dict(diagonals)

	#get keys for the top ten diagonal sums	sort(key=lambda x: x.count, reverse=True)
	regions.sort(key=lambda x: x.value, reverse=True)
	topRegions=regions[0:11]

	for diagNum in diagonals:
		for region in diagonals[diagNum]:
			if region not in topRegions:
				bestDiagonals[diagNum].remove(region)


	return bestDiagonals

def rescoreDiagonals0(blosum,bestDiagonals,hotspotRows):

	rescoredDiagonals=createDiagonalDictFrom(bestDiagonals)
	print rescoredDiagonals
	#iterate over diagonals and score the with blosum matrix
	for diag in bestDiagonals:
		for row in hotspotRows[diag]:
			col=row-diag
			print "row=", row, "col=", col, "diag=", diag
			#print "blosum(seq[", row,"], seqRef[", col,"]= blosum[", seq[row], seq[col],"]"
			rescoredDiagonals[diag]+=blosum[(seq[row], seqRef[col])]


	#remove diagonals with scores below a cutoff threshold
	bestRescoredDiagonals=dict(rescoredDiagonals)

	print rescoredDiagonals

	for diag in rescoredDiagonals:
		if rescoredDiagonals[diag]<cutoff:
			del bestRescoredDiagonals[diag]

	print
	print bestRescoredDiagonals
	return bestRescoredDiagonals

def rescoreDiagonals(seq, seqRef, blosum, bestDiagonals):
	rescoredDiagonals=createDiagonalDictFrom(bestDiagonals)
	#iterate over diagonals and score the with blosum matrix
	for diag in bestDiagonals:
		rescoredRegions=[]
		for region in bestDiagonals[diag]:
			region.value=0
			for row,col in region.hotspots:
				region.value+=blosum[(seq[row], seqRef[col])]
				#print (row,col), region.value,
			rescoredRegions.append(region)
			#print "//end reg"
		rescoredDiagonals[diag]=rescoredRegions
	#print
	#remove diagonals with scores below a cutoff threshold
	bestRescoredDiagonals=dict(rescoredDiagonals)

	#for diag in bestRescoredDiagonals:
		#print diag,
		#for reg in bestRescoredDiagonals[diag]:
			#print reg.hotspots, reg.value,
		#print


	for diag in bestDiagonals:
		for region in bestDiagonals[diag]:
			if region.value<cutoff:
				if len(rescoredDiagonals[diag])>1:
					bestRescoredDiagonals[diag].remove(region)
				else:
					del bestRescoredDiagonals[diag]

	#print "AFTER cutoff"
	#for diag in bestRescoredDiagonals:
		#print diag,
		#for reg in bestRescoredDiagonals[diag]:
			#print reg.hotspots, reg.value,
		#print
	return bestRescoredDiagonals

def listAllRegions(diagonalRegionsDict):
	regions=[]
	for diag in diagonalRegionsDict:
		#diagonalRegionsDict[diag].printIt()
		for region in diagonalRegionsDict[diag]:
			regions.append(region)
	print
	return regions

def createMatrixForDots0(diagonalSumDict, hotspotRows):
	"""
	just for debugging. drawing dot matrix is unnecessary since it takes too much memory to remember all values
	"""
	scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]
	for diag in diagonalSumDict:
		firstRow=0 if diag<=0 else diag
		for row in range(firstRow,rows-k+1):
			col=row-diag
			if(col<cols):
				#print "row=", row, "col=", col
				if(row in hotspotRows[diag]):
					scoreMatrix[row][col]=1
				else:
					scoreMatrix[row][col]=0
	return scoreMatrix

def createMatrixForDots(diagonalDict, seq, seqRef):
	"""
	just for debugging. drawing dot matrix is unnecessary since it takes too much memory to remember all values
	"""
	rows=len(seq)
	cols=len(seqRef)
	scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]

	for diag in diagonalDict:
		firstRow=0 if diag<=0 else diag
		for row in range(firstRow,rows-k+1):
			col=row-diag
			if(col<cols):
				#print "row=", row, "col=", col
				if (row,col) in diagonalDict[diag][0].hotspots:
					scoreMatrix[row][col]=1
				else:
					scoreMatrix[row][col]=0
	return scoreMatrix

def findRegions(diag, scoreMatrix):
	regionRows=[]
	regions=[]
	offset=diag
	firstRow=0 if diag<0 else offset
	print "Analizing diag=", diag
	sum=0
	print "rang=", rows-k+1
	for row in range(firstRow, rows-k+1):
		col=row-offset
		#print "rangeCol=",  col
		if(col<=cols-k+1):
			print "row=",row, "offset=",offset,
			if(scoreMatrix[row][col]!=0):
				sum+=scoreMatrix[row][col]
				regionRows.append(row)
				if(row==rows-k and sum>0):
					regions.append((sum, regionRows[:]))
			else:
				if(sum>0):
					regions.append((sum, regionRows[:]))
				sum=0
				regionRows[:]=[]
			print "col=",col, "val=",scoreMatrix[row][col],"sum=",sum


	return sorted(regions, key=lambda tup: tup[0], reverse=True)

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

def printDotMatrix(seq, seqRef, matrix):
	#print ref sequence's nukleotides
	rows=len(seq)
	cols=len(seqRef)
	print
	print" ",
	for n in range(cols):
		print seqRef[n],
	print
	#print sequence's nukleotides and "o" if a dot should be printed
	for i in range(0,rows):
		print seq[i],
		for j in range(0,cols):
			if matrix[i][j]==1:
				print "o",
			else:
				print " ",
		print

def printMatrix(seq, seqRef, matrix):
	#print ref sequence's nukleotides
	rows=len(seq)
	cols=len(seqRef)
	print
	print"        ",
	for n in range(cols):
		print seqRef[n], "  ",
	print
	#print sequence's nukleotides and "o" if a dot should be printed
	i=0
	for row in matrix:
		if i<len(seq):
			print seq[i],
		i+=1
		for col in row:
			if col!=0:
				print ('{0:>4}'.format(col)),
			else:
				print "    ",
		print

def findPathBetween(startPos, destRegion):
		path=[]
		dest=destRegion.hotspots[0]
		next=(startPos[0]+1, startPos[1]+1)
		while(next[0]<dest[0] and next[1]< dest[1]):
			path.append(next)
			next=(next[0]+1, next[1]+1)

		if (next!=dest):
			di=0
			dj=0
			if next[0]==dest[0]:
				next=(next[0]-1, next[1])
				dj=1
			elif next[1]==dest[1]:
				next=(next[0], next[1]-1)
				di=1

			while(next[0]!=dest[0] and next[1]!=dest[1]):
				path.append(next)
				next=(next[0]+di, next[1]+dj)

		return path

def findPathWithHoles(regionsPath):

	if len(regionsPath)==0:
		return None
	startRegion=regionsPath[0]
	path=startRegion.hotspots
	for region in regionsPath[1:]:
		path+=findPathBetween(path[-1], region)
		path+=region.hotspots

	return path

def calcE(seq, seqRef,blosum="", k=2):
	if len(blosum)==0:
		blosum=readBlosum("blosum.txt")

	rows=len(seq)
	cols=len(seqRef)

	# 1.identify common k-words between I (seq) and J(seqRef)
	#print "\n=======================================\n                STEP1\n=======================================\n"
	tuplesRef,tuplesRefDict=getTuplesList(seqRef)
	diagonalDict=calcDiagonalSums(seq, seqRef, tuplesRef,tuplesRefDict)

	matrix=createMatrixForDots(diagonalDict,seq, seqRef)#
	#printDotMatrix(seq, seqRef, matrix)#

	# 2a. Score diagonals with k-word matches and identify 10 best diagonals
	#print "\n=======================================\n                STEP2\n=======================================\n"
	bestTenDiagonals=scoreDiagonals(diagonalDict,seq, seqRef)


	# 3. Rescore initial regions with a substitution score matrix and get best 10 subregions
	#print "\n=======================================\n                STEP3\n=======================================\n"
	#rescoredDiagonals=rescoreDiagonals(blosum, betsTenDiagonals, hotspotRows)
	diagonalRegionsDict=rescoreDiagonals(seq, seqRef, blosum, bestTenDiagonals)

	if len(diagonalRegionsDict)>0:

		#4. Join initial regions using gaps, penalise for gaps
		#print "\n=======================================\n                STEP4\n=======================================\n"

		mygraph=graph.createGraph(listAllRegions(diagonalRegionsDict))
		if mygraph==None:
			return "ERROR","ERROR", 0

		path=[]
		value=0
		allPaths=[]
		for startNode in mygraph:
			#print "ooooooooooooooooooooooooooooooooooooooooo"
			for node in mygraph[startNode]:
				path, value=graph.findBestPath(mygraph, startNode, node[0])
				#print value,
				if path:
					allPaths.append([value,path])
					#for n in path:
						#print n, "->",
				#print

		allPaths.sort(reverse=True)
		bestPath=allPaths[0]
		# 5. Perform dynamic programming to find final alignments
		print "\n=======================================\n                STEP5\n=======================================\n"
		#list all cells on the path starting with first diagonal run and ending with the last one from the path found before
		rowColPath=findPathWithHoles(bestPath[1])
		if rowColPath==None:
			return "ERROR","ERROR", 0

		#print "Final path (before NW):"
		#print rowColPath

		matrix, seqAligned, seqRefAligned,score=FastaSW.SmithWaterman(seq, seqRef, rowColPath, k)


		print
		print seqAligned
		print seqRefAligned
		print score
		#printMatrix(seq, seqRef, matrix)

		return seqAligned, seqRefAligned, score

def readDb(fname):
	lines = open(fname, "rt").readlines()
	db=[]
	for line in lines:
		#get rid of CRLF
		db.append(line[:-1])
		#break
	return db

def fasta(seq, k=2, db=[], blosum={} ):
	if len(db)==0:
		db=readDb("applications/mbi/modules/db.txt")

	if len(blosum)==0:
		blosum=readBlosum("applications/mbi/modules/blosum.txt")

	answer=[]
	for ref in db:
		seqAligned, seqRefAligned, score=calcE(seq, ref, blosum, k)
		answer+=[seqAligned, seqRefAligned, score]

	return answer



"""
seq="ACTTGATAGCCGATTAGGAC"
seqRef="ATTGATTTAGTATATTATTAAATGTATATATTAATTCAATATTATTATTCTATTCATTTTTATTCATTTT"
calcE(seq,seqRef)

seqRef="CAAATTTATAATATATTAATCTATATATTAATTTAGAATTCTATTCTAATTCGAATTCAATTTTTAAATA"
calcE(seq,seqRef)

seqRef="TTCATATTCAATTAAAATTGAAATTTTTTCATTCGCGAGGAGCCGGATGAGAAGAAACTCTCATGTCCGG"
calcE(seq,seqRef)"""




""""




#matrix=needlemanWunsch(seq,seqRef, d, penalty)
#print_matrix(matrix)
#rightBottomCell=(len(seq), len(seqRef))
#seqAligned, seqRefAligned = traceback(matrix, rightBottomCell, d)
#print seqRefAligned
#print seqAligned

# 1.identify common k-words between I (seq) and J(seqRef)
# 2a. Score diagonals with k-word matches
diagonalSums, scoreMatrix=calcDiagonalSums()

#print dot mattrix
printDotMatrix(scoreMatrix)

print
printMatrix(scoreMatrix)

print "diagonal sums:\n"
print diagonalSums


# 2b. identify 10 best diagonals
filteredScoreMatrix,bestDiagonals=matrixWithTopDiagonals(diagonalSums, scoreMatrix)


# 3. Rescore initial regions with a substitution score matrix and get best 10 subregions
scoreMatrix, diagonalSums =rescoreRegions(scoreMatrix,blosum, bestDiagonals)
#Non-matching ends of the diagonal are trimmed?????
#filteredScoreMatrix=matrixWithTopDiagonals(diagonalSums, scoreMatrix)
printMatrix(scoreMatrix)
print


# 4. Join initial regions using gaps, penalise for gaps

# 5. Perform dynamic programming to find final alignments

	"""





"""
	The stages in the FASTA algorithm are as follows:
1. We specify an integer parameter called ktup (short for k respective tuples), and we look
for ktup-length matching substrings of the two strings. The standard recommended
ktup values are six for DNA sequence matching and two for protein sequence matching.
The matching ktup-length substrings are referred to as hot spots. Consecutive hot spots
are located along the dynamic programming matrix diagonals. This stage can be done
effciently by using a lookup table or a hash to store all the ktup-length substrings from
one string, and then search the table with the ktup-length substrings from the other
string.

2. In this stage we wish to find the 10 best diagonal runs of hot spots in the matrix. A
diagonal run is a sequence of nearby hot spots on the same diagonal (not necessarily
adjacent along the diagonal, i.e., spaces between these hot spots are allowed). A run
need not contain all the hot spots on its diagonal, and a diagonal may contain more
than one of the 10 best runs we find.
In order to evaluate the diagonal runs, FASTA gives each hot spot a positive score,
and the space between consecutive hot spots in a run is given a negative score that
decreases with the increasing distance. The score of the diagonal run is the sum of the
hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
runs under this evaluating scheme.

3. A diagonal run species a pair of aligned substrings. The alignment is composed of
matches (the hot spots) and mismatches (from the interspot regions), but it does not
contain any indels because it is derived from a single diagonal. We next evaluate
the runs using an amino acid (or nucleotide) substitution matrix, and pick the best
scoring run. The single best subalignment found in this stage is called init1. Apart
from computing init1, a ltration is performed and we discard of the diagonal runs
achieving relatively low scores.

4. Until now we essentially did not allow any indels in the subalignments. We now try to
combine \good" diagonal runs from close diagonals, thus achieving a subalignment with
indels allowed. We take \good" subalignments from the previous stage (subalignments
whose score is above some specied cuto) and attempt to combine them into a single
larger high-scoring alignment that allows some spaces. This can be done in the following
way:
We construct a directed weighted graph whose vertices are the subalignments found
in the previous stage, and the weight in each vertex is the score found in the previous
stage of the subalignment it represents. Next, we extend an edge from vertex u to
vertex v if the subalignment represented by v starts at a lower row than where the
3.3. BLAST - BASIC LOCAL ALIGNMENT SEARCH TOOL 3
subalignment represented by v ends. We give the edge a negative weight which depends
on the number of gaps that would be created by aligning according to subalignment v
followed by subalignment u. Essentially, FASTA then nds a maximum weight path in
this graph. The selected alignment species a single local alignment between the two
strings. The best alignment found in this stage is marked initn. As in the previous
stage, we discard alignments with relatively low score.

5. In this step FASTA computes an alternative local alignment score, in addition to initn. Recall that init1 denes a diagonal segment in the dynamic programming matrix. We
consider a narrow diagonal band in the matrix, centered along this segment. We observe
that it is highly likely that the best alignment path between the init1 substrings, lies
within the subtable dened by the band. We assume this is the case and compute
the optimal local alignment in this band, using the ordinary dynamic programming
algorithm. Assuming that the best local alignment is indeed within the dened band,
the local alignment algorithm essentially merges diagonal runs found in the previous
stages to achieve a local alignment which may contain indels. The band width is
dependent on the ktup choice. The best local alignment computed in this stage is
called opt.

6. In the last stage, the database sequences are ranked according to initn scores or opt
scores, and the full dynamic programming algorithm is used to align the query sequence
against each of the highest ranking result sequences.

"""
