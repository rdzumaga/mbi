import graph
import copy
import FastaSW

# Read BLOSUM table
def readBlosum(fname):
    """
    Read from file the Blosum table
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

def getTuplesList(str,k):
	"""
	Find k-length tuples of strings in provided DNA sequence
	"""
	tuples=set()
	tuplesDict={}
	
	for i in range(0, len(str)-k+1):
		tuple=str[i:i+k]
		
		if tuple not in tuples:
			tuples.add(tuple)
			tuplesDict[tuple]=[]
			
		tuplesDict[tuple].append(i)
	return tuples,tuplesDict

def createDiagonalDict(seq, seqRef,k):
	"""
	Create a dict with diagonal num as key and an array with one empty DiagonalRun object as value
	"""
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
		newDict[diag]=[]
	return newDict	
		
def calcDiagonalSums(seq, seqRef,tuplesRef, tuplesRefDict,k):
	diagonalDict=createDiagonalDict(seq, seqRef,k)
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

def createMatrixForDots(diagonalDict, seq, seqRef,k):
	"""
	Debug function creating a scoreMatrix for seq and seqRef. Drawing dot matrix is unnecessary since it takes too much memory to remember all values
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
		
def printDotMatrix(seq, seqRef, matrix):
	"""
	Debug function for printing Dot matrix
	"""
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

def listAllRegions(diagonalRegionsDict):
	regions=[]
	for diag in diagonalRegionsDict:
		#diagonalRegionsDict[diag].printIt() 
		for region in diagonalRegionsDict[diag]:
			regions.append(region)
	print
	return regions
	
def getDictWithTopRegions(diagonals):	
	regions=listAllRegions(diagonals)
	if len(regions)<=10:
		return diagonals
	
	#get updated diagonals dictionary
	bestDiagonals=dict(diagonals)
	
	#get keys for the top ten diagonal sums	sort(key=lambda x: x.count, reverse=True)		
	regions.sort(key=lambda x: x.value, reverse=True)	
	
	#if more multiple diagonalRuns have the same value, it is possible, that more than 10 diagonalRuns will be returned
	topRegions=regions[0:11]
	
	for diagNum in diagonals:
		for region in diagonals[diagNum]:
			if region not in topRegions:
				bestDiagonals[diagNum].remove(region)
	
	#remove diagonal entries with 0 diagonalRuns
	for diag in diagonals:
		if len(diagonals[diag])==0:
			del bestDiagonals[diag]

	return bestDiagonals
	
def scoreDiagonalRuns(diagonalDict, seq, seqRef,k, gapPenalty, reward):
	"""
	In order to evaluate each diagonal run, FASTA gives each hot spot a positive score,
	and the space between consecutive hot spots in a run is given a negative score that
	decreases with the increasing distance. The score of the diagonal run is the sum of the
	hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
	runs under this evaluating scheme. Each diagonal may contain more than 1 diagonal run
	"""
	
	rows=len(seq)
	cols=len(seqRef)
	newDiagonalDict=dict(diagonalDict)
	
	#group all consecutive hotspots into regions for each diagonal
	for diag in diagonalDict:
		regionNr=1
		regions=[]
		hspots=sorted(diagonalDict[diag][0].hotspots)
		region=graph.DiagonalRun(diag,str(diag)+"."+str(regionNr))
		newDiagonalDict[diag]=[]
		for i in range(len(hspots)):
			if(i>0):
				if hspots[i][0]!=hspots[i-1][0]+1:
					newDiagonalDict[diag].append(copy.copy(region))
					region=graph.DiagonalRun(diag,str(diag)+"."+str(regionNr))
					regionNr+=1
						
			region.add(hspots[i][0], hspots[i][1], reward)	
		newDiagonalDict[diag].append(copy.copy(region))
		
	#check if regions would have a better score when combined with others, inspite of gaps
	for diag in newDiagonalDict:
		for i in range( len(newDiagonalDict[diag]) -1 ):
			regA=newDiagonalDict[diag][i]
			regB=newDiagonalDict[diag][i+1]
			regAB=graph.DiagonalRun(diag, regA.name)
			
			for hotspot in regA.hotspots+regB.hotspots:
				regAB.add(hotspot[0], hotspot[1])
				
			pos=(regA.hotspots[-1][0]+1, regA.hotspots[-1][1]+1)
			#gapPenalty=gapOpenPenalty
			
			while (pos!=regB.hotspots[0]):
				regAB.add(pos[0],pos[1], gapPenalty)
				
				#make sure penalty stays a negativ value (in case there were too many consecutive gaps): for afinite penalty function
				#if gapPenalty<-1:
					#gapPenalty+=1
					
				pos=(pos[0]+1, pos[1]+1)
			
			if regAB.value>regA.value:
				#replace A and B with AB
				newDiagonalDict[diag][i]=regAB
				newDiagonalDict[diag].remove(regB)
					
	
			
	rescoredDiagonalsDict=getDictWithTopRegions(newDiagonalDict)
	

			
	return rescoredDiagonalsDict

	
def rescoreDiagonals(seq, seqRef, blosum, bestDiagonalsDict,k, cutoff):
	"""
	@retval (Dictionary with diagonal numbers as keys and DiagonalRun objects as items, init1 (best DiagonaRun))
	"""
	rescoredDiagonals=createDiagonalDictFrom(bestDiagonalsDict)	
	
	#iterate over diagonals and score the with blosum matrix
	for diag in bestDiagonalsDict:
		rescoredRegions=[]
		for region in bestDiagonalsDict[diag]:
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
	
	
	for diag in bestDiagonalsDict:
		for region in bestDiagonalsDict[diag]:
			if region.value<cutoff:
				if len(rescoredDiagonals[diag])>1:
					bestRescoredDiagonals[diag].remove(region)
				else:
					del bestRescoredDiagonals[diag]
			

	
	rescoredRegions=listAllRegions(bestRescoredDiagonals)
	print "len=", len(bestRescoredDiagonals)
	init1=max(rescoredRegions)

	return bestRescoredDiagonals, init1

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


def fastaScoreAlignment(seq, seqRef, k, gapPenalty=-10, rescoreCutoff=10, matchReward=20, blosum=""):
	"""
	
	@param	matchReward reward used when scoring initial regions. Must be positive
	"""
	rows=len(seq)
	cols=len(seqRef)
	blosum=readBlosum("blosum.txt")
	
	#1.identify common k-words between I (seq) and J(seqRef)
	tuplesRef,tuplesRefDict=getTuplesList(seqRef,k)
	
	#create a dict with diagonalNums being the keys and hotspots as values. To access one of the runs use diagonalDict[key][i], where i is the number of one of diagonal runs on this diagonal
	diagonalDict=calcDiagonalSums(seq, seqRef, tuplesRef,tuplesRefDict, k)
	if len(diagonalDict)==0:
		print "step1 error"
		return 0
		
	# 2 Score diagonals with k-word matches and identify 10 best diagonals
	#save the 10 best local regions, regardless of whether they are on the same of different diagonals.
	#print "\n=======================================\n                STEP2\n=======================================\n"
	bestTenDiagonalRunsDict=scoreDiagonalRuns(diagonalDict,seq, seqRef,k, gapPenalty, matchReward)

	# 3. Rescore best regions with a substitution score matrix,find init1 and dispose of diagonal runs with value below cutoff
	#print "\n=======================================\n                STEP3\n=======================================\n"
	diagonalRegionsDict, init1=rescoreDiagonals(seq, seqRef, blosum, bestTenDiagonalRunsDict,k, rescoreCutoff)

	if len(diagonalRegionsDict)==0:
		print "Step 3 error"
		return 0
	
	#4. Join initial regions using gaps (create a graph), penalise for gaps
	#print "\n======================================="                STEP4\n=======================================\n"
				
	mygraph=graph.createGraph(listAllRegions(diagonalRegionsDict), -1)
	
	if mygraph==None:
		regions=listAllRegions(diagonalRegionsDict)
		init_n=max(regions)
		
	path=[]
	value=0
	allPaths=[]
	for startNode in mygraph:
		for node in mygraph[startNode]:
			path, value=graph.findBestPath(mygraph, startNode, node[0])
			if path:
				allPaths.append([value,path])
				
	allPaths.sort(reverse=True)
	init_n=allPaths[0]
	
	
		
	# 5. Perform dynamic programming (Smith-Waterman)to find final alignments

	matrix, seqAligned, seqRefAligned,opt_score, alignedSeqRefStartIndex=FastaSW.SmithWaterman(seq, seqRef, init1.hotspots, k, gapPenalty)

	
	return init1.value, init_n[0], opt_score, seqAligned, seqRefAligned, alignedSeqRefStartIndex, len(seqRefAligned)
	
def readDb(fname):
	lines = open(fname, "rt").readlines()
	db=[]
	for line in lines:
		#get rid of CRLF
		db.append(line[:-1])
	return db
	
def fasta(seq, k=2, gapPenalty=-10, rescoreCutoff=10, matchReward=20, db=[], blosum="" ):
	if len(db)==0:
		db=readDb("db.txt")
	
	if len(blosum)==0:
		blosum=readBlosum("blosum.txt")
	
	answer=[]
	for i in range(len(db)):
		init1, init_n, opt, seqAligned, seqRefAligned, refStartIndex, length=fastaScoreAlignment(seq, db[i], k, gapPenalty, rescoreCutoff, matchReward, blosum)
		
		astring, idents, gaps, mismatched=createAlignmentString(seqAligned, seqRefAligned)
		answer.append([init1, init_n, opt, i, refStartIndex, length, mismatched])
		
		
	return answer
	