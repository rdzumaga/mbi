import copy
class ShortestPathResult(object):
	def __init__(self):
                self.d = {}
                self.parent = {}

				
def shortest_path(graph, s):
        """Single source shortest paths using DP on a DAG.
        Args:
        graph: weighted DAG.
        s: source
        """
        result = ShortestPathResult()
        result.d[s] = 0

        result.parent[s] = None
        for v in graph.itervertices():
                result=sp_dp(graph, v, result)
        return result

def sp_dp(graph, v, result):
        """Recursion on finding the shortest path to v.
         Args:
                 graph: weighted DAG.
                 v: a vertex in graph.
                 result: for memoization and keeping track of the result.
        """
        if v in result.d:
                return result.d[v]
        result.d[v] = float('inf')
        result.parent[v] = None
        for u in graph.inverse_neighbors(v): # Theta(indegree(v))
                new_distance = sp_dp(graph, u, result) + graph.weight(u, v)
                if new_distance < result.d[v]:
                        result.d[v] = new_distance
                        result.parent[v] = unic
        return result.d[v]


def shortest_path_bottomup(graph, s):
        """Bottom-up DP for finding single source shortest paths on a DAG.
        Args:
                graph: weighted DAG.
                s: source
        """
        order = topological_sort(graph)
        result = ShortestPathResult()
        for v in graph.itervertices():
                result.d[v] = float('inf')
                result.parent[v] = None
        result.d[s] = 0
        for v in order:
                for w in graph.neighbors(v):
                        new_distance = result.d[v] + graph.weight(v, w)
                        if result.d[w] > new_distance:
                                result.d[w] = new_distance
                                result.parent[w] = vars
        return result
		
def findShortestPath(graph, start, end, path=[], value=0):
	path=path+ [start]
	value+=start.value
	if(start==end):
		#print "START=END"
		return path, value
		
	if(start in graph):
		bestPath=None
		bestValue=-float('inf')
		best=None
		for node in graph[start]:
			if node not in path:
				#print "Looking for", start ,"->", end,
				#for p in path:
					#print p.name,
				#print value
				#value+=node[1]
				newPath, newValue=findShortestPath(graph, node[0], end, path, value)
				newValue+=node[1]
				if newPath:
					if not bestPath or bestValue<newValue:
						best=newPath
						bestValue=newValue
				
		#print "DOWN_RETURN"
		return best, bestValue
	return None, 0

def calcPathValue(path):
	val=0
	for node in path:
		val+=node.value
def findAllPaths(graph, start, end, path=[]):
	#print "Finding path for start=", start, " end=", end
	path=path+ [start]
	
	if(start==end):
		return [path]
		
	if(start in graph):
		paths=[]
		for node in graph[start]:
			if node not in path:
				newPaths=findAllPaths(graph, node[0], end, path)
				for newPath in newPaths:
					paths.append(newPath)
		return paths
	return []
	
class DiagonalRun:
	def __init__(self, diagonalNum, name=""):
		self.hotspots=[]
		self.value=0
		self.diag=diagonalNum
		self.name=name
	
	def add(self, i,j, val=2):
		if (i,j) not in self.hotspots:
			self.hotspots.append((i,j))
			self.value+=val
			
	def __str__(self):
		return self.name
		
	def printIt(self):
		print "diag=", self.name, self.diag, ":", self.value, self.hotspots

def createGraph(connections, subregions):
	graph={}
	
	print"!!!!!!!!!!!!!CREATING GRAPH!!!!!!!!!!!!!!!!"
	for con in connections:
		start=con[0]
		graph[start]=[]
	
	for con in connections:
		start=con[0]
		end=con[1]
		cost=con[2]
		#print start, end, cost
		#print "Start=", start.name
		graph[start].append( (end, cost) )
	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"	
	return graph
		

def distanceBetween(v, u):
	
	v_lastRow=v.hotspots[len(v.hotspots)-1][0]
	v_lastCol=v.hotspots[len(v.hotspots)-1][1]
	u_firstRow=u.hotspots[0][0]
	u_firstCol=u.hotspots[0][1]
	
	row_dist=u_firstRow-v_lastRow-1
	col_dist=u_firstCol-v_lastCol-1
	
	
	if row_dist<0 or col_dist<0:
		return -1
	
	#print"--------Comparing distance between-------------"
	#v.printIt()
	#u.printIt()
	#print 
	#print "row, col V:", v_lastRow, v_lastCol, "row, col U:", u_firstRow, u_firstCol
	#print "dist:", row_dist, col_dist
	#print "------------------------------------------------!"
	return max(row_dist, col_dist)
	
def printGraph(graph):
	print "graph={"
	for key in graph:
		print "\t", key, ": [",
		for dest in graph[key]:
			cost=dest[1]
			destNode=dest[0].name
			print (destNode, cost),
		print "]"
		
print "----------GRAPHS_-------------"

#create sample verices for graph
diagDict={}
diags=[]

diagNum=9
diag=DiagonalRun(diagNum,"b")
diag.add(10,1)
diag.add(11,2)
diag.add(12,3)
diag.add(13,4)
diag.add(14,5)
diag.add(15,6)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag


diagNum=4
diag=DiagonalRun(diagNum,"c")
diag.add(6,2)
diag.add(7,3)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diag=DiagonalRun(diagNum,"k")
diag.add(9,5)
diag.add(10,6)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=1
diag=DiagonalRun(diagNum, "f")
diag.add(11,10)
diag.add(12,11)
diag.add(13,12)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=0
diag=DiagonalRun(diagNum,"i")
diag.add(1,1)
diag.add(2,2)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=-1
diag=DiagonalRun(diagNum,"d")
diag.add(5,6)
diag.add(6,7)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=-2
diag=DiagonalRun(diagNum, "e")
diag.add(13,15)
diag.add(14,16)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=-3
diag=DiagonalRun(diagNum, "h")
diag.add(8,11)
diag.add(9,12)
diag.add(10,13)
diag.add(11,14)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=-5
diag=DiagonalRun(diagNum, "g")
diag.add(2,7)
diag.add(3,8)
diag.add(4,9)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diag=DiagonalRun(diagNum, "j")
diag.add(12,17)
diag.add(13,18)
diag.add(14,19)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag

diagNum=-7
diag=DiagonalRun(diagNum, "a")
diag.add(5,12)
diag.add(6,13)
diag.add(7,14)
diag.add(8,15)
diags.append(copy.copy(diag))
diagDict[diagNum]=diag


for d in diags:
	d.printIt()



#create graph
subregions=diags
connections=[]

for v in subregions:
	for u in subregions:
		if u!=v:
			dist=distanceBetween(v,u)
			if(dist>0):
				connections.append((v, u, -dist))
	
print "CONNECTIONS:"
for node in connections:
	print "(", node[0].name,",", node[1].name, ",", node[2], ")"
	
	
graph=createGraph(connections, subregions)

printGraph(graph)
print

paths=[]
#for startNode in graph:
	#print "\t", startNode
	#startNode.printIt()
	#for node in graph:
		#paths.append(findPath(graph, startNode, node))
		
paths2=findAllPaths(graph,connections[0][0], connections[4][1])
for path in paths2:
	for node in path:	
		print node,
	print

path=[]
value=0	
when=1
for startNode in graph:
	print "ooooooooooooooooooooooooooooooooooooooooo"
	for node in graph[startNode]:
		print"\t____________________"
		path, value=findShortestPath(graph, startNode, node[0])
		print value,
		if path:
			for n in path:	
				print n,
		print
		
	print
	if (when>100):
		break
	when+=1
		
