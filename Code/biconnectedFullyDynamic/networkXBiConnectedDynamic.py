import sys
import networkx as nx
from collections import OrderedDict
import time

class myGraph:
	totalVertices = -1
	totalEdges = -1
	edgeList = [] #consists of two values , u and v

def readInputs(G,inputFile):
	graphData = inputFile.readline().split(" ")
	G.totalVertices = int(graphData[0])
	G.totalEdges = int(graphData[1])
	graphData = inputFile.readlines()

	for currEdge in graphData:
		currEdge = currEdge.split(" ")
		l = (int(currEdge[0]),int(currEdge[1]))
		G.edgeList.append(l)

def readQuery(queries):
	queryFile = open(sys.argv[2],'r')
	totalQuery = queryFile.readline()
	queries = queryFile.readlines()	

def getCutVertices(G):
	return list(nx.articulation_points(G))

def getBridges(G):
	return list(nx.bridges(G))


def getBiconnectedComponents(currG):
	G = nx.Graph()
	G.add_edges_from(currG.edgeList)
	# print(G.number_of_nodes())
	# print(G.number_of_edges())
	return nx.biconnected_components(G)

def printComponents(Components):
	vertexList = []
	currComponentNumber = 1
	totalCount = 0
	for c in Components:
		totalCount += 1
		s = ""
		for v in c:
			s = s+str(v)+" "
		vertexList.append(s)

	vertexList.sort()
	print(len(vertexList))
	for item in vertexList:
		print(item)
	# print(totalCount)
	# print(len(vertexList))
	



def main():
	inputFile = open(sys.argv[1],'r')
	currG = myGraph() 
	readInputs(currG,inputFile) #reads data correctly
	G = nx.Graph()
	G.add_edges_from(currG.edgeList)

	

	start_time = time.time()
	cutVertices = getCutVertices(G)
	bridges = getBridges(G)
	end_time = time.time() - start_time;



	cutVertices.sort()

	# print(len(bridges))
	# print(len(cutVertices))
	# print(bridges)
	# print(cutVertices)
	print(nx.is_connected(G))
	print("Bridges")
	bridgesList = []
	for i in bridges:
		bridgesList.append(str(i[0])+" "+str(i[1]))
	bridgesList.sort()
	for i in bridgesList:
		print(i)
	print("Cut Vertices")
	# print(type(cutVertices))
	for i in cutVertices:
		print(i)
	print("")
	print("Time taken: "+str(end_time/1000)+" milliseconds")

	queries = open(sys.argv[2],'r').readlines()
	totalBatch = int(queries[0])
	currLine = 0
	# print(len(queries))
	while(currLine < len(queries)):
		currLine += 1
		if(currLine >= len(queries)):
			break
		# print(currLine)
		type = int(queries[currLine].split()[0])
		totalQuery = int(queries[currLine].split()[1])
		# print(str(type)+" "+str(totalQuery))
		for q in range(totalQuery):
			currLine += 1
			curr = queries[currLine].split()
			u = int(curr[0])
			v = int(curr[1])
			# print("Curr Query "+str(type)+" "+str(u)+" "+str(v))
			if type == 1:
				G.add_edge(u,v)
			else:
				G.remove_edge(u,v)
		
		start_time = time.time()
		cutVertices = getCutVertices(G)
		bridges = getBridges(G)
		cutVertices.sort()
		end_time = time.time() - start_time;
		print(nx.is_connected(G))
		print("Bridges")
		bridgesList = []
		for i in bridges:
			bridgesList.append(str(i[0])+" "+str(i[1]))
		bridgesList.sort()
		for i in bridgesList:
			print(i)
		print("Cut Vertices")
	
		for i in cutVertices:
			print(i)
		print("")
		print("Time taken: "+str(end_time/1000)+" milliseconds")
	
main()