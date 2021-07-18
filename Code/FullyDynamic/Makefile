CXX=g++
CXXFLAGS = -fopenmp -O3 -Wall
LINKFLAGS = -fopenmp -O3 -Wall


biconnected: main.o unweightedGraph.o spanningTree.o lca.o dynamicGraph.o bfs.o isCutVertex.o addEdge.o deleteEdge.o threads.h
	g++ -std=c++11 -fopenmp -o biconnected main.o unweightedGraph.o spanningTree.o lca.o dynamicGraph.o bfs.o isCutVertex.o addEdge.o deleteEdge.o

main.o: main.cpp threads.h
	g++ -std=c++11 -c -fopenmp -o main.o main.cpp

unweightedGraph.o: unweightedGraph.cpp threads.h
	g++ -std=c++11 -c -fopenmp unweightedGraph.cpp

dynamicGraph.o: dynamicGraph.cpp threads.h
	g++ -std=c++11 -c -fopenmp dynamicGraph.cpp

bfs.o: bfs.cpp bfs.h unweightedGraph.o threads.h
	g++ -std=c++11 -c -fopenmp bfs.cpp

spanningTree.o: spanningTree.cpp unweightedGraph.o threads.h
	g++ -std=c++11 -c -fopenmp spanningTree.cpp

lca.o: lca.cpp dynamicGraph.o unweightedGraph.o threads.h
	g++ -std=c++11 -c -fopenmp lca.cpp

isCutVertex.o: isCutVertex.cpp dynamicGraph.o threads.h
	g++ -std=c++11 -c -fopenmp isCutVertex.cpp

addEdge.o: addEdge.cpp lca.o unweightedGraph.o dynamicGraph.o threads.h
	g++ -std=c++11 -c -fopenmp addEdge.cpp

deleteEdge.o: deleteEdge.cpp lca.o unweightedGraph.o dynamicGraph.o addEdge.o threads.h
	g++ -std=c++11 -c -fopenmp deleteEdge.cpp

clean: 
	rm -f *.o biconnected
