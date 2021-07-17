CXX=g++
CXXFLAGS = -std=c++11 -fopenmp -O3 
LINKFLAGS = -std=c++11 -fopenmp -O3 

add: add.o unweightedGraph.o bfs.o isCutVertex.o 
	$(CXX) $(LINKFLAGS) -o add add.o unweightedGraph.o bfs.o isCutVertex.o
	
add.o: add.cpp threads.h unweightedGraph.o bfs.o isCutVertex.o
	$(CXX) $(CXXFLAGS) -c add.cpp

unweightedGraph.o: unweightedGraph.cpp unweightedGraph.h 
	$(CXX) $(CXXFLAGS) -c  unweightedGraph.cpp

bfs.o: bfs.cpp bfs.h unweightedGraph.o threads.h
	$(CXX) $(CXXFLAGS) -c bfs.cpp

isCutVertex.o: isCutVertex.cpp isCutVertex.h threads.h
	$(CXX) $(CXXFLAGS) -c isCutVertex.cpp

clean: 
	rm -f *.o add
