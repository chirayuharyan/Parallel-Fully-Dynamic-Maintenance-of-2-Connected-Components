CXX=g++
CXXFLAGS = -std=c++11 -fopenmp -O3
LINKFLAGS = -std=c++11 -fopenmp -O3

baseline: baseline.o unweightedGraph.o bfs.o isCutVertex.o
	$(CXX) $(LINKFLAGS) -o baseline baseline.o unweightedGraph.o bfs.o isCutVertex.o

baseline.o: baseline.cpp threads.h unweightedGraph.o bfs.o isCutVertex.o
	$(CXX) $(CXXFLAGS) -c baseline.cpp

unweightedGraph.o: unweightedGraph.cpp unweightedGraph.h 
	$(CXX) $(CXXFLAGS) -c unweightedGraph.cpp

bfs.o: bfs.cpp bfs.h unweightedGraph.o threads.h
	$(CXX) $(CXXFLAGS) -c bfs.cpp

isCutVertex.o: isCutVertex.cpp isCutVertex.h threads.h
	$(CXX) $(CXXFLAGS) -c isCutVertex.cpp


clean: 
	rm -f *.o baseline
