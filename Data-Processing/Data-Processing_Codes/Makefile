CXX=g++
CXXFLAGS = -std=c++11 -fopenmp -O3 -Wall
LINKFLAGS = -std=c++11 -fopenmp -O3 -Wall

process: dataProcessSNAP.o unweightedGraph.o bfs.o
	$(CXX) $(LINKFLAGS) -o process dataProcessSNAP.o unweightedGraph.o bfs.o

processYifanhu: dataProcessYifanhu.o unweightedGraph.o bfs.o
	$(CXX) $(LINKFLAGS) -o processYifanhu dataProcessYifanhu.o unweightedGraph.o bfs.o

dataProcessYifanhu.o: dataProcessYifanhu.cpp unweightedGraph.o bfs.o
	$(CXX) $(CXXFLAGS) -c dataProcessYifanhu.cpp
	
dataProcessSNAP.o: dataProcessSNAP.cpp unweightedGraph.o bfs.o
	$(CXX) $(CXXFLAGS) -c dataProcessSNAP.cpp

bfs.o: bfs.h bfs.cpp unweightedGraph.o
	$(CXX) $(CXXFLAGS) -c bfs.cpp

unweightedGraph.o: unweightedGraph.cpp unweightedGraph.h
	$(CXX) $(CXXFLAGS) -c unweightedGraph.cpp

clean: 
	rm -f *.o process processYifanhu
