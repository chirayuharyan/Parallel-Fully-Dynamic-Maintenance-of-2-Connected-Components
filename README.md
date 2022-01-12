# Parallel Fully Dynamic Maintenance of 2-Connected Components <mark>(Please download the report and view as sometimes the document doesn't render well on this platform.)</mark>

## Codes
This directory contains the implementation of our algorithms. 
The instructions to run these implementations is available in [README](Code/README.md) file inside the directory.

## Reports
The directory *Reports* includes an extended version of our paper. 

## Datasets

|    Dataset Name   | Link                                                                               | Vertices | Edges     |
|:-----------------:|------------------------------------------------------------------------------------|----------|-----------|
| amazon0302        | https://snap.stanford.edu/data/amazon0302.html | 262111   | 1234877   |
| web-Stanford      | https://snap.stanford.edu/data/web-Stanford.html | 281903   | 2312497   |
| web-Google        | https://snap.stanford.edu/data/web-Google.html | 875713   | 5105039   |
| roadNet-PA        | https://snap.stanford.edu/data/roadNet-PA.html | 1088092  | 1541898   |
| roadNet-CA        | https://snap.stanford.edu/data/roadNet-CA.html | 1965206  | 2766607   |
| cit-Patents       | https://snap.stanford.edu/data/cit-Patents.html | 3774768  | 16518948  |
| soc-LiveJournal1  | https://snap.stanford.edu/data/soc-LiveJournal1.html | 4847571  | 68993773  |
| europe_osm        | https://www.cise.ufl.edu/research/sparse/matrices/DIMACS10/europe_osm.html | 50912018 | 108109320 |
| com-orkut.ungraph | https://snap.stanford.edu/data/com-Orkut.html | 3072441  | 117185083 |
| webbase-2001        | https://www.cise.ufl.edu/research/sparse/MM/LAW/webbase-2001.tar.gz | 118142155   | 1019903190   |

##### Datasets Preprocessing
For our implementations we present the constraints and input format for the input graph:
###### Constraints:
1. The graph is undirected and connected.
2. The vertices are labelled from *0* to *n-1*, where *n* is the total number of vertices in graph.
###### Input Format:
- The first line contains two space-seperated integers *n* and *m* indicating total number of vertices and edges respectively.
- The following *m* lines contains the endpoints of edges space-seperated.
- As the graph is undirected, for an edge between vertex *x* and *y*, we provide two edges in the input namely *(x , y)* and *(y , x)*.

In the *Data Preprocessing* directory, we provide the scripts to convert the datasets and also the [instructions](Data-Processing/README.md) to run the scripts.

