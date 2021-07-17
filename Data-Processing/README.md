# Data-Preprocessing
For our implementations we present the constraints and format for the input graph file:
###### Constraints:
1. The graph is undirected and connected.
2. The vertices are labelled from *0* to *n-1*, where *n* is the total number of vertices in graph.
###### Input Format:
- The first line contains two space-seperated integers *n* and *m* indicating total number of vertices and edges respectively.
- The following *m* lines contains the endpoints of edges space-seperated.
- As the graph is undirected, for an edge between vertex *x* and *y*, we provide two edges in the input namely *(x , y)* and *(y , x)*.

#### Instructions to run the data-preprocessing code:
**Step 1: Extract the raw dataset and copy the graph file in *RawDataset* directory.**

Extract the files using *tar* or *gzip* whichever appropriate. Sometimes compressed files contains a folder, inside which there are graph files and some other files consisting of meta-data like coordinates. Move the graph file to the *RawDataset* directory and delete the rest of the files. 

**Step 2: Run preprocessing python script**

```python dataPreprocessing.py <RawDatasetFileName>```

eg: ```python dataPreprocessing.py webbase-2001.mtx```

After completion of the script a processed dataset file is created in the *Dataset* directory. In the above example, *webbase-2001.txt* file is created in Dataset directory. 

**Note:** These scripts have been tested on the graphs mentioned in the main [readme](../README.md). However, some of the raw dataset files could have a different format hence we recommend to compare new dataset's raw file format with the mentioned dataset files. 
