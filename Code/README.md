# Query Format
- The first line in the query file consists of total number of batches.
- Each batch is represented as follows:
    - The first line consists of two integers <type> <totalEdges> which indicates the batch type (1- Insert, 2- Delete) and total number of edges in the batch.
    - The following <totalEdges> lines contain edges that are to be inserted/deleted from the graph. 
- **Note:** While processing the batch, in case of insert batch, these edges should strictly not be part of the graph. Similarly, for the delete batch, the edge must be part of the graph. 

# Instructions to run the codes
**Input graph:** The graph should be connected and the format should be as mentioned in the main [README](../README.md) file.

## 1. Static

**Step 0: Compile the program**

In the directory *Static*, open a terminal and write the following command:

```make```

**Step 1: Run the program**

```./baseline <inputGraph> <totalNumberOfQueryFiles> <QueryFile1> <QueryFile2> ... output.txt```


## 2. Incremental Dynamic

**Step 0: Compile the program**

In the directory *IncrementalDynamic*, open a terminal and write the following command:

```make```

**Step 1: Run the program**

```./add <inputGraph> <totalNumberOfQueryFiles> <QueryFile1> <QueryFile2> ... output.txt```

**Note:** The query files should strictly not have any delete batches as this is a incremental dynamic program. 

## 3. Fully Dynamic

In this version, we provide a spanning tree of inputGraph as an additional input.   

**Step 0: Compile the program**

In the directory *FullyDynamic*, open a terminal and write the following command:

```make```

**Step 1: Run the program**

```./biconnected <inputGraph> <SpanningTree> <totalNumberOfQueryFiles> <QueryFile1> <QueryFile2> ... output.txt```
