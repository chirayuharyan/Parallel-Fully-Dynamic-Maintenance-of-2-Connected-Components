#ifndef ISCUTVERTEX_H
#define ISCUTVERTEX_H
#include<vector>
#include<utility>
#include"threads.h"


double markCutVertices(int** baseVertex1, int** baseVertex2, int* head, bool* unfinishedTraversal, bool *isLCA, bool* isBase, int* parent, bool *isCutVertex, int* componentNumber, int n, int root, bool* isSafe);

double markCutVerticesAffected(int** baseVertex1, int** baseVertex2, int* head, int** affectedLCA, int* headAffectedLCA, bool* unfinishedTraversal, bool *isLCA, bool* isBase, int* parent,int* offset, int* children, bool *isCutVertex, int* componentNumber, int n, int root, bool* isSafe);

// double markCutVerticesAffected(int**& baseVertex1, int**& baseVertex2, int*& head, int**& affectedLCA, int*& headAffectedLCA, bool*& unfinishedTraversal, bool *&isLCA, bool*& isBase, int*& parent,int*& offset, int*& children, bool*& isCutVertex, int*& componentNumber, int n, int root, bool*& isSafe);

#endif
