#ifndef ISCUTVERTEX_H
#define ISCUTVERTEX_H
#include<vector>
#include<utility>
#include"threads.h"


double markCutVertices(int** baseVertex1, int** baseVerter2, int* head, bool* isSafe, bool *isLCA, bool* isBase, int* parent, bool *isCutVertex, bool* isSafeOld, int* componentNumber, int n, int root);

#endif
