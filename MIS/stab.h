#ifndef _STAB_H_
#define _STAB_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <graph.h>

// Graph

typedef enum GT { graph, complement } graphType;

// Functions

void getStabGraph(adjMatrix * g);
void setStabGraph(adjMatrix * g);
int readGraphFile(FILE * graphFile, const graphType t, int * n, long long * m);
int mcrSort(const int n, int * R, int * antideg, int * maxantideg);
int coverSort(const int r, int * R, int * cover);
int readWeightFile(FILE * graphFile, int n, double * weight);
int mcrWeightedSort(const int n, int * R, double * weight, double * maxantiwdeg);
double coverWeightedSort(const int r, int * R, double * weight, double * cover);
int edgeProjSeparation(int nzcnt, int * vind, double * xval, int maxncut, int * cutsize, int ** cutvind, double * cutrhs);

#endif /*_STAB_H_*/
