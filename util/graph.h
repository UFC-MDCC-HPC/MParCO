/*
 * graph.h
 *
 *  
 *      Author: Corrêa
 */

#ifndef GRAPH_H_
#define GRAPH_H_

typedef struct AM {
	char * 	g;
	long	matrixsize;
	int		logrowsize;
} adjMatrix;

inline static int inline_ceillog2(int n) {
	int c = 0;
	int l = 0;
	while (n > 1) {
		c |= (n & 1);
		n >>= 1;
		l++;
	}
	l += c;
	return l;
}

static char * 	_GRAPH_H_g;
static long 	_GRAPH_H_matrixsize;
static int 		_GRAPH_H_logrowsize;

#define newAdjMatrixSize(n) (_GRAPH_H_matrixsize = ((long) 1 << (_GRAPH_H_logrowsize = inline_ceillog2(n)-3))*(n))
#define newAdjMatrix(n) _GRAPH_H_g = (char *) malloc(newAdjMatrixSize(n))
#define getAdjMatrix(a) a->g = _GRAPH_H_g; \
						a->matrixsize =	_GRAPH_H_matrixsize; \
						a->logrowsize = _GRAPH_H_logrowsize;
#define setAdjMatrix(a) _GRAPH_H_g = a->g; \
						_GRAPH_H_matrixsize = a->matrixsize; \
						_GRAPH_H_logrowsize = a->logrowsize;
#define cloneAdjMatrix(dest, orig) dest->g = (char *) malloc(orig->matrixsize); \
								   memcpy(dest->g, orig->g, orig->matrixsize); \
								   dest->matrixsize = orig->matrixsize; \
								   dest->logrowsize = orig->logrowsize;
#define freeAdjMatrix() free(_GRAPH_H_g)
#define hasGraph() (_GRAPH_H_g != NULL)
#define edgeIndex(i, j) (((long) i << _GRAPH_H_logrowsize)+(j >> 3))
#define edgeMask(j) (1 << (char) (j & 0x07))
#define addEdge(i, j) (_GRAPH_H_g[edgeIndex(i, j)] |= edgeMask(j))
#define delEdge(i, j) (_GRAPH_H_g[edgeIndex(i, j)] &= ~edgeMask(j))
#define addAllEdges() (memset(_GRAPH_H_g, 0xFF, _GRAPH_H_matrixsize))
#define delAllEdges() (memset(_GRAPH_H_g, 0, _GRAPH_H_matrixsize))
#define hasEdge(i, j) (_GRAPH_H_g[edgeIndex(i, j)] & edgeMask(j))


#endif /* GRAPH_H_ */
