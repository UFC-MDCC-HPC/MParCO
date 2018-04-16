/*
 * stab.c
 *
 *  Created on: Nov 10, 2010
 *      Author: Corrêa
 */

#include <mcr.h>
#include <stab.h>

#define MINVIOLATION 0.3
#define MAXITER 5
#define K 5
#define MAXGREEDYITER 10

void setStabGraph(adjMatrix * g) {
	setAdjMatrix(g);
}

void getStabGraph(adjMatrix * g) {
	getAdjMatrix(g);
}

int readGraphFile(FILE * graphFile, const graphType t, int * n, long long * m) {
	char               type  = ' ';
	char               linestr[100];
	char *             datastr;
	int                i, j;

	*n = 0;
	*m = 0;

	while (type != 'p') {
		type = fgetc(graphFile);
		if (type != EOF) {

			/* header */
			if (type == 'c') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr != NULL)
					printf("%s", linestr);
				else
					return -1;
			}

			/* Vertices */
			if (type == 'p') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr == NULL)
					return -1;

				datastr = strtok(linestr," ");

				datastr = strtok(NULL," ");
				*n = atoi(datastr);

				datastr = strtok(NULL," ");
				*m = atoll(datastr);
				if (t == complement)
					*m = ((((long long) (*n))*((long long) (*n)) - ((long long) *n)) >> ((long long) 1)) - *m;
			}
		}
	}

	////
	// Graph variables
	////

	printf("Graph with %d vertices and %lld edges.\n", *n, *m);
	newAdjMatrix(*n);

	long long nedges;
	if (t == graph) {
		nedges = 0;
	}
	else {
		addAllEdges();
		nedges = (*n)*(*n);
		for (i = 0; i < *n; i++) {
			delEdge(i, i);
			nedges--;
		}
		nedges = nedges >> 1;
	}

	type = fgetc(graphFile);
	while (type != EOF) {
		/* Edges */
		if (type == 'e') {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr == NULL)
				return -1;

			datastr = strtok(linestr," ");
			i = atoi(datastr) - 1;

			datastr = strtok(NULL," ");
			j = atoi(datastr) - 1;

			if (t == graph) {
				addEdge(i, j);
				addEdge(j, i);
				nedges++;
			}
			else {
				delEdge(i, j);
				delEdge(j, i);
				nedges--;
			}
		}
		else {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr != NULL)
				printf(" %s\n", linestr);
			else
				return -1;
		}
		type = fgetc(graphFile);
	}

	if (nedges != *m)
		printf("Ops! Reading %lld edges instead of %lld announced.\n", nedges, *m);

	return 0;
}

static int computeExantideg(const int * R, const int i, const int * antideg, const int k) {
	int d = 0;
	int j;
	for (j = k; j > 0; j--)
		if (!hasEdge(R[i], R[j]))
			d += antideg[R[j]];
	return d;
}

// if i is returned, then vertices R[0...i] induce a regular subgraph
int mcrSort(const int n, int * R, int * antideg, int * maxantideg) {
	// number of antineighbors
	memset(antideg, 0, n*sizeof(int));

	// maximum number of antineighbors (including itself)
	// used to complete vertex cover
	*maxantideg = 0;
	int iu, iv;
	for (iu = 0; iu < n; iu++) {
		antideg[iu]++;
		for (iv = iu+1; iv < n; iv++)
			if (iu != iv && !hasEdge(R[iu], R[iv])) {
				antideg[iu]++;
				antideg[iv]++;
			}
		if (antideg[iu] > *maxantideg)
			*maxantideg = antideg[iu];
	}

	// used to break ties
	int exantideg[n];
	// a negative number indicates that this parameter is not known
	memset(exantideg, -1, n*sizeof(int));

	// determine min antideg
	int Rmin[n];
	int nrmin = 1;
	Rmin[0] = 0;
	int antidegrmin = antideg[0];
	int i = n - 1;
	int j;
	for (j = 1; j <= i; j++)
		if (antideg[j] < antidegrmin) {
			nrmin = 1;
			Rmin[0] = j;
			antidegrmin = antideg[j];
		}
		else if (antideg[j] == antidegrmin)
			Rmin[nrmin++] = j;

	int rmin;
	while (i >= nrmin) {
		rmin = 0;
		if (nrmin > 1) {
			for (j = 1; j < nrmin; j++) {
				if (exantideg[Rmin[rmin]] < 0)
					exantideg[Rmin[rmin]]=computeExantideg(R, Rmin[rmin], antideg, i);
				if (exantideg[Rmin[j]] < 0)
					exantideg[Rmin[j]]=computeExantideg(R, Rmin[j], antideg, i);
				if (exantideg[Rmin[j]] < exantideg[Rmin[rmin]])
					rmin = j;
			}
		}

		rmin = Rmin[rmin];
		if (i != rmin) {
			int aux = R[i];
			R[i] = R[rmin];
			R[rmin] = aux;
			aux = antideg[i];
			antideg[i] = antideg[rmin];
			antideg[rmin] = aux;
			aux = exantideg[i];
			exantideg[i] = exantideg[rmin];
			exantideg[rmin] = aux;
		}

		for (j = 0; j < i; j++)
			if (!hasEdge(R[j], R[i])) {
				exantideg[j] -= antideg[i];
				antideg[j]--;
			}

		i--;
		Rmin[0] = 0;
		nrmin = 1;
		antidegrmin = antideg[0];
		for (j = 1; j <= i; j++)
			if (antideg[j] < antidegrmin) {
				nrmin = 1;
				Rmin[0] = j;
				antidegrmin = antideg[j];
			}
			else if (antideg[j] == antidegrmin)
				Rmin[nrmin++] = j;
	}

	return i;
}

// sort R according to a clique cover of the subgraph induced by R
// it returns the size of the clique cover and sets the cover array such that:
// for i in 0, ..., r-1, the entry cover[i] contains the size of a clique
// cover of the subgraph induced by R[0], ..., R[i]
int coverSort(const int r, int * R, int * cover) {
	int coversize = 0;
	int i, j;

	// linking the elements of R
	int next[r+1];
	// head of the candidates list
	int head = 0;
	for (i = head; i < r - 1; i++)
		next[i] = i + 1;
	next[i] = -1;
	// tail of the cliques list
	int tail = r;
	// head of the cliques list is r
	next[tail] = -1;

	for (i = head; i >= 0; i = head) {
		// current clique starts with size 1
		cover[i] = coversize + 1;
		// remove i from the candidates list
		head = next[i];
		// include i in the tail of cliques list
		next[i] = next[tail];
		next[tail] = i;
		tail = i;

		// search for new vertices to add to the current clique
		int antj = head, nj;
		for (j = head, nj = next[j]; j >= 0; j = nj, nj = next[j]) {
			int cc;
			int vaj = R[j];
			for (cc = i; cc >= 0 && hasEdge(R[cc], vaj); cc = next[cc]);
			if (cc == -1) {
				cover[j] = coversize + 1;
				// remove j from the candidates list
				if (j == head)
					head = next[j];
				else
					next[antj] = next[j];
				// include j in the tail of cliques list
				next[j] = next[tail];
				next[tail] = j;
				tail = j;
			}
			else
				antj = j;
		}
		coversize++;
	}
	// sort R array according to the list
	int laux[r];
	for (i = next[r], j = 0; i >= 0; i = next[i], j++)
		laux[j] = R[i];
	memcpy(R, laux, r*sizeof(int));
	int iaux[r];
	for (i = next[r], j = 0; i >= 0; i = next[i], j++)
		iaux[j] = cover[i];
	memcpy(cover, iaux, r*sizeof(int));

	return coversize;
}

int readWeightFile(FILE * graphFile, int n, double * weight) {
	char               type  = ' ';
	char               linestr[100];
	char *             datastr;
	int                i, j;

	while (type != 'w') {
		type = fgetc(graphFile);
		if (type != EOF) {

			/* header */
			if (type == 'c') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr != NULL)
					printf("%s", linestr);
				else
					return -1;
			}
		}
	}

	////
	// Vertex weights
	////

	for (i = 0; i < n && type != EOF; i++) {
		/* Edges */
		if (type == 'w') {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr == NULL)
				return -1;

			datastr = strtok(linestr," ");
			weight[i] = atof(linestr);
		}
		else {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr != NULL)
				printf(" %s\n", linestr);
			else
				return -1;
		}
		type = fgetc(graphFile);
	}

	return 0;
}

static double computeWeightedExantideg(const int * R, const double * weight, const int i, const double * antiwdeg, const int k) {
	double d = 0.0;
	int j;
	for (j = k; j > 0; j--)
		if (!hasEdge(R[i], R[j]))
			d += antiwdeg[R[j]];
	return d;
}

int mcrWeightedSort(const int n, int * R, double * weight, double * maxantiwdeg) {
	// number of antineighbors
	double antiwdeg[n];
	memset(antiwdeg, 0, n*sizeof(double));

	// maximum number of antineighbors (including itself)
	// used to complete vertex cover
	*maxantiwdeg = 0.0;
	int iu, iv;
	for (iu = 0; iu < n; iu++) {
		antiwdeg[iu] += weight[R[iu]];
		for (iv = iu+1; iv < n; iv++)
			if (iu != iv && !hasEdge(R[iu], R[iv])) {
				antiwdeg[iu] += weight[R[iv]];
				antiwdeg[iv] += weight[R[iu]];
			}
		if (antiwdeg[iu] > *maxantiwdeg)
			*maxantiwdeg = antiwdeg[iu];
	}

	// used to break ties
	double exantiwdeg[n];
	// a negative number indicates that this parameter is not known
	for (iu = 0; iu < n; iu++)
		exantiwdeg[iu] = -1.0;

	// determine min antideg
	int Rmin[n];
	int nrmin = 1;
	Rmin[0] = 0;
	double antiwdegrmin = antiwdeg[0];
	int i = n - 1;
	int j;
	for (j = 1; j <= i; j++)
		if (antiwdeg[j] < antiwdegrmin) {
			nrmin = 1;
			Rmin[0] = j;
			antiwdegrmin = antiwdeg[j];
		}
		else if (antiwdeg[j] == antiwdegrmin)
			Rmin[nrmin++] = j;

	int rmin;
	while (i >= nrmin) {
		rmin = 0;
		if (nrmin > 1) {
			for (j = 1; j < nrmin; j++) {
				if (exantiwdeg[Rmin[rmin]] < 0)
					exantiwdeg[Rmin[rmin]]=computeWeightedExantideg(R, weight, Rmin[rmin], antiwdeg, i);
				if (exantiwdeg[Rmin[j]] < 0)
					exantiwdeg[Rmin[j]]=computeWeightedExantideg(R, weight, Rmin[j], antiwdeg, i);
				if (exantiwdeg[Rmin[j]] < exantiwdeg[Rmin[rmin]])
					rmin = j;
			}
		}

		rmin = Rmin[rmin];
		if (i != rmin) {
			int iaux = R[i];
			R[i] = R[rmin];
			R[rmin] = iaux;
			double daux = antiwdeg[i];
			antiwdeg[i] = antiwdeg[rmin];
			antiwdeg[rmin] = daux;
			daux = exantiwdeg[i];
			exantiwdeg[i] = exantiwdeg[rmin];
			exantiwdeg[rmin] = daux;
		}

		for (j = 0; j < i; j++)
			if (!hasEdge(R[j], R[i])) {
				exantiwdeg[j] -= antiwdeg[i];
				antiwdeg[j] -= weight[R[i]];
			}

		i--;
		Rmin[0] = 0;
		nrmin = 1;
		antiwdegrmin = antiwdeg[0];
		for (j = 1; j <= i; j++)
			if (antiwdeg[j] < antiwdegrmin) {
				nrmin = 1;
				Rmin[0] = j;
				antiwdegrmin = antiwdeg[j];
			}
			else if (antiwdeg[j] == antiwdegrmin)
				Rmin[nrmin++] = j;
	}

	return i;
}

// R is an array of vertex indices. Weight of vertex R[i] is stored in weight[R[i]].
// sort R according to a clique cover of the subgraph induced by R
// it returns the weight of the clique cover and sets the cover array such that:
// for i in 0, ..., r-1, the entry cover[i] contains the weight of a clique
// cover of the subgraph induced by R[0], ..., R[i]
double coverWeightedSort(const int r, int * R, double * weight, double * cover) {
	double coversize = 0.0;
	int i, j;

	// linking the elements of R
	int next[r+1];
	// head of the candidates list
	int head = 0;
	for (i = head; i < r - 1; i++)
		next[i] = i + 1;
	next[i] = -1;
	// tail of the cliques list
	int tail = r;
	// head of the cliques list is r
	next[tail] = -1;

	int anti;
	for (i = head; i >= 0; i = head) {
		// remove i from the candidates list
		head = next[i];
		// include i in the tail of cliques list
		next[i] = next[tail];
		next[tail] = i;
		anti = tail;
		tail = i;

		// search for new vertices to add to the current clique
		int antj = head, nj;
		for (j = head, nj = next[j]; j >= 0; j = nj, nj = next[j]) {
			int cc;
			int vaj = R[j];
			for (cc = i; cc >= 0 && hasEdge(R[cc], vaj); cc = next[cc]);
			if (cc == -1) {
				// remove j from the candidates list
				if (j == head)
					head = next[j];
				else
					next[antj] = next[j];
				if (weight[R[j]] <= weight[R[i]]) {
					// include j in the tail of cliques list
					next[j] = next[tail];
					next[tail] = j;
					tail = j;
				}
				else {
					// include j in the head of current clique
					next[anti] = j;
					next[j] = i;
					i = j;
				}
			}
			else
				antj = j;
		}
		// maximum weight of current clique is given by its first vertex
		coversize += weight[R[i]];
		for (j = i; j >= 0; j = next[j])
			cover[j] = coversize;
	}
	// sort R array according to the list
	int laux[r];
	for (i = next[r], j = 0; i >= 0; i = next[i], j++)
		laux[j] = R[i];
	memcpy(R, laux, r*sizeof(int));
	double daux[r];
	for (i = next[r], j = 0; i >= 0; i = next[i], j++)
		daux[j] = cover[i];
	memcpy(cover, daux, r*sizeof(double));

	return coversize;
}

typedef struct E {
	int 		u;
	double 		xu;
	int 		v;
	double 		xv;
	struct E * 	next;
} edge;

typedef struct EP {
	int v;
	struct EP * next;
} endpoint;

typedef struct P {
	int   		u;
	int   		v;
	edge * 		falseedge;
	endpoint * 	removededge;
	endpoint *	commoneig;
	int			inset;
	struct P * 	next;
} projection;

typedef struct EJ {
	int 			v;
	projection * 	p;
	struct EJ *		puw;
	struct EJ *		pvz;
	struct EJ * 	next;
} eprojected;

static eprojected 	nulleprojected[1];
static endpoint 	nullendpoint[1];
static edge 		nulledge[1];

static void restoreGraph(int nfe, int * vind, eprojected * falseEdge) {
	int iu;
	for (iu = 0; iu < nfe; iu++) {
		eprojected * ep = falseEdge[iu].next;
		while (ep != nulleprojected) {
			endpoint * ed = ep->p->removededge;
			while (ed != nullendpoint) {
				vind[ep->p->u] < ed->v ? addEdge(vind[ep->p->u], ed->v) : addEdge(ed->v, vind[ep->p->u]);
				endpoint * fed = ed;
				ed = ed->next;
				free(fed);
			}
			ep->p->removededge = nullendpoint;

			delEdge(vind[iu], ep->v);
			eprojected * fep = ep;
			ep = ep->next;
			free(fep);
		}
		falseEdge[iu].next = nulleprojected;
	}
}

static int antiProjection(int *q, int * Q, int * vind, projection * p) {
	int cq = *q;
	int i;
//	printf("(%d,%d) is a projected edge\n", vind[p->u], vind[p->v]);
//	for (i = 0; i < *q && Q[i] != p->u; i++);
//	if (i == *q)
		Q[cq++] = p->u;
//	else
//		printf("POOOINNN-u\n");
//	for (i = 0; i < *q && Q[i] != p->v; i++);
//	if (i == *q)
		Q[cq++] = p->v;
//	else
//		printf("POOOINNN-v\n");
	endpoint * common = p->commoneig;
	while (common != nullendpoint) {
//		for (i = 0; i < *q && Q[i] != common->v; i++);
//		if (i == *q) {
//			printf("(%d) is a common neig vertex (i=%d q=%d cq=%d)\n", vind[common->v], i, *q, cq);
			Q[cq++] = common->v;
//		}
//		else
//			printf("POOOINNN-nuv\n");
		common = common->next;
	}

	*q = cq;

	return 0;
}

static eprojected * getFalseEdge(int js, int jg, int * vind, eprojected * falseEdge) {
	int w = vind[js];
	int z = vind[jg];

	eprojected * e = falseEdge[jg].next;
	if (hasEdge(z, w)) {
		while (e->v != w) {
			e = e->next;
			if (e == nulleprojected) {
				printf("Error in false edge list: %d is not in the %d's list.\n", w, z);
				exit(0);
			}
		}
	}
	else if (e != nulleprojected) {
		printf("UAI (%d,%d)! Já foi false edge?\n", w, z);
		while (e != nulleprojected && e->v != w) {
			e = e->next;
		}
	}

	return e;
}

static int recAntiProjection(eprojected * e, int *q, int * Q, int * alpha, int * vind, eprojected * falseEdge) {
	if (!e->p->inset) {
//		printf("false edges of projection of (%d,%d):", vind[e->p->u], vind[e->p->v]);
		//		edge * fe = e->p->falseedge;
		//		while (fe != nulledge) {
		//			printf(" (%d,%d)", fe->u, fe->v);
		//			fe = fe->next;
		//		}
		//		printf("\n");

		e->p->inset = 1;
		int cq = *q;
		antiProjection(q, Q, vind, e->p);
		int dq = *q;
		(*alpha)++;

		int i, j, jg, js, u, v, w, z;
		for (j = cq; j < dq; j++)
			for (i = 0; i < j; i++) {
				u = vind[Q[i]];
				v = vind[Q[j]];
				if (u > v) {
					w = v;
					z = u;
					jg = Q[i];
					js = Q[j];
				}
				else {
					z = v;
					w = u;
					jg = Q[j];
					js = Q[i];
				}
				if (hasEdge(z, w)) {
					e = getFalseEdge(js, jg, vind, falseEdge);
					if (e != nulleprojected) {
						//				printf("(%d,%d) is a false edge", z, w);
						//				fflush(NULL);
						//				printf(" of (%d,%d)\n", vind[e->p->u], vind[e->p->v]);
						//				fflush(NULL);
						recAntiProjection(e, q, Q, alpha, vind, falseEdge);
					}
					else {
						printf("UAI (%d,%d)! Nem edge, nem false edge?\n", w, z);
						fflush(NULL);
						exit(0);
					}
				}
				//		else
				//			printf("(%d,%d) is an edge\n", z, w);
			}
	}
//	else
//		printf("(%d,%d) is a projected edge, but considered already\n", vind[e->p->u], vind[e->p->v]);
	if (e->puw != nulleprojected)
		recAntiProjection(e->puw, q, Q, alpha, vind, falseEdge);
	if (e->pvz != nulleprojected)
		recAntiProjection(e->pvz, q, Q, alpha, vind, falseEdge);

	return 0;
}

static int lifting(int q, int * Q, int * alpha, int * vind, eprojected * falseEdge) {

	*alpha = 1;
	int jg, js, jw, jz, w, z, cq = q;
	for (jw = 0; jw < q; jw++)
		for (jz = jw + 1; jz < q; jz++) {
			int jg, js;
			if (vind[Q[jw]] > vind[Q[jz]]) {
				w = vind[Q[jz]];
				z = vind[Q[jw]];
				jg = jw;
				js = jz;
			}
			else {
				z = vind[Q[jz]];
				w = vind[Q[jw]];
				jg = jz;
				js = jw;
			}
			if (hasEdge(z, w)) {
				eprojected * e = getFalseEdge(Q[js], Q[jg], vind, falseEdge);
				if (e != nulleprojected) {
					recAntiProjection(e, &cq, Q, alpha, vind, falseEdge);
				}
				else {
					printf("UAI (%d,%d)! Nem edge, nem false edge?\n", w, z);
					fflush(NULL);
					exit(0);
				}
			}
		}

	return cq;
}

static int * 	greedy = NULL;
static double * greedyd = NULL;

static int cmpxd(const void * a, const void * b) {
	double diff = greedyd[*(int *)b] - greedyd[*(int *)a];
	return diff < 0 ? -1 : 1;
}

static inline void swap(int a, int b, int * vind, double * xval) {
	int iaux;
	double daux;

	iaux = vind[a];
	vind[a] = vind[b];
	vind[b] = iaux;
	daux = xval[a];
	xval[a] = xval[b];
	xval[b] = daux;
}

// vind and xval: subgraph
// elements of vind and xval may be rearranged
// cutvind and cutrhs: cuts generated, not more than maxncut
// return: number of cuts returned
int edgeProjSeparation(int nzcnt, int * vind, double * xval, int maxncut, int * cutsize, int ** cutvind, double * cutrhs) {
	int ncut = 0;
	int j;

	int riu, riv, iu, iv;
	eprojected * falseEdge = NULL;
	projection projlist[1];

	// edge projection for each subgraph connecting edge
	int step;
	int removed[nzcnt];
	int firstiu = 0;
	int * Q = cutvind[ncut];
	for (step = 0; step < MAXITER && firstiu >= 0 && ncut < maxncut; step++) {
//		printf("STEP %d out of %d, size %d\n", step, MAXITER, nzcnt);
		memset(removed, 0, nzcnt*sizeof(int));
		int nremainv = nzcnt;
		projlist->next = projlist;
		int h = 0;
//		swap(firstiu, firstiu != 0 ? 0 : drand48()*nzcnt, vind, xval);
		int attemptsu;
//		printf("STEP=%d FIRSTIU=%d\n", step, (firstiu+1) % nzcnt);
		for (riu = (firstiu+1) % nzcnt, firstiu = -1, attemptsu = 0; attemptsu < nzcnt && nremainv > 0 && ncut < maxncut; (riu = (riu+1) % nzcnt), attemptsu++) {
			int foundh = 0;
			if (!removed[riu]) {
				int attemptsv;
				for (riv = ((riu+1) % nzcnt), attemptsv = 0; attemptsv + attemptsu < nzcnt && ncut < maxncut; (riv = (riv+1) % nzcnt), attemptsv++)
					if (!removed[riv] && (hasEdge(vind[riu], vind[riv]) || hasEdge(vind[riv], vind[riu])) && xval[riu] + xval[riv] > 0.6) {
						iu = riu;
						iv = riv;

//						printf("edge (%d,%d) to project with weights: %6.3f and %6.3f\n", vind[iu], vind[iv], xval[iu], xval[iv]);
						if (falseEdge == NULL) {
							falseEdge = (eprojected *) calloc(nzcnt, sizeof(eprojected));

							for (j = 0; j < nzcnt; j++)
								falseEdge[j].next = nulleprojected;
						}

						projection * p = (projection *) malloc(sizeof(projection));
						p->next = projlist->next;
						projlist->next = p;
						p->falseedge = nulledge;
						p->removededge = nullendpoint;
						p->commoneig = nullendpoint;
						p->u = iu;
						p->v = iv;

						int u = vind[iu];
						int v = vind[iv];
						removed[iu] = 1;
						removed[iv] = 1;
						nremainv -= 2;

						// neighborhoods and degrees
						int nu = 0, nv = 0;
						int neig[nzcnt << 1];
						int * neigu = neig;
						for (j = 0; j < nzcnt; j++)
							if (!removed[j] && (hasEdge(u, vind[j]) || hasEdge(vind[j], u)))
								neigu[nu++] = j;
						int * neigv = &neigu[nu];
						for (j = 0; j < nzcnt; j++)
							if (!removed[j] && (hasEdge(v, vind[j]) || hasEdge(vind[j], v)))
								neigv[nv++] = j;

						// u has to have the smallest degree
						if (nv < nu) {
							int aux;
							aux = iu;
							iu = iv;
							iv = aux;
							aux = u;
							u = v;
							v = aux;
							aux = p->u;
							p->u = p->v;
							p->v = aux;
							aux = nu;
							nu = nv;
							nv = aux;
							int * neigaux = neigu;
							neigu = neigv;
							neigv = neigaux;
						}

						// maximal clique in the neighborhood of u
						int jw, w, z;
						double dgu[nu];
						memset(dgu, 0, nu*sizeof(int));
						for (j = 0; j < nu; j++) {
							dgu[j] = xval[neigu[j]];
							z = vind[neigu[j]];
							for (jw = 0; jw < nu; jw++) {
								w = vind[neigu[jw]];
								if (hasEdge(z, w) || hasEdge(w, z))
									dgu[j]+=xval[neigu[j]];
							}
						}
						double edgu[nu];
						for (j = 0; j < nu; j++) {
							z = vind[neigu[j]];
							edgu[j] = 1;
							for (jw = 0; jw < nu; jw++) {
								w = vind[neigu[jw]];
								if (hasEdge(z, w) || hasEdge(w, z))
									edgu[j]+=dgu[jw];
							}
						}
						for (j = 0; j < nu; j++) {
							int jmax=j;
							for (jw = j+1; jw < nu; jw++) {
								if (dgu[jw] > dgu[jmax] || (dgu[jw] == dgu[jmax] && edgu[jw] > edgu[jmax]))
									jmax = jw;
							}
							if (jmax != j) {
								int aa;
								aa = neigu[j];
								neigu[j] = neigu[jmax];
								neigu[jmax] = aa;
								double daa = dgu[j];
								dgu[j] = dgu[jmax];
								dgu[jmax] = daa;
							}
						}

						double xq = 0.0;
						int qu = 0;
						for (j = 0; j < nu; j++) {
							z = vind[neigu[j]];
							for (jw = 0; jw < qu; jw++) {
								w = vind[Q[jw]];
								if (!hasEdge(z, w) && !hasEdge(w, z))
									break;
							}
							if (jw == qu) {
								Q[qu++] = neigu[j];
								xq += xval[neigu[j]];
							}
							else {
								// be careful: a false edge may be deleted
								if ((u < z && hasEdge(u, z)) || (z < u && hasEdge(z, u))) {
									endpoint * removed = (endpoint *) malloc(sizeof(endpoint));
									removed->next = p->removededge;
									p->removededge = removed;
									removed->v = z;
									neigu[j--] = neigu[--nu];

									if (u < z)
										delEdge(u, z);
									else
										delEdge(z, u);
								}
								else {
									// a false edge is deleted
									delEdge(u, z);
									delEdge(z, u);
//									printf("remove false edge: (%d,%d)\n", u,z);
//									printf("removed, but false edges of (%d):", u < z ? z : u);
//									eprojected * e = falseEdge[u < z ? neigu[j] : iu].next;
//									while (e != nulleprojected) {
//										printf(" (%d)", e->v);
//										e = e->next;
//									}
//									printf("\n");
								}
							}
						}
						// projection: removing common neighborhood
						for (j = 0; j < qu; j++)
							if (hasEdge(vind[Q[j]], v) || hasEdge(v, vind[Q[j]])) {
								endpoint * common = (endpoint *) malloc(sizeof(endpoint));
								common->next = p->commoneig;
								p->commoneig = common;
								common->v = Q[j];
								removed[Q[j]] = 1;
								nremainv--;
								xq -= xval[Q[j]];
								Q[j--] = Q[--qu];
							}

						// projection: including false edges
						h++;
						int qv = qu;
						for (j = 0; j < nv; j++) {
							z = vind[neigv[j]];
							if (!hasEdge(z, u) && !hasEdge(u, z)) {
								for (jw = 0; jw < qu; ++jw) {
									w = vind[Q[jw]];
									// w is a u's neighbor which has been removed if v's neighbor as well
									if (!hasEdge(w, z) && !hasEdge(z, w)) {
										edge * false = (edge *) malloc(sizeof(edge));
										false->next = p->falseedge;
										p->falseedge = false;
										eprojected * eproj = (eprojected *) malloc(sizeof(eprojected));
										eproj->p = p;
										if (u > w) {
											eproj->puw = falseEdge[iu].next;
											while (eproj->puw != nulleprojected && eproj->puw->v != w)
												eproj->puw = eproj->puw->next;
										}
										else {
											eproj->puw = falseEdge[Q[jw]].next;
											while (eproj->puw != nulleprojected && eproj->puw->v != u)
												eproj->puw = eproj->puw->next;
										}
										if (v > z) {
											eproj->pvz = falseEdge[iv].next;
											while (eproj->pvz != nulleprojected && eproj->pvz->v != z)
												eproj->pvz = eproj->pvz->next;
										}
										else {
											eproj->pvz = falseEdge[neigv[j]].next;
											while (eproj->pvz != nulleprojected && eproj->pvz->v != v)
												eproj->pvz = eproj->pvz->next;
										}

										if (z < w) {
											addEdge(w, z);
											false->u = z;
											false->xu = xval[neigv[j]];
											false->v = w;
											false->xv = xval[Q[jw]];
											eproj->v = z;
											eproj->next = falseEdge[Q[jw]].next;
											falseEdge[Q[jw]].next = eproj;
										}
										else {
											addEdge(z, w);
											false->u = w;
											false->xu = xval[Q[jw]];
											false->v = z;
											false->xv = xval[neigv[j]];
											eproj->v = w;
											eproj->next = falseEdge[neigv[j]].next;
											falseEdge[neigv[j]].next = eproj;
										}
//										printf("included in false edges of (%d):", w < z ? z : w);
//										eprojected * e = falseEdge[w < z ? neigv[j] : Q[jw]].next;
//										while (e != nulleprojected) {
//											printf(" (%d)", e->v);
//											e = e->next;
//										}
//										printf("\n");
									}
								}
								if ((h % K) || foundh) {
									for (; jw < qv; ++jw) {
										w = vind[Q[jw]];
										if (!hasEdge(z, w) && !hasEdge(w, z))
											break;
									}
									if (jw == qv) {
										Q[qv++] = neigv[j];
										xq += xval[neigv[j]];
									}
								}
							}
						}
						int q;
						if (!(h % K) && !foundh) {
//								printf("ENUMERATE CLIQUES, h=%d step=%d.\n", h, step);
								if (!greedy) {
									greedy = (int *) calloc(nzcnt, sizeof(int));
									greedyd = (double *) calloc(nzcnt, sizeof(double));
								}
								int sgy = 0, l;
								for (j = 0; j < nzcnt; j++)
									if (!removed[j]) {
										greedy[sgy++] = j;
										greedyd[j] = xval[j];
										for (l = 0; l < nzcnt; l++)
											if (!removed[l] && (hasEdge(vind[j], vind[l]) || hasEdge(vind[l], vind[j])))
												greedyd[j] += xval[l];
									}
								if (sgy > 0) {
									qsort(greedy, sgy, sizeof(int), cmpxd);
									int gy = 0;
									double xq1, lxq = 0.0;
									int a, b;
									do {
										q = 0;
										xq = 0.0;
										for (j = 0; j < sgy; j++) {
											for (l = 0; l < q && (hasEdge(vind[Q[l]], vind[greedy[j]]) || hasEdge(vind[greedy[j]], vind[Q[l]])); l++);
											if (l == q) {
												Q[q++] = greedy[j];
												xq += xval[greedy[j]];
											}
										}
										xq1 = xq;

										if ((xq < lxq - 1.0e-10 || xq > lxq + 1.0e-10) && xq > 1.0 + 1.0e-10) {
											if (firstiu < 0)
												firstiu = riu;

											projection * p = projlist->next;
											while (p != projlist) {
												p->inset = 0;
												p = p->next;
											}

											int alpha;
											int cq = lifting(q, Q, &alpha, vind, falseEdge);
											for (j = q; j < cq; j++)
												xq += xval[Q[j]];
//											printf("new value %8.4f    alpha = %d\n", xq, alpha);
											if (xq > alpha + MINVIOLATION) {
												cutrhs[ncut] = alpha;
												cutsize[ncut++] = cq;
												foundh = 1;
												if (ncut < maxncut)
													Q = cutvind[ncut];
											}
										}

										if (gy < MAXGREEDYITER && ncut < maxncut) {
											if (xq1 <= lxq) {
												int aux = greedy[a];
												greedy[a] = greedy[b];
												greedy[b] = aux;
											}
											else
												lxq = xq1;
											do {
												a = drand48()*sgy;
												b = drand48()*sgy;
											} while (a == b);
											int aux = greedy[a];
											greedy[a] = greedy[b];
											greedy[b] = aux;
										}
									} while (gy++ < MAXGREEDYITER && ncut < maxncut);
								}
							}
						else {
							if (!(h % K))
								foundh = 0;

							q = qv;
							for (j = 0; j < nzcnt; j++) {
								z = vind[j];
								if (!removed[j] && !hasEdge(z, u) && !hasEdge(u, z) && !hasEdge(z, v) && !hasEdge(v, z)) {
									for (jw = 0; jw < q; ++jw) {
										w = vind[Q[jw]];
										if (!hasEdge(z, w) && !hasEdge(w, z))
											break;
									}
									if (jw == q) {
										Q[q++] = j;
										xq += xval[j];
									}
								}
							}

							if (xq > 1.0 + 1.0e-10) {
								if (firstiu < 0)
									firstiu = riu;

								projection * p = projlist->next;
								while (p != projlist) {
									p->inset = 0;
									p = p->next;
								}

								int alpha;
								int cq = lifting(q, Q, &alpha, vind, falseEdge);
								for (j = q; j < cq; j++)
									xq += xval[Q[j]];
	//							printf("new value %8.4f    alpha = %d\n", xq, alpha);
								if (xq > alpha + MINVIOLATION) {
									cutrhs[ncut] = alpha;
									cutsize[ncut++] = cq;
									foundh = 1;
									if (ncut < maxncut)
										Q = cutvind[ncut];
								}
							}
						}
						attemptsv = nzcnt;
					}
			}
		}

		if (falseEdge != NULL) {
			restoreGraph(nzcnt, vind, falseEdge);
			projection * p = projlist->next;
			while (p != projlist) {
				projection * fp = p;
				p = p->next;
				free(fp);
			}
		}
	}

	if (falseEdge != NULL) {
		free(falseEdge); // ainda falta otimizar toda a memória reservada
	}

	if (greedy) {
		free(greedy);
		greedy = NULL;
		free(greedyd);
		greedyd = NULL;
	}

	return ncut;
}
