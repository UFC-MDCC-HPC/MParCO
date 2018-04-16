/*
 * repository.c
 *
 *  Created on: Nov 17, 2010
 *      Author: correa
 */
#include <repository.h>
#include <stdio.h>

typedef unsigned char byte;

// node of a repository's chained list
typedef struct repE {
	void        * elem;
	struct repE * next;
} repoNode;

// a repository
typedef struct repo {
	repoNode * head;     // of its chained list
	repoNode * buf;      // of the buffers' chained list
	size_t     elemsz;   // size of each element
	size_t     step;     // number of elements allocated when this repository is empty
} repository;

// repository of repoNode
static repository nodeRepo = { NULL, NULL, sizeof(repoNode), -1 };
static int repocount = 0;

// creates a new repository
void * newRepo(size_t elemsz, size_t step) {
	repository * r = (repository *) malloc(sizeof(repository));
	r->head = NULL;
	r->buf = NULL;
	r->elemsz = elemsz;
	r->step = step;

	repocount++;

	return (void *) r;
}

// creates new repoNodes, assigning specified elements to them
static inline repoNode * increaseNodeRepo(size_t nelem, size_t elemsz, byte * base) {
	if (nelem <= 0)
		return NULL;

	repoNode * b = (repoNode *) malloc((nelem + 1)*sizeof(repoNode));
	b->next = nodeRepo.buf;
	nodeRepo.buf = b;
	int j;
	size_t displacement = 0;
	for (j = 1; j < nelem; j++) {
		b[j].next = &b[j + 1];
		b[j].elem = base+displacement;
		displacement += elemsz;
	}
	b[j].next = NULL;
	b[j].elem = base+displacement;

	return &b[1];
}

void * getFromRepo(void * r) {
	repository * repo = (repository *) r;
	repoNode * l;
	int i;
	if (repo->head == NULL) {
		byte * u = (byte *) malloc(repo->step*repo->elemsz);
		// get an empty nodeRepo to store this new buffer u
		if (nodeRepo.head == NULL)
			nodeRepo.head = increaseNodeRepo(repo->step+1, 0, NULL);
		// store new buffer
		l = nodeRepo.head;
		nodeRepo.head = nodeRepo.head->next;
		l->elem = u;
		l->next = repo->buf;
		repo->buf = l;

		// get empty repoNodes to store new elements from u
		i = 0;
		size_t displac = 0;
		l = NULL;
		repo->head = nodeRepo.head;
		while (nodeRepo.head != NULL && i < repo->step) {
			l = nodeRepo.head;
			nodeRepo.head = nodeRepo.head->next;
			l->elem = u+displac;
			displac += repo->elemsz;
			i++;
		}
		if (i < repo->step) {
			if (l != NULL)
				l->next = increaseNodeRepo(repo->step - i, repo->elemsz, u+displac);
			else
				repo->head = increaseNodeRepo(repo->step - i, repo->elemsz, u+displac);
		}
		else
			l->next = NULL;
	}
	l = repo->head;
	repo->head = l->next;
	// repoNode becomes empty
	l->next = nodeRepo.head;
	nodeRepo.head = l;

	return l->elem;
}

void release(void * r, void * b) {
	repository * repo = (repository *) r;
	// get an empty repoNode
	if (nodeRepo.head == NULL)
		nodeRepo.head = increaseNodeRepo(repo->step, 0, NULL);
	repoNode * e = nodeRepo.head;
	nodeRepo.head = nodeRepo.head->next;

	e->elem = b;
	e->next = repo->head;
	repo->head = e;
}

void freeRepo(void * r) {
	repository * repo = (repository *) r;
	repoNode * b;
	while (repo->buf != NULL) {
		b = repo->buf;
		repo->buf = repo->buf->next;
		free(b->elem);
		b->next = nodeRepo.head;
		nodeRepo.head = b;
	}
	repo->head = NULL;

	repocount--;
}

// returns number of active repositories
int freeBaseRepo() {
	if (repocount > 0)
		return repocount;

	repoNode * b;
	while (nodeRepo.buf != NULL) {
		b = nodeRepo.buf;
		nodeRepo.buf = nodeRepo.buf->next;
		free(b);
	}
	nodeRepo.head = NULL;

	return 0;
}
