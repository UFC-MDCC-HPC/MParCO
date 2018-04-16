/*
 * repository.h
 *
 *  Created on: Nov 17, 2010
 *      Author: Corrêa
 */

#ifndef REPOSITORY_H_
#define REPOSITORY_H_
#include <stdlib.h>

// new repository initialization
// elemsz: size, in bytes, of each element
// step: number of elements allocated whenever repository is empty
// return: a pointer to be used with repository functions
void * newRepo(size_t elemsz, size_t step);

// get an element from a specified repository
void * getFromRepo(void * r);

// release a specified element using the specified repository
void release(void * r, void * b);

// free memory space used by a specified repository
void freeRepo(void * r);

// free memory space used to manage repositories
int freeBaseRepo();

#endif /* REPOSITORY_H_ */
