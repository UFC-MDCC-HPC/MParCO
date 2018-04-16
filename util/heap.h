/*
 * heap.h
 *
 *  Created on: Feb 17, 2011
 *      Author: Corrêa
 */
 
#ifndef HEAP_H_
#define HEAP_H_

// Partially sort the specified array of <tt>nmemb</tt> elements of size <tt>size</tt> in a heap manner.  
// The <tt>base</tt> argument points to the start of the array.
// The  contents  of the array are partially sorted in ascending order according to a comparison function pointed to by <tt>compar</tt>, 
// which is called with two arguments that point to the objects being compared.
// The comparison function must return an integer less than, equal to, or greater than zero if the first argument is considered  to  be  respectively
// less than, equal to, or greater than the second.  If two members compare as equal, their order in the sorted array is undefined.
void heapify(void *base, size_t nmemb, size_t size, int(*compar)(const void *, const void *));

// Add a new element <tt>x</tt> to a specified array, which is assumed to be partially sorted in a heap manner.
// It is also assumed that the specified array has enough place for the new element.
// The remaining arguments are similar to the ones of <tt>heapify</tt>.
void heapoffer(void *x, void *base, size_t nmemb, size_t size, int(*compar)(const void *, const void *));

void heappoll(void *base, size_t nmemb, size_t size, int(*compar)(const void *, const void *));

#endif /*HEAP_H_*/
