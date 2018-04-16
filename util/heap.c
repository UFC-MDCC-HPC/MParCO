/*
 * heap.c
 *
 *  Created on: Feb 17, 2011
 *      Author: Corrêa
 */
 
#include <stdio.h>
#include <string.h>
#include <heap.h>

typedef unsigned char byte;
static int ls;

void heapoffer(void *x, void *b, size_t nmemb, size_t size, int(*compar)(const void *, const void *)) {
	int is = nmemb*size;
	byte * base = (byte *) b;
	memcpy(base+is, x, size);
	int f = nmemb+1;
	int fs = f & 0x01 ? (is >> 1) - size: (is >> 1) - (size >> 1);
	f = f >> 1;
	while (is > 0 && compar(x, base+fs) < 0) {
		memcpy(base+is, base+fs, size);
		is = fs;
		fs = f & 0x01 ? (is >> 1) - size: (is >> 1) - (size >> 1);
		f = f >> 1;
	}
	memcpy(base+is, x, size);
}

static void down(int is, void *aux, byte *base, size_t nmemb, size_t size, int(*compar)(const void *, const void *)) {
	int os, ss;
	int pls = nmemb & 0x01 ? (ls >> 1) - size: (ls >> 1) - (size >> 1);
	do {
		os = (is << 1) + size;
		ss = os <= ls && (compar(base+os, base+is) < 0) ? os : is;
		os += size;
		ss = os <= ls && (compar(base+os, base+ss) < 0) ? os : ss;
		if (ss != is) {
			memcpy(aux, base+is, size);
			memcpy(base+is, base+ss, size);
			memcpy(base+ss, aux, size);
			is = ss;
		}
		else
			is = pls;
	} while (is < pls);
}

void heappoll(void *b, size_t nmemb, size_t size, int(*compar)(const void *, const void *)) {
	byte * base = (byte *) b;
	ls = (nmemb-2)*size;
	memcpy(base, base+ls+size, size);
	down(0, base+ls+size, base, nmemb-1, size, compar);
}

void heapify(void *b, size_t nmemb, size_t size, int(*compar)(const void *, const void *)) {
	char aux[size];
	ls = (nmemb-1)*size;
	int is = nmemb & 0x01 ? (ls >> 1) - size: (ls >> 1) - (size >> 1);
	for (; is >= 0; is -= size)
		down(is, aux, (byte *) b, nmemb, size, compar);
}
