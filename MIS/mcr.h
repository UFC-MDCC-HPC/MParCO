/*
 * mcr.h
 *
 *  
 *  Author: Allberson/Corrêa
 */

#ifndef MCR_H_
#define MCR_H_
#include <graph.h>

int getStabSet(adjMatrix * g, int n, int * vset, int qq, double ldens, double udens);
double getWeightedStabSet(adjMatrix * g, int nvset, int * vset, double * w, double qq, double ldens, double udens);

#endif /* MCR_H_ */
