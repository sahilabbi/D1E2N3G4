#ifndef __ISOMORPHISM__

#define __ISOMORPHISM__
#include <stdbool.h>

#include "graph.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


//Used to avoid rounding errors
bool float_equals(double a, double b);

gsl_matrix * graph_to_matrix(graph_t * g);

//A permutation set is a set of indices from
//the first graph (leftSet) to the
//second graph (rightSet)
//size gives the number of indicies
typedef struct
{
    int * leftSet;
    int * rightSet;
    int size;
} PermutationSet;

typedef struct
{
    PermutationSet * permSets;
    int num_perm_sets;
} PartialPermutation;

//combines both permutations and stores the result in acc.  perm is destroyed in the process
void combinePartialPermutations(PartialPermutation * acc, PartialPermutation * perm);

PartialPermutation * getPermutationFromVectors(gsl_vector * left, gsl_vector * right);

PartialPermutation * isomorphism(graph_t * g1, graph_t * g2);

bool is_isomorphic(graph_t * g1, graph_t * g2);

#endif
