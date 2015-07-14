#ifndef __ISOMORPHISM__

#define __ISOMORPHISM__
#include <stdbool.h>

#include "graph.h"

typedef struct
{
    float values[N][N];
} matrix;

matrix * graph_to_matrix(graph_t * g);

float * eigenvalues(matrix * m);
float * matrix_vector_mult(matrix * g, float * v);

//for the power iteration method
float getEigenvalue(matrix * m);
float remove_eigenvalue(matrix * m, float * eigenvector, float eigenvalue);

bool accurate_isomorphism(graph_t * g1, graph_t * g2);


#endif
