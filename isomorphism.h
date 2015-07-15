#ifndef __ISOMORPHISM__

#define __ISOMORPHISM__
#include <stdbool.h>

#include "graph.h"

typedef struct
{
    float * values;
    //all matrices are square, so only one size is needed
    int size;
} matrix;

matrix * create_matrix(int size);


typedef struct
{
    float * values;
    int size;
} vector;

vector * create_vector(int size);

matrix * graph_to_matrix(graph_t * g);

vector * eigenvalues(matrix * m);
vector * matrix_vector_mult(matrix * g, vector * v);

//for the power iteration method
float getEigenvalue(matrix * m);
float remove_eigenvalue(matrix * m, vector * eigenvector, float eigenvalue);

bool accurate_isomorphism(graph_t * g1, graph_t * g2);

void print_matrix(matrix * m);
void print_vector(vector * v);
#endif
