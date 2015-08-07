#ifndef __ISOMORPHISM__

#define __ISOMORPHISM__
#include <stdbool.h>

#include "graph.h"


//Used to avoid rounding errors
bool float_equals(float a, float b);

typedef struct
{
    float * values;
    //all matrices are square, so only one size is needed
    int size;
} matrix;

matrix * create_matrix(int size);
void delete_matrix(matrix * m);

typedef struct
{
    float * values;
    int size;
} vector;

vector * create_vector(int size);
void delete_vector(vector * v);

typedef struct
{
    vector ** vectors;
    int num_vectors;
    //the size of the vectors in the space
    int size;
} vector_space;
vector_space * create_vector_space(int num_vectors, int size);
void delete_vector_space(vector_space * vspace);

matrix * graph_to_matrix(graph_t * g);

matrix * matrix_transpose(matrix * m);

vector * eigenvalues(matrix * m);
vector * matrix_vector_mult(matrix * g, vector * v);

float dot_product(vector * v1, vector * v2);
//multiplies the vector by a scalar to make it a unit vector
void normailize(vector * v);
//projects the vector v onto the space spanned by the vectors in vspace
//vspacesize gives the number of vectors in vspace
//normalized states whether the vectors in vspace have been normalized
vector * project_vector(vector * v, vector_space * vspace,  bool normailized);
//used for when the eigenspace has dimension greater that 1
vector * angles_to_standard_basis(vector_space * vspace,  bool normalized);

matrix * reduced_echelon_form(matrix * m);

matrix * householder(matrix * m);

typedef struct
{
    float * diagonal;
    matrix * eigenvectors;
} diagnolization;

diagnolization * diagnolize(matrix * m);


bool accurate_isomorphism(graph_t * g1, graph_t * g2);

void print_matrix(matrix * m);
void print_vector(vector * v);
#endif
