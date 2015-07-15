#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

matrix * create_matrix(int size){
    matrix * ret = (matrix *) malloc(sizeof(matrix));
    ret->size = size;
    ret->values = (float *) malloc(size * size * sizeof(float));
    memset(ret->values, 0, size * size * sizeof(float));
    return ret;
}

vector * create_vector(int size){
    vector * ret = (vector *) malloc(sizeof(vector));
    ret->size = size;
    ret->values = (float *) malloc(size * sizeof(float));
    memset(ret->values, 0, size * sizeof(float));
    return ret;
}

vector * matrix_vector_mult(matrix * m, vector * v){
    if(m->size != v->size){
	printf("Error: wrong sizes for matrix-vector multiplication\n");
	return 0;
    }
    vector * ret = create_vector(m->size);
    int i,j;
    for(i = 0; i < m->size; i++){
	int sum = 0;
	for(j = 0; j < m->size; j++){
	    sum += m->values[i*m->size + j] * v->values[j];
	}
	ret->values[i] = sum;
    }
    return ret;
}

matrix * matrix_matrix_mult(matrix * m1, matrix * m2){
    if(m1->size != m2->size){
	printf("Error: wrong sizes for matrix-matrix multiplication\n");
	return 0;
    }
    int size = m1->size;
    matrix * ret = create_matrix(size);
    int i,j,k;
    for(i = 0; i < size; i++){
	for(j = 0; j < size; j++){
	    for(k = 0; k < size; k++){
		ret->values[i * ret->size + j] += m1->values[i * m1->size + k] *
		                                  m2->values[k * m2->size + j];
	    }
	}
    }
    return ret;
}



bool accurate_isomorphism(graph_t * g1, graph_t * g2){
}

void print_matrix(matrix * m){
    int i,j;
    for(i = 0; i < m->size; i++){
	for(j = 0; j < m->size; j++){
	    printf("%f ", m->values[i * m->size + j]);
 	}
	printf("\n");
    }
}

void print_vector(vector * v){
    int i;
    for(i = 0; i < v->size; i++){
	printf("%f ", v->values[i]);
    }
    printf("\n");
}

#ifdef MATRIX_TEST

int main(){

	matrix * a = create_matrix(2);
	matrix * b = create_matrix(2);

	a->values[0] = 1; a->values[1] = 0;
	a->values[2] = 2; a->values[3] = 1;

	b->values[0] = 5; b->values[1] = 4;
	b->values[2] = -5; b->values[3] = 1;

	printf("AB:\n");
	print_matrix(matrix_matrix_mult(a,b));
	printf("\nBA:\n");
	print_matrix(matrix_matrix_mult(b,a));


    return 0;
}

#endif

