#include "isomorphism.h"

int * matrix_vector_mult(matrix * m, float * v){
    float * ret = (float *) malloc(N * sizoef(float));
    int i,j;
    for(i = 0; i < N; i++){
	int sum = 0;
	for(j = 0; j < N; j++){
	    sum += matrix->values[i][j] * v[j];
	}
	ret[i] = sum;
    }
    return ret;
}



bool accurate_isomorphism(graph_t * g1, graph_t * g2){
}
