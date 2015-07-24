#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.001

bool float_equals(float a, float b){
    float diff = a - b;
    if(diff < 0) diff = -diff;
    return diff < EPSILON;
}

matrix * create_matrix(int size){
    matrix * ret = (matrix *) malloc(sizeof(matrix));
    ret->size = size;
    ret->values = (float *) malloc(size * size * sizeof(float));
    memset(ret->values, 0, size * size * sizeof(float));
    return ret;
}

void delete_matrix(matrix * m){
    free(m->values);
    free(m);
}

vector * create_vector(int size){
    vector * ret = (vector *) malloc(sizeof(vector));
    ret->size = size;
    ret->values = (float *) malloc(size * sizeof(float));
    memset(ret->values, 0, size * sizeof(float));
    return ret;
}

void delete_vector(vector * v){
    free(v->values);
    free(v);
}

vector_space * create_vector_space(int num_vectors, int size){
    vector_space * ret = (vector_space *) malloc(sizeof(vector_space));
    ret->vectors = (vector **) malloc(num_vectors * sizeof(vector *));
    int i;
    for(i = 0; i < num_vectors; i++) ret->vectors[i] = create_vector(size);
    ret->size = size;
    return ret;
}

matrix * matrix_transpose(matrix * m){
    matrix * ret = create_matrix(m->size);
    int i,j;
    for(i = 0; i < ret->size; i++){
	for(j = 0; j < ret->size; j++){
	    ret->values[j * ret->size + i] = m->values[i * m->size + j];
	}
    }
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

float dot_product(vector * v1, vector * v2){
    if(v1->size != v2->size){
	printf("Error: invalid sizes in dot product\n");
	return 0;
    }
    int i;
    float ret = 0;
    for(i = 0; i < v1->size; i++) ret += v1->values[i] * v2->values[i];
    return ret;
    
}

void normalize(vector * v){
    float length = sqrt(dot_product(v,v));
    float correction = 1 / length;
    int i;
    for(i = 0; i < v->size; i++) v->values[i] *= correction;
}

vector * project_vector(vector * v, vector_space * vspace, bool normalized){
    if(vspace->size != v->size){
	printf("Error: projecting a vector onto a space of incorrect size\n");
	return 0;
    }
    vector * ret = create_vector(v->size);
    int i,j;
    for(i = 0; i < vspace->num_vectors; i++){
	//the length of the projection vector
	float length = dot_product(v,vspace->vectors[i]);
	//if the vector is not a unit vector, correct by dividing out the length of the vector
	if(!normalized) length /= dot_product(vspace->vectors[i],vspace->vectors[i]);
	for(j = 0; j < ret->size; j++) ret->values[j] += vspace->vectors[i]->values[j] * length;
    }
    return ret;
}

vector * angles_to_standard_basis(vector_space * vspace,  bool normalized){
    int i;
    vector * ret = create_vector(vspace->size);
    vector * basis_vector = create_vector(vspace->size);
    for(i = 0; i < ret->size; i++){
	if(i > 0) basis_vector->values[i-1] = 0;
	basis_vector->values[i] = 1;
	vector * proj = project_vector(basis_vector, vspace,  normalized);
	normalize(proj);
	ret->values[i] = dot_product(basis_vector, proj);
	delete_vector(proj);
    }
    delete_vector(basis_vector);
    return ret;
}

//Functions to apply row operations
static void row_swap(matrix * m, int row1, int row2){
    int i;
    float temp;
    for(i = 0; i < m->size; i++){
	temp = m->values[row1 * m->size + i];
	m->values[row1 * m->size + i] = m->values[row2 * m->size + i];
	m->values[row2 * m->size + i] = temp;
    }
}

static void row_mult(matrix * m, int row, float multiplier){
    int i;
    for(i = 0; i < m->size; i++) m->values[row * m->size + i] *= multiplier;
}

static void row_mult_and_add(matrix * m, int row_from, int row_to, float multiplier){
    int i;
    for(i = 0; i < m->size; i++){
	m->values[row_to * m->size + i] += m->values[row_from * m->size + i] * multiplier;
    }
}

#ifndef ABS
#define ABS(A) ((A) < 0 ? -(A) : (A))
#endif

matrix * row_reduced_echelon_form(matrix * m){
    matrix * ret = create_matrix(m->size);
    int i;
    for(i = 0; i < ret->size * ret->size; i++) ret->values[i] = m->values[i];
    //Throughout the process, anything below the current row and
    //column will be in row reduced echelon form
    int row = 0;
    int column = 0;
    while(row < ret->size && column < ret->size){
	//get the largest element in column
	float max_value = 0;
	int max_value_index = 0;
	for(i = row; i < m->size; i++){
	    float value_temp = ret->values[i * m->size + column];
	    if(value_temp < 0) value_temp = -value_temp;
	    if(value_temp > max_value){
		max_value = value_temp;
		max_value_index = i;
	    }
	}
	//if the columnn was all 0s, move on to the next columnn and retry
	if(float_equals(max_value, 0)){
	    column++;
	    continue;
	}
	//swap the row with the largest leading entry to just below the already
	//reduced rows
	if(row != max_value_index) row_swap(ret, row, max_value_index);

	//Turn the leading entry in the row into a 1
	row_mult(ret, row, 1 / ret->values[row * ret->size + column]);
	
	//remove all non-zero entries in the column
	for(i = 0; i < m->size; i++){
	    //Don't zero out the row itself
	    if(i == row) continue;
	    //and don't zero out anything that's already zeroed out
	    if(float_equals(ret->values[i * ret->size + column], 0)) continue;
	    row_mult_and_add(ret, row, i, -ret->values[i * ret->size + column]);
	}
	row++;
	column++;
	/*printf("Matrix after zeroing out:\n");
	print_matrix(ret);
	printf("\n");*/
    }
    return ret;
}

#ifndef SIGN
#define SIGN(A) ((A) < 0 ? -1 : 1)
#endif

matrix * householder(matrix * m){
    matrix * ret = create_matrix(m->size);
    int i,j,k;
    matrix * p = create_matrix(m->size);
    for(i = 0; i < ret->size; i++)
	for(j = 0; j < ret->size; j++)
	    ret->values[i * ret->size + j] = m->values[i * m->size + j];
    for(i = 0; i < ret->size - 1; i++){
	float alpha = -SIGN(ret->values[(i+1) * ret->size + i]);
	float sum = 0;
	for(j = i+1; j < ret->size; j++){
	    sum += ret->values[j * ret->size + i] *
		   ret->values[j * ret->size + i];
	}
	//printf("sum: %f\n", sum);
	alpha *= sqrt(sum);
	//printf("alpha: %f\n", alpha);
	float r = 0.5 * (alpha * alpha - ret->values[(i+1) * ret->size + i] * alpha);
	r = sqrt(r);
	r = 1 / (2*r);
	//printf("r: %f\n", r);
	vector * v = create_vector(ret->size);
	v->values[i+1] = ret->values[(i+1) * ret->size + i] - alpha;
	v->values[i+1] *= r;
	for(j = i+2; j < v->size; j++)
	    v->values[j] = ret->values[j * ret->size + i] * r;
	float * pos = p->values;
	for(j = 0; j < p->size; j++){
	    for(k = 0; k < p->size; k++){
		*pos = -2 * v->values[k] * v->values[j];
		if(j == k) *pos += 1;
		pos++;
	    }
	}
	//printf("Householder matrix:\n");
	//print_matrix(p);
	//printf("\n");
	matrix * temp1 = matrix_matrix_mult(p, ret);
	matrix * temp2 = matrix_matrix_mult(temp1, p);
	delete_matrix(ret);
	delete_matrix(temp1);
	ret = temp2;
	//printf("A_k:\n");
	//print_matrix(ret);
	//printf("\n");

    }
    return ret;
}


static void delete_diagnolization(diagnolization * d){
    delete_matrix(d->diagnol);
    delete_matrix(d->eigenvectors);
    free(d);
}

#define SMALL_MATRIX_BOUNDARY 4



static diagnolization * small_diagnolization(matrix * m){
    if(m->size > SMALL_MATRIX_BOUNDARY){
	printf("Error: Matrix of size %d inputed to small_diagnolization.  The maximum size should be %d", m->size, SMALL_MATRIX_BOUNDARY);
	return 0;
    }
    int i,j;
    matrix * A = create_matrix(m->size);
    int k = 0;
    while(k++ < 100){
	for(i = 0; i < m->size; i++)
	    for(j = 0; j < m->size; j++)
		A->values[i * A->size + j] = m->values[i * m->size + j];
	vector ** columns = (vector **) malloc(m->size * sizeof(vector *));
	for(i = 0; i < m->size; i++){
	    columns[i] = create_vector(m->size);
	    for(j = 0; j < columns[i]->size; j++)
		columns[i]->values[i] = m->values[j * m->size * i];
	}
	for(i = 0; i < m->size; i++) normalize(columns[i]);
	for(i = 1; i < m->size; i++){
	    vector * diff = project_vector(columns[i], columns, i, true);
	    for(j = 0; j < m->size; j++) columns[i]->values[j] -= diff->values[j];
	    delete_vector(diff);
	    normalize(columns[i]);
	}
	matrix * Q = create_matrix(m->size);
	for(i = 0; i < Q->size; i++)
	    for(j = 0; j < Q->size; j++)
		Q->values[i * Q->size + j] = columns[j]->values[i];
	matrix * Qtranspose = matrix_transpose(Q);
	matrix * R = matrix_matrix_mult(Qtranspose, A);
	matrix * temp = matrix_matrix_mult(R, Q);
	delete_matrix(A);
	A = temp;
    }
    
   
}

diagnolization * diagnolize(matrix * m){
    if(m->size <= SMALL_MATRIX_BOUNDARY) return small_diagnolization(m);
    diagnolization * ret = malloc(sizeof(diagnolization));
    ret->diagnol = create_matrix(m->size);
    ret->eigenvectors = create_matrix(m->size);
    //A tridiagnol matrix which is similar to m
    matrix * tridiagnol = householder(m);
    matrix * T1 = create_matrix(m->size / 2);
    matrix * T2 = create_matrix(m->size - T1->size);
    
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

    float asdf[9] = {1,   2 ,  -1,
		     2  , 3 ,  -1,
		     -2 ,  0 ,  -3 };
    matrix * m = create_matrix(3);
    int i;
    for(i = 0; i < m->size * m->size; i++) m->values[i] = asdf[i];
    printf("Matrix:\n");
    print_matrix(m);
    printf("\n");

    matrix * rref = row_reduced_echelon_form(m);
    printf("rref:\n");
    print_matrix(rref);

    
    return 0;
}
    
#endif

