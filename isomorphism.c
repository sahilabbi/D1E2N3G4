#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.00001

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
    ret->num_vectors = num_vectors;
    int i;
    for(i = 0; i < num_vectors; i++) ret->vectors[i] = create_vector(size);
    ret->size = size;
    return ret;
}

void delete_vector_space(vector_space * vspace){
    int i;
    for(i = 0; i < vspace->num_vectors; i++) delete_vector(vspace->vectors[i]);
    free(vspace->vectors);
    free(vspace);
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
    //printf("Vector:\n");
    //print_vector(v);
    for(i = 0; i < v->size; i++) v->values[i] *= correction;
    //printf("Normalized\n");
    //print_vector(v);
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
    for(i = 0; i < ret->size * ret->size; i++){
	if(float_equals(ret->values[i], 0))
	    ret->values[i] = 0;
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
	//printf("corrr: %f\n", r);
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
	printf("Householder matrix:\n");
	print_matrix(p);
	printf("\n");
	matrix * temp1 = matrix_matrix_mult(p, ret);
	matrix * temp2 = matrix_matrix_mult(temp1, p);
	delete_matrix(ret);
	delete_matrix(temp1);
	ret = temp2;
	printf("A_k:\n");
	print_matrix(ret);
	printf("\n");

    }
    return ret;
}


static void delete_diagnolization(diagnolization * d){
    free(d->diagonal);
    delete_matrix(d->eigenvectors);
    free(d);
}

#define SMALL_MATRIX_BOUNDARY 4

static void orthagonalize(vector_space * vspace){
    int i,j;
    int temp = vspace->num_vectors;
    normalize(vspace->vectors[0]);
    for(i = 1; i < temp; i++){
	vspace->num_vectors = i;
	vector * diff = project_vector(vspace->vectors[i], vspace, true);
	for(j = 0; j < vspace->size; j++){
	    vspace->vectors[i]->values[j] -= diff->values[j];
	}
	delete_vector(diff);
	normalize(vspace->vectors[i]);
    }
    vspace->num_vectors = temp;
}

static vector_space * get_nullspace(matrix * m){
    printf("Matrix:\n");
    print_matrix(m);
    matrix * rref = row_reduced_echelon_form(m);
    printf("RREF:\n");
    print_matrix(rref);
    //represents index of the first pivot in each row
    //puts m->size if there is no pivot
    int * pivots = (int *) malloc(m->size * sizeof(int));
    //same as pivots, but reversing index/values
    int i,j;
    for(i = 0; i < m->size; i++){
	pivots[i] = m->size;
	for(j = i; j < m->size; j++){
	    if(float_equals(rref->values[i * rref->size + j], 1)){
		pivots[i] = j;
		break;
	    }
	}
    }
    
    int dimension = m->size;
    for(i = 0; i < m->size; i++){
	if(pivots[i] == m->size){
	    dimension = m->size - i;
	    break;
	}
    }
    printf("Nullspace dimension: %d\n", dimension);
    vector_space * ret = create_vector_space(dimension, m->size);
    //contains the columns which do not have pivots in them
    int * pivot_holes = (int *) malloc(dimension * sizeof(int));
    //contains the row in which the pivot for the given column appears
    int * pivot_reverse = (int *) malloc(m->size * sizeof(int));
    for(i = 0; i < m->size; i++) pivot_reverse[i] = m->size;
    int last = 0;
    int k = 0;
    for(i = 0; i < m->size; i++){
	if(pivots[i] > last + 1){
	    for(j = last + 1; j < pivots[i]; j++) pivot_holes[k++] = j;
	}
	pivot_reverse[pivots[i]] = i;
    }
    for(i = 0; i < ret->num_vectors; i++){
	ret->vectors[i] = create_vector(ret->size);
	ret->vectors[i]->values[pivot_holes[i]] = 1;
	//iterate across the values in columnn i and balance
	//them off with values in pivot columns to make the
	//matrix-vector product 0
	for(j = 0; j < rref->size; j++){
	    ret->vectors[i]->values[pivot_reverse[j]] = -rref->values[j * rref->size + i];
	}
    }
    free(pivots);
    free(pivot_holes);
    free(pivot_reverse);
    delete_matrix(rref);
    return ret;
}

static diagnolization * small_diagnolization(matrix * m){
    if(m->size > SMALL_MATRIX_BOUNDARY){
	printf("Error: Matrix of size %d inputed to small_diagnolization.  The maximum size should be %d", m->size, SMALL_MATRIX_BOUNDARY);
    }
    int i,j;
    matrix * A = householder(m);
    bool done = false;
    while(!done){
	vector_space * columns = create_vector_space(m->size, m->size);
	for(i = 0; i < m->size; i++){
	    for(j = 0; j < columns->size; j++)
		columns->vectors[i]->values[j] = A->values[j * A->size + i];
	}
	orthagonalize(columns);
	//printf("columns->num_vectors: %d\n", columns->num_vectors);
	
       	matrix * Q = create_matrix(m->size);
	for(i = 0; i < Q->size; i++)
	    for(j = 0; j < Q->size; j++)
		Q->values[i * Q->size + j] = columns->vectors[j]->values[i];
	matrix * Qtranspose = matrix_transpose(Q);
	matrix * R = matrix_matrix_mult(Qtranspose, A);
	matrix * temp = matrix_matrix_mult(R, Q);
	/*printf("A:\n")
	print_matrix(A);
	printf("\nQ:\n");
	print_matrix(Q);
	printf("\nR:\n");
	print_matrix(R);
	printf("\n\n");*/
	delete_matrix(A);
	delete_matrix(Qtranspose);
	delete_matrix(R);
	delete_matrix(Q);
	A = temp;
	done = true;
	for(i = 1; i < A->size; i++){
	    for(j = 0; j < i; j++){
		if(!float_equals(A->values[i * A->size + j], 0)) done = false;
	    }
	}
    }
    diagnolization * ret = (diagnolization *) malloc(sizeof(diagnolization));
    ret->diagonal = (float *) malloc(m->size * sizeof(float));
    for(i = 0; i < m->size; i++){
	ret->diagonal[i] = A->values[i * A->size + i];
    }
    //Getting the Eigenvectors
    ret->eigenvectors = create_matrix(m->size);
    matrix * tempMat = create_matrix(m->size);
    int k;
    for(i = 0; i < m->size; i++){
	for(j = 0; j < m->size * m->size; j++) tempMat->values[j] = m->values[j];
	//The null space of tempMat (= m - lambda * I) is the eigenspace
	//corresponding to lambda
	for(j = 0; j < tempMat->size; j++){
	    tempMat->values[j * tempMat->size + j] -= ret->diagonal[i];
	}
	vector_space * eigenspace = get_nullspace(tempMat);
	int multiplicity = 0;
	for(j = i; j < m->size; j++){
	    if(ret->diagonal[j] == ret->diagonal[i]) multiplicity++;
	    else break;
	}
	for(j = i; j - i < multiplicity && j - i < eigenspace->num_vectors; j++){
	    for(k = 0; k < m->size; k++){
		ret->eigenvectors->values[k * ret->eigenvectors->size + j]
		    = eigenspace->vectors[j]->values[k];
	    }
	}
	delete_vector_space(eigenspace);
    }
    delete_matrix(tempMat);
    return ret;
}

diagnolization * diagnolize(matrix * m){
    if(m->size <= SMALL_MATRIX_BOUNDARY) return small_diagnolization(m);
    diagnolization * ret = malloc(sizeof(diagnolization));
    ret->diagonal = (float *) malloc(m->size * sizeof(float));
    ret->eigenvectors = create_matrix(m->size);
    //A tridiagnol matrix which is similar to m
    matrix * tridiagonal = householder(m);
    matrix * T1 = create_matrix(m->size / 2);
    matrix * T2 = create_matrix(m->size - T1->size);
    
    delete_matrix(tridiagonal);
    delete_matrix(T1);
    delete_matrix(T2);
    return NULL;
    
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

    float asdf[16] = {4,1,-2,2,
		     1,2,0,1,
		      -2,0,3,-1,
		     2,1,-1,-1};
    matrix * m = create_matrix(4);
    int i,j;
    for(i = 0; i < m->size * m->size; i++) m->values[i] = asdf[i];
    //    printf("Matrix:\n");
    //    print_matrix(m);
    //    printf("\n");

    vector_space * vspace = create_vector_space(2,3);
    for(i = 0; i < 2; i++)
	for(j = 0; j < 3; j++)
	    vspace->vectors[i]->values[j] = 1;
    vspace->vectors[0]->values[0] = 0;

    vector * angles = angles_to_standard_basis(vspace, false);

    printf("Angles:\n");
    print_vector(angles);

    //    diagnolization * h = householder(m);

    //    printf("Householder:\n");
    //    print_matrix(h);

    
    return 0;
}
    
#endif

