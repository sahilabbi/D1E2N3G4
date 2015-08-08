#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.00001

bool float_equals(double a, double b){
    double diff = a - b;
    if(diff < 0) diff = -diff;
    return diff < EPSILON;
}

int float_compare(double a, double b){
    double diff = a - b;
    if(diff > EPSILON) return 1;
    if(diff < -EPSILON) return -1;
    return 0;
}

typedef struct
{
    int index;
    double value;
} indexValuePair;

static int compareTuple(const void * tup1, const void * tup2){
    indexValuePair * pair1 = (indexValuePair *) tup1;
    indexValuePair * pair2 = (indexValuePair *) tup2;

    double val1 = pair1->value;
    double val2 = pair2->value;

    return float_compare(val1, val2);
}

PartialPermutation * getPermutationFromVectors(gsl_vector * left, gsl_vector * right){
    if(left->size != right->size){
	printf("Error: incorrectly sized vectors in getPermutationFromVectors\n");
	return NULL;
    }
    indexValuePair * leftPairs = (indexValuePair *) malloc(left->size * sizeof(indexValuePair));
    indexValuePair * rightPairs = (indexValuePair *) malloc(right->size * sizeof(indexValuePair));
    int i;
    for(i = 0; (size_t) i < left->size; i++){
	leftPairs[i].index = i;
	leftPairs[i].value = gsl_vector_get(left, i);

	rightPairs[i].index = i;
	rightPairs[i].value = gsl_vector_get(right, i);
    }

    qsort(leftPairs, left->size, sizeof(indexValuePair), compareTuple);
    qsort(rightPairs, right->size, sizeof(indexValuePair), compareTuple);

    indexValuePair * currLeftPos = leftPairs;
    indexValuePair * currRightPos = rightPairs;

    PartialPermutation * ret = (PartialPermutation *) malloc(sizeof(PartialPermutation));
    ret->permSets = (PermutationSet *) malloc(left->size * sizeof(PermutationSet));
    ret->num_perm_sets = 0;
    i = 0;
    while((size_t) i < left->size){
	//If the values don't agree, return false
	if(!float_equals(currLeftPos->value, currRightPos->value)) return NULL;
	double currVal = currLeftPos->value;
	PermutationSet * currPermSet = &ret->permSets[ret->num_perm_sets];
	currPermSet->size = 0;
	currPermSet->leftSet = (int *) malloc(left->size * sizeof(int));
	currPermSet->rightSet = (int *) malloc(left->size * sizeof(int));
	while((size_t) i < left->size &&
	      float_equals(currLeftPos->value, currVal)){
	    if(!float_equals(currLeftPos->value, currRightPos->value)) return NULL;
	    currPermSet->leftSet[currPermSet->size] = currLeftPos->index;
	    currPermSet->rightSet[currPermSet->size] = currRightPos->index;
	    currPermSet->size++;
	    currLeftPos++;
	    currRightPos++;
	    i++;
	}
	realloc(currPermSet->leftSet, currPermSet->size * sizeof(int));
	ret->num_perm_sets++;
    }
    realloc(ret, ret->num_perm_sets);

    return ret;
}

static void print_perm_set(PermutationSet * permSet){
    printf("[");
    int i;
    for(i = 0; i < permSet->size; i++){
	printf("%d", permSet->leftSet[i]);
	if(i + 1 < permSet->size) printf(",");
    }
    printf("] -> [");
    for(i = 0; i < permSet->size; i++){
	printf("%d", permSet->rightSet[i]);
	if(i + 1 < permSet->size) printf(",");
    }
    printf("]");
}

static void print_permutation(PartialPermutation * perm){
    int i;
    for(i = 0; i < perm->num_perm_sets; i++){
	print_perm_set(&perm->permSets[i]);
	printf("\n");
    }
}
    

#ifdef MATRIX_TEST

int main(){
    double left[7] = {1.7, 2.2, 0.45, 1.2, -2.3, -1.4, 2.2};
    double right[7] = {2.2, 0.45, -2.3, 2.2, 1.7, -1.4, 1.2};

    gsl_vector * vec1 = gsl_vector_alloc(7);
    gsl_vector * vec2 = gsl_vector_alloc(7);

    int i;
    for(i = 0; i < 7; i++){
	gsl_vector_set(vec1, i, left[i]);
	gsl_vector_set(vec2, i, right[i]);
    }

    print_permutation(getPermutationFromVectors(vec1, vec2));
    
    return 0;
}
    
#endif

