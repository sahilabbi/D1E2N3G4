#include "isomorphism.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.0001

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

void delete_permutation(PartialPermutation * perm){
    int i;
    for(i = 0; i < perm->num_perm_sets; i++){
	free(perm->permSets[i].leftSet);
	free(perm->permSets[i].rightSet);
    }
    free(perm->permSets);
    free(perm);
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

//Splits up permSet into several permutation sets based on the contents of perm
static PartialPermutation * splitPermGroup(PermutationSet * permSet, PartialPermutation * perm){
    //The value stored in the tuple is the position of the index in perm
    indexValuePair * leftIndices  = (indexValuePair*) malloc(permSet->size * sizeof(PermutationSet));
    indexValuePair * rightIndices = (indexValuePair*) malloc(permSet->size * sizeof(PermutationSet));

    
    int i,j;

    int permSize = 0;
    for(i = 0; i < perm->num_perm_sets; i++) permSize += perm->permSets[i].size;
    //leftPermIndexMap[n] contains the index of the perm set in perm which contains n
    int * leftPermIndexMap  = (int *) malloc(permSize * sizeof(int));
    int * rightPermIndexMap = (int *) malloc(permSize * sizeof(int));
    for(i = 0; i < perm->num_perm_sets; i++){
	for(j = 0; j < perm->permSets[i].size; j++){
	    leftPermIndexMap[perm->permSets[i].leftSet[j]] = i;
	    rightPermIndexMap[perm->permSets[i].rightSet[j]] = i;
	}
    }
    
    for(i = 0; i < permSet->size; i++){
	leftIndices[i].index = permSet->leftSet[i];
	leftIndices[i].value = leftPermIndexMap[leftIndices[i].index];
	rightIndices[i].index = permSet->rightSet[i];
	rightIndices[i].value = rightPermIndexMap[rightIndices[i].index];
    }

    qsort(leftIndices, permSet->size, sizeof(indexValuePair), compareTuple);
    qsort(rightIndices, permSet->size, sizeof(indexValuePair), compareTuple);

    free(leftPermIndexMap);
    free(rightPermIndexMap);

    for(i = 0; i < permSet->size; i++)
	if(!float_equals(leftIndices[i].value, rightIndices[i].value))
	    return NULL;
    
    PartialPermutation * ret = (PartialPermutation*) malloc(sizeof(PartialPermutation));
    ret->num_perm_sets = permSet->size;
    ret->permSets = (PermutationSet *) malloc(permSet->size * sizeof(PermutationSet));
    ret->permSets[0].size = 0;
    ret->permSets[0].leftSet = (int *) malloc(permSet->size * sizeof(int));
    ret->permSets[0].rightSet = (int *) malloc(permSet->size * sizeof(int));

    int currPermSet = 0;
    //currPos is where we are in the current perm set
    for(i = 0; i < permSet->size; i++){
	if(i != 0 && !float_equals(leftIndices[i].value, leftIndices[i-1].value)){
	    ret->permSets[currPermSet].leftSet = realloc(ret->permSets[currPermSet].leftSet,
					       ret->permSets[currPermSet].size * sizeof(int));
	    ret->permSets[currPermSet].rightSet = realloc(ret->permSets[currPermSet].rightSet,
						ret->permSets[currPermSet].size * sizeof(int));
	    currPermSet++;
	    ret->permSets[currPermSet].size = 0;
	    ret->permSets[currPermSet].leftSet = (int *) malloc(permSet->size * sizeof(int));
	    ret->permSets[currPermSet].rightSet = (int *) malloc(permSet->size * sizeof(int));
	}
	ret->permSets[currPermSet].leftSet[ret->permSets[currPermSet].size]
	                                                           = leftIndices[i].index;
	ret->permSets[currPermSet].rightSet[ret->permSets[currPermSet].size]
	                                                           = rightIndices[i].index;
	ret->permSets[currPermSet].size++;
    }
    ret->permSets[currPermSet].leftSet = realloc(ret->permSets[currPermSet].leftSet,
				       ret->permSets[currPermSet].size * sizeof(int));
    ret->permSets[currPermSet].rightSet = realloc(ret->permSets[currPermSet].rightSet,
					ret->permSets[currPermSet].size * sizeof(int));
    ret->num_perm_sets = currPermSet + 1;
    ret->permSets = realloc(ret->permSets, ret->num_perm_sets * sizeof(PermutationSet));

    free(leftIndices);
    free(rightIndices);

    return ret;
}

PartialPermutation *combinePartialPermutations(PartialPermutation * acc, PartialPermutation * perm){
    if(acc == NULL || perm == NULL) {
	return NULL;
    }
    int i;
    
    int accSize, permSize;
    accSize = 0;
    permSize = 0;
    for(i = 0; i < acc->num_perm_sets; i++) accSize += acc->permSets[i].size;
    for(i = 0; i < perm->num_perm_sets; i++) permSize += perm->permSets[i].size;

    if(accSize != permSize){
	return NULL;
    }
    
    PartialPermutation * ret = (PartialPermutation *) malloc(sizeof(PartialPermutation));
    ret->num_perm_sets = 0;
    ret->permSets = (PermutationSet *) malloc(accSize * sizeof(PermutationSet));
    PartialPermutation * temp;
    for(i = 0; i < perm->num_perm_sets; i++){
	temp = splitPermGroup(&perm->permSets[i], acc);
	memcpy(&ret->permSets[ret->num_perm_sets],
	       temp->permSets,
	       temp->num_perm_sets * sizeof(PermutationSet));
	ret->num_perm_sets += temp->num_perm_sets;
	free(temp->permSets);
	free(temp);
    }
    ret->permSets = realloc(ret->permSets, ret->num_perm_sets * sizeof(PermutationSet));

    return ret;
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

    /*    for(i = 0; (size_t) i < left->size; i++){
	printf("%d, %lf\t%d, %lf\n", leftPairs[i].index, leftPairs[i].value,
	       rightPairs[i].index, rightPairs[i].value);
    }
    printf("\n");*/

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

void print_permutation(PartialPermutation * perm){
    if(perm == NULL){
	printf("NULL\n");
	return;
    }
    int i;
    for(i = 0; i < perm->num_perm_sets; i++){
	print_perm_set(&perm->permSets[i]);
	printf("\n");
    }
}
    

#ifdef MATRIX_TEST

int main(){
    double left[7] = {1.7, 2.2, 2.2, -1.4, -2.3, -1.4, 2.2};
    double right[7] = {2.2, 2.2, -2.3, -1.4, 1.7, -1.4, 2.2};
    double left2[7] = {-2.8, 1.31, 1.31, 1.31, 1.31, 0.92, -0.21};
    double right2[7] = {1.31, 1.31, 1.31, 0.92, -2.8, 1.31, -0.21};

    gsl_vector * vec1 = gsl_vector_alloc(7);
    gsl_vector * vec2 = gsl_vector_alloc(7);
    gsl_vector * vec3 = gsl_vector_alloc(7);
    gsl_vector * vec4 = gsl_vector_alloc(7);

    int i;
    for(i = 0; i < 7; i++){
	gsl_vector_set(vec1, i, left[i]);
	gsl_vector_set(vec2, i, right[i]);
	gsl_vector_set(vec3, i, left2[i]);
	gsl_vector_set(vec4, i, right2[i]);
    }


    PartialPermutation * thePerm = getPermutationFromVectors(vec1, vec2);
    PartialPermutation * thePerm2 = getPermutationFromVectors(vec3, vec4);

    print_permutation(thePerm);
    printf("\n");
    print_permutation(thePerm2);
    printf("\n");

    PartialPermutation * combinded = combinePartialPermutations(thePerm, thePerm2);
    print_permutation(combinded);
    
    return 0;
}
    
#endif


