#include "isomorphism.h"
#include "graph_set.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>

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

gsl_matrix * graph_to_matrix(graph_t * g){
    gsl_matrix * ret = gsl_matrix_alloc(NUM_NODES, NUM_NODES);
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	for(j = 0; j < NUM_NODES; j++){
	    gsl_matrix_set(ret, i, j, g->adj[i][j]);
	}
    }
    return ret;
}


void print_vector(gsl_vector * v){
    int i;
    for(i = 0; (size_t) i < v->size; i++){
	printf("%lf ", gsl_vector_get(v, i));
    }
    fflush(stdout);
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

static void vector_normalize(gsl_vector * v){
    double length = 0;
    gsl_blas_ddot(v, v, &length);
    length = 1 / sqrt(length);
    gsl_vector_scale(v, length);
}

gsl_vector * project_vector(gsl_vector * v, gsl_vector ** vspace, int num_vectors){
    int i;
    gsl_vector * ret = gsl_vector_calloc(v->size);
    gsl_vector * temp = gsl_vector_alloc(v->size);
    for(i = 0; i < num_vectors; i++){
	gsl_vector_memcpy(temp, vspace[i]);
	double scale = 0;
	gsl_blas_ddot(vspace[i], v, &scale);
	gsl_vector_scale(temp, scale);
	gsl_vector_add(ret, temp);
    }
    gsl_vector_free(temp);
    return ret;
}

//takes the absolute value of each value in v
static void vector_abs(gsl_vector * v){
    int i;
    for(i = 0; (size_t) i < v->size; i++){
	double val = gsl_vector_get(v, i);
	if(val < 0) val = -val;
	gsl_vector_set(v, i, val);
    }
}

//For M_PI
#include <math.h>

gsl_vector * reduce_vector_space(gsl_vector ** vectors, int num_vectors){
    if(num_vectors == 1){
	gsl_vector * ret = gsl_vector_alloc((*vectors)->size);
	gsl_vector_memcpy(ret, vectors[0]);
	vector_abs(ret);
	return ret;
    }
    gsl_vector * standard_basis_vector = gsl_vector_calloc((*vectors)->size);
    gsl_vector * ret = gsl_vector_alloc((*vectors)->size);
    int i;
    for(i = 0; (size_t) i < (*vectors)->size; i++){
	if(i > 0) gsl_vector_set(standard_basis_vector, i-1, 0);
	gsl_vector_set(standard_basis_vector, i, 1);
	gsl_vector * proj = project_vector(standard_basis_vector, vectors, num_vectors);
	double projLength = 0;
	gsl_blas_ddot(proj, proj, &projLength);
	double angle = 0;
	if(float_equals(projLength, 0)){
	    //If the basis vector is orthogonal to the space
	    angle = M_PI / 2;
	}
	else {
	vector_normalize(proj);
	gsl_blas_ddot(standard_basis_vector, proj, &angle);
	}
	gsl_vector_set(ret, i, angle);

	gsl_vector_free(proj);
    }
    gsl_vector_free(standard_basis_vector);
    return ret;
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

    for(i = 0; i < permSet->size; i++){
	if(!float_equals(leftIndices[i].value, rightIndices[i].value)){
	    free(leftIndices);
	    free(rightIndices);
	    return NULL;
	}
    }
    
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
	if(temp == NULL){
	    delete_permutation(ret);
	    return NULL;
	}
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
	if(!float_equals(currLeftPos->value, currRightPos->value)){
	    delete_permutation(ret);
	    free(leftPairs);
	    free(rightPairs);
	    return NULL;
	}
	double currVal = currLeftPos->value;
	PermutationSet * currPermSet = &ret->permSets[ret->num_perm_sets];
	currPermSet->size = 0;
	currPermSet->leftSet = (int *) malloc(left->size * sizeof(int));
	currPermSet->rightSet = (int *) malloc(left->size * sizeof(int));
	while((size_t) i < left->size &&
	      float_equals(currLeftPos->value, currVal)){
	    if(!float_equals(currLeftPos->value, currRightPos->value)){
		free(currPermSet->leftSet);
		free(currPermSet->rightSet);
		delete_permutation(ret);
		free(leftPairs);
		free(rightPairs);
		return NULL;
	    }
	    currPermSet->leftSet[currPermSet->size] = currLeftPos->index;
	    currPermSet->rightSet[currPermSet->size] = currRightPos->index;
	    currPermSet->size++;
	    currLeftPos++;
	    currRightPos++;
	    i++;
	}
	currPermSet->leftSet = realloc(currPermSet->leftSet, currPermSet->size * sizeof(int));
	currPermSet->rightSet = realloc(currPermSet->rightSet, currPermSet->size * sizeof(int));
	ret->num_perm_sets++;
    }
    ret->permSets = realloc(ret->permSets, ret->num_perm_sets * sizeof(PermutationSet));

    free(leftPairs);
    free(rightPairs);

    return ret;
}


static void sort_diagonalization(gsl_vector * eigval, gsl_matrix * eigvec){
    gsl_permutation * p = gsl_permutation_alloc(eigval->size);
    gsl_sort_vector_index(p, eigval);
    int i;
    for(i = 0; (size_t) i < eigval->size; i++){
	gsl_vector_view theRow = gsl_matrix_row(eigvec, i);
	gsl_permute_vector(p, &(theRow.vector));
    }
    gsl_permute_vector(p, eigval);
    gsl_permutation_free(p);
}

void print_matrix(gsl_matrix * mat){
    int i,j;
    for(i = 0; (size_t) i < mat->size1; i++){
	for(j = 0; (size_t) j < mat->size2; j++){
	    printf("%lf ", gsl_matrix_get(mat, i, j));
	}
	printf("\n");
    }
    fflush(stdout);
}

#include <gsl/gsl_eigen.h>

PartialPermutation * isomorphism(graph_t * g1, graph_t * g2){
    gsl_matrix * mat1 = graph_to_matrix(g1);
    gsl_matrix * mat2 = graph_to_matrix(g2);

    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc(NUM_NODES);
    
    gsl_vector * eigval1 = gsl_vector_alloc(NUM_NODES);
    gsl_vector * eigval2 = gsl_vector_alloc(NUM_NODES);

    gsl_matrix * eigvec1 = gsl_matrix_alloc(NUM_NODES,NUM_NODES);
    gsl_matrix * eigvec2 = gsl_matrix_alloc(NUM_NODES,NUM_NODES);

    gsl_eigen_symmv(mat1, eigval1, eigvec1, workspace);
    gsl_eigen_symmv(mat2, eigval2, eigvec2, workspace);

    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(mat1);
    gsl_matrix_free(mat2);

    sort_diagonalization(eigval1, eigvec1);
    sort_diagonalization(eigval2, eigvec2);

    /*    printf("eigval1:\n");
    print_vector(eigval1);
    printf("\neigvec1:\n");
    print_matrix(eigvec1);
    
    printf("eigval2:\n");
    print_vector(eigval2);
    printf("\neigvec2:\n");
    print_matrix(eigvec2); */

    int i;
    for(i = 0; i < NUM_NODES; i++){
	if(!float_equals(gsl_vector_get(eigval1, i),
			 gsl_vector_get(eigval2, i))) return NULL;
    }
    //Since the two vectors are equal, we no longer need the second one
    gsl_vector_free(eigval2);

    gsl_vector ** columns1 = (gsl_vector **) malloc(NUM_NODES * sizeof(gsl_vector *));
    gsl_vector ** columns2 = (gsl_vector **) malloc(NUM_NODES * sizeof(gsl_vector *));
    for(i = 0; i < NUM_NODES; i++){
	gsl_vector theCol = gsl_matrix_column(eigvec1, i).vector;
	columns1[i] = gsl_vector_alloc(NUM_NODES);
	gsl_vector_memcpy(columns1[i], &theCol);
	theCol = gsl_matrix_column(eigvec2, i).vector;
	columns2[i] = gsl_vector_alloc(NUM_NODES);
	gsl_vector_memcpy(columns2[i], &theCol);
    }

    //Deleting the diagonalization matrices
    gsl_matrix_free(eigvec1);
    gsl_matrix_free(eigvec2);
    
    /*    printf("columns1:\n");
    for(i = 0; i < NUM_NODES; i++){
	print_vector(columns1[i]);
	printf("\n");
    }
    printf("columns2:\n");
    for(i = 0; i < NUM_NODES; i++){
	print_vector(columns2[i]);
	printf("\n");
	} */

    gsl_vector ** reducedColumns1 = (gsl_vector **) malloc(NUM_NODES * sizeof(gsl_vector *));
    gsl_vector ** reducedColumns2 = (gsl_vector **) malloc(NUM_NODES * sizeof(gsl_vector *));
    int numSpaces = 0;
    i = 0;
    while(i < NUM_NODES){
	int currSpaceSize = 0;
	do {
	    currSpaceSize++;
	    i++;
	}
	while(i < NUM_NODES && float_equals(gsl_vector_get(eigval1, i),
					    gsl_vector_get(eigval1, i-1)));
	reducedColumns1[numSpaces] = reduce_vector_space(columns1 + i - currSpaceSize,
							 currSpaceSize);
	reducedColumns2[numSpaces] = reduce_vector_space(columns2 + i - currSpaceSize,
							 currSpaceSize);
	numSpaces++;
    }
    reducedColumns1 = (gsl_vector **) realloc(reducedColumns1, numSpaces * sizeof(gsl_vector *));
    reducedColumns2 = (gsl_vector **) realloc(reducedColumns2, numSpaces * sizeof(gsl_vector *));

    for(i = 0; i < NUM_NODES; i++){
	gsl_vector_free(columns1[i]);
	gsl_vector_free(columns2[i]);
    }
    free(columns1);
    free(columns2);

    gsl_vector_free(eigval1);

    PartialPermutation ** perms
	= (PartialPermutation **) malloc(numSpaces * sizeof(PartialPermutation *));
    for(i = 0; i < numSpaces; i++)
	perms[i] = getPermutationFromVectors(reducedColumns1[i], reducedColumns2[i]);

    for(i = 0; i < numSpaces; i++) {
	gsl_vector_free(reducedColumns1[i]);
	gsl_vector_free(reducedColumns2[i]);
    }
    free(reducedColumns1);
    free(reducedColumns2);

    PartialPermutation * ret = (PartialPermutation *) malloc(sizeof(PartialPermutation));
    if(perms[0] != NULL){
	ret->permSets = (PermutationSet *) malloc(perms[0]->num_perm_sets * sizeof(PermutationSet));
	ret->num_perm_sets = perms[0]->num_perm_sets;
	for(i = 0; i < ret->num_perm_sets; i++){
	    ret->permSets[i].size = perms[0]->permSets[i].size;
	    ret->permSets[i].leftSet  = (int *) malloc(ret->permSets[i].size * sizeof(int));
	    ret->permSets[i].rightSet = (int *) malloc(ret->permSets[i].size * sizeof(int));
	    memcpy(ret->permSets[i].leftSet,
		   perms[0]->permSets[i].leftSet,
		   ret->permSets[i].size * sizeof(int));
	    memcpy(ret->permSets[i].rightSet,
		   perms[0]->permSets[i].rightSet,
		   ret->permSets[i].size * sizeof(int));
	}
    }
    else {
	free(ret);
	ret = NULL;
    }
    if(ret != NULL) {
	for(i = 1; i < numSpaces; i++){
	    PartialPermutation * temp = combinePartialPermutations(ret, perms[i]);
	    delete_permutation(ret);
	    if(temp == NULL){
		ret = NULL;
		break;
	    }
	    ret = temp;
	}
    }

    for(i = 0; i < numSpaces; i++){
	if(perms[i] != NULL) delete_permutation(perms[i]);
    }
    free(perms);
    
    return ret;
}

bool is_isomorphic(graph_t * g1, graph_t * g2){
    gsl_permutation * iso = getIsomorphism(g1, g2);
    if(iso == NULL) return false;
    gsl_permutation_free(iso);
    return true;
}

static bool isComplete(PartialPermutation * perm){
    return perm->num_perm_sets == NUM_NODES;
}

gsl_permutation * getIsomorphism(graph_t * g1, graph_t * g2){
    PartialPermutation * iso = isomorphism(g1,g2);
    if(iso == NULL) return NULL;
    int i,j;
    
    //used to make sure the two corresponding vertecies are mapped to each other
    //the first perm group is v -> whichever vertex in g2 corresonds
    //and the last one is all the remaining elements
    PartialPermutation * extraRestriction = malloc(sizeof(PartialPermutation));
    extraRestriction->num_perm_sets = 2;
    extraRestriction->permSets = malloc(3 * sizeof(PermutationSet));
    extraRestriction->permSets[0].leftSet = malloc(sizeof(int));
    extraRestriction->permSets[0].rightSet = malloc(sizeof(int));
    extraRestriction->permSets[0].size = 1;
    extraRestriction->permSets[1].leftSet = malloc((NUM_NODES - 1) * sizeof(int));
    extraRestriction->permSets[1].rightSet = malloc((NUM_NODES - 1) * sizeof(int));
    extraRestriction->permSets[1].size = NUM_NODES - 1;

    /*printf("Started with inital permutation:\n");
    print_permutation(iso);
    printf("\n");*/

    bool done = false;
    while(!isComplete(iso)){
	//Select a vertex whose other isomorphism are not determined
	int v = 0;
	PermutationSet * currPermSet = 0;
	for(i = 0; i < iso->num_perm_sets; i++){
	    if(iso->permSets[i].size > 1){
		currPermSet = iso->permSets + i;
		v = currPermSet->leftSet[0];
		break;
	    }
	}
	
	//Remove or add two links to v
	//theLink1 contains a link from withing the same perm set as v
	int theLink1 = currPermSet->leftSet[1];
	int link1Value = g1->adj[v][theLink1];
	int oppositeLink1Value = link1Value == 0 ? 1 : 0;
	
	//theLink2Set is the index of the permutation set from which
	//the second link (theLink2) comes from
	int theLink2Set = 0;
	if(i == 0) theLink2Set++;
	//Incase there is only one perm set, so theLink2 must come from
	//the same set as v
	int theLink2;
	if(iso->num_perm_sets == 1){
	    theLink2Set = 0;
	    theLink2 = iso->permSets[theLink2Set].leftSet[2];
	}
	else {
	    theLink2 = iso->permSets[theLink2Set].leftSet[0];
	}
	int link2Value = g1->adj[v][theLink2];
	int oppositeLink2Value = link2Value == 0 ? 1 : 0;



	extraRestriction->permSets[0].leftSet[0] = v;
	int k = 0;
	for(i = 0; i < NUM_NODES; i++){
	    if(i == v) continue;
	    extraRestriction->permSets[1].leftSet[k++] = i;
	}

	//Go through every vertex in g2 and try to get an isomorphism
	done = false;
	for(i = 0; i < currPermSet->size; i++){
	    int u = currPermSet->rightSet[i];
	    for(j = 0; j < currPermSet->size; j++){
		if(j == i) continue;
		int testLink1 = currPermSet->rightSet[j];
		//We're testing the link between u and testLink1
		//to match with the link between v and theLink
		if(g2->adj[u][testLink1] != link1Value) continue;
		//Remove both the links we're comparing
		g1->adj[v][theLink1] = oppositeLink1Value;
		g1->adj[theLink1][v] = oppositeLink1Value;
		g2->adj[u][testLink1] = oppositeLink1Value;
		g2->adj[testLink1][u] = oppositeLink1Value;
		//Get the new isomorphism
		PartialPermutation * newIso1 = isomorphism(g1, g2);
		//Reset the links back to what they were
		g1->adj[v][theLink1] = link1Value;
		g1->adj[theLink1][v] = link1Value;
		g2->adj[u][testLink1] = link1Value;
		g2->adj[testLink1][u] = link1Value;
		//If these links didn't work out, try again
		if(newIso1 == NULL) continue;
		for(k = 0; k < iso->permSets[theLink2Set].size; k++){
		    int testLink2 = iso->permSets[theLink2Set].rightSet[k];
		    if(g2->adj[u][testLink2] != link2Value) continue;
		    //Remove the links we're comparing
		    g1->adj[v][theLink2] = oppositeLink2Value;
		    g1->adj[theLink2][v] = oppositeLink2Value;
		    g2->adj[u][testLink2] = oppositeLink2Value;
		    g2->adj[testLink2][u] = oppositeLink2Value;
		    //Get the new ismoorphism
		    PartialPermutation * newIso2 = isomorphism(g1, g2);
		    //Reset the links back to what they were
		    g1->adj[v][theLink2] = link2Value;
		    g1->adj[theLink2][v] = link2Value;
		    g2->adj[u][testLink2] = link2Value;
		    g2->adj[testLink2][u] = link2Value;

		    if(newIso2 == NULL) continue;

		    //We now know that an isomorphism carries v to u

		    //Add info from the two ismorphisms (newIso1 and newIso2)
		    PartialPermutation * temp1 = combinePartialPermutations(iso, newIso1);
		    PartialPermutation * temp2 = combinePartialPermutations(temp1, newIso2);
		    //Prepare extraRestriction, which ensures [v] -> [u] is in the final permutation
		    extraRestriction->permSets[0].rightSet[0] = u;
		    int l = 0;
		    for(k = 0; k < NUM_NODES; k++){
			if(k == u) continue;
			extraRestriction->permSets[1].rightSet[l++] = k;
		    }
		    PartialPermutation * temp3 = combinePartialPermutations(temp2,
									    extraRestriction);
		    //Sometimes the permutations don't combine properly so we
		    //cannot assume that temp1 and temp2 aren't NULL
		    if(temp1 != NULL) delete_permutation(temp1);
		    if(temp2 != NULL) delete_permutation(temp2);
		    delete_permutation(newIso1);
		    delete_permutation(newIso2);
		    if(temp3 != NULL){
			delete_permutation(iso);
			/*printf("Found that an isomorphism carries %d to %d\n", v, u);
			printf("New partial permutation:\n");
			print_permutation(temp3);
			printf("\n");*/
			iso = temp3;
			done = true;
			break;
		    }
		}
		if(done) break;
	    }
	    if(done) break;	    
	}


	
	//If none of the vertecies in g2 correspond, there is no ismorphism
	//log the fact that the inital test was incorrect
	if(!done){
	    delete_permutation(iso);
	    printf("Could not find isomorphism between two graphs\n");
	    printf("g1:\n");
	    print_graph(g1);
	    printf("\ng2:\n");
	    print_graph(g2);
	    printf("\n");
	    return NULL;
	}
    }
    
    delete_permutation(extraRestriction);
    //Now the partial permutation is a complete permutation
    gsl_permutation * leftPerm = gsl_permutation_alloc(NUM_NODES);
    gsl_permutation * rightPerm = gsl_permutation_alloc(NUM_NODES);
    for(i = 0; i < NUM_NODES; i++){
	leftPerm->data[i] = iso->permSets[i].leftSet[0];
	rightPerm->data[i] = iso->permSets[i].rightSet[0];
    }
    gsl_permutation * inverse = gsl_permutation_alloc(NUM_NODES);
    gsl_permutation_inverse(inverse, leftPerm);
    gsl_permutation * ret = gsl_permutation_alloc(NUM_NODES);

    gsl_permutation_mul(ret, inverse, rightPerm);

    gsl_permutation_free(leftPerm);
    gsl_permutation_free(rightPerm);
    gsl_permutation_free(inverse);

    return ret;
}

static void print_perm_group(PermutationSet * permSet){
    int i;
    printf("[");
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
	print_perm_group(perm->permSets + i);
	printf("\n");
    }
}

#ifdef MATRIX_TEST

int main(){

    int adj1[8][8] = {{0,1,1,0,1,0,0,0},
		      {1,0,0,1,1,1,1,0},
		      {1,0,0,1,1,1,0,1},
		      {0,1,1,0,0,1,0,0},
		      {1,1,1,0,0,0,1,0},
		      {0,1,1,1,0,0,0,1},
		      {0,1,0,0,1,0,0,1},
		      {0,0,1,0,0,1,1,0}};
    int adj2[8][8] = {{0,1,0,1,0,0,0,1},
		      {1,0,1,0,0,0,1,1},
		      {0,1,0,1,0,1,0,0},
		      {1,0,1,0,1,1,0,0},
		      {0,0,0,1,0,1,0,1},
		      {0,0,1,1,1,0,1,1},
		      {0,1,0,0,0,1,0,1},
		      {1,1,0,0,1,1,1,0}};
    int iso2[6][6] ={{0, 0, 0, 0, 1, 0},
		     {0, 0, 0, 1, 1, 1},
		     {0, 0, 0, 0, 0, 1},
		     {0, 1, 0, 0, 0, 0},
		     {1, 1, 0, 0, 0, 1},
		     {0, 1, 1, 0, 1, 0}};
	
    int iso1[6][6] = {{0, 0, 0, 0, 0, 1},
		      {0, 0, 0, 0, 1, 0},
		      {0, 0, 0, 1, 1, 1},
		      {0, 0, 1, 0, 0, 0},
		      {0, 1, 1, 0, 0, 1},
		      {1, 0, 1, 0, 1, 0}};

    int test1[8][8] = {{0,0,0,0,0,1,0,1},
		       {0,0,1,1,0,0,0,0},
		       {0,1,0,0,0,0,0,0},
		       {0,1,0,0,0,1,1,0},
		       {0,0,0,0,0,0,1,1},
		       {1,0,0,1,0,0,0,0},
		       {0,0,0,1,1,0,0,0},
		       {1,0,0,0,1,0,0,0}};
	    
    int test2[8][8] = {{0,0,0,0,0,0,0,1},
		       {0,0,0,0,0,1,1,0},
		       {0,0,0,0,1,0,1,0},
		       {0,0,0,0,1,1,0,0},
		       {0,0,1,1,0,0,0,0},
		       {0,1,0,1,0,0,0,0},
		       {0,1,1,0,0,0,0,1},
		       {1,0,0,0,0,0,1,0}};

    graph_t * g1 = malloc(sizeof(graph_t));
    graph_t * g2 = malloc(sizeof(graph_t));

    memcpy(&g1->adj[0][0], test1, sizeof(test1));
    memcpy(&g2->adj[0][0], test2, sizeof(test2));

    gsl_permutation * theIso = getIsomorphism(g1, g2);
    int i;
    if(theIso == NULL) printf("NULL\n");
    else {
	for(i = 0; i < NUM_NODES; i++){
	    printf("%d ", theIso->data[i]);
	}
	printf("\n");
	gsl_permutation_free(theIso);
    }

    PartialPermutation * perm = isomorphism(g1, g2);
    printf("Partial perm:\n");
    print_permutation(perm);

    printf("\nIso test: %d\n", is_isomorphic(g1, g2));
	
    free(g1);
    free(g2);
    
    return 0;
}
    
#endif

#ifdef ISO_TEST

int main(){
    graph_t * g = (graph_t *) malloc(sizeof(graph_t));
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	for(j = 0; j < NUM_NODES; j++){
	    g->adj[i][j] = 0;
	}
    }
    graph_set * gs = graph_set_alloc();
    bool done = false;
    while(!done){
	if(!check_isomorphism(gs, g)){
	    insert_graph(gs, g);
	}
	done = increment_graph(g);
    }

    int numGraphs = 0;
    for(i = 0; i < gs->size; i++){
	numGraphs += gs->iso_groups[i].num_graphs;
    }
    printf("Total graphs with N=%d: %d\n", NUM_NODES, numGraphs);
    //    printf("Final Graph: \n");

    //print_graph(g);

    //    printf("\nThe Graphs:\n");
    //    print_graph_set(gs);
}

#endif

