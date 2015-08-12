#include "graph_set.h"
#include "isomorphism.h"
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

graph_set * graph_set_alloc(){
    graph_set * ret = malloc(sizeof(graph_set));
    ret->eigenvalues = NULL;
    ret->graphs = NULL;
    ret->graph_set_sizes = NULL;
    ret->size = 0;
    return ret;
}

void delete_graph_set(graph_set * gs){
    int i,j;
    for(i = 0; (size_t) i < gs->size; i++){
	gsl_vector_free(gs->eigenvalues[i]);
	for(j = 0; (size_t) j < gs->graph_set_sizes[i]; j++){
	    free(gs->graphs[i][j]);
	}
	free(gs->graphs[i]);
    }
    free(gs->graph_set_sizes);
    free(gs->eigenvalues);
    free(gs->graphs);
    free(gs);
}

static int vector_compare(gsl_vector * v1, gsl_vector * v2){
    size_t i;
    for(i = 0; i < v1->size; i++){
	int cmp = float_compare(gsl_vector_get(v1, i), gsl_vector_get(v2, i));
	if(cmp != 0) return cmp;
    }
    return 0;
}

static int find_eigenvalues_iter(graph_set * gs, gsl_vector * eigenval, int start, int end){
    if(end - start < 2){
	if(vector_compare(gs->eigenvalues[start], eigenval) == 0)
	    return start;
	else if(vector_compare(gs->eigenvalues[end], eigenval) == 0)
	    return end;
	else return -(start + 1);
    }
    int mid = (start + end)/2;
    int cmp = vector_compare(gs->eigenvalues[mid], eigenval);
    if(cmp == 0) return mid;
    if(cmp == 1) //search [start, mid]
	return find_eigenvalues_iter(gs, eigenval, start, mid);
    else //search [mid, end]
	return find_eigenvalues_iter(gs, eigenval, mid, end);
}

int find_eigenvalues(graph_set * gs, gsl_vector * eigenval){
    if(gs->size == 0) return -1;
    return find_eigenvalues_iter(gs, eigenval, 0, gs->size - 1);
}

#include <gsl/gsl_sort_vector.h>


void insert_graph(graph_set * gs, graph_t * g){
    gsl_vector * eigval = gsl_vector_alloc(NUM_NODES);
    gsl_eigen_symm_workspace * workspace = gsl_eigen_symm_alloc(NUM_NODES);

    gsl_matrix * mat = graph_to_matrix(g);

    gsl_eigen_symm(mat, eigval, workspace);
    gsl_eigen_symm_free(workspace);
    gsl_matrix_free(mat);

    gsl_sort_vector(eigval);

    int index = find_eigenvalues(gs, eigval);
    //The eigenvalues are already in the set
    if(index >= 0){
	gsl_vector_free(eigval);
	gs->graph_set_sizes[index]++;
	gs->graphs[index] = (graph_t **) realloc(gs->graphs[index],
						 gs->graph_set_sizes[index] * sizeof(graph_t *));
	gs->graphs[index][gs->graph_set_sizes[index] - 1] = (graph_t *) malloc(sizeof(graph_t));
	memcpy(&gs->graphs[index][gs->graph_set_sizes[index] - 1]->adj[0][0],
	       &g->adj[0][0], sizeof(g->adj));
	       
    }
    else { //The eigenvalues must be added to the set
	index = -index - 1;
	gs->size++;
	//Reallocating memory
	gs->eigenvalues = (gsl_vector **) realloc(gs->eigenvalues, gs->size * sizeof(gsl_vector *));
	gs->graphs = (graph_t ***) realloc(gs->graphs, gs->size * sizeof(graph_t **));
	gs->graph_set_sizes = (size_t *) realloc(gs->graph_set_sizes, gs->size * sizeof(size_t));
	//Shifting values down
	int i;
	for(i = gs->size - 1; i > index; i--){
	    gs->eigenvalues[i] = gs->eigenvalues[i-1];
	    gs->graphs[i] = gs->graphs[i-1];
	    gs->graph_set_sizes[i] = gs->graph_set_sizes[i-1];
	}
	//Inserting the graph
	gs->eigenvalues[index] = eigval;
	gs->graphs[index] = (graph_t **) malloc(sizeof(graph_t *));
	gs->graphs[index][0] = (graph_t *) malloc(sizeof(graph_t));
	memcpy(&gs->graphs[index][0]->adj[0][0], &g->adj[0][0], sizeof(g->adj));
	gs->graph_set_sizes[index] = 1;
    }
}


#ifdef GRAPH_SET_TEST

static void print_graphs(gsl_vector * eigenvalues, graph_t ** graphs, size_t graph_set_size){
    printf("Eigenvalues:\n");
    print_vector(eigenvalues);
    printf("\nGraphs:\n");
    int i;
    for(i = 0; i < graph_set_size; i++){
	print_graph(graphs[i]);
	printf("\n");
    }
}

static void print_graph_set(graph_set * gs){
    int i;
    for(i = 0; i < gs->size; i++){
	print_graphs(gs->eigenvalues[i], gs->graphs[i], gs->graph_set_sizes[i]);
	printf("\n\n");
    }
}

int main(){
    graph_set * gs = graph_set_alloc();
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
    graph_t * g1 = (graph_t *) malloc(sizeof(graph_t));
    graph_t * g2 = (graph_t *) malloc(sizeof(graph_t));
    graph_t * g3 = (graph_t *) malloc(sizeof(graph_t));

    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	for(j = 0; j < NUM_NODES; j++){
	    if(i == j) continue;
	    g3->adj[i][j] = 1;
	}
    }

    memcpy(&g1->adj[0][0], &adj1[0][0], sizeof(adj1));
    memcpy(&g2->adj[0][0], &adj2[0][0], sizeof(adj2));

    insert_graph(gs, g3);
    insert_graph(gs, g1);
    insert_graph(gs, g2);

    print_graph_set(gs);
}

#endif
