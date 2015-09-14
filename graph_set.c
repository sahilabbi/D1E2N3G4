#include "graph_set.h"
#include "isomorphism.h"
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

graph_set * graph_set_alloc(){
    graph_set * ret = malloc(sizeof(graph_set));
    ret->iso_groups = NULL;
    ret->size = 0;
    ret->amount_alloced = 0;
    return ret;
}


void delete_graph_set(graph_set * gs){
    int i,j;
    //printf("Deleting graph set %d.\ngs->size: %d\n", gs, gs->size);
    //print_graph_set(gs);
    for(i = 0; (size_t) i < gs->size; i++){
	//printf("Deleting group %d: %d\n", i, gs->iso_groups + i);
	isospectral_group * currGroup = gs->iso_groups + i;
	gsl_vector_free(currGroup->eigenvalues);
	for(j = 0; (size_t) j < currGroup->num_graphs; j++)
	    //printf("\tDeleting graph %d\n", j);
	    free(currGroup->graphs[j]);
    }
    //printf("Done deleting groups\n");
    if(gs->iso_groups != NULL) free(gs->iso_groups);
    free(gs);
}

int vector_compare(const void * a, const void * b){
    gsl_vector * v1 = (gsl_vector *) a;
    gsl_vector * v2 = (gsl_vector *) b;
    size_t i;
    for(i = 0; i < v1->size; i++){
	int cmp = float_compare(gsl_vector_get(v1, i), gsl_vector_get(v2, i));
	if(cmp != 0) return cmp;
    }
    return 0;
}

int isospectral_group_compare(const void * a, const void * b){
    isospectral_group * gr1 = (isospectral_group *) a;
    isospectral_group * gr2 = (isospectral_group *) b;
    return vector_compare(gr1->eigenvalues, gr2->eigenvalues);
}

isospectral_group * find_eigenvalues(graph_set * gs, gsl_vector * eigenval){
    if(gs->size == 0) return NULL;
    //printf("Got here!\n");
    /*printf("gs: %d\n", gs);
    printf("gs->size: %d\n", gs->size);
    printf("gs->iso_groups: %d\n", gs->iso_groups);
    printf("gs->iso_groups->eigenvalues: %d\n", gs->iso_groups->eigenvalues);
    print_vector(gs->iso_groups->eigenvalues);
    printf("\n\n");*/
    isospectral_group key = {.eigenvalues = eigenval,
			     .graphs = NULL,
			     .num_graphs = 0};
    isospectral_group * ret = (isospectral_group *) bsearch(&key,
							    gs->iso_groups,
							    gs->size,
							    sizeof(isospectral_group),
							    isospectral_group_compare);
    return ret;
}

#include <gsl/gsl_sort_vector.h>


void insert_graph(graph_set * gs, graph_t * g){
    gsl_vector * eigval = gsl_vector_alloc(NUM_NODES);
    gsl_eigen_symm_workspace * workspace = gsl_eigen_symm_alloc(NUM_NODES);

    gsl_matrix * mat = graph_to_matrix(g);

    //printf("Eigenvalue size: %d\nMatrix size: %dx%d\n", eigval->size,
//	   mat->size1, mat->size2);

    gsl_eigen_symm(mat, eigval, workspace);
    gsl_eigen_symm_free(workspace);
    gsl_matrix_free(mat);

    gsl_sort_vector(eigval);

    /*   printf("Inserting Vector:\n");
    print_vector(eigval);
    printf("\n"); */

    isospectral_group * index = find_eigenvalues(gs, eigval);
    if(index == NULL){ //The eigenvalues need to be added to the set
	//Create the new group to be added to gs
	isospectral_group * newGroup = malloc(sizeof(isospectral_group));
	newGroup->eigenvalues = eigval;
	newGroup->graphs = (comp_graph **) malloc(sizeof(comp_graph *));
	newGroup->graphs[0] = (comp_graph *) malloc(sizeof(comp_graph));
	newGroup->num_graphs = 1;
	Compress_graph(g, newGroup->graphs[0], 1);

	//allocate the extra memory for the new set
	gs->size++;
	if(gs->size > gs->amount_alloced){
	    if(gs->amount_alloced == 0) gs->amount_alloced = 1;
	    else gs->amount_alloced *= 2;
	    gs->iso_groups = realloc(gs->iso_groups,
				     gs->amount_alloced * sizeof(isospectral_group));
	}

	//Shifting the groups to the right until we reach the spot for the new group
	int i;
	for(i = gs->size - 1; i > 0; i--){
	    //If the next group to be shifted is less than newGroup, insert newGroup at index i
	    if(isospectral_group_compare(&gs->iso_groups[i-1], newGroup) == -1){
		memcpy(&gs->iso_groups[i], newGroup, sizeof(isospectral_group));
		break;
	    }
	    else
		memcpy(&gs->iso_groups[i], &gs->iso_groups[i - 1], sizeof(isospectral_group));
	}
	//If newGroup is smaller than every element in the group
	if(i == 0) memcpy(gs->iso_groups, newGroup, sizeof(isospectral_group));

	//removes only the top-level elements so the elements in gs are not deleted
	free(newGroup);
    }
    else { //The eigenvalues are already in the set
	//reallocate the memory to allow for the new graph
	index->num_graphs++;
	index->graphs = (comp_graph **) realloc(index->graphs, index->num_graphs * sizeof(comp_graph *));
	index->graphs[index->num_graphs - 1] = (comp_graph *) malloc(sizeof(graph_t));

	//Inserting the graph
	Compress_graph(g, index->graphs[index->num_graphs - 1], 1);

	//Since the vector was already in the set we can free the memory
	gsl_vector_free(eigval);
    }
}

bool check_isomorphism(graph_set * gs, graph_t * g){
    gsl_vector * eigval = gsl_vector_alloc(NUM_NODES);
    gsl_eigen_symm_workspace * workspace = gsl_eigen_symm_alloc(NUM_NODES);

    gsl_matrix * mat = graph_to_matrix(g);

    gsl_eigen_symm(mat, eigval, workspace);
    gsl_eigen_symm_free(workspace);
    gsl_matrix_free(mat);

    gsl_sort_vector(eigval);

    isospectral_group * index = find_eigenvalues(gs, eigval);
    gsl_vector_free(eigval);
    //If the eigenvalues are not in the set, the graph is not ismorphic
    if(index == NULL) return false;
    graph_t decompressed_graph;
    int i;
    for(i = 0; (size_t) i < index->num_graphs; i++){
	Compress_graph(&decompressed_graph, index->graphs[i], 0);
	if(is_isomorphic(&decompressed_graph, g)) return true;
    }
    return false;
}




static void print_isospectral_group(isospectral_group * gr){
    printf("Eigenvalues:\n");
    print_vector(gr->eigenvalues);
    printf("\nGraphs:\n");
    int i;
    graph_t decompressed_graph;
    for(i = 0; (size_t) i < gr->num_graphs; i++){
	Compress_graph(&decompressed_graph, gr->graphs[i], 0);
	print_graph(&decompressed_graph);
	printf("\n");
    }
}

void print_graph_set(graph_set * gs){
    int i;
    for(i = 0; (size_t) i < gs->size; i++){
	print_isospectral_group(&gs->iso_groups[i]);
	printf("\n\n");
    }
}

//switches the order of the bits
//See Hacker's Delight p. 101
static unsigned char flip_byte(unsigned char c){
    char x = c;
    x = (x & 0x55) << 1 | (x & 0xAA) >> 1;
    x = (x & 0x33) << 2 | (x & 0xCC) >> 2;
    x = (x & 0x0F) << 4 | (x & 0xF0) >> 4;
    return x;
}

void print_graph_compressed(comp_graph * cg){
    int i;
    for(i = 0; i < COMPRESS_SIZE; i++){
	printf("%02x", flip_byte(cg->comp[i]));
    }
}

static void print_isospectral_group_compressed(isospectral_group * gr){
    int i;
    for(i = 0; (size_t) i < gr->num_graphs; i++){
	print_graph_compressed(gr->graphs[i]);
	printf("\n");
    }
}

void print_graph_set_compressed(graph_set * gs){
    int i;
    for(i = 0; (size_t) i < gs->size; i++){
	print_isospectral_group_compressed(&gs->iso_groups[i]);
    }
}

#ifdef GRAPH_SET_TEST

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
	    if(i == j) g3->adj[i][j] = 0;
	    else g3->adj[i][j] = 1;
	}
    }

    memcpy(&g1->adj[0][0], &adj1[0][0], sizeof(adj1));
    memcpy(&g2->adj[0][0], &adj2[0][0], sizeof(adj2));

    insert_graph(gs, g3);
    insert_graph(gs, g1);
    insert_graph(gs, g2);

    print_graph_set(gs);

    printf("\nIsomorphism check: %d, %d\n", check_isomorphism(gs, g2),
	   check_isomorphism(gs, g3));
}

#endif
