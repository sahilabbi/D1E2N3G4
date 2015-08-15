#include "graph.h"
#include <stdlib.h>
#include <string.h>

void copy_graph(graph_t * dst, graph_t * src){
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	for(j = 0; j < NUM_NODES; j++){
	    dst->adj[i][j] = src->adj[i][j];
	}
    }
}

//gives whether a node is one of the initial edges.
//these edges are not changed to avoid isomorphisms,
//as any graph obtained by altering them would be isomporphic
//to some other graph  obtained without altering them
static bool is_protected(int row, int column){
    //If ISO_TEST is defined then we are cycling through
    //all the graphs, so only the diagonal is protected
#ifndef ISO_TEST
    if(row == 0) return true;
#endif
    if(row == column) return true;
#ifdef ISO_TEST
    return false;
#else
    if(row - column != 1 && column - row != 1) return false;
    int x = row > column ? row : column;
    return x % 2 == NODE_DEGREE % 2;
#endif
}

//moves the cursor right to left on the row
//if the middle of the graph has been reached,
//then the cursor jumps up the the right edge of
//the row above
static void increment_cursor(int * i, int * j){
    *j -= 1;
    if(is_protected(*i,*j)){
	*i -= 1;
	*j = NUM_NODES-1;
	if(*i < 0) return;
    }
}

bool is_regular(int * degree){
    int i;
    for(i = 0; i < NUM_NODES; i++){
	if(degree[i] != NODE_DEGREE) return false;
    }
    return true;
}

#include <stdio.h>

bool increment_graph(graph_t * graph){
    int degree[NUM_NODES];
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	int sum = 0;
	for(j = 0; j < NUM_NODES; j++){
	    sum += graph->adj[i][j];
	}
	degree[i] = sum;
    }
    int iStart = NUM_NODES-2;
    int jStart = NUM_NODES-1;
    //incase the first edge selected is one of the initial edges
    if(is_protected(iStart,jStart)) increment_cursor(&iStart,&jStart);
    do {
	i = iStart;
	j = jStart;
	while(graph->adj[i][j] == 1){
	    graph->adj[i][j] = 0;
	    graph->adj[j][i] = 0;
	    degree[i]--;
	    degree[j]--;
	    increment_cursor(&i,&j);
	    if(i < 0) return true;
	}
	graph->adj[i][j] = 1;
	graph->adj[j][i] = 1;
	degree[i]++;
	degree[j]++;
	//	printf("(i,j) = (%d,%d)\n", i, j);
	//	print_graph(graph);
	//	printf("degree:\n");
	//	int k;
	//	for(k = 0; k < NUM_NODES; k++) printf("%d ", degree[k]);
	//	printf("\nis_regular(degree): %d\n\n", is_regular(&degree[0]));
    } while(!is_regular(&degree[0]));
    return false;
}

graph_t * initial_graph(){
    graph_t * ret = (graph_t *) malloc(sizeof(graph_t));
    memset(&ret->adj[0][0], 0, sizeof(ret->adj));
    int i;
    for(i = 1; i <= NODE_DEGREE; i++){
	ret->adj[0][i] = 1;
	ret->adj[i][0] = 1;
    }

    for(i = NODE_DEGREE+1; i < NODE_DEGREE; i += 2){
	ret->adj[i][i+1] = 1;
	ret->adj[i+1][i] = 1;
    }
	
    return ret;
}

#include <stdio.h>
void print_graph(graph_t * g){
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	for(j = 0; j < NUM_NODES; j++){
	    printf("%d ",g->adj[i][j]);
	}
	printf("\n");
    }
}

#ifdef INCREMENT_TEST

int main(){
    graph_t * g = initial_graph();
    increment_graph(g);
    print_graph(g);
    printf("\n\n\n");
    int i;
    bool done = false;
    for(; !done; ){
	done = increment_graph(g);
	print_graph(g);
	printf("\n");

    }



    printf("is_protected(%d,%d)=%d\n",NODE_DEGREE+1,NODE_DEGREE+2,
	   is_protected(NODE_DEGREE+1,NODE_DEGREE+2));
    return 0;
}


#endif
