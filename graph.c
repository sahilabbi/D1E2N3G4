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
    if(row == 0 || column == 0) return true;
#endif
    //Check if we're on the diagonal
    if(row == column) return true;
#ifdef ISO_TEST
    return false;
#else
    //The scattered nodes along the diagnoal don't exist in the
    //first K rows and columns, so if we're there then we don't
    //need to worry about them
    if(row <= NODE_DEGREE || column <= NODE_DEGREE) return false;
    //Checking for the scattered nodes along the diagnoal
    if(row - column != 1 && column - row != 1) return false;
    int x = row > column ? row : column;
    return x % 2 == NODE_DEGREE % 2;
#endif
}

bool is_regular(int * degree){
    int i;
    for(i = 0; i < NUM_NODES; i++){
	if(degree[i] != NODE_DEGREE) return false;
    }
    return true;
}

#include <stdio.h>

/*//increments the graph by one, starting the cursor at (istart, jStart)
static bool increment_graph_single_from(graph_t * graph, int * degree, int iStart, int jStart){
	int i = iStart;
	int j = jStart;
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
	//printf("(i,j) = (%d,%d)\n", i, j);
	//print_graph(graph);
	//printf("degree:\n");
	//int k;
	//for(k = 0; k < NUM_NODES; k++) printf("%d ", degree[k]);
	//printf("\n");
	//printf("\nis_regular(degree): %d\n\n", is_regular(&degree[0]));

	return false;
	}*/

//Gives the number of non-protected elements to
//the right of the diagnoal in the row
//Used in increment_graph()
static int lengthOfRowRemaining(int row){
    int ret = NUM_NODES - row - 1;
    if(ret <= 0) return 0;
    if(is_protected(row, row + 1)) ret--;
    return ret;
}

static bool lexographicaly_increment(graph_t * g, int * degree, int row){
    //printf("In:\n");
    //print_graph(g);
    
    int length = lengthOfRowRemaining(row);
    int * nums = &g->adj[row][NUM_NODES] - length;
    int i;
    //start and end will denote the bounds of the rightmost
    //grouping of 1s
    int start = -1;
    int end = -1;
    for(i = length - 1; i >= 0; i--){
	if(nums[i] == 1 && start == -1) start = i;
	if(nums[i] == 0 && start > -1){
	    end = i;
	    break;
	}
    }
    //If this is the last set
    if(end == -1) return true;

    //Moving start and end to where they are in the full graph
    start = start + NUM_NODES - length;
    end = end + NUM_NODES - length;

    //for(i = 0; i < NUM_NODES; i++) printf("%d ", degree[i]);
    //printf("\n");
    //printf("(Start, End): (%d, %d)\n", start, end);

    //move the leftmost 1 to the left, reset all the others to the right
    g->adj[row][end] = 1;
    g->adj[end][row] = 1;
    g->adj[row][end + 1] = 0;
    g->adj[end + 1][row] = 0;

    degree[end]++;
    degree[end + 1]--;

    int k = NUM_NODES - 1;
    for(i = start; i > end + 1; i--){
	if(g->adj[row][i] == 1) degree[i]--;
	g->adj[row][i] = 0;
	g->adj[i][row] = 0;

	if(g->adj[row][k] == 0) degree[k]++;
	g->adj[row][k] = 1;
	g->adj[k][row] = 1;

	k--;
    }

    //printf("Out:\n");
    //print_graph(g);
    //for(i = 0; i < NUM_NODES; i++) printf("%d ", degree[i]);
    //printf("\n");
    return false;
    
}

//Starts from a regular graph and gets the next regular graph
bool increment_graph(graph_t * graph){

    bool test = false;
    int degree[NUM_NODES];
    int i,j;
    for(i = 0; i < NUM_NODES; i++){
	int sum = 0;
	for(j = 0; j < NUM_NODES; j++){
	    sum += graph->adj[i][j];
	}
	degree[i] = sum;
    }
    int currNode = NUM_NODES - 2;
    if(lengthOfRowRemaining(currNode) == 0) currNode--;
    do {

	//PART 1: increment the current row
	
	//Get the next permutation of the 1s
	bool layer_done = lexographicaly_increment(graph, degree, currNode);

	if(test){
	    printf("Entering loop at row %d:\n", currNode);
	    print_graph(graph);
	    printf("\n");
	}
	//PART 2: if this row is bad, move up
	
	//If we're at the last permutation for this layer, move up
	//Also moves up if the current node is not of the correct degree,
	//which should only happen under first_graph()
	if(layer_done || degree[currNode] != NODE_DEGREE){
	    //If we've reached the top, there are no more graphs so we return true
	    if(currNode <= 1){
		//printf("Reached top row.  Returning true\n\n");
		//print_graph(graph);
		return true;
	    }
	    //Zero out this layer
	    for(j = currNode; j < NUM_NODES; j++){
		if(is_protected(currNode, j)) continue;
		if(graph->adj[currNode][j] == 1){
		    degree[currNode]--;
		    degree[j]--;
		}
		graph->adj[currNode][j] = 0;
		graph->adj[j][currNode] = 0;
	    }
	    currNode--;
	    continue;
	}

	//PART 3: if this row is good, move down
	
	//If we've found a possibly valid permutation, move down
	//if there we've given a node too many commections, don't move down
	bool needRestart = false;
	for(i = 0; i < NUM_NODES; i++)
	    if(degree[i] > NODE_DEGREE){
		//printf("Restarting loop due to node %d having %d connections\n\n", i, degree[i]);
		needRestart = true;
		break;
	    }
	if(!needRestart) {
	    //If the next row doesn't have enough places to give it NODE_DEGREE connections, try again
	    if(lengthOfRowRemaining(currNode + 1) < NODE_DEGREE - degree[currNode + 1]){
		//printf("Restarting loop since row %d needs %d empty slots but only has %d\n\n", currNode + 1,
		//       NODE_DEGREE - degree[currNode + 1], lengthOfRowRemaining(currNode + 1));
		continue;
	    }
	    
	    //If we've finished, break
	    if(is_regular(&degree[0])) break;
	    
	    //otherwise, move down
	    currNode++;
	    //printf("Moving down to row %d\n", currNode);
	    //printf("Adding %d 1s to the row\n", NODE_DEGREE - degree[currNode]);
	    //Fill in the first 
	    for(i = NUM_NODES - 1; degree[currNode] < NODE_DEGREE; i--){
		if(is_protected(currNode, i)){
		    //The code shouldn't get here due to the check above, which makes
		    //sure we have enough empty solts
		    printf("Error in increment_graph(): check lengthOfRowRemaining\n");
		    return false;
		}
		graph->adj[currNode][i] = 1;
		graph->adj[i][currNode] = 1;
		degree[currNode]++;
		degree[i]++;
	    }
	}
    } while(true);
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

    for(i = NODE_DEGREE+1; i < NUM_NODES - 1; i += 2){
	ret->adj[i][i+1] = 1;
	ret->adj[i+1][i] = 1;
    }

    return ret;
}

graph_t * first_graph(){
    graph_t * ret = initial_graph();
    int i, j;
    //This pattern will force increment_graph to
    //quickly get to an early state
    for(i = 1; i < NUM_NODES; i++){
	for(j = NUM_NODES - 1; j > NUM_NODES - NODE_DEGREE; j--){
	    if(is_protected(i,j)) break;
	    ret->adj[i][j] = 1;
	    ret->adj[j][i] = 1;
	}
    }
    /*printf("Initial graph:\n");
    print_graph(ret);
    printf("\n\n");*/
    increment_graph(ret);
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
    graph_t * g = first_graph();
    printf("Initial graph:\n");
    print_graph(g);
    printf("\n");
    graph_t * copy = malloc(sizeof(graph_t));
    bool done = false;
    for(; !done; ){
	copy_graph(copy, g);
	done = increment_graph(g);
	print_graph(g);
	printf("\n");

    }

    printf("Copy:\n");
    print_graph(copy);

    //increment_graph(copy, true);
    
    printf("is_protected(%d,%d)=%d\n",NODE_DEGREE+1,NODE_DEGREE+2,
    is_protected(NODE_DEGREE+1,NODE_DEGREE+2));
    
    return 0;
}


#endif
