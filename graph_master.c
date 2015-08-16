#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "graph.h"
#include "distance.h"
#include "comp_graph.h"
#include "graph_set.h"



static comp_graph     cgraph_mast[MAX_GRAPH_ARR];


static int num_procs, graph_arr_size, n_graph_files;

int lowest_diam = L; // defined in Makefile
unsigned int lowest_dist = L * NUM_NODES * NUM_NODES + 5;


void
print_graph_array()
{

    printf("Array size %d\n", graph_arr_size);
}

void
del_graph_files()
{
    int i;
    for(i=0; i < n_graph_files; i++){
	char    fpath[100];
	sprintf(fpath, "%s%d.%s", FPATH, i, "ori");
	unlink(fpath);
    }
    n_graph_files=0;
}

int
save_graphs_file()
{
    char    fpath[100];
    sprintf(fpath, "%s%d.%s", FPATH, n_graph_files, "ori");
    printf("Saving graph files to %s...\n", fpath);
    FILE* fp = fopen(fpath, "w");

    if(fp){

        fwrite((char* )cgraph_mast, sizeof(comp_graph),  graph_arr_size, fp);
        fclose(fp);
	n_graph_files++;
       	graph_arr_size = 0;
	return 1;
    }
    return 0;
}

#include <string.h>


void
add_graph(graph_set ** all_graphs, compress_graph* p)
{
    //printf("Entering add_graph\n");
    graph_t decompressed_graph;
    comp_graph comp;
    memcpy(&comp.comp[0], &p->comp[0], sizeof(p->comp));
    Compress_graph(&decompressed_graph, &comp, 0);
    /*printf("Adding graph:\n");
    print_graph(&decompressed_graph);
    printf("\n");*/
	if(p->diameter < lowest_diam || (p->diameter == lowest_diam && p->sum_dist < lowest_dist)){
	    //printf("Found better graph: Diameter=%d, Distance Sum=%d\n", p->diameter, p->sum_dist);
	    delete_graph_set(*all_graphs);
	    *all_graphs = graph_set_alloc();
	    insert_graph(*all_graphs, &decompressed_graph);
	    /*del_graph_files();
       		graph_arr_size = 1;
       		cgraph_mast[0] = *((comp_graph*)p);
        	cgraph_mast[0] = *((comp_graph*)p);*/
       		lowest_diam = p->diameter;
		lowest_dist= p->sum_dist;
	} else if(p->diameter == lowest_diam && p->sum_dist == lowest_dist){
	    /*cgraph_mast[graph_arr_size ++] = *((comp_graph*)p);
		if(graph_arr_size >= MAX_GRAPH_ARR)
		save_graphs_file();*/
	    //printf("Found optimal graph, checking isomorphism\n");
	    if(!check_isomorphism(*all_graphs, &decompressed_graph)){
		//printf("Graph found not isomorphic.  Inserting...\n");
		insert_graph(*all_graphs, &decompressed_graph);
		//printf("Done.\n");
	    }
    	} else {
		// discard
	}
}



// generates graph, sends to slaves, waits for diameter and average distance
void
generate_graphs()
{
    int	no_more_graphs=0;
    int 	ierr;
    int	recv_count=0;
	
    int ng=0;
    compress_graph	cg;
    MPI_Status status;


    graph_t * g = initial_graph();

    graph_set ** all_graphs = malloc(sizeof(graph_set *));
    *all_graphs = graph_set_alloc();

    printf("Starting best average: %d\n", lowest_dist);

    int i;
    for(; ; ){
	if(no_more_graphs){
	    printf("Printing out graphs...\n");
	    //print_graph_set(*all_graphs);
	    //printf("\n");
	    print_graph_set_compressed(*all_graphs);
	    break;
	}
	/*if(ng >= MAX_GRAPH_ARR){
	  save_graphs_file();
	  ng=0;
	  }*/
	int num_sent = 1;
	for(i=1; i < num_procs; i++){
	    ng++;
	    if(true == increment_graph(g)){
		printf("incr graphs returned false %d\n", graph_arr_size);
		no_more_graphs=1;
		break;
	    }
	    Compress_graph(g, (comp_graph*)&cg, 1);
	    //printf("gen graphs: bef send\n");
	    ierr = MPI_Send(&cg, sizeof(compress_graph), MPI_CHAR, 
			    i, COMP_DIST_TAG, MPI_COMM_WORLD);
	    num_sent++;
	    //printf("gen graphs: aft send\n");
	}
	for(i=1; i<num_sent; i++){
	    // printf("BEF-mpi-recv %d\n", i);
	    ierr = MPI_Recv (&cg, sizeof(compress_graph), MPI_CHAR, 
			     MPI_ANY_SOURCE, MPI_ANY_TAG /*COMP_DIST_TAG*/,
			     MPI_COMM_WORLD, &status);
	    if(0 == (++recv_count) % 1000)
		printf("%d AFT-mpi-recv %d Diam: %d Avg-Dist: %d %d\n",
		       ng, i, cg.diameter, cg.sum_dist, recv_count);
	    //printf("Received graph, Diameter %d, Sum dists %d\n", cg.diameter, cg.sum_dist);
	    if(cg.diameter < lowest_diam || cg.diameter==lowest_diam){
		add_graph(all_graphs, &cg);
					  
	    } else {
		// disregard
	    }
	}
    }
    for(i = 1; i < num_procs; i++){
	ierr = MPI_Send(NULL, 0, MPI_CHAR, i, END_TASK_TAG, MPI_COMM_WORLD);
    }
}

// receives graph , computes diameter/distance and sends back to graph generator
void
process_distance_calcs()
{
        int ierr;
	compress_graph cg;
	graph_t	g;
      	MPI_Status status;

	for( ; ; ){
		// printf("process-dist: bef recv\n");
		ierr = MPI_Recv (&cg, sizeof(compress_graph), MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//If the program is over, end the program
		if(status.MPI_TAG == END_TASK_TAG) break;
		// printf("process-dist: aft recv\n");
		Compress_graph(&g, (comp_graph*)&cg, 0);
		graph_info *info = getinfo((int*)g.adj, NUM_NODES);
		//printf("Recieved Graph:\n");
		//print_graph(&g);
		cg.diameter = info->diameter;
		cg.sum_dist = info->sum_of_distances;
		//printf("Diameter: %d Average: %d\n", (int)cg.diameter, cg.sum_dist);
                ierr = MPI_Send(&cg, sizeof(compress_graph), MPI_CHAR, 0, COMP_DIST_TAG, MPI_COMM_WORLD);
	}

}



int 
main(int argc, char **argv)
{
      int ierr, my_id; // i; recv_dist;


      ierr = MPI_Init(&argc, &argv);

      /* find out MY process ID, and how many processes were started. */

      ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
      ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

      printf("Graph Master! I'm process %d out of %d processes\n",
        		my_id, num_procs);


	if(my_id == 0){
		generate_graphs();
	} else {
		process_distance_calcs();

	}

	printf("Process %d checking out\n", my_id);



     ierr = MPI_Finalize();

     return 0;
}
