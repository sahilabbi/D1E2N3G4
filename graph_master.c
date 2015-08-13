#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "graph.h"
#include "distance.h"
#include "comp_graph.h"



static comp_graph	cgraph_mast[MAX_GRAPH_ARR];


static int num_procs, graph_arr_size, n_graph_files;

int lowest_diam = L; // defined in Makefile
float lowest_dist = 99999.0;

#define EPSILON 0.0001

// copied for convenience,  original copy in graph.c
static bool float_equals(double a, double b){
    double diff = a - b;
    if(diff < 0) diff = -diff;
    return diff < EPSILON;
}

static int float_compare(double a, double b){
    double diff = a - b;
    if(diff > EPSILON) return 1;
    if(diff < -EPSILON) return -1;
    return 0;
}



void
print_graph_array()
{

    printf("Array size %d\n", graph_arr_size);
}

void
del_graph_files()
{
	for(int i=0; i < n_graph_files; i++){
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


void
add_graph(compress_graph* p)
{
	if(p->diameter < lowest_diam || (p->diameter == lowest_diam && p->avg_dist < lowest_dist)){
		del_graph_files();
       		graph_arr_size = 1;
       		cgraph_mast[0] = *((comp_graph*)p);
        	cgraph_mast[0] = *((comp_graph*)p);
       		lowest_diam = p->diameter;
		lowest_dist= p->avg_dist;
	} else if(p->diameter == lowest_diam && float_equals(p->avg_dist, lowest_dist)){
        	cgraph_mast[graph_arr_size ++] = *((comp_graph*)p);
		if(graph_arr_size >= MAX_GRAPH_ARR)
			save_graphs_file();
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

	int i;
	for(; ; ){
	  if(no_more_graphs){
		save_graphs_file();
		break;
	  }
	  if(ng >= MAX_GRAPH_ARR){
		save_graphs_file();
		ng=0;
	  }
          for(i=1; i < num_procs; i++){
		ng++;
		if(true == increment_graph(g)){
			printf("incr graphs returned false %d\n", graph_arr_size);
			no_more_graphs=1;
			break;
		}
		Compress_graph(g, (comp_graph*)&cg, 1);
		// printf("gen graphs: bef send\n");
                ierr = MPI_Send(&cg, sizeof(compress_graph), MPI_CHAR, 
						i, COMP_DIST_TAG, MPI_COMM_WORLD);
		// printf("gen graphs: aft send\n");
          }
          for(i=1; i<num_procs; i++){
		  // printf("BEF-mpi-recv %d\n", i);
                  ierr = MPI_Recv (&cg, sizeof(compress_graph), MPI_CHAR, 
                     MPI_ANY_SOURCE, MPI_ANY_TAG /*COMP_DIST_TAG*/, MPI_COMM_WORLD, &status);
		  if(0 == (++recv_count) % 1000)
		 	printf("%d AFT-mpi-recv %d Diam: %d Avg-Dist: %f %d\n", ng, i, cg.diameter, cg.avg_dist, recv_count);
		  if(cg.diameter < lowest_diam || cg.diameter==lowest_diam){
			add_graph(&cg);
					  
		  } else {
				// disregard
		  }
	  }
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
		ierr = MPI_Recv (&cg, sizeof(compress_graph), MPI_CHAR, MPI_ANY_SOURCE, COMP_DIST_TAG, MPI_COMM_WORLD, &status);
		// printf("process-dist: aft recv\n");
		Compress_graph(&g, (comp_graph*)&cg, 0);
		graph_info *info = getinfo((int*)g.adj, N);
		// print_graph(&g);
		cg.diameter = info->diameter;
		cg.avg_dist = (float)info->sum_of_distances / 56.0f;
		// printf("Diameter: %d Average: %f\n", (int)cg.diameter, cg.avg_dist);
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



     ierr = MPI_Finalize();
      
     return 0;
}
