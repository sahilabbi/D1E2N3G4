//
//  ISO_checker.c
//  
//

#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include "graph.h"
#include "distance.h"
#include "comp_graph.h"
#include "isomorphism.h"

void ISOmain();

int num_ISOprocs;

static int
iso_check(graph_t *pg1, graph_t *pg2)
{
    return is_isomorphic(pg1, pg2);
}

int
ISO_Checker(comp_graph* p1, comp_graph* p2)
{
	graph_t g1, g2;

	Compress_graph(&g1, p1, 0);
	Compress_graph(&g2, p2, 0);

	return iso_check(&g1, &g2);	// To be replaced by actual iso_checker
}

static void
process_iso_check_requests()
{
        MPI_Status status;
        int arr_length;
        int ierr, val;
	static comp_graph* parr = NULL;


	if(!parr)
		parr = (comp_graph*)malloc(num_ISOprocs*2*sizeof(comp_graph));
	if(!parr){
		printf("Error: process_iso_check: parr is NULL\n");
		return;
	}
        
        ierr = MPI_Recv ((char*)parr, num_ISOprocs*2*sizeof(comp_graph), MPI_CHAR,
                         MPI_ANY_SOURCE, MPI_ANY_TAG /*COMP_ISO_TAG*/, MPI_COMM_WORLD, &status);
        
        MPI_Get_count(&status, MPI_CHAR, &arr_length);
        arr_length /= (2 * sizeof(comp_graph));
        
        // if any pair in the array is isomorphic, we discard the array
        for(int i = 0; i < arr_length; i++){
            
            //call ISO Checker, if they are iso, return -1, else return index or first element in first pair.
            // val = ISO_Checker(&rarr[i][0], &rarr[i][1]); // to be implemented
            val = ISO_Checker(&parr[2*i], &parr[2*i + 1]);

            if(val){
                val = -1;	//if isomorphism, return 1
                ierr = MPI_Send(&val, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                break;
            }

        }
        
        //if no isomorphism, return 1
        val = 1;
        ierr = MPI_Send(&val, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
}


int main(argc,argv)
int argc;
char *argv[];
{
    

    int numprocs, myid;
    
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    // printf("I'm process %d out of %d processes\n", myid, numprocs);
    num_ISOprocs = numprocs - 1;
    
    
    
    //very first node, this will send out array of graph pairs to each available node to check for isomorphism in any pair

    if(myid == 0){
        
        ISOmain();
        
    }
    
    if(myid){	// receive array of graph pairs to check for isomorphism in any pair
	for(;;)
		process_iso_check_requests();
        
        
    }
    
    MPI_Finalize();
    
}
