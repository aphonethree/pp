#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
int main(int argc,char* argv[])
{
    int *min_local;
    int rank_number,cpu_number,tag;
    int *min_global,local_min=5;
    min_global = (int *)malloc(8*sizeof(int));
    min_local = (int *)malloc(8*sizeof(int));
    int *vector_swap_forward = (int*)malloc(1*8*sizeof(int));
     int *read_buf_back = (int*)malloc(8*sizeof(int));
     MPI_Request request;
    MPI_Status status;

    MPI_DATATYPE COLUMN;
    MPI_Type_contiguous(8, MPI_INT,&COLUMN);
    MPI_Type_commit(&COLUMN);


      MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_number);
    MPI_Comm_size(MPI_COMM_WORLD, &cpu_number);
    for(int i=0;i<8;i++)
       vector_swap_forward[i]=i;
    if(rank_number==1)
        local_min=10;
    if(rank_number==2)
	local_min=20;

     MPI_Gather(&local_min,1,MPI_INT,min_global,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank_number==0){
   	 for(int i=0;i<cpu_number;i++)
        	printf("%d  \n",min_global[i]);
	}
    if(rank_number==1){
        for(int i=0;i<8)
            printf("%d  ",vector_swap_forward[i]);
        MPI_Isend(vector_swap_forward,1,COLUMN,2,tag,MPI_COMM_WORLD,&request);
    }
    if(rank_number==2){
        MPI_Irecv(read_buf_back,1,COLUMN,1,tag,MPI_COMM_WORLD,&request);
        for(int i=0;i<8)
            printf("%d  ",read_buf_back[i]);
    }
        MPI_Finalize();
        return 0;
}
