#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
in tmain()
{
    int *min_local;
    int rank_number,cpu_number;
    int *min_global,local_min=10;
    min_global = (int *)malloc(8*sizeof(int));
    min_local = (int *)malloc(8*sizeof(int));
      MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_number);
    MPI_Comm_size(MPI_COMM_WORLD, &cpu_number);
    for(int i=0;i<8;i++)
        min_local[i] = i;
    if(rank_number==5)
        min_local[rank_number]=10;

     MPI_Gather(&local_min,1,MPI_INT,min_global,1,MPI_INT,0,MPI_COMM_WORLD);

    for(int i=0;i<cpu_number;i++)
        printf("%d  ",min_global[i]);
        MPI_Finalize();
        return 0;
}