#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#ifndef W
#define W 20                                    // Width
#endif
int main(int argc, char **argv) {
  int L = atoi(argv[1]);                        // Length
  int iteration = atoi(argv[2]);                // Iteration
  srand(atoi(argv[3]));                         // Seed
  float d = (float) random() / RAND_MAX * 0.2;  // Diffusivity
  int *temp = malloc(L*W*sizeof(int));          // Current temperature
  int *next = malloc(L*W*sizeof(int));          // Next time step
  int rank_number,cpu_number;
  int *vector_swap_forward = (int*)malloc(1*W*sizeof(int));
  int *vector_swap_backward = (int*)malloc(1*W*sizeof(int));
  int *read_buf_front = (int*)malloc(W*sizeof(int));
  int *read_buf_back = (int*)malloc(W*sizeof(int));
    MPI_Datatype COLUMN;
   


  int local_l,tag=1,min=0,read_buff_min,read_buff_balance,*global_min,*global_balance;
  MPI_Request request;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_number);
  MPI_Comm_size(MPI_COMM_WORLD, &cpu_number);
   MPI_Type_contiguous(W, MPI_INT,&COLUMN);
    MPI_Type_commit(&COLUMN);
  global_min = (int*) malloc(cpu_number*sizeof(int));
  global_balance = (int*) malloc(cpu_number*sizeof(int));
  /*for(int i=0;i<cpu_number;i++)
    finish_number[i]=0;*/
  local_l = L/cpu_number;
   
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < W; j++) {
      temp[i*W+j] = random()>>3;
    }
  }
  
  int local_count = 0, balance = 0, local_min;
  while (iteration--) {       // Compute with up, left, right, down points
    balance = 1;

    if(cpu_number==1){
        local_count++;
        for (int i = local_l*rank_number; i < local_l*rank_number+local_l; i++) {
            for (int j = 0; j < W; j++) {
                    float t = temp[i*W+j] / d;
                    t += temp[i*W+j] * -4;
                    t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
                    t += temp[(i + 1 >= L ? i : i + 1)*W+j];
                    t += temp[i*W+(j - 1 <  0 ? 0 : j - 1)];
                    t += temp[i*W+(j + 1 >= W ? j : j + 1)];
                    t *= d;
                    next[i*W+j] = t ;
                    if (next[i*W+j] != temp[i*W+j]) {
                        balance = 0;
                    }
                    if(t<local_min)
                        local_min = t;
            }
        }
        if (balance) {
            break;
        }
        int *tmp = temp;
        temp = next;
        next = tmp;
    }else{
         for (int i = local_l*rank_number; i < local_l*rank_number+local_l; i++) {
            for (int j = 0; j < W; j++) {
                float t = temp[i*W+j] / d;
                t += temp[i*W+j] * -4;
                t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
                t += temp[(i + 1 >= L ? i : i + 1)*W+j];
                t += temp[i*W+(j - 1 <  0 ? 0 : j - 1)];
                t += temp[i*W+(j + 1 >= W ? j : j + 1)];
                t *= d;
                next[i*W+j] = t ;
                if (next[i*W+j] != temp[i*W+j]) {
                    balance = 0;
                }
                
                if(rank_number==0 && i==(local_l*rank_number+local_l-1)){

                    vector_swap_backward[j]= next[i*W+j];

                }else if(rank_number>0 && rank_number<(cpu_number-1)){

                    if(i==local_l*rank_number)
                        vector_swap_forward[j]= next[i*W+j]; 
                    else if(i==(local_l*rank_number+local_l-1))
                        vector_swap_backward[j]= next[i*W+j]; 

                }else if(rank_number==(cpu_number-1)&&i==local_l*rank_number ){
                    vector_swap_forward[j]= next[i*W+j]; 
                }          
            }
        }
        
        if(rank_number==0){
            MPI_Isend(vector_swap_backward,1,COLUMN,rank_number+1,tag,MPI_COMM_WORLD,&request);

            MPI_Recv(read_buf_back,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&status);
            for(int i=0;i<W;i++)
                next[(local_l*rank_number+local_l)*W+i] = read_buf_back[i];
            
            local_count++;
        }else if(rank_number>0 && rank_number<cpu_number-1){
            
            MPI_Isend(vector_swap_forward,1,COLUMN,rank_number-1,tag,MPI_COMM_WORLD,&request);
            MPI_Isend(vector_swap_backward,1,COLUMN,rank_number+1,tag,MPI_COMM_WORLD,&request);
        
            MPI_Recv(read_buf_back,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&status);
            MPI_Recv(read_buf_front,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&status);
           
            for(int i=0;i<W;i++){
                next[(local_l*rank_number-1)*W+i] = read_buf_front[i];
                next[(local_l*rank_number+local_l)*W+i] = read_buf_back[i];
            }
        }else if(rank_number==(cpu_number-1)){
            
            MPI_Isend(vector_swap_forward,1,COLUMN,rank_number-1,tag,MPI_COMM_WORLD,&request);
            
            MPI_Recv(read_buf_front,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&status);
            
            for(int i=0;i<W;i++)
                next[(local_l*rank_number-1)*W+i] = read_buf_front[i];
        }
        
        MPI_Allgather(&balance,1,MPI_INT,global_balance,1,MPI_INT,MPI_COMM_WORLD); 
        int flag = 0;
        for(int i=0;i<cpu_number;i++)
            if(global_balance[i]==0)
                flag=1;

        if(!flag){
            break;
        }

        int *tmp = temp;
        temp = next;
        next = tmp;
        }
    }


    local_min = temp[local_l*rank_number*W];
        for(int i = local_l*rank_number; i < local_l*rank_number+local_l; i++){
            for(int j=0;j<W;j++){
                if(local_min>temp[i*W+j])
                    local_min = temp[i*W+j];
            }
        }   
    MPI_Allgather(&local_min,1,MPI_INT,global_min,1,MPI_INT,MPI_COMM_WORLD); 

        local_min = global_min[0];
        for(int i=0;i<cpu_number;i++)
            if(global_min[i]<local_min)
                local_min = global_min[i];
  
  
  
    if(rank_number==0)
        printf("Size: %d*%d, Iteration: %d, Min Temp: %d\n", L, W, local_count, local_min);
    
    MPI_Finalize();
    return 0;
}
/*
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#ifndef W
#define W 20                                    // Width
#endif
int main(int argc, char **argv) {
  int L = atoi(argv[1]);                        // Length
  int iteration = atoi(argv[2]);                // Iteration
  srand(atoi(argv[3]));                         // Seed
  float d = (float) random() / RAND_MAX * 0.2;  // Diffusivity
  int *temp = malloc(L*W*sizeof(int));          // Current temperature
  int *next = malloc(L*W*sizeof(int));          // Next time step

  
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < W; j++) {
      temp[i*W+j] = random()>>3;
    }
  }
  int count = 0, balance = 0;
  while (iteration--) {     // Compute with up, left, right, down points
    balance = 1;
    count++;
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < W; j++) {
        float t = temp[i*W+j] / d;
        t += temp[i*W+j] * -4;
        t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
        t += temp[(i + 1 >= L ? i : i + 1)*W+j];
        t += temp[i*W+(j - 1 <  0 ? 0 : j - 1)];
        t += temp[i*W+(j + 1 >= W ? j : j + 1)];
        t *= d;
        next[i*W+j] = t ;
        if (next[i*W+j] != temp[i*W+j]) {
          balance = 0;
        }
      }
    }
    if (balance) {
      break;
    }
    int *tmp = temp;
    temp = next;
    next = tmp;
  }
  int min = temp[0];
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < W; j++) {
      if (temp[i*W+j] < min) {
        min = temp[i*W+j];
      }
    }
  }
  printf("Size: %d*%d, Iteration: %d, Min Temp: %d\n", L, W, count, min);
  return 0;
}*/
