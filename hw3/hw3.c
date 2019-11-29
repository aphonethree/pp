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
  int rank_number,cpu_number,int_buff;
  int *vector_swap_forward = (int*)malloc(1*W*sizeof(int));
  int *vector_swap_backward = (int*)malloc(1*W*sizeof(int));
  int *read_buf_front = (int*)malloc(W*sizeof(int));
  int *read_buf_back = (int*)malloc(W*sizeof(int));
  int local_l,tag=0,min=0,read_buff_min,finish_number=0,read_buff_balance;
  MPI_Request request;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_number);
  MPI_Comm_size(MPI_COMM_WORLD, &cpu_number);
  
  local_l = L/cpu_number;
   
  for (int i = local_l*rank_number; i < local_l*rank_number+local_l; i++) {
    for (int j = 0; j < W; j++) {
      temp[i*W+j] = random()>>3;
    }
  }
  
  int local_count = 0, balance = 0, local_min=temp[0];
  while (iteration>local_count) {       // Compute with up, left, right, down points
    balance = 1;
    local_count++;

    if(cpu_number==1){
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
                    if (next[i*W+j] != temp[i*W+j]) 
                        balance = 0;
                    if(t<local_min)
                        local_min = t;
                    if(rank_number==0){
                        vector_swap_backward[j]= next[i*W+j];
                    }else if(rank_number>0 && rank_number<cpu_number){
                        if(i==local_l*rank_number)
                            vector_swap_forward[j]= next[i*W+j]; 
                        else if(i==(local_l*rank_number+local_l-1))
                            vector_swap_backward[j]= next[i*W+j]; 
                    }else if(rank_number==cpu_number){
                            vector_swap_forward[j]= next[i*W+j]; 
                    }          
            }
        }
         if (balance) {
                if(rank_number>0){
                     MPI_Isend(&balance,1,MPI_INT,0,tag,MPI_COMM_WORLD,&request);
                    break;
                }
            }
        if(rank_number==0){
            MPI_Isend(vector_swap_backward,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&request);
            MPI_Irecv(read_buf_back,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&request);
            for(int i=0;i<W;i++)
                next[(local_l*rank_number+local_l)*W+i] = read_buf_back[i];
            for(int i=1;i<cpu_number;i++){
                MPI_Irecv(&read_buff_min,1,MPI_INT,i,tag,MPI_COMM_WORLD,&request);
                MPI_Irecv(&read_buff_balance,1,MPI_INT,i,tag,MPI_COMM_WORLD,&request);
                
                if(read_buff_min<local_min)
                    local_min = read_buff_min;
                if(read_buff_balance==1)
                    finish_number++;
                if(finish_number==cpu_number && balance)
                    break;
                 
            }
           MPI_Wait(&request,&status);
        }else if(rank_number>0 && rank_number<cpu_number){
            MPI_Isend(vector_swap_forward,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&request);
            MPI_Isend(vector_swap_backward,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&request);
            MPI_Isend(&local_min,1,MPI_INT,0,tag,MPI_COMM_WORLD,&request);
            MPI_Irecv(read_buf_back,W,MPI_INT,rank_number+1,tag,MPI_COMM_WORLD,&request);
            MPI_Irecv(read_buf_front,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&request);
           
            for(int i=0;i<W;i++){
                next[(local_l*rank_number-1)*W+i] = read_buf_front[i];
                next[(local_l*rank_number+local_l)*W+i] = read_buf_back[i];
            }
        }else if(rank_number==cpu_number){
            MPI_Isend(vector_swap_forward,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&request);
            MPI_Isend(&local_min,1,MPI_INT,0,tag,MPI_COMM_WORLD,&request);
            MPI_Irecv(read_buf_front,W,MPI_INT,rank_number-1,tag,MPI_COMM_WORLD,&request);
            
             for(int i=0;i<W;i++)
                 next[(local_l*rank_number-1)*W+i] = read_buf_front[i];
        }
        int *tmp = temp;
        temp = next;
        next = tmp;
    }
  }

  
 
  
  
  MPI_Finalize();
  if(rank_number==0)
    printf("Size: %d*%d, Iteration: %d, Min Temp: %d\n", L, W, local_count, local_min);
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