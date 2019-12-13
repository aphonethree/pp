/**********************************************************************
 * DESCRIPTION:
 *   Serial Concurrent Wave Equation - C Version
 *   This program implements the concurrent wave equation
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXPOINTS 1000000
#define MAXSTEPS 1000000
#define MINPOINTS 20
#define PI 3.14159265

void check_param(void);
void init_line(void);
void update (void);
void printfinal (void);




/**********************************************************************
 *	Checks input values from parameters
 *********************************************************************/
__host__ void check_param(int *tpoints,int *nsteps)
{
   char tchar[20];

   /* check number of points, number of iterations */
   while ((*tpoints < MINPOINTS) || (*tpoints > MAXPOINTS)) {
      printf("Enter number of points along vibrating string [%d-%d]: "
           ,MINPOINTS, MAXPOINTS);
      scanf("%s", tchar);
      *tpoints = atoi(tchar);
      if ((*tpoints < MINPOINTS) || (*tpoints > MAXPOINTS))
         printf("Invalid. Please enter value between %d and %d\n", 
                 MINPOINTS, MAXPOINTS);
   }
   while ((*nsteps < 1) || (*nsteps > MAXSTEPS)) {
      printf("Enter number of time steps [1-%d]: ", MAXSTEPS);
      scanf("%s", tchar);
      *nsteps = atoi(tchar);
      if ((*nsteps < 1) || (*nsteps > MAXSTEPS))
         printf("Invalid. Please enter value between 1 and %d\n", MAXSTEPS);
   }

   printf("Using points = %d, steps = %d\n", *tpoints, *nsteps);

}

/**********************************************************************
 *      Calculate new values using wave equation
 *********************************************************************/

__device__ float do_math(float currentvalue,float oldval)
{
   float dtime, c, dx, tau, sqtau;
  
   dtime = 0.3;
   c = 1.0;
   dx = 1.0;
   tau = (c * dtime / dx);
   sqtau = tau * tau;
   return ((2.0 * currentvalue) - oldval + (sqtau *  (-2.0)*currentvalue));

}



/**********************************************************************
 *     Initialize points on line
 *********************************************************************/
__global__ void init_line(float *values,int tpoint,int nsteps)
{
   float x, fac;
   float currentval,oldval,newval;
   int indx = blockIdx.x * blockDim.x + threadIdx.x+1;

   
   /* Calculate initial values based on sine curve */
   fac = 2.0 * PI;
   
   x = (indx -1.0)/(float)(tpoint-1.0);
   currentval = sin(fac * x);
   oldval= currentval;
   
   /*
        k = 0.0; 
        tmp = tpoints - 1;
        for (j = 1; j <= tpoints; j++) {
            x = k/tmp;
            values[j] = sin (fac * x);
            k = k + 1.0;
        } 
   */

    /* Initialize old values array */
    /*
        for (i = 1; i <= tpoints; i++) 
            oldval[i] = values[i];
    */
    /**********************************************************************
    *     Update all values along line a specified number of times
    *********************************************************************/
    int i;
    #pragma unroll 1024
    for ( i = 1; i<=nsteps; i++) {
        if (indx ==0 || indx == tpoint)
            currentval =0.0;
        else{
            newval =  do_math(currentval,oldval);
            oldval = currentval;
            currentval = newval;
        }
    }
    values[indx] = currentval;


}





/**********************************************************************
 *     Print final results
 *********************************************************************/
__host__ void printfinal(float values[],int tpoints)
{
   int i;
  
   for (i = 1; i <= tpoints; i++) {
      printf("%6.4f ", values[i]);
      if (i%10 == 0)
         printf("\n");
   }
}

/**********************************************************************
 *	Main program
 *********************************************************************/
__host__ int main(int argc, char *argv[])
{
    float *values,*final_result;

    int nsteps,                 	/* number of time steps */
        tpoints;    	     		/* total points along string */ 
    int blocknumber,datasize;
   
	sscanf(argv[1],"%d",&tpoints);
	sscanf(argv[2],"%d",&nsteps);
    blocknumber = tpoints/1024+1;
    datasize = (blocknumber * 1024 + 1) * sizeof(float);
	check_param(&tpoints,&nsteps);
    cudaMalloc( (void**)&values,  datasize);
    final_result =(float *) malloc(datasize);
	
    printf("Initializing points on the line...\n");
    printf("Updating all points for all time steps...\n");
	init_line<<<blocknumber,1024>>>(values,tpoints,nsteps);


	printf("Printing final results...\n");
    cudaMemcpy( final_result, values, datasize, cudaMemcpyDeviceToHost ); 
	printfinal(final_result,tpoints);
	printf("\nDone.\n\n");
	cudaFree(values);
	return 0;
}