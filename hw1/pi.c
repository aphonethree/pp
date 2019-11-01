#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
pthread_mutex_t mutex;
unsigned long long number_in_circle;
void *local_count(void *local_toss)
{
    unsigned long long times,local_sum=0;
    //unsigned long long toss_times = (unsigned long long) local_toss;
    unsigned int random_seed = time(NULL);
    double x,y;
    for (times = 0; times < (unsigned long long) local_toss; times++) {
        x = rand_r(&random_seed)/(RAND_MAX+1.0);
	y = rand_r(&random_seed)/(RAND_MAX+1.0);
        //x = drand48();
        //y = drand48();
        if (x*x + y*y <= 1)
            local_sum++;
    }
    pthread_mutex_lock(&mutex);
    number_in_circle += local_sum;
    pthread_mutex_unlock(&mutex);
    return NULL;
}
int main(int argc, char **argv)
{
    long i;
    double pi_estimate, distance_squared, x, y;
    unsigned long long  number_of_cpu, number_of_tosses, toss;
    pthread_t *thread_array;

    if ( argc < 2) {
        exit(-1);
    }
    number_of_cpu = atoi(argv[1]);
    number_of_tosses = atoi(argv[2]);
    toss = number_of_tosses / number_of_cpu;
    thread_array = (pthread_t *)malloc(sizeof(pthread_t)*number_of_cpu);
    if (( number_of_cpu < 1) || ( number_of_tosses < 0)) {
        exit(-1);
    }

    number_in_circle = 0;
    pthread_mutex_init(&mutex, NULL);
    for(i=0;i<number_of_cpu;i++)
        pthread_create(&thread_array[i], NULL,local_count, (void*)toss);
    for(i=0;i<number_of_cpu;i++)
        pthread_join(thread_array[i],NULL);

    printf("%f\n",4*number_in_circle/((double) number_of_tosses));
    pthread_mutex_destroy(&mutex);
    free(thread_array);
    return 0;
}
