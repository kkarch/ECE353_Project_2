/*  Purpose: Calculate definite integral using trapezoidal rule.
 *
 * Input:   a, b, n, num_threads
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids, with num_threads.
 *
 * Compile: gcc -o trap trap.c -O3 -std=c99 -Wall -lpthread -lm
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 * Author: Naga Kandasamy
 * Date modified: February 21, 2020
 *
 */

#define _REENTRANT /* Make sure the library functions are MT (muti-thread) safe */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* Function prototype for the thread routines */
void *my_func (void *);

/* Structure used to pass arguments to the worker threads */
typedef struct args_for_thread_t {
    int num_threads;        /* Number of worker threads */
    int tid;                /* Thread ID */
    int arg1;               /* First argument */
    float arg2;             /* Second argument */
    double *partial_integral;    /* Third argument */
} ARGS_FOR_THREAD; 

double compute_using_pthreads (float, float, int, float, int);
double compute_gold (float, float, int, float);

int 
main (int argc, char **argv) 
{
    if (argc < 5) {
        printf ("Usage: %s lower-limit upper-limit num-trapezoids num-threads\n", argv[0]);
        printf ("lower-limit: The lower limit for the integral\n");
        printf ("upper-limit: The upper limit for the integral\n");
        printf ("num-trapezoids: Number of trapeziods used to approximate the area under the curve\n");
        printf ("num-threads: Number of threads to use in the calculation\n");
        exit (EXIT_FAILURE);
    }

  float a = atoi (argv[1]); /* Lower limit */
	float b = atof (argv[2]); /* Upper limit */
	float n = atof (argv[3]); /* Number of trapezoids */

	float h = (b - a)/(float) n; /* Base of each trapezoid */  
	printf ("The base of the trapezoid is %f\n", h);

    struct timeval start, stop;	
	gettimeofday (&start, NULL);
	double reference = compute_gold (a, b, n, h);
    gettimeofday (&stop, NULL);
    printf ("Reference solution computed using single-threaded version = %f\n", reference);
    printf ("Execution time = %fs\n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
	printf ("\n");
	/* Write this function to complete the trapezoidal rule using pthreads. */
    int num_threads = atoi (argv[4]); /* Number of threads */
    gettimeofday (&start, NULL);
	double pthread_result = compute_using_pthreads (a, b, n, h, num_threads);
	gettimeofday (&stop, NULL);
    printf ("Solution computed using %d threads = %f\n", num_threads, pthread_result);
    printf ("Execution time = %fs\n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
	printf ("\n");
    exit (EXIT_SUCCESS);
} 

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Defines the integrand
 * Input args:  x
 * Output: sqrt((1 + x^2)/(1 + x^4))

 */
float 
f (float x) 
{
    return sqrt ((1 + x*x)/(1 + x*x*x*x));
}

/*------------------------------------------------------------------
 * Function:    compute_gold
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids using a single-threaded version
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
double 
compute_gold (float a, float b, int n, float h) 
{
   double integral;
   int k;

   integral = (f(a) + f(b))/2.0;

   for (k = 1; k <= n-1; k++){
     integral += f(a+k*h);
     //printf("Single thread computed f(%f)\n",a+k*h);
     }
   
   integral = integral*h;

   return integral;
}  

/* FIXME: Complete this function to perform the trapezoidal rule using pthreads. */
double 
compute_using_pthreads (float a, float b, int n, float h, int num_threads)
{
	double integral = (f(a) + f(b))/2.0;
    //printf("Thread 0 is done with partial integral %f\n", integral);
 
    //added from vector_dot_product_v1
    pthread_t *thread_id;                       /* Data structure to store the thread IDs */
    pthread_attr_t attributes;                  /* Thread attributes */
    pthread_attr_init(&attributes);             /* Initialize thread attributes to default values */
    int i;
    int trap_count=1;
    int close_trap=1;

    /* Allocate memory on heap for data structures and create worker threads */
    thread_id = (pthread_t *) malloc (sizeof (pthread_t) * n);
    double *partial_integral = (double *) malloc (sizeof (double) * n);
    ARGS_FOR_THREAD *args_for_thread = (ARGS_FOR_THREAD *) malloc (sizeof (ARGS_FOR_THREAD) * n);

    for (i = 1; i <= n-1; i++) {
        args_for_thread[i].num_threads = num_threads;
        args_for_thread[i].tid = i; 
        args_for_thread[i].arg1 = a; 
        args_for_thread[i].arg2 = h; 
        args_for_thread[i].partial_integral = partial_integral; 
    }

    while (trap_count<=n-1){
        //for (i = 1; i <= num_threads; i++){//only run 4 threads at a time
        while(trap_count-close_trap<=4){
            if (trap_count<=n-1){
                pthread_create (&thread_id[trap_count], &attributes, my_func, (void *) &args_for_thread[trap_count]);
                //printf("Thread #%d created.\n",trap_count);
                trap_count++;
                }
        }
                        
        /* Wait for workers to finish */
        //for (i = 1; i <= num_threads; i++){
            if (close_trap<=n-1){
                pthread_join (thread_id[close_trap], NULL);
                //printf("Current trap count is %d\n", close_trap);
                close_trap++;
            }
       // }
            
        /* Accumulate partial integrals computed by worker threads */ 
    }
    for (i = 1; i <= n-1; i++){
            integral += partial_integral[i];
        }
        //printf("Group of %d threads complete.\n", num_threads);

    integral=integral*h;
    //printf("Thread %d is done with partial integral %f\n", n, integral);
    return integral;
}

/* Function that will be executed by all the worker threads */
void *
my_func (void *this_arg)
{
    ARGS_FOR_THREAD *args_for_me = (ARGS_FOR_THREAD *) this_arg; /* Typecast argument passed to function to appropriate type */

    /* Compute partial sum that this thread is responsible for */
    double partial_integral = 0.0; 
	  partial_integral=f(args_for_me->arg1+(args_for_me->tid*args_for_me->arg2));

    /* Store partial integrals into the partial_sum array */
    args_for_me->partial_integral[args_for_me->tid] = (float) partial_integral;
    //printf ("Thread %d is done with partial integral= %f using f(%f)\n", args_for_me->tid, partial_integral, (args_for_me->arg1+(args_for_me->tid*args_for_me->arg2)));
    pthread_exit ((void *) 0);
}



