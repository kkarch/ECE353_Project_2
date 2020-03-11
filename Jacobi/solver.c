/* Code for the Jacobi equation solver. 
 * Author: Naga Kandasamy
 * Date created: April 19, 2019
 * Date modified: February 21, 2020
 *
 * Compile as follows:
 * gcc -o solver solver.c solver_gold.c -O3 -Wall -std=gnu99 -lm -lpthread
 *
 * If you wish to see debug info, add the -D DEBUG option when compiling the code.
 */



#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "grid.h" 

/* Shared data structure used by the threads */
typedef struct args_for_thread_t {
    int tid;                          /* The thread ID */
    int num_threads;                  /* Number of worker threads */
    int offset;                       /* Starting offset for thread within the vectors */
    int chunk_size;                   /* Chunk size */
    double *diff;                      /* Location of the shared variable sum */
    int *numels;
    int *ret_numiter;
    pthread_mutex_t *mutex_for_diff;   /* Location of the lock variable protecting sum */
    grid_t *ref_grid;                      /* Location of the grid */
    grid_t *update_grid                      /* Location of the grid */
} ARGS_FOR_THREAD;

pthread_barrier_t barrier; 




extern int compute_gold (grid_t *);
int compute_using_pthreads_jacobi (grid_t *, int);
void *Jacobi(void *);
void compute_grid_differences(grid_t *, grid_t *);
grid_t *create_grid (int, float, float);
grid_t *copy_grid (grid_t *);
void print_grid (grid_t *);
void print_stats (grid_t *);
double grid_mse (grid_t *, grid_t *);



int 
main (int argc, char **argv)
{	
	if (argc < 5) {
        printf ("Usage: %s grid-dimension num-threads min-temp max-temp\n", argv[0]);
        printf ("grid-dimension: The dimension of the grid\n");
        printf ("num-threads: Number of threads\n"); 
        printf ("min-temp, max-temp: Heat applied to the north side of the plate is uniformly distributed between min-temp and max-temp\n");
        exit (EXIT_FAILURE);
    }
    
    /* Parse command-line arguments. */
    int dim = atoi (argv[1]);
    int num_threads = atoi (argv[2]);
    float min_temp = atof (argv[3]);
    float max_temp = atof (argv[4]);
    
    /* Generate the grids and populate them with initial conditions. */
 	grid_t *grid_1 = create_grid (dim, min_temp, max_temp);
    /* Grid 2 should have the same initial conditions as Grid 1. */
    grid_t *grid_2 = copy_grid (grid_1); 

	/* Compute the reference solution using the single-threaded version. */
	printf ("\nUsing the single threaded version to solve the grid\n");
	int num_iter = compute_gold (grid_1);
	printf ("Convergence achieved after %d iterations\n", num_iter);
    /* Print key statistics for the converged values. */
	printf ("Printing statistics for the interior grid points\n");
    print_stats (grid_1);
#ifdef DEBUG
    print_grid (grid_1);
#endif
	
	/* Use pthreads to solve the equation using the jacobi method. */
	printf ("\nUsing pthreads to solve the grid using the jacobi method\n");
	num_iter = compute_using_pthreads_jacobi (grid_2, num_threads);
	printf ("Convergence achieved after %d iterations\n", num_iter);			
    printf ("Printing statistics for the interior grid points\n");
	print_stats (grid_2);
#ifdef DEBUG
    print_grid (grid_2);
#endif
    
    /* Compute grid differences. */
    double mse = grid_mse (grid_1, grid_2);
    printf ("MSE between the two grids: %f\n", mse);

	/* Free up the grid data structures. */
	free ((void *) grid_1->element);	
	free ((void *) grid_1); 
	free ((void *) grid_2->element);	
	free ((void *) grid_2);

	exit (EXIT_SUCCESS);
}

/* FIXME: Edit this function to use the jacobi method of solving the equation. The final result should be placed in the grid data structure. */
int 
compute_using_pthreads_jacobi (grid_t *grid, int num_threads)
{	
    int bar_ret;
    grid_t *grid_3 = copy_grid (grid);


    pthread_t *tid = (pthread_t *) malloc (sizeof (pthread_t) * num_threads); /* Data structure to store the thread IDs */
    if (tid == NULL) {
        perror ("malloc");
        exit (EXIT_FAILURE);
    }

    pthread_attr_t attributes;                  /* Thread attributes */
    pthread_mutex_t mutex_for_diff;              /* Lock for the shared variable sum */

    bar_ret = pthread_barrier_init(&barrier, NULL, num_threads);
    if (bar_ret != 0){
        printf("Barrier Init Failure\n");
        exit(EXIT_FAILURE);
        }

    
    pthread_attr_init (&attributes);            /* Initialize the thread attributes to the default values */
    pthread_mutex_init (&mutex_for_diff, NULL);  /* Initialize the mutex */


    /* Allocate memory on the heap for the required data structures and create the worker threads */
    int i;
    double diff = 0;
    int numels = 0;
    int ret_numiter = 0; 
    ARGS_FOR_THREAD **args_for_thread;
    args_for_thread = malloc (sizeof (ARGS_FOR_THREAD) * num_threads);
    for (i = 0; i < num_threads; i++){
        args_for_thread[i] = (ARGS_FOR_THREAD *) malloc (sizeof (ARGS_FOR_THREAD));
        args_for_thread[i]->tid = i; 
        args_for_thread[i]->num_threads = num_threads;
        //args_for_thread[i]->offset = i * chunk_size; 
        //args_for_thread[i]->chunk_size = chunk_size; 
        args_for_thread[i]->diff = &diff;
        args_for_thread[i]->numels = &numels;
        args_for_thread[i]->ret_numiter = &ret_numiter;
        args_for_thread[i]->mutex_for_diff = &mutex_for_diff;
        args_for_thread[i]->ref_grid = grid;
        args_for_thread[i]->update_grid=grid_3;
    }

    for (i = 0; i < num_threads; i++)
        pthread_create (&tid[i], &attributes, Jacobi, (void *) args_for_thread[i]);
					 
    /* Wait for the workers to finish */
    for(i = 0; i < num_threads-1; i++)
        pthread_join (tid[i], NULL);
        //printf("Joined\n");
  
    pthread_barrier_destroy(&barrier);

    /* Free data structures */
    for(i = 0; i < num_threads; i++)
        free ((void *) args_for_thread[i]);

    return ret_numiter;
}


void *
Jacobi(void *args)
{
    int debug = 0;
    ARGS_FOR_THREAD *targs = (ARGS_FOR_THREAD *) args;

    int num_iter = 0;
	int done = 0;
    int i, j;
	float old, new; 
    float eps = 1e-2; /* Convergence criteria. */
    int num_elements;
    double t_diff;
    double avg; 

    grid_t *ref_grid = targs->ref_grid;
    grid_t *update_grid = targs->update_grid;
    grid_t *t_grid;
    
    if (targs->tid==1 && debug == 1)
    {
        printf("ref_grid %d\n",ref_grid->dim);
        printf("update_grid %d\n",update_grid->dim);

        for ( int i = 0; i < 5; i++)
        {
            printf("Thread %d at iteration: %d\n",targs->tid,i);
            pthread_mutex_lock(targs->mutex_for_diff);
            *(targs->diff) += 1;
            pthread_mutex_unlock(targs->mutex_for_diff);

            pthread_barrier_wait(&barrier);
        }
    }
    
	while(!done) { /* While we have not converged yet. */
        t_diff = 0.0;
        num_elements = 0;
        //printf("Thread %d at iteration: %d\n",targs->tid,num_iter);

        for (i = targs->tid+1; i < (ref_grid->dim - 1); i+=targs->num_threads) {
            for (j = 1; j < (ref_grid->dim - 1); j++) {
                old = ref_grid->element[i * ref_grid->dim + j]; /* Store old value of grid point. */
                /* Apply the update rule. */	
                new = 0.25 * (ref_grid->element[(i - 1) * ref_grid->dim + j] +\
                              ref_grid->element[(i + 1) * ref_grid->dim + j] +\
                              ref_grid->element[i * ref_grid->dim + (j + 1)] +\
                              ref_grid->element[i * ref_grid->dim + (j - 1)]);

                update_grid->element[i * update_grid->dim + j] = new; /* Update the grid-point value. */
                t_diff = t_diff + fabs(new - old); /* Calculate the difference in values. */
                num_elements++;
            }
        }
        pthread_mutex_lock(targs->mutex_for_diff);
        *(targs->diff) += t_diff;
        *(targs->numels) += num_elements;
        pthread_mutex_unlock(targs->mutex_for_diff);
        //printf ("Iteration %d. DIFF: %f.\n", num_iter, diff);
        
        pthread_barrier_wait(&barrier);
        num_iter++;

        avg = *(targs->diff)/ *(targs->numels);

        pthread_barrier_wait(&barrier);

        if (targs->tid == 0)
        {
            printf ("Iteration %d. DIFF: %f.\n", num_iter, avg);
            *(targs->ret_numiter) = num_iter;
            *(targs->diff) = 0;
            *(targs->numels) = 0;
        }
        
            
        if (avg < eps) 
            done = 1;
        
        t_grid = ref_grid;
        ref_grid = update_grid;
        update_grid = t_grid;
        
    }
    
    //printf("Thread %d exiting!\n",targs->tid);
    pthread_exit ((void *)0);
}

/* Create a grid with the specified initial conditions. */
grid_t * 
create_grid (int dim, float min, float max)
{
    grid_t *grid = (grid_t *) malloc (sizeof (grid_t));
    if (grid == NULL)
        return NULL;

    grid->dim = dim;
	printf("Creating a grid of dimension %d x %d\n", grid->dim, grid->dim);
	grid->element = (float *) malloc (sizeof (float) * grid->dim * grid->dim);
    if (grid->element == NULL)
        return NULL;

    int i, j;
	for (i = 0; i < grid->dim; i++) {
		for (j = 0; j < grid->dim; j++) {
            grid->element[i * grid->dim + j] = 0.0; 			
		}
    }

    /* Initialize the north side, that is row 0, with temperature values. */ 
    srand ((unsigned) time (NULL));
	float val;		
    for (j = 1; j < (grid->dim - 1); j++) {
        val =  min + (max - min) * rand ()/(float)RAND_MAX;
        grid->element[j] = val; 	
    }

    return grid;
}

/* Creates a new grid and copies over the contents of an existing grid into it. */
grid_t *
copy_grid (grid_t *grid) 
{
    grid_t *new_grid = (grid_t *) malloc (sizeof (grid_t));
    if (new_grid == NULL)
        return NULL;

    new_grid->dim = grid->dim;
	new_grid->element = (float *) malloc (sizeof (float) * new_grid->dim * new_grid->dim);
    if (new_grid->element == NULL)
        return NULL;

    int i, j;
	for (i = 0; i < new_grid->dim; i++) {
		for (j = 0; j < new_grid->dim; j++) {
            new_grid->element[i * new_grid->dim + j] = grid->element[i * new_grid->dim + j] ; 			
		}
    }

    return new_grid;
}

/* This function prints the grid on the screen. */
void 
print_grid (grid_t *grid)
{
    int i, j;
    for (i = 0; i < grid->dim; i++) {
        for (j = 0; j < grid->dim; j++) {
            printf ("%f\t", grid->element[i * grid->dim + j]);
        }
        printf ("\n");
    }
    printf ("\n");
}


/* Print out statistics for the converged values of the interior grid points, including min, max, and average. */
void 
print_stats (grid_t *grid)
{
    float min = INFINITY;
    float max = 0.0;
    double sum = 0.0;
    int num_elem = 0;
    int i, j;

    for (i = 1; i < (grid->dim - 1); i++) {
        for (j = 1; j < (grid->dim - 1); j++) {
            sum += grid->element[i * grid->dim + j];

            if (grid->element[i * grid->dim + j] > max) 
                max = grid->element[i * grid->dim + j];

             if(grid->element[i * grid->dim + j] < min) 
                min = grid->element[i * grid->dim + j];
             
             num_elem++;
        }
    }
                    
    printf("AVG: %f\n", sum/num_elem);
	printf("MIN: %f\n", min);
	printf("MAX: %f\n", max);
	printf("\n");
}

/* Calculate the mean squared error between elements of two grids. */
double
grid_mse (grid_t *grid_1, grid_t *grid_2)
{
    double mse = 0.0;
    int num_elem = grid_1->dim * grid_1->dim;
    int i;

    for (i = 0; i < num_elem; i++) 
        mse += (grid_1->element[i] - grid_2->element[i]) * (grid_1->element[i] - grid_2->element[i]);
                   
    return mse/num_elem; 
}



		

