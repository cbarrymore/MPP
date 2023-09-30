#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "../include/io.h"

extern double aplicar_ga(const double *, int, int, int, int *);

static double mseconds() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec*1000 + t.tv_usec/1000;
}

int main(int argc, char **argv)
{
//	Check Number of Input Args
	if(argc < 3) {
		fprintf(stderr,"Ayuda:\n"); 
		fprintf(stderr,"  ./programa n nGen tamPob\n");
		return(EXIT_FAILURE);
	}
	
	int n = atoi(argv[1]);
	int n_gen = atoi(argv[2]);
	int tam_pob = atoi(argv[3]);
	
//	Generate matrix D with distance values among elements
	double *d = generar_matriz_distancias(n);
	
	#ifdef DEBUG
//		print_matrix(d, n);
	#endif
	
//	Allocate memory for output data
	int *sol = (int *) malloc(n * sizeof(int));
	
	#ifdef TIME
		double ti = mseconds();
	#endif
	
//	Call Metaheuristic
	double value = aplicar_ga(d, n, n_gen, tam_pob, sol);
	
	#ifdef TIME
		double tf = mseconds();
		printf("Execution Time: %.2lf sec\n", (tf - ti)/1000);
	#endif
	
        // Free Allocated Memory	
	

	#ifdef DEBUG
		print_solution(n, sol, value);
	#endif
	free(sol);
	free(d);
	return(EXIT_SUCCESS);
}
