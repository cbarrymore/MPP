#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include "../include/ga.h"

#define PRINT 1

#pragma omp threadprivate(seed)

int seed;

int aleatorio(int n) {
	return (rand_r(&seed) % n);  // genera un numero aleatorio entre 0 y n-1
}

int search_element(int *array, int end, int element)
{
	int i=0;
	int found=0;
	
	// comprueba que un elemento no está incluido en el individuo (el cual no admite enteros repetidos)
	while((i < end) && ! found) {
		if(array[i] == element) {
			found = 1;
		}
		i++;
	}
        return found;
}

int find_element(int *array, int end, int element)
{
        int pos = 0;
	for(int i = 0; i < end; i++) {
             if(array[i] == element) {
                 pos = i;
                 break;
             }
        }
        return pos; // Posición del elemento encontrado
}

int *crear_individuo(int n)
{
        // El primer elemento del individuo siempre será el 0, por ejemplo.
	int i=1, value;
	int *individuo = (int *) malloc(n * sizeof(int));
	
	// inicializa array de elementos
	memset(individuo, 0, n * sizeof(int));
	
	while(i < n) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!search_element(individuo, i, value)) {
			individuo[i] = value;  // lo incluimos
			i++;
		}
	}
	return individuo;
}

int comp_fitness(const void *a, const void *b) {
	/* qsort pasa un puntero al elemento que está ordenando */
	return (*(Individuo **)a)->fitness - (*(Individuo **)b)->fitness;
}

double aplicar_ga(const double *d, int n, int n_gen, int tam_pob, float m_rate, int *sol)
{
	int i, g, mutation_start;
	int fitness_anterior = __INT32_MAX__;
	int criterio = 1;
	

	seed = getpid() + (int)time(NULL);
	// crea poblacion inicial (array de individuos)
	Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
	assert(poblacion);
	
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		poblacion[i]->array_int = crear_individuo(n);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n);
	}
	
	// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
	
	// evoluciona la poblacion durante un numero de generaciones

	for(g = 0; g < n_gen && criterio; g++)
	{
		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		for(i = 0; i < (tam_pob/2) - 1; i += 2) {
			cruzar(poblacion[i], poblacion[i+1], poblacion[tam_pob/2 + i], poblacion[tam_pob/2 + i + 1], n);
		}
		
		// por ejemplo, inicia la mutacion a partir de 1/4 de la poblacion.
                // puede haber otras opciones pero dejar el primer individuo sin modificar siempre
		mutation_start = tam_pob/4;
		
		// muta 3/4 partes de la poblacion
		for(i = mutation_start; i < tam_pob; i++) {
			mutar(poblacion[i], n, m_rate);
		}
		
		// recalcula el fitness del individuo
		for(i = 0; i	 < tam_pob; i++) {
			fitness(d, poblacion[i], n);
		}
		
		// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
		qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
		
		if (PRINT) {
			printf("Generacion %d - ", g);
			printf("Fitness = %.0lf\n", (poblacion[0]->fitness));
		}
		if(poblacion[0]->fitness <= (fitness_anterior - fitness_anterior*0.01))
			fitness_anterior = poblacion[0]->fitness;
		else
			criterio=0;

	}
	
	memmove(sol, poblacion[0]->array_int, n*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0]->fitness);
	
	// se libera la memoria reservada
	for(i = 0; i < tam_pob; i++) {
		free(poblacion[i]->array_int);
		free(poblacion[i]);
	} 
	free(poblacion); //Antes solo se liberaba el puntero poblacion, y no tambien el contenido del array
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n)
{
	// Elegir un punto (o puntos) de corte aleatorio a partir del que se realiza el intercambio de los genes. 
	
	// Entonces, por ejemplo, los primeros genes del padre1 van al hijo1, y los primeros del padre2 al hijo2.
        // Se debe evitar en cada paso la introduccion de duplicados en los hijos
	// Y los restantes genes de cada hijo son del otro padre, respectivamente.

        // Otra opción podría ser factibilizar a posteriori, despues de generar los descendientes: eliminar posibles 
        // repetidos de ambos hijos. Si encuentro algún elemento repetido en el hijo, lo cambio por otro que no este el array

	// Establecemos punto de corte aleatorio
	srand(time(NULL));
    int punto_de_corte = rand() % (n - 1) + 1; // Rango válido: 1 a n-1

    // Copiar los primeros genes del padre1 a hijo1 y del padre2 a hijo2
    for (int i = 0; i < punto_de_corte; i++) {
        hijo1->array_int[i] = padre1->array_int[i];
        hijo2->array_int[i] = padre2->array_int[i];
    }


	// Inicializar un arry auxiliar para rastrear los genes utilizados 
	int * genes_usados = (int*) malloc(n*sizeof(int));
	memset(genes_usados,0,n*sizeof(int));
	for (int j = 0; j < punto_de_corte; j++) {
			genes_usados[hijo1->array_int[j]] = 1;
    }
	// Llenar el resto de hijo1 con genes de padre2 y rastrear su uso
    for (int i = punto_de_corte; i < n; i++) {
        int gen_padre2 = padre2->array_int[i];
       
		if(!genes_usados[gen_padre2]) {
			hijo1->array_int[i] = gen_padre2;
			genes_usados[gen_padre2] =1;
		}
		else{
			for (int j = 0; j < n; j++) {
                if (!genes_usados[j]) {
                    hijo1->array_int[i] = j;
                    genes_usados[j] = 1;
                    break;
                }
			}
		}
	}

	memset(genes_usados,0,n*sizeof(int));
	// Llenar el resto de hijo2 con genes de padre1 y rastrear su uso
	for (int i = punto_de_corte; i < n; i++) {
        int gen_padre1 = padre1->array_int[i];

        // Verificar si el gen ya se utilizó en hijo1
        int gen_repetido = 0;
        for (int j = 0; j < punto_de_corte; j++) {
            if (hijo2->array_int[j] == gen_padre1) {
                gen_repetido = 1;
                break;
            }
        }
		if(!gen_repetido) {
			hijo2->array_int[i] = gen_padre1;
			genes_usados[gen_padre1] =1;
		}
		else{
			for (int j = 0; j < n; j++) {
                if (!genes_usados[j]) {
                    hijo2->array_int[i] = j;
                    genes_usados[j] = 1;
                    break;
                }
			}
		}
	}

}

void invertir(int *a, int k)
{
        int t;
	// Uno por uno invierte los elementos de a[0..k-1]
}

void mutar(Individuo *actual, int n, float m_rate)
{
		int i;
		int j;
		srand(time(NULL) + getpid());
		if ((double) rand_r(&seed) / ((double) RAND_MAX) < m_rate){ //Solo se muta cuando número aleatorio es menor que la tasa de mutación
			//Se establecen aleatoriamente las posiciones i y j
			j =rand_r(&seed) % n + 1;  // generamos un número aleatorio entre 1 y n
			i = rand_r(&seed) % (j-1) + 1; // 0 a j-2 <+1> 1 a j-1 generamos un número aleatorio entre 1 y j-1
			int swap_aux;
			//Invertimos en bucle el orden de los elementos en el rango [i,j, aumentando i y reduciendo j 
			//hasta que i == j, indicando que ya se han invertido todos los elementos.
			for (int i = i; i<j; i++){
				swap_aux = actual->array_int[i];
				actual->array_int[i] = actual->array_int[j];
				actual->array_int[j] = swap_aux;
				j--;
				
			}
		}
	// Implementación recomendada (aunque puede ser cualquier otra que se considere adecuada para este problema): 
	// Reverse Sequence Mutation (RSM), donde elegimos una secuencia S limitada por dos posiciones i, j
        // elegidas aleatoriamente con i<j, e i>0 para no modificar nunca el 1er elemento. El orden de los elementos en 
	// esta secuencia será invertido, por ejemplo con i=1, j=4: (1,2,3,4,5,6) --> (1,5,4,3,2,6). 
	// i=1 , j = 5 (1,2,3,4,5,6) --> 
	// 1,6,3,4,5,2  i=2 j=4
	// 1,6,5,4,3,2  i=3 j=3
	// 1,6,5,4,3,2 i=4, j=2
	// 1 

        // Usar la variable m_rate para establecer la intensidad (iteraciones) de la mutación, teniendo en cuenta que
	// si el valor es demasiado pequeño la convergencia es muy pequeña y si es demasiado puede diverger.
}

double distancia_ij(const double *d, int i, int j, int n)
{
	double dist = 0.0;
	
	// Devuelve la distancia entre dos elementos de la matriz 'd'
	return dist;
}

void fitness(const double *d, Individuo *individuo, int n)
{
	// Determina la calidad del individuo calculando la suma de la distancia entre cada par de ciudades consecutivas en el array
	double suma = 0;
    int d_i, d_j;
    int * ciudades = individuo->array_int;
	omp_set_num_threads(8);
    #pragma omp parallel private(d_i, d_j) shared(suma	)
    {
        #pragma omp for
        for (int i = 0; i < n; i++) {
			//Obtenemos los pares de ciudades consecutivas y accedemos a la posición del array que contiene su distancia,
			//sumandola a la variable suma.
            d_i = ciudades[i - 1];
            d_j = ciudades[i];
			#pragma omp critical
			{
				suma += d[d_i*n + d_j];
			}
        }

        
        	
    }
	//Como por simplicidad en la codificación no se incluye el nodo final por coincidir con el origen,
	// sumamos fuera del bucle la distancia entre la última ciudad y la primera, la cual es el origen.
	d_i = ciudades[n - 1];
	d_j = ciudades[0];
	suma += d[d_i * n + d_j];
    
	//establecemos el valor de fitness del individuo al de suma.
    individuo->fitness = suma;
}
