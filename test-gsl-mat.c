#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spmatrix.h>
#define DIM1 (300000)
#define DIM2 (300000)

int main() {
  gsl_spmatrix* A = gsl_spmatrix_alloc(DIM1, DIM2);
  size_t i, j;

  for( i = 0; i < 10000000; i++ ) {
    gsl_spmatrix_set(A, rand()%DIM1, rand()%DIM2, rand());
  }
  printf( "Done\n" );
}
