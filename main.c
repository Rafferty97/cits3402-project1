#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "site.h"

int main(int argc, char *argv[])
{
  float time = 0;
  for (int i = 0; i < 1000; i++) {
    percolation_results results;
    site_percolation(256, 0.6, atoi(argv[1]), &results);
    time += results.time_taken;
    if ((i - 1) % 100 == 0) printf("%i iterations...\n", i - 1);
  }
  printf("256 x 256, 1000 iterations\n%i threads: %f seconds", atoi(argv[1]), time);
}
