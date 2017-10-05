#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "site.h"

int main(int argc, char *argv[])
{
  int size = 1024;
  int iter = 1000;
  float dfs = 0, time = 0;
  for (int i = 0; i < iter; i++) {
    percolation_results results;
    site_percolation(size, 0.6, atoi(argv[1]), &results);
    dfs += results.dfs_time;
    time += results.time_taken;
    // if ((i + 1) % 10 == 0) printf("%i iterations...\n", i + 1);
  }
  if (argc > 2) printf("Type\tSize\tIter's\tThreads\tTime\n");
  printf("SITE\t%i\t%i\t%i\t%f\n", size, iter, atoi(argv[1]), time);
}
