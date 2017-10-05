#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "site.h"

int main(int argc, char *argv[])
{
  int size = 1024;
  int iter = 100;
  float dfs = 0, time = 0;
  for (int i = 0; i < iter; i++) {
    percolation_results results;
    site_percolation(size, 0.6, atoi(argv[1]), &results);
    dfs += results.dfs_time;
    time += results.time_taken;
    if ((i + 1) % 10 == 0) printf("%i iterations...\n", i + 1);
  }
  printf("%i x %i, %i iterations\n%i threads: %f seconds\n", size, size, iter, atoi(argv[1]), time);
  printf("DFS time: %f seconds\n\n", dfs);
}
