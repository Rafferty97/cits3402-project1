#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "site.h"

int main(int argc, char *argv[])
{
  int size = 256;
  int iter = 250;
  int threads = 4;
  float p = atof(argv[1]);
  float dfs = 0, tt = 0;
  int perc = 0;
  for (int i = 0; i < iter; i++) {
    percolation_results results;
    site_percolation(size, p, threads, &results, time(NULL) + (i * 321));
    dfs += results.dfs_time;
    tt += results.time_taken;
    if (results.percolates) perc++;
    // if ((i + 1) % 10 == 0) printf("%i iterations...\n", i + 1);
    // printf("%s\t%i\n", results.percolates ? "Yes" : "No", results.largest_cluster);
  }
  if (argc > 2) printf("Type\tSize\tp\tIter's\tThreads\tPerc\tTime\n");
  printf("SITE\t%i\t%f\t%i\t%i\t%i\t%f\n", size, p, iter, threads, perc * 4, tt);
}
