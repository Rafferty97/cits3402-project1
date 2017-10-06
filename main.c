#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>

#include "results.h"
#include "site.h"
#include "bond.h"

void random_percolation(char type, int size, float p, int iter, int threads, bool l)
{
  unsigned long seed = time(NULL);
  percolation_results results;
  char mode = 0;
  if (tolower(type) == 'b') mode += 1;
  if (threads > 1) mode += 2;
  float total_time = 0, avg_largest_cluster = 0;
  int perc_count = 0, n = 0;
  if (l) printf("Percolating...");
  for (int i=0; i<iter; i++) {
    switch (mode) {
      case 0:
      site_percolation(size, p, seed, &results); break;
      case 1:
      bond_percolation(size, p, seed, &results); break;
      case 2:
      site_percolation_parallel(size, p, seed, &results, threads); break;
      case 3:
      bond_percolation_parallel(size, p, seed, &results, threads); break;
    }
    if (results.percolates) perc_count++;
    total_time += results.time_taken;
    avg_largest_cluster += results.largest_cluster;
    if (l) {
      n += size * size;
      if (n < 5000) break;
      n = 0;
      printf("\33[2K\rPercolating... %i / %i iterations.", i + 1, iter);
      fflush(stdout);
    }
    seed += 1923 * i;
  }
  if (l) printf("\n\n");
  avg_largest_cluster /= iter;
  if (l) {
    printf(
      "%s percolation, %i x %i grid, p = %f\n%i iterations, %i thread%s.\n\n",
      (tolower(type) == 'b' ? "Bond" : "Site"), size, size, p, iter, threads,
      (threads == 1 ? "" : "s")
    );
    printf(
      "Total time: %f seconds\nPercolating grids: %i / %i\nAvg. largest cluster size: %f\n\n",
      total_time, perc_count, iter, avg_largest_cluster
    );
  } else {
    printf(
      "%s\t%i\t%f\t%i\t%i\t%f\t%i\t%f\n", (tolower(type) == 'b' ? "bond" : "site"),
      size, p, iter, threads, total_time, perc_count, avg_largest_cluster
    );
  }
}

int main(int argc, char *argv[])
{
  random_percolation('b', 1024, 0.5, 100, 1, false);
  random_percolation('b', 1024, 0.5, 100, 2, false);
  random_percolation('b', 1024, 0.5, 100, 4, false);
  random_percolation('s', 1024, 0.5, 100, 1, false);
  random_percolation('s', 1024, 0.5, 100, 2, false);
  random_percolation('s', 1024, 0.5, 100, 4, false);
}
