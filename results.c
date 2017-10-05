#include "stdio.h"

#include "results.h"

void print_results(percolation_results results)
{
  printf(
    "Type: %s\nSize: %i x %i\nSeed probability: %f\nPercolates: %s\nLargest cluster: %i\nThreads: %i\nTime taken: %fs\n\n",
    results.type,
    results.lattice_size,
    results.lattice_size,
    results.p,
    results.percolates ? "Yes" : "No",
    results.largest_cluster,
    results.threads,
    results.time_taken
  );
}

void print_results_batch(percolation_batch_results results)
{
  printf(
    "Type: %s\nSize: %i x %i\nSeed probability: %f\nTrials: %i\nPercolates: %i / %i\nAvg Largest cluster: %f\nThreads: %i\nAvg Time taken: %fs\n\n",
    results.type,
    results.lattice_size,
    results.lattice_size,
    results.p,
    results.trials,
    results.num_percolates,
    results.trials,
    results.avg_largest_cluster,
    results.threads,
    results.avg_time_taken
  );
}
