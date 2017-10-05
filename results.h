#include "stdbool.h"

typedef struct {
  char *type;
  int lattice_size;
  float p;
  int threads;
  bool percolates;
  int largest_cluster;
  float time_taken;
} percolation_results;

typedef struct {
  char *type;
  int lattice_size;
  float p;
  int threads;
  int trials;
  int num_percolates;
  float avg_largest_cluster;
  float avg_time_taken;
} percolation_batch_results;

void print_results(percolation_results results);
void print_results_batch(percolation_batch_results results);
