#ifndef __RESULTS_H__
#define __RESULTS_H__

#include <stdbool.h>

typedef struct {
  char *type;
  int lattice_size;
  float p;
  int threads;
  bool percolates;
  int largest_cluster;
  float time_taken;
} percolation_results;

#endif