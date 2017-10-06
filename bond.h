#ifndef __BOND_H__
#define __BOND_H__

#include "results.h"

struct bond_grid {
  short *bond;         // 2D array specifying if a lattice point has a bond to right and down
  int *cluster;        // 2D array specifying which cluster a site belongs to (0 if none or unknown)
  int *cluster_size;   // The number of sites in each cluster (index 0 = cluster 1)
  int cluster_cap;     // The current capacity of the cluster_size array
  int sx;              // Horizontal size of the grid
  int sy;              // Vertical size of the grid
  int c_c;             // The first valid cluster ID (prevents overlap between subgrids)
  int c_n;             // The gap between each valid cluster ID (prevents overlap between subgrids)
  int max_cluster;     // The largest cluster ID
  int largest_cluster; // The size of the largest cluster in the grid
};

void bond_percolation(int size, float p, unsigned long seed, percolation_results *results);

void bond_percolation_parallel(int size, float p, unsigned long seed, percolation_results *results, int threads);

#endif