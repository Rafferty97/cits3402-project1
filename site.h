#ifndef __SITE_H__
#define __SITE_H__

#include "results.h"

struct site_grid {
  bool *site;          // 2D array specifying if a lattice point (site) is occupied
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

void site_percolation(int size, float p, unsigned long seed, percolation_results *results);

void site_percolation_parallel(int size, float p, unsigned long seed, percolation_results *results, int threads);

void load_site_grid(bool *sites);

#endif