#ifndef __CLUSTERS_H__
#define __CLUSTERS_H__

#include <stdbool.h>

#include "site.h"
#include "bond.h"

typedef struct cluster cluster;

struct cluster {
  int id;
  cluster *parent;
  short ox;
  short oy;
  short rank;
  int size;
};

typedef struct {
  cluster **slots;       // The slots in the hash table
  int slots_num;         // The numbers of slots in the hash table
  cluster *cluster_pool; // Cluster pool
  int cluster_cap;       // Capacity of the cluster pool
  int cluster_num;       // Clusters in the cluster pool
  int largest_cluster;   // The size of the largest cluster
  int **cluster_sizes;   // An array of cluster_size arrays, for each subgrid
  int num_grids;         // The number of subgrids
  bool percolates;       // Whether or not the grid percolates
} cluster_info;

cluster_info create_clusters_from_site_grids(struct site_grid *grids, int num_grids);

cluster_info create_clusters_from_bond_grids(struct bond_grid *grids, int num_grids);

void merge_clusters(cluster_info *ci, int c1, int c2, short ox, short oy);

void free_cluster_info(cluster_info *ci);

#endif