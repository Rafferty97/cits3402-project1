#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include "clusters.h"

cluster_info create_clusters_from_site_grids(struct site_grid *grids, int num_grids)
{
  cluster_info ci;
  ci.num_grids = num_grids;
  int largest_cluster = 0, max_cluster = 0;
  ci.cluster_sizes = malloc(num_grids * sizeof(int*));
  for (int i = 0; i < num_grids; i++) {
    if (grids[i].max_cluster > max_cluster) {
      max_cluster = grids[i].max_cluster;
    }
    if (grids[i].largest_cluster > largest_cluster) {
      largest_cluster = grids[i].largest_cluster;
    }
    ci.cluster_sizes[i] = grids[i].cluster_size;
  }
  ci.largest_cluster = largest_cluster;
  ci.slots_num = max_cluster;
  ci.slots = calloc(ci.slots_num, sizeof(cluster*));
  ci.percolates = false;
  ci.cluster_cap = 1024;
  ci.cluster_pool = malloc(ci.cluster_cap * sizeof(cluster));
  ci.cluster_num = 0;
  return ci;
}

cluster *get_cluster(cluster_info *ci, int id)
{
  cluster **c = &ci->slots[id - 1];
  if (*c == NULL) {
    if (ci->cluster_num >= ci->cluster_cap) {
      ci->cluster_cap *= 2;
      ci->cluster_pool = realloc(ci->cluster_pool, ci->cluster_cap * sizeof(cluster));
    }
    *c = &ci->cluster_pool[ci->cluster_num];
    ci->cluster_num++;
    (*c)->id = id;
    (*c)->parent = NULL;
    (*c)->rank = 0;
    (*c)->size = ci->cluster_sizes[(id - 1) % ci->num_grids][(id - 1) / ci->num_grids];
  }
  return *c;
}

void merge_clusters(cluster_info *ci, int c1, int c2, short ox, short oy)
{
  cluster *clust1 = get_cluster(ci, c1);
  cluster *clust2 = get_cluster(ci, c2);
  while (clust1->parent != NULL) {
    ox -= clust1->ox;
    oy -= clust1->oy;
    clust1 = clust1->parent;
  }
  while (clust2->parent != NULL) {
    ox += clust2->ox;
    oy += clust2->oy;
    clust2 = clust2->parent;
  }
  if (clust1 == clust2) {
    if (ox != 0 || oy != 0) ci->percolates = true;
  } else {
    clust1->parent = clust2;
    clust1->ox = ox;
    clust1->oy = oy;
    clust2->size += clust1->size;
    clust1->size = 0;
    if (clust2->size > ci->largest_cluster) {
      ci->largest_cluster = clust2->size;
    }
  }
}

void free_cluster_info(cluster_info *ci)
{
  free(ci->slots);
  free(ci->cluster_pool);
  free(ci->cluster_sizes);
}