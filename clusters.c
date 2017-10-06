#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include "clusters.h"

cluster_info create_clusters_from_site_grids(struct site_grid *grids, int num_grids)
{
  cluster_info ci;
  ci.num_grids = num_grids;
  int largest_cluster = 0, max_cluster = 0, slots = 0;
  ci.cluster_sizes = malloc(num_grids * sizeof(int*));
  for (int i = 0; i < num_grids; i++) {
    slots += 8 * (grids[i].sx + grids[i].sy);
    if (grids[i].max_cluster > max_cluster) {
      max_cluster = grids[i].max_cluster;
    }
    if (grids[i].largest_cluster > largest_cluster) {
      largest_cluster = grids[i].largest_cluster;
    }
    ci.cluster_sizes[i] = grids[i].cluster_size;
  }
  ci.largest_cluster = largest_cluster;
  ci.slots_num = slots;
  ci.slots = calloc(slots, sizeof(cluster));
  ci.percolates = false;
  return ci;
}

int hash_key(int key, int slots)
{
  return key % slots;
}

cluster_info create_clusters_from_bond_grids(struct bond_grid *grids, int num_grids)
{
  cluster_info ci;
  ci.num_grids = num_grids;
  int largest_cluster = 0, max_cluster = 0, slots = 0;
  ci.cluster_sizes = malloc(num_grids * sizeof(int*));
  for (int i = 0; i < num_grids; i++) {
    slots += 8 * (grids[i].sx + grids[i].sy);
    if (grids[i].max_cluster > max_cluster) {
      max_cluster = grids[i].max_cluster;
    }
    if (grids[i].largest_cluster > largest_cluster) {
      largest_cluster = grids[i].largest_cluster;
    }
    ci.cluster_sizes[i] = grids[i].cluster_size;
  }
  ci.largest_cluster = largest_cluster;
  ci.slots_num = slots;
  ci.slots = calloc(slots, sizeof(cluster));
  ci.percolates = false;
  return ci;
}

cluster *get_cluster(cluster_info *ci, int id)
{
  cluster *c = &ci->slots[hash_key(id, ci->slots_num)];
  cluster *end = ci->slots + ci->slots_num;
  while (c->id != 0 && c->id != id) {
    c++;
    if (c == end) c = ci->slots;
  }
  if (c->id == 0) {
    c->id = id;
    c->parent = NULL;
    c->rank = 0;
    c->size = ci->cluster_sizes[(id - 1) % ci->num_grids][(id - 1) / ci->num_grids];
  }
  return c;
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
  free(ci->cluster_sizes);
}