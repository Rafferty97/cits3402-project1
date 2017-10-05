#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

#include "site.h"

#define PERCOLATION_TYPE "SITE"

typedef struct {
  bool *site;          // 2D array specifying if a lattice point (site) is occupied
  int *cluster;        // 2D array specifying which cluster a site belongs to (0 if none or unknown)
  int *cluster_size;   // The number of sites in each cluster (index 0 = cluster 1)
  int cluster_cap;     // The current capacity of the cluster_size array
  int sx;              // Horizontal size of the grid
  int sy;              // Vertical size of the grid
  int c_c;             // The first valid cluster ID (prevents overlap between subgrids)
  int c_n;             // The gap between each valid cluster ID (prevents overlap between subgrids)
  int max_cluster;     // The largest cluster ID
} grid;

typedef struct {
  int *cluster_size;   // The number of sites in each cluster (index 0 = cluster 1)
  int *cluster_alias;  // The smallest connected cluster ID (defines an equivalence relation)
  int max_cluster;     // The largest cluster ID (and length of the arrays)
} cluster_info;

grid create_grid(int sx, int sy, int c_c, int c_n)
{
  grid g;
  g.site = malloc(sx * sy * sizeof(bool));
  g.cluster = calloc(sx * sy, sizeof(int));
  g.cluster_size = malloc(4 * sx * sizeof(int));
  g.cluster_cap = 4 * sx;
  g.sx = sx;
  g.sy = sy;
  g.c_c = c_c;
  g.c_n = c_n;
  return g;
}

void seed_grid(grid g, float p)
{
  for (int i = 0; i < g.sx * g.sy; i++) {
    g.site[i] = ((float)rand() / (float)RAND_MAX) < p;
  }
}

void grid_do_dfs(grid *g)
{
  int cluster = g->c_c, clusters = 0;
  int *stack = malloc(2 * g->sx * g->sy * sizeof(int));
  for (int x=0; x<g->sx; x++) {
    for (int y=0; y<g->sy; y++) {
      int ind = (y * g->sx) + x;
      if (!g->site[ind] || g->cluster[ind] != 0) continue;
      int sp = 0;
      g->cluster[ind] = cluster;
      stack[sp++] = x;
      stack[sp++] = y;
      int clust_sz = 1;
      while (sp > 0) {
        int cy = stack[--sp];
        int cx = stack[--sp];
        if (cx > 0) {
          int left = (cy * g->sx) + cx - 1;
          if (g->site[left] && g->cluster[left] == 0) {
            g->cluster[left] = cluster;
            clust_sz++;
            stack[sp++] = cx - 1;
            stack[sp++] = cy;
          }
        }
        if (cx < g->sx - 1) {
          int right = (cy * g->sx) + cx + 1;
          if (g->site[right] && g->cluster[right] == 0) {
            g->cluster[right] = cluster;
            clust_sz++;
            stack[sp++] = cx + 1;
            stack[sp++] = cy;
          }
        }
        if (cy > 0) {
          int top = ((cy - 1) * g->sx) + cx;
          if (g->site[top] && g->cluster[top] == 0) {
            g->cluster[top] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy - 1;
          }
        }
        if (cy < g->sy - 1) {
          int bottom = ((cy + 1) * g->sx) + cx;
          if (g->site[bottom] && g->cluster[bottom] == 0) {
            g->cluster[bottom] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy + 1;
          }
        }
      }
      if (clusters >= g->cluster_cap) {
        g->cluster_cap += 4 * g->sx;
        g->cluster_size = realloc(g->cluster_size, g->cluster_cap * sizeof(int));
      }
      g->cluster_size[clusters++] = clust_sz;
      cluster += g->c_n;
    }
  }
  free(stack);
  g->max_cluster = cluster - g->c_n;
}

cluster_info create_cluster_info(grid *grids, int num_grids)
{
  cluster_info c;
  int max = 0;
  for (int i = 0; i < num_grids; i++) {
    if (grids[i].max_cluster > max) max = grids[i].max_cluster;
  }
  c.cluster_size = calloc(max, sizeof(int));
  for (int i = 0; i < num_grids; i++) {
    int k = 0;
    for (int j = grids[i].c_c; j <= grids[i].max_cluster; j += grids[i].c_n) {
      c.cluster_size[j - 1] = grids[i].cluster_size[k++];
    }
  }
  c.cluster_alias = calloc(max, sizeof(int));
  c.max_cluster = max;
  return c;
}

void merge_grids(grid *g1, grid *g2, cluster_info *clusters)
{
  int *sizes = clusters->cluster_size;
  int *alias = clusters->cluster_alias;
  int sx = g1->sx;
  for (int x = 0; x < sx; x++) {
    int ind = ((g1->sy - 1) * sx) + x;
    int ind2 = x;
    if (!g1->site[ind] || !g2->site[ind2]) continue;
    int c1 = g1->cluster[ind];
    int c2 = g2->cluster[ind2];
    while (alias[c1 - 1] != 0) c1 = alias[c1 - 1];
    while (alias[c2 - 1] != 0) c2 = alias[c2 - 1];
    if (c1 == c2) continue;
    if (c1 > c2) {
      int t = c2;
      c2 = c1;
      c1 = t;
    }
    alias[c2 - 1] = c1;
    sizes[c1 - 1] += sizes[c2 - 1];
    sizes[c2 - 1] = 0;
  }
}

bool grids_do_percolate(grid *grids, int num_grids, cluster_info c)
{
  int *alias = c.cluster_alias;
  int sx = grids[0].sx;
  grid *first_g = &grids[0];
  grid *last_g = &grids[num_grids - 1];
  for (int x = 0; x < sx; x++) {
    int f_ind = x;
    int l_ind = ((last_g->sy - 1) * sx) + x;
    if (!first_g->site[f_ind] || !last_g->site[l_ind]) continue;
    int c1 = first_g->cluster[f_ind];
    int c2 = last_g->cluster[l_ind];
    while (alias[c1 - 1] != 0) c1 = alias[c1 - 1];
    while (alias[c2 - 1] != 0) c2 = alias[c2 - 1];
    if (c1 == c2) return true;
  }
  for (int i = 0; i < num_grids; i++) {
    grid *g = &grids[i];
    for (int y = 0; y < g->sy; y++) {
      int ind = y * sx;
      int ind2 = ((y + 1) * sx) - 1;
      if (!g->site[ind] || !g->site[ind2]) continue;
      int c1 = g->cluster[ind];
      int c2 = g->cluster[ind2];
      while (alias[c1 - 1] != 0) c1 = alias[c1 - 1];
      while (alias[c2 - 1] != 0) c2 = alias[c2 - 1];
      if (c1 == c2) return true;
    }
  }
  return false;
}

int largest_cluster(cluster_info c)
{
  int max = 0;
  for (int i = 0; i < c.max_cluster; i++) {
    if (c.cluster_size[i] > max) max = c.cluster_size[i];
  }
  return max;
}

void free_grid(grid g)
{
  free(g.site);
  free(g.cluster);
  free(g.cluster_size);
}

void site_percolation(int size, float p, int threads, percolation_results *results)
{
  // Start timing
  struct timeval start, end;
  gettimeofday(&start, NULL);
  // Simulation parameters
  results->type = PERCOLATION_TYPE;
  results->lattice_size = size;
  results->p = p;
  results->threads = threads;
  // Divide the grid
  int subgrid_h = size / threads;
  grid *grids = malloc(threads * sizeof(grid));
  omp_set_num_threads(threads);
  #pragma omp parallel
  {
    int i = omp_get_thread_num();
    int hei = (i == threads - 1 ? size - (i * subgrid_h) : subgrid_h);
    grids[i] = create_grid(size, hei, i + 1, threads);
    seed_grid(grids[i], p);
    grid_do_dfs(&grids[i]);
  }
  // Merge cluster stats
  cluster_info clusters = create_cluster_info(grids, threads);
  // Merge clusters
  for (int i = 0; i < threads - 1; i++) {
    merge_grids(&grids[i], &grids[i + 1], &clusters);
  }
  // Calculate results
  results->percolates = grids_do_percolate(grids, threads, clusters);
  results->largest_cluster = largest_cluster(clusters);
  // Free grid memory
  for (int i = 0; i < threads; i++) free_grid(grids[i]);
  // Stop timing
  gettimeofday(&end, NULL);
  float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  results->time_taken = delta;
}