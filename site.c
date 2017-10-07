#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

#include "site.h"
#include "clusters.h"

#define PERCOLATION_TYPE "SITE"

typedef struct site_grid grid;

static bool *site_grid_buffer = NULL;

static grid create_grid(int sx, int sy, int c_c, int c_n)
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

static void seed_grid(grid g, float p, unsigned long seed)
{
  unsigned int rand_buffer = seed;
  for (int i = 0; i < g.sx * g.sy; i++) {
    g.site[i] = ((float)rand_r(&rand_buffer) / (float)RAND_MAX) < p;
  }
}

static void grid_do_dfs(grid *g)
{
  int cluster = g->c_c, clusters = 0;
  int *stack = malloc(2 * g->sx * g->sy * sizeof(int));
  int largest_cluster = 0;
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
      if (clust_sz > largest_cluster) {
        largest_cluster = clust_sz;
      }
      cluster += g->c_n;
    }
  }
  free(stack);
  g->max_cluster = cluster - g->c_n;
  g->largest_cluster = largest_cluster;
}

static void merge_grids_horiz(grid *g1, grid *g2, cluster_info *ci, short oy)
{
  int sx = g1->sx;
  int last_row = (g1->sy - 1) * g1->sx;
  for (int x = 0; x < sx; x++) {
    int ind = last_row + x;
    int ind2 = x;
    if (!g1->site[ind] || !g2->site[ind2]) continue;
    int c1 = g1->cluster[ind];
    int c2 = g2->cluster[ind2];
    merge_clusters(ci, c1, c2, 0, oy);
  }
}

static void merge_grids_vert(grid *g1, grid *g2, cluster_info *ci, short ox)
{
  int sx1 = g1->sx;
  int sx2 = g2->sx;
  int sy = g1->sy;
  int last_col = sx1 - 1;
  for (int y = 0; y < sy; y++) {
    int ind = last_col + (y * sx1);
    int ind2 = y * sx2;
    if (!g1->site[ind] || !g2->site[ind2]) continue;
    int c1 = g1->cluster[ind];
    int c2 = g2->cluster[ind2];
    merge_clusters(ci, c1, c2, ox, 0);
  }
}

static void free_grid(grid g)
{
  free(g.site);
  free(g.cluster);
  free(g.cluster_size);
}

void site_percolation(int size, float p, unsigned long seed, percolation_results *results)
{
  // Start timing
  struct timeval start, end;
  gettimeofday(&start, NULL);
  // Simulation parameters
  results->type = PERCOLATION_TYPE;
  results->lattice_size = size;
  results->p = p;
  results->threads = 1;
  // Create grid
  grid g = create_grid(size, size, 1, 1);
  // Seeg grid
  if (site_grid_buffer != NULL) {
    memcpy(g.site, site_grid_buffer, size * size * sizeof(bool));
  } else {
    seed_grid(g, p, seed);
  }
  /* printf("\n");
  for (int y=0; y<size; y++) {
    for (int x=0; x<size; x++) {
      printf("%s", g.site[x + (y * size)] ? "X" : " ");
    }
    printf("\n");
  } */
  // Perform DFS
  grid_do_dfs(&g);
  // Merge clusters
  cluster_info clusters = create_clusters_from_site_grids(&g, 1);
  merge_grids_horiz(&g, &g, &clusters, 1);
  merge_grids_vert(&g, &g, &clusters, 1);
  // Calculate results
  results->percolates = clusters.percolates;
  results->largest_cluster = clusters.largest_cluster;
  // Free memory
  free_grid(g);
  free_cluster_info(&clusters);
  // Stop timing
  gettimeofday(&end, NULL);
  float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  results->time_taken = delta;
}

void site_percolation_parallel(int size, float p, unsigned long seed, percolation_results *results, int threads)
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
    grid g = create_grid(size, hei, i + 1, threads);
    if (site_grid_buffer != NULL) {
      int offset = i * subgrid_h * size;
      memcpy(g.site, site_grid_buffer + offset, size * hei * sizeof(bool));
    } else {
      seed_grid(g, p, seed + (i * 1961));
    }
    grid_do_dfs(&g);
    grids[i] = g;
  }
  // Merge clusters
  cluster_info clusters = create_clusters_from_site_grids(grids, threads);
  for (int i = 1; i < threads; i++) {
    merge_grids_horiz(&grids[i - 1], &grids[i], &clusters, 0);
    merge_grids_vert(&grids[i], &grids[i], &clusters, 1);
  }
  merge_grids_horiz(&grids[threads - 1], &grids[0], &clusters, 1);
  merge_grids_vert(&grids[0], &grids[0], &clusters, 1);
  // Calculate results
  results->percolates = clusters.percolates;
  results->largest_cluster = clusters.largest_cluster;
  // Free memory
  for (int i = 0; i < threads; i++) free_grid(grids[i]);
  free_cluster_info(&clusters);
  // Stop timing
  gettimeofday(&end, NULL);
  float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  results->time_taken = delta;
}

void load_site_grid(bool *sites)
{
  site_grid_buffer = sites;
}