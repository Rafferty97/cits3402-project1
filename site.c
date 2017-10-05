#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

#include "site.h"

#define PERCOLATION_TYPE "Site percolation"

typedef struct {
  bool *site;          // 2D array specifying if a lattice point (site) is occupied
  int *cluster;        // 2D array specifying which cluster a site belongs to (0 if none or unknown)
  int *cluster_sizes;  // The number of sites in each cluster (index 0 = cluster 1)
  int stride;          // If cell (x, y) is at index i, cell (x, y+1) is at index i+stride
  int sx;              // Horizontal size of the grid
  int sy;              // Vertical size of the grid
  int c_c;             // The first valid cluster ID (prevents overlap between subgrids)
  int c_n;             // The gap between each valid cluster ID (prevents overlap between subgrids)
  int *cluster_alias;  // The lowest cluster ID that connects to each cluster (defines an equivalence relation)
  int max_cluster;     // The largest cluster ID
} grid;

void grid_print(grid g)
{
  printf("Size: %i x %i\n", g.sx, g.sy);
  printf("c_c: %i, c_n: %i\n\n", g.c_c, g.c_n);
  for (int y=0; y<g.sy; y++) {
    for (int x=0; x<g.sx; x++) {
      printf("%02d ", g.cluster[x + y*g.stride]);
    }
    printf("\n\n");
  }
  printf("Clusters (%i): ", (*g.max_cluster));
  for (int i=0; i<(*g.max_cluster); i++) {
    printf("%i ", g.cluster_sizes[i]);
  }
  printf("\nAliases: ");
  for (int i=0; i<(*g.max_cluster); i++) {
    printf("%i ", g.cluster_alias[i]);
  }
}

grid create_grid(int size)
{
  grid g;
  g.site = malloc(size * size * sizeof(bool));
  g.cluster = malloc(size * size * sizeof(int));
  g.cluster_sizes = NULL;
  g.sx = size;
  g.sy = size;
  g.stride = size;
  g.c_c = 1;
  g.c_n = 1;
  g.max_cluster = malloc(sizeof(int));
  return g;
}

grid create_subgrid(int sx, int sy, int c_c, int c_n)
{
  grid g;
  g.site = malloc(sx * sy * sizeof(bool));
  g.cluster = calloc(sx * sy, sizeof(int));
  g.cluster_sizes = calloc(sx * sy, sizeof(int));
  g.sx = sx;
  g.sy = sy;
  g.stride = sx;
  g.c_c = c_c;
  g.c_n = c_n;
  g.max_cluster = malloc(sizeof(int));
  return g;
}

void seed_grid(grid g, float p)
{
  for (int i = 0; i < g.sx * g.sy; i++) {
    g.site[i] = ((float)rand() / (float)RAND_MAX) < p;
  }
}

/* grid create_subgrid(grid g, int x, int y, int w, int h, int c_c, int c_n)
{
  grid subgrid;
  int offset = (y * g.stride) + x;
  subgrid.site = g.site + offset;
  subgrid.cluster = g.cluster + offset;
  subgrid.cluster_sizes = g.cluster_sizes;
  subgrid.stride = g.stride;
  subgrid.sx = w;
  subgrid.sy = h;
  subgrid.c_n = g.c_n * c_n;
  subgrid.c_c = g.c_c + (g.c_n * c_c);
  subgrid.max_cluster = g.max_cluster;
  return subgrid;
} */

void merge_subgrid(grid *g, grid sg, int x, int y)
{
  for (int yy = 0; yy < sg.sy; yy++) {
    memcpy(
      g->site + (x + (g->stride * (y + yy))),
      sg.site + (sg.stride * yy),
      sg.sx * sizeof(bool)
    );
    memcpy(
      g->cluster + (x + (g->stride * (y + yy))),
      sg.cluster + (sg.stride * yy),
      sg.sx * sizeof(int)
    );
  }
  if ((*g->max_cluster) < (*sg.max_cluster)) {
    (*g->max_cluster) = (*sg.max_cluster);
  }
  g->cluster_sizes = realloc(g->cluster_sizes, (*g->max_cluster) * sizeof(int));
  for (int i = sg.c_c; i <= (*sg.max_cluster); i += sg.c_n) {
    g->cluster_sizes[i - 1] = sg.cluster_sizes[i - 1];
  }
}

void grid_do_dfs(grid g)
{
  int cluster = g.c_c;
  int *stack = malloc(2 * g.sx * g.sy * sizeof(int));
  for (int x=0; x<g.sx; x++) {
    for (int y=0; y<g.sy; y++) {
      int ind = (y * g.stride) + x;
      if (!g.site[ind] || g.cluster[ind] != 0) continue;
      int sp = 0;
      g.cluster[ind] = cluster;
      stack[sp++] = x;
      stack[sp++] = y;
      int clust_sz = 1;
      while (sp > 0) {
        int cy = stack[--sp];
        int cx = stack[--sp];
        if (cx > 0) {
          int left = (cy * g.stride) + cx - 1;
          if (g.site[left] && g.cluster[left] == 0) {
            g.cluster[left] = cluster;
            clust_sz++;
            stack[sp++] = cx - 1;
            stack[sp++] = cy;
          }
        }
        if (cx < g.sx - 1) {
          int right = (cy * g.stride) + cx + 1;
          if (g.site[right] && g.cluster[right] == 0) {
            g.cluster[right] = cluster;
            clust_sz++;
            stack[sp++] = cx + 1;
            stack[sp++] = cy;
          }
        }
        if (cy > 0) {
          int top = ((cy - 1) * g.stride) + cx;
          if (g.site[top] && g.cluster[top] == 0) {
            g.cluster[top] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy - 1;
          }
        }
        if (cy < g.sy - 1) {
          int bottom = ((cy + 1) * g.stride) + cx;
          if (g.site[bottom] && g.cluster[bottom] == 0) {
            g.cluster[bottom] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy + 1;
          }
        }
      }
      g.cluster_sizes[cluster - 1] = clust_sz;
      cluster += g.c_n;
    }
  }
  free(stack);
  int max_cluster = cluster - g.c_n;
  (*g.max_cluster) = max_cluster;
}

void grid_init_aliases(grid *g)
{
  int max = (*g->max_cluster);
  g->cluster_alias = malloc(max * sizeof(int));
  for (int i = 0; i<max; i++) {
    g->cluster_alias[i] = i + 1;
  }
}

void grid_merge_clusters(grid g, int pos, bool vert)
{
  int *alias = g.cluster_alias;
  int max = (*g.max_cluster);
  int adj = vert ? 1 : g.stride;
  int len = vert ? g.sy : g.sx;
  for (int i=0; i<len; i++) {
    int ind = vert ? (i * g.stride) + pos : (pos * g.stride) + i;
    int ind2 = ind + adj;
    if (!g.site[ind] || !g.site[ind2]) continue;
    int c1 = alias[g.cluster[ind] - 1], c2 = alias[g.cluster[ind2] - 1];
    if (c1 == c2) continue;
    if (c1 > c2) {
      int t = c2;
      c2 = c1;
      c1 = t;
    }
    for (int i=c2 - 1; i<max; i++) {
      if (alias[i] == c2) {
        alias[i] = c1;
      }
    }
    g.cluster_sizes[c1 - 1] += g.cluster_sizes[c2 - 1];
    g.cluster_sizes[c2 - 1] = 0;
  }
}

bool grid_does_percolate(grid g, int threads)
{
  int *alias = g.cluster_alias;
  for (int x=0; x<g.sx; x++) {
    int ind = x;
    int ind2 = ((g.sy - 1) * g.stride) + x;
    if (!g.site[ind] || !g.site[ind2]) continue;
    if (alias[g.cluster[ind] - 1] == alias[g.cluster[ind2] - 1]) {
      return true;
    }
  }
  for (int y=0; y<g.sy; y++) {
    int ind = y * g.stride;
    int ind2 = (y * g.stride) + (g.sx - 1);
    if (!g.site[ind] || !g.site[ind2]) continue;
    if (alias[g.cluster[ind] - 1] == alias[g.cluster[ind2] - 1]) {
      return true;
    }
  }
  return false;
}

int grid_largest_cluster(grid g)
{
  int max = 0;
  for (int i = 0; i < (*g.max_cluster); i++) {
    if (g.cluster_sizes[i] > max) max = g.cluster_sizes[i];
  }
  return max;
}

void free_grid(grid g)
{
  free(g.site);
  free(g.cluster);
  free(g.cluster_sizes);
  free(g.cluster_alias);
  free(g.max_cluster);
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
  // Create the grid
  grid g = create_grid(size);
  // Divide the grid
  int subgrid_w = size / threads;
  // Perform depth-first search
  omp_set_num_threads(threads);
  #pragma omp parallel
  {
    int i = omp_get_thread_num();
    int wid = (i == threads - 1 ? size - (i * subgrid_w) : subgrid_w);
    grid sg = create_subgrid(wid, size, i + 1, threads);
    seed_grid(sg, p);
    grid_do_dfs(sg);
    merge_subgrid(&g, sg, i * subgrid_w, 0);
  }
  // Create cluster alias array
  grid_init_aliases(&g);
  // Merge clusters
  for (int i = 1; i < threads; i++) {
    grid_merge_clusters(g, (i * subgrid_w) - 1, true);
  }
  // Calculate results
  results->percolates = grid_does_percolate(g, threads);
  results->largest_cluster = grid_largest_cluster(g);
  // Free grid memory
  free_grid(g);
  // Stop timing
  gettimeofday(&end, NULL);
  float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  results->time_taken = delta;
}

void site_percolation_batch(int size, float p, int threads, int trials, percolation_batch_results *results)
{
  srand(time(NULL));
  results->type = PERCOLATION_TYPE;
  results->lattice_size = size;
  results->p = p;
  results->threads = threads;
  results->trials = trials;
  results->num_percolates = 0;
  results->avg_largest_cluster = 0;
  results->avg_time_taken = 0;
  printf("Performing %i trials...\n", trials);
  for (int i = 0; i < trials; i++) {
    percolation_results trial;
    site_percolation(size, p, threads, &trial);
    if (trial.percolates) {
      results->num_percolates++;
    }
    results->avg_largest_cluster += trial.largest_cluster;
    results->avg_time_taken += trial.time_taken;
    if ((i + 1) % 100 == 0) printf("%i trials done.\n", i + 1);
  }
  printf("Finished.\n\n");
  results->avg_largest_cluster /= trials;
  results->avg_time_taken /= trials;
}
