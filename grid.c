#include "header.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define SEED 1000

grid alloc_grid(char t, int sx, int sy)
{
  grid g;
  g.t = t;
  g.sx = sx;
  g.sy = sy;
  g.grid = malloc(sx * sy * sizeof(char));
  g.cluster = calloc(sx * sy, sizeof(int));
  g.n_cluster = 0;
  g.cluster_size = NULL;
  g.max_cluster = 0;
  return g;
}

void free_grid(grid g)
{
  free(g.grid);
  free(g.cluster);
  free(g.cluster_size);
}

void seed_grid(grid g, float p)
{
  int d = (g.t == 'b' ? 2 : 1);
  unsigned int rand_buffer = SEED;
  for (int i = 0; i < d * (grid_ox + (grid_oy * grid_tx)); i++) rand_r(&rand_buffer);
  for (int i = 0; i < g.sx * g.sy; i++) {
    char j = 0;
    if (((float)rand_r(&rand_buffer) / (float)RAND_MAX) < p) {
      j += 1;
    }
    if (g.t == 'b') {
      if (((float)rand_r(&rand_buffer) / (float)RAND_MAX) < p) {
        j += 2;
      }
    }
    g.grid[i] = j;
    if ((i + 1) % g.sx == 0) {
      for (int j = 0; j < d * (grid_tx - g.sx); j++) rand_r(&rand_buffer);
    }
  }
}

void print_grid_bond(FILE *f, grid g)
{
  for (int y = 0; y < g.sy; y++) {
    for (int x = 0; x < g.sx; x++) {
      fprintf(f, "O%s", g.grid[x + (y * g.sx)] & 1 ? "-" : " ");
    }
    fprintf(f, "\n");
    for (int x = 0; x < g.sx; x++) {
      fprintf(f, "%s ", g.grid[x + (y * g.sx)] & 2 ? "|" : " ");
    }
    fprintf(f, "\n");
  }
}

void print_grid(FILE *f, grid g, bool clust)
{
  if (g.t == 'b' && !clust) {
    print_grid_bond(f, g);
    return;
  }
  for (int y = 0; y < g.sy; y++) {
    for (int x = 0; x < g.sx; x++) {
      if (!clust) {
        fprintf(f, "%i", g.grid[x + (y * g.sx)]);
      } else {
        fprintf(f, "%i", g.cluster[x + (y * g.sx)]);
      }
    }
    fprintf(f, "\n");
  }
}

void grid_do_dfs(grid *g)
{
  int gsx = g->sx, gsy = g->sy;
  int cluster = 1;
  int *stack = malloc(2 * gsx * gsy * sizeof(int));
  int largest_cluster = 0;
  int cluster_cap = 0;
  for (int x=0; x<gsx; x++) {
    for (int y=0; y<gsy; y++) {
      int ind = (y * gsx) + x;
      if (g->cluster[ind] != 0) continue;
      if (g->t == 's' && g->grid[ind] == 0) continue;
      int sp = 0;
      g->cluster[ind] = cluster;
      stack[sp++] = x;
      stack[sp++] = y;
      int clust_sz = 1;
      while (sp > 0) {
        int cy = stack[--sp];
        int cx = stack[--sp];
        int node = (cy * gsx) + cx;
        if (cx > 0) {
          int left = (cy * gsx) + cx - 1;
          if (g->cluster[left] == 0 && (g->t == 'b' ? g->grid[left] & 1 : g->grid[left])) {
            g->cluster[left] = cluster;
            clust_sz++;
            stack[sp++] = cx - 1;
            stack[sp++] = cy;
          }
        }
        if (cx < gsx - 1) {
          int right = (cy * gsx) + cx + 1;
          if (g->cluster[right] == 0 && (g->t == 'b' ? g->grid[node] & 1 : g->grid[right])) {
            g->cluster[right] = cluster;
            clust_sz++;
            stack[sp++] = cx + 1;
            stack[sp++] = cy;
          }
        }
        if (cy > 0) {
          int top = ((cy - 1) * gsx) + cx;
          if (g->cluster[top] == 0 && (g->t == 'b' ? g->grid[top] & 2 : g->grid[top])) {
            g->cluster[top] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy - 1;
          }
        }
        if (cy < gsy - 1) {
          int bottom = ((cy + 1) * gsx) + cx;
          if (g->cluster[bottom] == 0 && (g->t == 'b' ? g->grid[node] & 2 : g->grid[bottom])) {
            g->cluster[bottom] = cluster;
            clust_sz++;
            stack[sp++] = cx;
            stack[sp++] = cy + 1;
          }
        }
      }
      if (cluster >= cluster_cap) {
        cluster_cap += 4 * gsx;
        g->cluster_size = realloc(g->cluster_size, cluster_cap * sizeof(int));
      }
      g->cluster_size[cluster - 1] = clust_sz;
      if (clust_sz > largest_cluster) {
        largest_cluster = clust_sz;
      }
      cluster++;
    }
  }
  free(stack);
  g->n_cluster = cluster - 1;
  g->max_cluster = largest_cluster;
}
