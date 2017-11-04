#include "header.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

#define INT_HBIT ~(~0u >> 1)

outline alloc_outline(int sx, int sy, int n_cluster)
{
  outline o;
  o.sx = sx;
  o.sy = sy;
  o.max_cluster = 0;
  o.n_cluster = n_cluster;
  o.buffer = malloc( (n_cluster + sx + sx + sy + sy) * sizeof(int) );
  o.top = o.buffer;
  o.bottom = o.buffer + sx;
  o.left = o.buffer + (2 * sx);
  o.right = o.buffer + (2 * sx) + sy;
  o.cluster_size = o.buffer + (2 * sx) + (2 * sy);
  o.percolates = false;
  return o;
}

void free_outline(outline o)
{
  free(o.buffer);
}

void send_outline(outline o, int dest)
{
  int info[4];
  info[0] = o.sx;
  info[1] = o.sy;
  info[2] = o.n_cluster;
  info[3] = o.max_cluster;
  MPI_Send(info, 4, MPI_INT, dest, TAG_OUTLINE_INFO, MPI_COMM_WORLD);
  int buff_size = o.n_cluster + (2 * o.sx) + (2 * o.sy);
  MPI_Send(o.buffer, buff_size, MPI_INT, dest, TAG_OUTLINE_DATA, MPI_COMM_WORLD);
}

outline recv_outline(int src)
{
  MPI_Status status;
  int info[4];
  MPI_Recv(info, 4, MPI_INT, src, TAG_OUTLINE_INFO, MPI_COMM_WORLD, &status);
  outline o = alloc_outline(info[0], info[1], info[2]);
  o.max_cluster = info[3];
  int buff_size = o.n_cluster + (2 * o.sx) + (2 * o.sy);
  MPI_Recv(o.buffer, buff_size, MPI_INT, src, TAG_OUTLINE_DATA, MPI_COMM_WORLD, &status);
  return o;
}

static int find_clust_id(int oc, int os, int *parent, int *c, int *sizes, int *sizes2)
{
  if (oc == -1) return -1;
  oc += os;
  while ((parent[oc] != 0) && ((parent[oc] & INT_HBIT) == 0)) {
    oc = parent[oc];
  }
  if (parent[oc] == 0) {
    parent[oc] = INT_HBIT | *c;
    sizes2[*c] = sizes[oc];
    (*c)++;
  }
  return parent[oc] & ~INT_HBIT;
}

outline merge_outlines_horiz(outline o1, outline o2)
{
  if (o1.sy != o2.sy) {
    fprintf(stderr, "Outlines must have same height to be horizontally merged.\n");
  }
  int c1 = o1.n_cluster, c2 = o2.n_cluster;
  outline o3 = alloc_outline(o1.sx + o2.sx, o1.sy, o1.n_cluster + o2.n_cluster);
  int *parent = calloc(c1 + c2, sizeof(int));
  int *sizes = malloc((c1 + c2) * sizeof(int));
  memcpy(sizes, o1.cluster_size, c1 * sizeof(int));
  memcpy(sizes + c1, o2.cluster_size, c2 * sizeof(int));
  int mc = (o1.max_cluster > o2.max_cluster) ? o1.max_cluster : o2.max_cluster;
  // Merge the edge
  for (int y = 0; y < o1.sy; y++) {
    int oc1 = o1.right[y];
    if (oc1 == -1) continue;
    while (parent[oc1] != 0) oc1 = parent[oc1];
    int oc2 = o2.left[y];
    if (oc2 == -1) continue;
    oc2 += c1;
    while (parent[oc2] != 0) oc2 = parent[oc2];
    if (oc1 == oc2) continue;
    parent[oc1] = oc2;
    sizes[oc2] += sizes[oc1];
    if (sizes[oc2] > mc) mc = sizes[oc2];
  }
  // Reallocate cluster numbers
  int c = 0;
  // Left edge
  for (int y = 0; y < o1.sy; y++) {
    o3.left[y] = find_clust_id(o1.left[y], 0, parent, &c, sizes, o3.cluster_size);
  }
  // Right edge
  for (int y = 0; y < o2.sy; y++) {
    o3.right[y] = find_clust_id(o2.right[y], c1, parent, &c, sizes, o3.cluster_size);
  }
  // Top edge
  for (int x = 0; x < o1.sx; x++) {
    o3.top[x] = find_clust_id(o1.top[x], 0, parent, &c, sizes, o3.cluster_size);
  }
  for (int x = 0; x < o2.sx; x++) {
    o3.top[o1.sx + x] = find_clust_id(o2.top[x], c1, parent, &c, sizes, o3.cluster_size);
  }
  // Bottom edge
  for (int x = 0; x < o1.sx; x++) {
    o3.bottom[x] = find_clust_id(o1.bottom[x], 0, parent, &c, sizes, o3.cluster_size);
  }
  for (int x = 0; x < o2.sx; x++) {
    o3.bottom[o1.sx + x] = find_clust_id(o2.bottom[x], c1, parent, &c, sizes, o3.cluster_size);
  }
  free(parent);
  free(sizes);
  o3.n_cluster = c;
  o3.max_cluster = mc;
  return o3;
}

outline merge_outlines_vert(outline o1, outline o2)
{
  if (o1.sx != o2.sx) {
    fprintf(stderr, "Outlines must have same width to be vertically merged.\n");
  }
  int c1 = o1.n_cluster, c2 = o2.n_cluster;
  outline o3 = alloc_outline(o1.sx, o1.sy + o2.sy, o1.n_cluster + o2.n_cluster);
  int *parent = calloc(c1 + c2, sizeof(int));
  int *sizes = malloc((c1 + c2) * sizeof(int));
  memcpy(sizes, o1.cluster_size, c1 * sizeof(int));
  memcpy(sizes + c1, o2.cluster_size, c2 * sizeof(int));
  int mc = (o1.max_cluster > o2.max_cluster) ? o1.max_cluster : o2.max_cluster;
  // Merge the edge
  for (int x = 0; x < o1.sx; x++) {
    int oc1 = o1.bottom[x];
    if (oc1 == -1) continue;
    while (parent[oc1] != 0) oc1 = parent[oc1];
    int oc2 = o2.top[x];
    if (oc2 == -1) continue;
    oc2 += c1;
    while (parent[oc2] != 0) oc2 = parent[oc2];
    if (oc1 == oc2) continue;
    parent[oc1] = oc2;
    sizes[oc2] += sizes[oc1];
    if (sizes[oc2] > mc) mc = sizes[oc2];
  }
  // Reallocate cluster numbers
  int c = 0;
  // Left edge
  for (int y = 0; y < o1.sy; y++) {
    o3.left[y] = find_clust_id(o1.left[y], 0, parent, &c, sizes, o3.cluster_size);
  }
  for (int y = 0; y < o2.sy; y++) {
    o3.left[o1.sy + y] = find_clust_id(o2.left[y], c1, parent, &c, sizes, o3.cluster_size);
  }
  // Right edge
  for (int y = 0; y < o1.sy; y++) {
    o3.right[y] = find_clust_id(o1.right[y], 0, parent, &c, sizes, o3.cluster_size);
  }
  for (int y = 0; y < o2.sy; y++) {
    o3.right[o1.sy + y] = find_clust_id(o2.right[y], c1, parent, &c, sizes, o3.cluster_size);
  }
  // Top edge
  for (int x = 0; x < o1.sx; x++) {
    o3.top[x] = find_clust_id(o1.top[x], 0, parent, &c, sizes, o3.cluster_size);
  }
  // Bottom edge
  for (int x = 0; x < o2.sx; x++) {
    o3.bottom[x] = find_clust_id(o2.bottom[x], c1, parent, &c, sizes, o3.cluster_size);
  }
  free(parent);
  free(sizes);
  o3.n_cluster = c;
  o3.max_cluster = mc;
  return o3;
}

outline outline_from_grid(grid g)
{
  outline o = alloc_outline(g.sx, g.sy, g.n_cluster);
  o.max_cluster = g.max_cluster;
  memcpy(o.cluster_size, g.cluster_size, o.n_cluster * sizeof(int));
  bool bond;
  for (int x = 0; x < o.sx; x++) {
    int it = x;
    int ib = x + ((o.sy - 1) * o.sx);
    bond = g.t == 's' ? (g.grid[it] != 0) : true;
    o.top[x] = bond ? g.cluster[it] - 1 : -1;
    bond = g.t == 's' ? (g.grid[ib] != 0) : (g.grid[ib] & 2);
    o.bottom[x] = bond ? g.cluster[ib] - 1 : -1;
  }
  for (int y = 0; y < o.sy; y++) {
    int il = y * o.sx;
    int ir = ((y + 1) * o.sx) - 1;
    bond = g.t == 's' ? (g.grid[il] != 0) : true;
    o.left[y] = bond ? g.cluster[il] - 1 : -1;
    bond = g.t == 's' ? (g.grid[ir] != 0) : (g.grid[ir] & 1);
    o.right[y] = bond ? g.cluster[ir] - 1 : -1;
  }
  return o;
}

void percolate_outline(outline *o)
{
  int *parent = calloc(o->n_cluster, sizeof(int));
  char *ox = malloc(o->n_cluster * sizeof(char));
  char *oy = malloc(o->n_cluster * sizeof(char));
  for (int i = 0; i < o->sx + o->sy; i++) {
    int c1, c2;
    char nox, noy;
    if (i < o->sx) {
      c1 = o->top[i];
      c2 = o->bottom[i];
      nox = 0;
      noy = 1;
    } else {
      c1 = o->left[i - o->sx];
      c2 = o->right[i - o->sx];
      nox = 1;
      noy = 0;
    }
    if (c1 == -1 || c2 == -1) continue;
    while (parent[c1] - 1 != -1) {
      nox += ox[c1]; noy += oy[c1];
      c1 = parent[c1] - 1;
    }
    while (parent[c2] - 1 != -1) {
      nox -= ox[c2]; noy -= oy[c2];
      c2 = parent[c2] - 1;
    }
    if (c1 == c2) {
      if (nox != 0 || noy != 0) o->percolates = true;
    } else {
      parent[c2] = c1 + 1;
      ox[c2] = nox;
      oy[c2] = noy;
      int clust_sz = o->cluster_size[c1] + o->cluster_size[c2];
      if (clust_sz > o->max_cluster) {
        o->max_cluster = clust_sz;
      }
      o->cluster_size[c1] = clust_sz;
    }
  }
  free(parent);
  free(ox);
  free(oy);
}

void print_outline(FILE *f, outline o)
{
  fprintf(f, " ");
  for (int x=0; x<o.sx; x++) {
    if (o.top[x] == -1) fprintf(f, "-"); else fprintf(f, "%i", o.top[x]);
  }
  fprintf(f, "\n");
  char *gap = malloc(o.sx + 1);
  for (int x=0; x<o.sx; x++) gap[x] = ' ';
  gap[o.sx] = 0;
  for (int y=0; y<o.sy; y++) {
    if (o.left[y] == -1 && o.right[y] == -1)
      fprintf(f, "-%s-", gap);
    if (o.left[y] != -1 && o.right[y] == -1)
      fprintf(f, "%i%s-", o.left[y], gap);
    if (o.left[y] == -1 && o.right[y] != -1)
      fprintf(f, "-%s%i", gap, o.right[y]);
    if (o.left[y] != -1 && o.right[y] != -1)
      fprintf(f, "%i%s%i", o.left[y], gap, o.right[y]);
    fprintf(f, "\n");
  }
  fprintf(f, " ");
  for (int x=0; x<o.sx; x++) {
    if (o.bottom[x] == -1) fprintf(f, "-"); else fprintf(f, "%i", o.bottom[x]);
  }
  fprintf(f, "\n");
  free(gap);
}

void print_outline_consistent(FILE *f, outline o)
{
  if (o.top[0] != o.left[0]) fprintf(stderr, "Top-left corner inconsistent.\n");
  if (o.top[o.sx - 1] != o.right[0]) fprintf(stderr, "Top-right corner inconsistent.\n");
  if (o.bottom[0] != o.left[o.sy - 1]) fprintf(stderr, "Bottom-left corner inconsistent.\n");
  if (o.bottom[o.sx - 1] != o.right[o.sy - 1]) fprintf(stderr, "Bottom-right corner inconsistent.\n");
  for (int x=0; x<o.sx; x++) {
    if (o.top[x] == -1) fprintf(f, "-"); else fprintf(f, "%i", o.top[x]);
  }
  fprintf(f, "\n");
  char *gap = malloc(o.sx - 1);
  for (int x=0; x<o.sx-2; x++) gap[x] = ' ';
  gap[o.sx-2] = 0;
  for (int y=1; y<o.sy-1; y++) {
    if (o.left[y] == -1 && o.right[y] == -1)
      fprintf(f, "-%s-", gap);
    if (o.left[y] != -1 && o.right[y] == -1)
      fprintf(f, "%i%s-", o.left[y], gap);
    if (o.left[y] == -1 && o.right[y] != -1)
      fprintf(f, "-%s%i", gap, o.right[y]);
    if (o.left[y] != -1 && o.right[y] != -1)
      fprintf(f, "%i%s%i", o.left[y], gap, o.right[y]);
    fprintf(f, "\n");
  }
  for (int x=0; x<o.sx; x++) {
    if (o.bottom[x] == -1) fprintf(f, "-"); else fprintf(f, "%i", o.bottom[x]);
  }
  fprintf(f, "\n");
  free(gap);
}
