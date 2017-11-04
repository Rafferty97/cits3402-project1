#include <stdbool.h>
#include <stdio.h>

/* Constants */

#define TAG_OUTLINE_INFO 1
#define TAG_OUTLINE_DATA 2

/* Global variables */

int grid_tx;
int grid_ox;
int grid_oy;

/* Data structures */

typedef struct {
  char t;
  int sx;
  int sy;
  char *grid;
  int *cluster;
  int n_cluster;
  int *cluster_size;
  int max_cluster;
} grid;

typedef struct {
  int sx;
  int sy;
  int n_cluster;
  int max_cluster;
  int *cluster_size;
  int *top;
  int *left;
  int *bottom;
  int *right;
  int *buffer;
} outline;

/* grid.c */

grid alloc_grid(char t, int sx, int sy);

void free_grid(grid g);

void seed_grid(grid g, float p);

void print_grid(FILE *f, grid g, bool cluster);

void grid_do_dfs(grid *g);

/* outline.c */

outline alloc_outline(int sx, int sy, int n_cluster);

void free_outline(outline o);

void send_outline(outline o, int dest);

outline recv_outline(int src);

outline merge_outlines_horiz(outline o1, outline o2);

outline merge_outlines_vert(outline o1, outline o2);

outline outline_from_grid(grid g);

void print_outline(FILE *f, outline o);

/* percolate.c */

outline percolate_mpi(char t, int sx, int sy, float p, int nodes, int rank);
