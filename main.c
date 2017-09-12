#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <time.h>

bool *gen_grid(int size, float p)
{
  srand(time(NULL));
  bool *grid = malloc(size * size * sizeof(bool));
  for (int i=0; i<size*size; i++) {
    grid[i] = ((float)rand() / (float)RAND_MAX) < p;
  }
  return grid;
}

int get_edge_ind(int size, int side, int n) {
  // side: 0=top, 1=bottom, 2=left, 3=right
  switch (side) {
    case 0: return n;
    case 1: return size * (size - 1) + n;
    case 2: return n * size;
    case 3: return (size - 1) + n * size;
  }
  return 0;
}

void find_clusters(bool *grid, int sz, int *edges)
{
  int *grid_c = calloc(sz * sz, sizeof(int));
  int *stack = malloc(sz * sz * sizeof(int));
  int cluster = 1;
  for (int side=0; side<4; side++) {
    for (int n=0; n<sz; n++) {
      int start_point = get_edge_ind(sz, side, n);
      if (grid_c[start_point] != 0) continue;
      int ss = 1;
      stack[0] = start_point;
      while (ss > 0) {
        ss--;
        int curr = stack[ss];
        if (grid_c[curr] != 0) continue;
        if (!grid[curr]) continue;
        grid_c[curr] = cluster;
        if (curr >= sz && grid[curr - sz]) {
          stack[ss++] = curr - sz;
        }
        if (curr < sz * (sz - 1) && grid[curr + sz]) {
          stack[ss++] = curr + sz;
        }
        if (curr % sz > 0 && grid[curr - 1]) {
          stack[ss++] = curr - 1;
        }
        if (curr % sz < sz - 1 && grid[curr + 1]) {
          stack[ss++] = curr + 1;
        }
      }
      cluster++;
    }
  }
  for (int side=0; side<4; side++) {
    for (int n=0; n<sz; n++) {
      int ind = get_edge_ind(sz, side, n);
      edges[side * sz + n] = grid_c[ind];
    }
  }
  free(grid_c);
  free(stack);
  /* for (int i=0; i<sz; i++) {
    for (int j=0; j<sz; j++) {
      printf("%2d ", grid_c[i*sz + j]);
    }
    printf("\n");
  } */
}

int main(int argc, char *argv[])
{
  // Get parameters
  int size = atoi(argv[1]);
  float p = atof(argv[2]);
  printf("Percolating a %i x %i grid with p=%f\n", size, size, p);

  // Generate the grid
  bool *grid = gen_grid(size, p);

  // Find the clusters
  int *edges = malloc(4 * size * sizeof(int));
  find_clusters(grid, size, edges);

  /* Print grid
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      if (grid[i*size + j]) printf("*"); else printf("O");
    }
    printf("\n");
  } */

  /* Print edge clusters
  for (int i=0; i<size *4; i++) {
    printf("%3d ", edges[i]);
    if ((i + 1) % size == 0) printf("\n");
  } */

  return 0;
}