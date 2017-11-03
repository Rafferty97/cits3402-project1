#include "header.h"
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

#define SEED 1000

outline null_outline;

outline percolate(char t, int sx, int sy, float p, int rank)
{
  char filename[256];
  sprintf(filename, "rank%i", rank);
  FILE *f = fopen(filename, "w");
  fprintf(f, "Percolating a %i x %i %c grid, p=%f.\n\n", sx, sy, t, p);
  grid g = alloc_grid(t, sx, sy);
  seed_grid(g, p, SEED + 763 * rank);
  print_grid(f, g, false);
  fprintf(f, "\n");
  grid_do_dfs(&g);
  print_grid(f, g, true);
  fprintf(f, "\n");
  outline o = outline_from_grid(g);
  print_outline(f, o);
  fprintf(f, "\n");
  fclose(f);
  return o;
}

outline percolate_mpi(char t, int sx, int sy, float p, int nodes, int rank, int rank_start)
{
  if (nodes > 1) {
    // Split the lattice in half
    bool horiz = sx > sy;
    int dim = horiz ? sx : sy;
    int nodes_l = nodes / 2;
    int nodes_r = nodes - nodes_l;
    int dim_l = (dim * nodes_l) / nodes;
    int dim_r = dim - dim_l;
    // Determine if this node is in the left or right partition
    if (rank < nodes_l) {
      int sx_l = horiz ? dim_l : sx;
      int sy_l = horiz ? sy : dim_l;
      outline l = percolate_mpi(t, sx_l, sy_l, p, nodes_l, rank, rank_start);
      if (rank != 0) return l;
      outline r = recv_outline(rank_start + nodes_l);
      printf("Recieved grid %ix%i from %i to %i\n", r.sx, r.sy, rank_start + nodes_l, rank_start);
      if (horiz) {
        return merge_outlines_horiz(l, r);
      } else {
        return merge_outlines_vert(l, r);
      }
    } else {
      rank -= nodes_l;
      rank_start += nodes_l;
      int sy_r = horiz ? sy : dim_r;
      int sx_r = horiz ? dim_r : sx;
      outline r = percolate_mpi(t, sx_r, sy_r, p, nodes_r, rank, rank_start);
      if (rank == 0) {
        printf("Sending grid %ix%i from %i to %i\n", r.sx, r.sy, rank_start, rank_start - nodes_l);
        send_outline(r, rank_start - nodes_l);
      }
      return r;
    }
  } else {
    return percolate(t, sx, sy, p, rank_start + rank);
  }
}
