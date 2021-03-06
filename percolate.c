#include "header.h"
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <sys/time.h>
#include <inttypes.h>

double t_alloc, t_seed, t_dfs, t_outline, t_stitch, t_percolate;

double lap()
{
  static double start = 0;
  double end = MPI_Wtime();
  double diff = end - start;
  start = end;
  return diff;
}

outline percolate(char t, int sx, int sy, float p)
{
  lap();
  grid g = alloc_grid(t, sx, sy);
  t_alloc = lap();
  seed_grid(g, p);
  t_seed = lap();
  grid_do_dfs(&g);
  t_dfs = lap();
  outline o = outline_from_grid(g);
  t_outline = lap();
  return o;
}

outline percolate_mpi_rec(char t, int sx, int sy, float p, int nodes, int rank, int rank_start)
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
      outline l = percolate_mpi_rec(t, sx_l, sy_l, p, nodes_l, rank, rank_start);
      if (rank != 0) return l;
      outline r = recv_outline(rank_start + nodes_l);
      lap();
      outline comb;
      if (horiz) {
        comb = merge_outlines_horiz(l, r);
      } else {
        comb = merge_outlines_vert(l, r);
      }
      t_stitch += lap();
      return comb;
    } else {
      rank -= nodes_l;
      rank_start += nodes_l;
      int sx_r = horiz ? dim_r : sx;
      int sy_r = horiz ? sy : dim_r;
      if (horiz) grid_ox += dim_l; else grid_oy += dim_l;
      outline r = percolate_mpi_rec(t, sx_r, sy_r, p, nodes_r, rank, rank_start);
      if (rank == 0) {
        send_outline(r, rank_start - nodes_l);
      }
      return r;
    }
  } else {
    t_stitch = 0;
    return percolate(t, sx, sy, p);
  }
}

outline percolate_mpi(char t, int sx, int sy, float p, int nodes, int rank)
{
  grid_tx = sx;
  grid_ox = 0;
  grid_oy = 0;
  outline out = percolate_mpi_rec(t, sx, sy, p, nodes, rank, 0);
  lap();
  if (rank == 0) {
    percolate_outline(&out);
  }
  t_percolate = lap();
  char fn[256];
  sprintf(fn, "rank%i", mpi_rank);
  FILE *f = fopen(fn, "w");
  fprintf(f, "alloc: %f\nseed: %f\ndfs: %f\noutline: %f\nstitch: %f\npercolate: %f\n", t_alloc, t_seed, t_dfs, t_outline, t_stitch, t_percolate);
  fclose(f);
  return out;
}
