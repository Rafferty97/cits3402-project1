#include "header.h"
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <sys/time.h>

struct timeval time_start, time_now;

outline percolate(char t, int sx, int sy, float p)
{
  // char filename[256];
  // sprintf(filename, "grid%i,%i-%ix%i", grid_ox, grid_oy, sx, sy);
  // FILE *f = fopen(filename, "w");
  // fprintf(stdout, "Percolating a %i x %i %c grid at %i, %i.\n", sx, sy, t, grid_ox, grid_oy);
  grid g = alloc_grid(t, sx, sy);
  seed_grid(g, p);
  gettimeofday(&time_now, NULL);
  printf("%i seeded: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
  // printf("%i, %i:\n", grid_ox, grid_oy);
  // print_grid(stdout, g, false);
  // printf("\n");
  grid_do_dfs(&g);
  gettimeofday(&time_now, NULL);
  printf("%i DFS done: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
  // print_grid(f, g, true);
  // fprintf(f, "\n");
  outline o = outline_from_grid(g);
  gettimeofday(&time_now, NULL);
  printf("%i outlined: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
  // print_outline(f, o);
  // fprintf(f, "\n");
  // fclose(f);
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
      gettimeofday(&time_now, NULL);
      printf("%i recv'd: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
      // printf("Recieved grid %ix%i from %i to %i\n", r.sx, r.sy, rank_start + nodes_l, rank_start);
      if (horiz) {
        return merge_outlines_horiz(l, r);
      } else {
        return merge_outlines_vert(l, r);
      }
    } else {
      rank -= nodes_l;
      rank_start += nodes_l;
      int sx_r = horiz ? dim_r : sx;
      int sy_r = horiz ? sy : dim_r;
      if (horiz) grid_ox += dim_l; else grid_oy += dim_l;
      outline r = percolate_mpi_rec(t, sx_r, sy_r, p, nodes_r, rank, rank_start);
      if (rank == 0) {
        // printf("Sending grid %ix%i from %i to %i\n", r.sx, r.sy, rank_start, rank_start - nodes_l);
        send_outline(r, rank_start - nodes_l);
        gettimeofday(&time_now, NULL);
        printf("%i sent: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
      }
      return r;
    }
  } else {
    return percolate(t, sx, sy, p);
  }
}

outline percolate_mpi(char t, int sx, int sy, float p, int nodes, int rank)
{
  gettimeofday(&time_start, NULL);
  printf("%i start: %li %i\n", mpi_rank, time_start.tv_sec, time_start.tv_usec);
  grid_tx = sx;
  grid_ox = 0;
  grid_oy = 0;
  outline out = percolate_mpi_rec(t, sx, sy, p, nodes, rank, 0);
  if (rank == 0) {
    gettimeofday(&time_now, NULL);
    printf("%i recv'd: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
    percolate_outline(&out);
    gettimeofday(&time_now, NULL);
    printf("%i done: %li\n", mpi_rank, 1000000 * (time_now.tv_sec - time_start.tv_sec) + (time_now.tv_usec - time_start.tv_usec));
  }
  return out;
}
