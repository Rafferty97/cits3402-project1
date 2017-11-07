#include "header.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

void man()
{
  printf("Usage:\n");
  printf("  percolate <infile>\n");
  printf("  percolate (s|b) <size> <p> [<iterations>]\n\n");
}

int main(int argc, char *argv[])
{
  double total, start, end, lap;
  total = 0;
  lap = 0;
  // Initialise MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  // Parse command line arguments
  if (argc < 4 || argc > 5) {
    if (mpi_rank > 0) return EXIT_FAILURE;
    fprintf(stderr, "Wrong number of arguments.\n\n");
    man();
    return EXIT_FAILURE;
  }
  char type = argv[1][0];
  int size = atoi(argv[2]);
  float p = atof(argv[3]);
  int iter = 1;
  if (argc >= 5) iter = atoi(argv[4]);
  if (size <= 0 || p < 0 || p > 1 || iter < 1) {
    if (mpi_rank > 0) return EXIT_FAILURE;
    fprintf(stderr, "Incorrect arguments supplied.\n\n");
    man();
    return EXIT_FAILURE;
  }
  // Run the percolation and print results
  if (mpi_rank == 0) {
    printf("Percolating with %i nodes.\n", mpi_size);
    printf("Percolating... 0 / %i iterations.", iter);
    fflush(stdout);
  }
  for (int i = 0; i < iter; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    percolate_mpi(type, size, size, p, mpi_size, mpi_rank);
    end = MPI_Wtime();
    total += end - start;
    lap -= end - start;
    if (lap < 0) {
      lap = 0.2;
      printf("\33[2K\rPercolating... %i / %i iterations.", i + 1, iter);
      fflush(stdout);
    }
  }
  if (mpi_rank == 0) {
    printf("\nParameters: %c %i %f\nIterations: %i\nTotal time: %f\n\n", type, size, p, iter, total);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
