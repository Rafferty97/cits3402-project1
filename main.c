#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>

#include "results.h"
#include "site.h"
#include "bond.h"

void percolate(char type, int size, void *data, float p, int iter, int threads, bool l)
{
  unsigned long seed = time(NULL);
  percolation_results results;
  char mode = 0;
  if (tolower(type) == 'b') mode += 1;
  if (threads > 1) mode += 2;
  float total_time = 0, avg_largest_cluster = 0;
  int perc_count = 0, n = 0;
  if (l) printf("Percolating...");
  if (data != NULL) {
    if (tolower(type) == 'b') {
      load_bond_grid(data);
    } else {
      load_site_grid(data);
    }
  }
  for (int i=0; i<iter; i++) {
    switch (mode) {
      case 0:
      site_percolation(size, p, seed, &results); break;
      case 1:
      bond_percolation(size, p, seed, &results); break;
      case 2:
      site_percolation_parallel(size, p, seed, &results, threads); break;
      case 3:
      bond_percolation_parallel(size, p, seed, &results, threads); break;
    }
    if (results.percolates) perc_count++;
    total_time += results.time_taken;
    avg_largest_cluster += results.largest_cluster;
    if (l) {
      n += size * size;
      if (n < 5000) break;
      n = 0;
      printf("\33[2K\rPercolating... %i / %i iterations.", i + 1, iter);
      fflush(stdout);
    }
    seed += 1923 * i;
  }
  if (l) printf("\n\n");
  avg_largest_cluster /= iter;
  if (l) {
    printf(
      "%s percolation, %i x %i grid, p = %f\n%i iteration%s, %i thread%s.\n\n",
      (tolower(type) == 'b' ? "Bond" : "Site"), size, size, p, iter,
      (iter == 1 ? "" : "s"), threads, (threads == 1 ? "" : "s")
    );
    printf(
      "Total time: %f seconds\nPercolating grids: %i / %i\nAvg. largest cluster size: %f\n\n",
      total_time, perc_count, iter, avg_largest_cluster
    );
  } else {
    printf(
      "%s\t%i\t%f\t%i\t%i\t%f\t%i\t%f\n", (tolower(type) == 'b' ? "bond" : "site"),
      size, p, iter, threads, total_time, perc_count, avg_largest_cluster
    );
  }
}

void interactive_mode()
{
  char *line = NULL;
  size_t linesz = 0;
  char *delims = " \t\n\r";
  while (getline(&line, &linesz, stdin) > 0) {
    char *token = strtok(line, delims);
    if (!token) break;
    if (strcmp(token, "h") == 0) {
      printf("Type\tSize\tp\tIter's\tThreads\tTotal Time\tPerc Count\tAvg. Largest Cluster\n");
      continue;
    }
    if (strcmp(token, "s") && strcmp(token, "b")) {
      printf("Percolation type must be s (site) or b (bond).\n");
      exit(EXIT_FAILURE);
    }
    char type = token[0];
    token = strtok(NULL, delims);
    int size = atoi(token);
    if (size <= 0) {
      printf("Lattice size must be a positive integer.\n");
      exit(EXIT_FAILURE);
    }
    token = strtok(NULL, delims);
    char *buffer = NULL;
    float p = 0;
    int threads = 1, iter = 1;
    if (strcmp(token, "d") == 0) {
      token = strtok(NULL, delims);
      if (token) {
        threads = atoi(token);
        if (threads == 0) threads = 1;
      }
      if (getline(&line, &linesz, stdin) <= 0) {
        printf("No lattice data provided.\n");
        exit(EXIT_FAILURE);
      }
      buffer = malloc(size * size);
      char *c = line;
      for (int i = 0; i < size * size; i++) {
        if (*c == '\0') break;
        buffer[i] = *c++ - '0';
      }
    } else {
      p = atof(token);
      token = strtok(NULL, delims);
      if (token) {
        threads = atoi(token);
        if (threads == 0) threads = 1;
        token = strtok(NULL, delims);
        if (token) {
          iter = atoi(token);
          if (iter == 0) iter = 1;
        }
      }
    }
    percolate(type, size, buffer, p, iter, threads, false);
    free(buffer);
  }
}

void man()
{
  printf("No manual as of yet.\n\n");
}

int main(int argc, char *argv[])
{
  if (argc == 1) {
    interactive_mode();
    return EXIT_SUCCESS;
  }
  if (argc == 2) {
    if (strcmp(argv[1], "man") == 0) {
      man();
      return EXIT_SUCCESS;
    }
  }
  if (argc >= 4) {
    char type = argv[1][0];
    int size = atoi(argv[2]);
    float p = atof(argv[3]);
    int threads = 1;
    int iter = 1;
    if (argc >= 5) threads = atoi(argv[4]);
    if (argc >= 6) iter = atoi(argv[5]);
    if (argc >= 7) {
      printf("Wrong number of arguments.\n\n");
      man();
      return EXIT_FAILURE;
    }
    if (size <= 0 || p < 0 || p > 1 || threads < 1 || iter < 1) {
      printf("Incorrect arguments supplied.\n\n");
      man();
      return EXIT_FAILURE;
    }
    percolate(type, size, NULL, p, iter, threads, true);
    return EXIT_SUCCESS;
  }
  printf("Wrong number of arguments.\n\n");
  man();
  return EXIT_FAILURE;
}