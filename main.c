#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "site.h"

int main(int argc, char *argv[])
{
  percolation_results results;
  site_percolation(4096, 0.6, atoi(argv[1]), &results);
  print_results(results);
}
