#include <stdio.h>

#include "Configuration.h"


/**
 * Entry point.
 */
int main(int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf(stderr, "Usage: ./SSCA2 <SCALE>\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Configuring the Benchmark...\n");

  Configuration config;
  configure(atoi(argv[1]), &config);

  fprintf(stderr, "Scalable Data Generation...\n");
  // TODO.
}
