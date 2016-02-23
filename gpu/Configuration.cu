#include "Configuration.h"

void configure(int scale, Configuration *config)
{
  // Apply the Binary Scaling Heuristic.
  config->scale = scale;
  config->n = 1 << scale;
  config->m = 8 * config->n;
  config->maxIntWeight = 1 << scale;
  config->subGraphPathLength = 3;

  if (scale < 10) {
    config->k4Approx = scale;
  } else {
    config->k4Approx = 10;
  }

  // Set R-Mat params.
  config->a = 0.55;
  config->b = 0.1;
  config->c = config->b;
  config->d = 0.25;
}
