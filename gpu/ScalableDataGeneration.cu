#include "sprng.h"

#include "ScalableDataGeneration.h"


/**
 * Port of the SSCA2v2.2 genScalData implementation. 
 */
void generateScalableData(Configuration *config, int *permV, GraphSDG *tuples)
{
  int i, j, u, v, step, tmpVal;
  double av, bv, cv, dv, p, s, var;

  tuples->n = config->n;
  tuples->m = config->m;

  // sprng seed.
  int seed = 2387;
  // Initialize a random number stream.
  int *stream = init_sprng(0, 0, 1, seed, SPRNG_DEFAULT);

  // Generate edges.
  for (i = 0; i < tuples->m; ++i)
  {
    u = 1;
    v = 1;
    step = tuples->n / 2;

    av = config->a;
    bv = config->b;
    cv = config->c;
    dv = config->d;

    p = sprng(stream);

    if (p < av)
    {
      // Do nothing.
    }
    else if ((p >= av) && (p < av + bv))
    {
      v += step;
    }
    else if ((p >= av + bv) && (p < av + bv + cv))
    {
      u += step;
    }
    else
    {
      u += step;
      v += step;
    }

    for (j = 1; j < config->scale; ++j)
    {
      step = step / 2;

      // Vary a, b, c, and d up to 10%.
      var = 0.1;
      av *= 0.95 + var * sprng(stream);
      bv *= 0.95 + var * sprng(stream);
      cv *= 0.95 + var * sprng(stream);
      dv *= 0.95 + var * sprng(stream);

      s = av + bv + cv + dv;
      av = av / s;
      bv = bv / s;
      cv = cv / s;
      dv = dv / s;

      // Choose partition.
      p = sprng(stream);

      if (p < av)
      {
        // Do nothing.
      }
      else if ((p >= av) && (p < av + bv))
      {
        v += step;
      }
      else if ((p >= av + bv) && (p < av + bv + cv))
      {
        u += step;
      }
      else
      {
        u += step;
        v += step;
      }
    }

    tuples->startVertex[i] = u - 1;
    tuples->endVertex[i] = v - 1;
  }

  // Generate vertex id permutations.
  for (i = 0; i < tuples->n; ++i)
  {
    permV[i] = i;
  }


  for (i = 0; i < tuples->n; ++i)
  {
    j = tuples->n * sprng(stream);

    if (i != j)
    {
      tmpVal = permV[i];
      permV[i] = permV[j];
      permV[j] = tmpVal;
    }
  }

  // Permute the vertices.
  for (i = 0; i < tuples->m; ++i)
  {
    tuples->startVertex[i] = permV[tuples->startVertex[i]];
    tuples->endVertex[i] = permV[tuples->endVertex[i]];
  }

  // Generate edge weights.
  for (i = 0; i < tuples->m; ++i)
  {
    tuples->weight[i] = 1 + config->maxIntWeight * sprng(stream);
  }
}

void generatePermutation(int n, int *permutation)
{
  int seed = 2387;
  int tid = 0;
  int nthreads = 1;

  // Initialize RNG stream.
  int *stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);

  for (int i = 0; i < n; ++i)
  {
    permutation[i] = i;
  }

  for (int i = 0; i < n; ++i)
  {
    int j = n * sprng(stream);

    if (i != j)
    {
      int k = permutation[i];
      permutation[i] = permutation[j];
      permutation[j] = k;
    }
  }
}

void printTuples(FILE *stream, GraphSDG *tuples)
{
  fprintf(stream, "Tuples:\n");

  for (int i = 0; i < tuples->m; ++i)
  {
    fprintf(
        stream,
        "%d -> %d : %d\n",
        tuples->startVertex[i],
        tuples->endVertex[i],
        tuples->weight[i]);
  }
}
