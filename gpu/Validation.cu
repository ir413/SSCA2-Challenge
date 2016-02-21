#include <stdio.h>

#include "Validation.h"


/***************************************************************
 * Validation LN and PHJK (execution times on 4-threads i7-4790)
 * This is a table of reference outputs, but since the algorithm
 * is non-deterministic they are only representative.
 * You will need to add to this table if you run for larger
 * SCALE.
 **************************************************************/
bool isValid(Configuration *config, double *bc)
{
  int i = 0;
  double maxBC = bc[0], minBC = bc[0], sum = 0., avg = 0.;
  long int intMaxBC = 0, intMinBC = 0, intAvg = 0;

  #define NResults 50
  long int prestoredResults[NResults][3] = { 
  /* SCALE=1 */     { 0, 0, 0 }, 
  /* SCALE=2 */     { 0, 2500, 1250 }, 
  /* SCALE=3 */     { 0, 30999, 6249 }, 
  /* SCALE=4 */     { 0, 48028, 14812 }, 
  /* SCALE=5 */     { 0, 133238, 44000 }, 
  /* SCALE=6 */     { 0, 368702, 115578 }, 
  /* SCALE=7 */     { 0, 1790304, 257476 }, 
  /* SCALE=8 */     { 0, 5326002, 547863 }, 
  /* SCALE=9 */     { 0, 21598987, 1237150 }, 
  /* SCALE=10 */     { 0, 40382717, 2415442 }, 
  /* SCALE=11 */     { 0, 77220870, 2828273 }, 
  /* SCALE=12 */     { 0, 155192044, 3002415}, 
  /* SCALE=13 */     { 0, 237379503, 3098723}, 
  /* SCALE=14 */     { 0, 375125133, 3250859 }, 
  /* SCALE=15 */     { 0, 701218499, 3277281 }, 
  /* SCALE=16 */     { 0, 1158382011, 3338237 }, 
  /* SCALE=17 */     { 0, 1791767763, 3408419 }, 
  /* SCALE=18 */     { 0, 3341777600, 3491043 }, 
  /* SCALE=19 */     { 0, 5511190429, 3522221  }, /* ca.100 seconds */
  /* SCALE=19 */     /* { 0, 5321333625, 3523534  }, another run */
  /* SCALE=20 */     { 0, 11269213138, 3577907 }, /* ca.200 seconds */
  /* SCALE=20 */     /* { 0, 10994335868, 3553423 }, another run */
  /* SCALE=21 */     { 0, 15593566622, 3556022 } /* 2000 seconds */
  };
  
  fprintf(stderr, "N=%d, SCALE=%d\n", config->n, config->scale); 
  for (i = 0; i < config->n; i++)
  {
      if (bc[i] > maxBC)
      {
        maxBC = bc[i];
      }
      if (bc[i] < minBC)
      {
        minBC = bc[i];
      }
      sum += bc[i];
  }    
  avg = sum / config->n; 

  // It is easier to convert them in integer to do the comparison
  // and store the array in prestoredResults.
  intMaxBC = maxBC * 1000;
  intMinBC = minBC * 1000;
  intAvg = avg * 1000;

  fprintf(stderr, "min=%ld, max=%ld, avg=%ld\n", intMinBC, intMaxBC, intAvg); 
  fprintf(stderr, "This result should match:\n"); 
  fprintf(
      stderr,
      "min=%ld, max=%ld, avg=%ld\n", 
      prestoredResults[config->scale - 1][0],
      prestoredResults[config->scale - 1][1],
      prestoredResults[config->scale - 1][2]); 
  fprintf(stderr, "%ld, %ld, %ld\n", intMinBC, intMaxBC, intAvg); 

  /* PHJK: the betweenness centrality algorithm is slightly non-deterministic
   * so we check that the min, max and avg are within reasonable bounds.
   * This might need adjusting if you manage very large problems.
   */
  if (labs(intMinBC - prestoredResults[config->scale - 1][0]) <= 1 + (intMinBC / 10) && 
      labs(intMaxBC - prestoredResults[config->scale - 1][1]) <= 1 + (intMaxBC / 5) && 
      labs(intAvg   - prestoredResults[config->scale - 1][2]) <= 1 + (intAvg / 20))
  {
    return true; 
  } else {
    return false; 
  }
}
