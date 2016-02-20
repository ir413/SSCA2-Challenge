#include <sys/time.h>

#include "Timer.h"


double getSeconds()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (double) (tp.tv_sec + ((1e-6) * tp.tv_usec));
}
