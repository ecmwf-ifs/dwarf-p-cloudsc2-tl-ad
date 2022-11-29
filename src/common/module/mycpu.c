#define _GNU_SOURCE
#if defined(__APPLE__)
static int sched_getcpu() { return 0; }
#else
#include <sched.h>
#endif

/*
 * Find the core the thread belongs to
 */

int mycpu_ ()
{
  /* int sched_getcpu(void); */
  int cpu;
  cpu = sched_getcpu();
  return cpu;
}
int mycpu() { return mycpu_(); }
