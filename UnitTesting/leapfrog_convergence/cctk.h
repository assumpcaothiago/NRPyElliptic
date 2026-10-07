/* Minimal standalone Cactus interface; numerical routines are not stubbed. */
#ifndef TEST_CCTK_H
#define TEST_CCTK_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef double CCTK_REAL;
typedef int CCTK_INT;
typedef struct { int unused; } cGH;
#define CCTK_ATTRIBUTE_UNUSED __attribute__((unused))
#define CCTK_EQUALS(a, b) (strcmp((a), (b)) == 0)
#define CCTK_INFO(...) ((void)0)
#define CCTK_VINFO(...) ((void)0)
/* Fatal Cactus errors terminate the process. _Exit also avoids leak reports
 * for allocations deliberately left behind by that fatal production path. */
#define CCTK_VERROR(...) do { fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); fflush(stderr); _Exit(2); } while (0)
#define CCTK_ERROR(message) CCTK_VERROR("%s", message)
#endif
