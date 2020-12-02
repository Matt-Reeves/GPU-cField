#ifndef HAS_GRID_H
#define HAS_GRID_H

typedef struct ivec {
  int x;
  int y;
  int z;
} ivec;

typedef struct vec {
  double x;
  double y;
  double z;
} vec;

typedef struct grid {
  ivec n;
  vec delta;
  vec offset;
  vec dk;
} grid;

#endif
 
