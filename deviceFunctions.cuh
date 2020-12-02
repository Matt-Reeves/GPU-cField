#ifndef HAS_DEVICE_FUNCTIONS_H
#define HAS_DEVICE_FUNCTIONS_H

#include <getCUDA.h>
__device__ double2 mulrc (double a, double2 b);
__device__ double modsq (double2 a);
__device__ double2 addc (double2 a, double2 b);
__device__ double2 mulc (double2 a, double2 b);
__device__ double absdiffc(double2 a, double2 b);

__forceinline__ __device__ double2 mulrc (double a, double2 b)
{
  double2 c;
  c.x = a*b.x;
  c.y = a*b.y;
  return c;
}

__forceinline__ __device__ double modsq (double2 a)
{
  double b;
  b = a.x*a.x + a.y*a.y;
  return b;
}

__forceinline__ __device__ double2 addc (double2 a, double2 b)
{
  double2 c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return c;
}

__forceinline__ __device__ double2 mulc (double2 a, double2 b)
{
  double2 c;
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.x * b.y + a.y * b.x;
  return c;
}

__forceinline__ __device__ double absdiffc(double2 a, double2 b)
{
  double c;
  c = hypot(a.x-b.x, a.y-b.y);
  return c;
}


#endif

