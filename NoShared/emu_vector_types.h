// types as in CUDA

struct int2 {
  int x, y;
};

struct float2 {
  float x, y;
};

struct dim3 {
  unsigned x, y, z;
  dim3(unsigned a = 1, unsigned b = 1, unsigned c = 1) { x = a, y = b, z = c; }
};
