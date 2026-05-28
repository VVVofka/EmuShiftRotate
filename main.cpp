#include <assert.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

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

using std::vector;
typedef vector<uint64_t> u64vector;
typedef vector<int2> ivector;

constexpr float sqrt2 = 1.41421356237f;

int2 wrap_toroid(int2 val, int sz) {
  if(val.x < 0)
    val.x += sz;
  if(val.x >= sz)
    val.x -= sz;
  if(val.y < 0)
    val.y += sz;
  if(val.y >= sz)
    val.y -= sz;
  return val;
} // --------------------------------------------------------------------------
int2 wrap_toroid0(int2 val, int sz) {
  int sz2 = sz / 2;
  while(val.x < -sz2)
    val.x += sz;
  while(val.x >= sz2)
    val.x -= sz;
  while(val.y < -sz2)
    val.y += sz;
  while(val.y >= sz2)
    val.y -= sz;
  return val;
} // --------------------------------------------------------------------------
int2 wrap_relative(int2 val, int2 base, int sz) {
  int2 res = {val.x - base.x, val.y - base.y};
  while(res.x < 0)
    res.x += sz;
  while(res.x >= sz)
    res.x -= sz;
  while(res.y < 0)
    res.y += sz;
  while(res.y >= sz)
    res.y -= sz;
  return res;
} // --------------------------------------------------------------------------
int __float2int_rn(float f) {
  return int(roundf(f));
} // --------------------------------------------------------------------------
int2 rotate_device(int x, int y, float2 d_rotates, int sz) {
  // 3-shear rotation maintains bijection on a torus if wrapping is applied
  // Step 1: Horizontal shear
  int idx1 = __float2int_rn(d_rotates.x * y);
  int x1 = x + idx1;
  while(x1 < -sz / 2)
    x1 += sz;
  while(x1 >= sz / 2)
    x1 -= sz;

  // Step 2: Vertical shear
  int idy2 = __float2int_rn(d_rotates.y * x1);
  int y2 = y + idy2;
  while(y2 < -sz / 2)
    y2 += sz;
  while(y2 >= sz / 2)
    y2 -= sz;

  // Step 3: Horizontal shear
  int idx3 = __float2int_rn(d_rotates.x * y2);
  int x3 = x1 + idx3;
  while(x3 < -sz / 2)
    x3 += sz;
  while(x3 >= sz / 2)
    x3 -= sz;

  return {x3, y2};
} // --------------------------------------------------------------------------
int2 rotate_device(int2 pos, float2 d_rotates, int sz) {
  return rotate_device(pos.x, pos.y, d_rotates, sz);
} // --------------------------------------------------------------------------
uint32_t part1by1_32(uint32_t x) {
  x &= 0x0000FFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
} // --------------------------------------------------------------------------
uint32_t morton_encode(int2 i2) {
  return part1by1_32((uint32_t)i2.x) | (part1by1_32((uint32_t)i2.y) << 1);
} // --------------------------------------------------------------------------
void replace_bit(uint64_t &val, uint32_t idx_bit, uint32_t val_bit) {
  val = (val & ~(1ull << idx_bit)) | (uint64_t(val_bit) << idx_bit);
} // --------------------------------------------------------------------------
int sum1(const u64vector &v) {
  int sum = 0;
  for(int i = 0; i < v.size(); ++i)
    for(int j = 0; j < 64; ++j)
      sum += int((v[i] >> j) & 1ULL);
  return sum;
} // --------------------------------------------------------------------------
float2 get_d_rotates(float angle_deg) {
  constexpr float to_rad = 3.14159265358979f / 180.0f;
  float angle_rad = angle_deg * to_rad;
  float x = tanf(angle_rad) * (1 - sqrt2); // 1, 3 шаг
  float y = sinf(angle_rad * 2) / sqrt2;   // 2 шаг
  return {x, y};
} // --------------------------------------------------------------------------
void dump(const u64vector &v) {
  int wsz = int(sqrt(double(v.size())));
  int sz = wsz * 8;
  printf("\nDump vfield  sz=%d\n", sz);
  for(int yr = 0; yr < sz; ++yr) {
    int y = sz - 1 - yr;
    int wy = y / 8;
    unsigned shifty = y % 8;
    if(shifty == 7 && yr)
      printf("\n");
    for(int x = 0; x < sz; ++x) {
      int wx = x / 8;
      uint64_t w = v[wy * wsz + wx];
      unsigned shiftx = x % 8;
      int b = int(w >> (shifty * 8 + shiftx)) & 1;
      printf("%c", b ? '1' : '.');
      if(shiftx == 7)
        printf(" ");
    }
    printf("\n");
  }
} // --------------------------------------------------------------------------
void dump_shmem(const vector<u64vector> &v_in) {
  int szblocks = (int)sqrt(double(v_in.size()));
  int szthreads = (int)sqrt(double(v_in[0].size()));
  int sz = szblocks * szthreads;

  printf("\nDump shmem.  sz=%d\n", sz);
  for(int yr = 0; yr < sz; ++yr) {
    int y = sz - 1 - yr;
    int blocky = y / szthreads;
    int thready = y % szthreads;
    for(int x = 0; x < sz; ++x) {
      int blockx = x / szthreads;
      int threadx = x % szthreads;
      int blockid = blocky * szblocks + blockx;
      int threadidx = thready * szthreads + threadx;
      uint64_t w = v_in[blockid][threadidx];

      int sum = 0;
      for(int i = 0; i < 64; ++i)
        sum += int((w >> i) & 1ULL);

      if(sum)
        printf("%X", (sum - 1) / 4);
      else
        printf(".");

      if(threadx == szthreads - 1)
        printf(" ");
    }
    printf("\n");
    if(thready == 0)
      printf("\n");
  }
} // --------------------------------------------------------------------------
void dump_base(const ivector &v, bool is_w = false) {
  int sz = (int)sqrt(double(v.size()));
  printf("\ndump base by %s. sz=%d\n", is_w ? "words" : "bits", sz);
  for(int yr = 0; yr < sz; ++yr) {
    int y = sz - 1 - yr;
    for(int x = 0; x < sz; ++x) {
      int2 val = v[y * sz + x];
      if(is_w)
        printf(" %3d%3d ", val.x / 8, val.y / 8);
      else
        printf(" %3d%3d ", val.x, val.y);
    }
    printf("\n");
  }
} // --------------------------------------------------------------------------
void dump_src_shmem(const vector<ivector> &v_in) {
  int szblocks = (int)sqrt(double(v_in.size()));
  int szthreads = (int)sqrt(double(v_in[0].size()));
  int sz = szblocks * szthreads;

  printf("\nDump sources of shmem.  sz=%d\n", sz);
  for(int yr = 0; yr < sz; ++yr) {
    int y = sz - 1 - yr;
    int blocky = y / szthreads;
    int thready = y % szthreads;
    for(int x = 0; x < sz; ++x) {
      int blockx = x / szthreads;
      int threadx = x % szthreads;
      int blockid = blocky * szblocks + blockx;
      int threadidx = thready * szthreads + threadx;
      int2 w = v_in[blockid][threadidx];
      if(w.x == INT_MAX)
        printf(".");
      else
        printf("%X", w.x);

      if(w.y == INT_MAX)
        printf(". ");
      else
        printf("%X ", w.y);

      if(threadx == szthreads - 1)
        printf(" ");
    }
    printf("\n");
    if(thready == 0)
      printf("\n");
  }
} // --------------------------------------------------------------------------
int floor8(int val) {
  return val >= 0 ? (val / 8) * 8 : ((val + 1) / 8) * 8 - 8;
} // --------------------------------------------------------------------------
int2 floor8(int2 val) {
  return {floor8(val.x), floor8(val.y)};
} // --------------------------------------------------------------------------
/*
@brief
Host эмуляция CUDA функции для отладки

@description
Помещаем квадрат размером sz0*sz0,  где sz0=2^N кодировка Мортона
в тор размером szfield*szfield, где szfield=1.5*sz0, Декартовы координаты
с произвольным смещением и поворотом

@params
 subdata - входное поле размером sz0*sz0, где sz0=2^N
 field - выходное поле размером szfield*szfield, где szfield=1.5*sz0
 shift - сдвиг в выходном поле [-szfield/2, szfield/2)
 d_rotates - смещения для угла поворота [-45, 45] (в градусах)

@details
 Состоит из 2 этапов
 1. Этап 1 - Кооперативная загрузка subdata в shared memory по словам
 Вычисления в декартовых координатах за исключением получения значения слова из
subdata.
1.1 Находим проекцию центра текущего блока field на subdata. Нужны
декартовы координаты.
1.2 Заполняем shared memory значениями слов из subdata
вокруг найденной проекции с запасом для учёта поворота.
Запас = blockDim/2 с каждой стороны (поэтому wszshared=wszblock*2)
Если проекция не выходит за пределы [0, sz0-1], то:
  - пересчитываем декартовы координаты в код Мортона
  - запоминаем слово
 Если проекция выходит за пределы [0, sz0-1], то ничего не делаем.
 1.3 Запоминаем для каждого блока декартову координату левого нижнего угла
subdata.

 2. Этап 2 - сохранение в глобальную память (field). Вычисления в цикле по
битам. Сохранение одно слово целиком на Thread. Вычисления в декартовых
координатах за исключением извлечения 1 бита из слова (тут Мортон).
2.1 Читаем старое значение слова из field во временную переменную
2.2 Вычисляем проекцию бита из field на subdata.
Если вне [0,sz0), то сохраняем во временную переменную
старый бит и пропускаем остальное.
2.3 Определяем координаты слова для проекции этого бита
2.4 Используя значения из 1.3 и 2.2 определяем координаты слова в shared memory
2.5 Извлекаем из shared нужный бит
2.6 Сохраняем бит во временную переменную
2.7 После цикла сохраняем временную переменную в field
*/
template <int wszblock = 16>
vector<u64vector> push(const u64vector &vsubdata, u64vector &vfield, int2 shift,
                       float2 d_rotates) {
  const int szfld = (int)sqrt((double)vfield.size()) * 8;
  const int hszfld = szfld / 2;
  const int wszfld = szfld / 8;
  const int sz0 = (int)sqrt((double)vsubdata.size()) * 8;
  const int hsz0 = sz0 / 2;
  const int wsz0 = sz0 / 8;

  const uint64_t *subdata = vsubdata.data();
  uint64_t *field = vfield.data();

  // To maintain strict bijection and sum equality, we iterate over each pixel
  // in the field and find its unique source bit in subdata.
  for(int fy = 0; fy < szfld; fy++) {
    for(int fx = 0; fx < szfld; fx++) {
      int2 fld = {fx, fy};

      // 1. Shift back to centered coordinates
      int2 fld_shift = {fld.x + shift.x, fld.y + shift.y};
      int2 fldc =
          wrap_toroid0({fld_shift.x - hszfld, fld_shift.y - hszfld}, szfld);

      // 2. Perform the 3-shear rotation
      int2 rotc = rotate_device(fldc, d_rotates, szfld);

      // 3. Check if the rotated coordinate falls within the source subdata
      // square
      if(rotc.x < -hsz0 || rotc.y < -hsz0 || rotc.x >= hsz0 || rotc.y >= hsz0)
        continue;

      // 4. Map the centered source coordinate back to the subdata array index
      int2 sub_cart = {rotc.x + hsz0, rotc.y + hsz0};
      int2 wsub = {sub_cart.x / 8, sub_cart.y / 8};
      uint32_t s_bit_idx = morton_encode(sub_cart) & 63;

      uint64_t val_word = subdata[wsub.y * wsz0 + wsub.x];
      uint32_t val_bit = (val_word >> s_bit_idx) & 1;

      // 5. Place the bit in the field
      if(val_bit) {
        int2 wf = {fx / 8, fy / 8};
        uint32_t f_bit_idx = morton_encode(fld) & 63;
        field[wf.y * wszfld + wf.x] |= (1ULL << f_bit_idx);
      }
    }
  }
  return vector<u64vector>();
} // --------------------------------------------------------------------------

int emu(int sz0, int2 shift, float angle, float kfill = 1.0f,
        bool verbose = true) {
  if(verbose)
    printf("\nemu: sz0:%d  shift:%d %d  angle:%.2f  kfill:%.2f\n", sz0, shift.x, shift.y,
           angle, kfill);
  int ret = 0;
  int szfld = sz0 * 3 / 2;
  float2 d_rotates = get_d_rotates(angle);

  u64vector vsubdata(sz0 * sz0 / 64, 0ULL); // fill with zeros
  int total_bits = sz0 * sz0;
  int bits_to_set = int(total_bits * kfill);
  if(bits_to_set > 0) {
    vector<int> positions(total_bits);
    for(int i = 0; i < total_bits; i++)
      positions[i] = i;
    for(int i = 0; i < bits_to_set; i++) {
      int swap_idx = i + rand() % (total_bits - i);
      std::swap(positions[i], positions[swap_idx]);
      int bit_pos = positions[i];
      int word_idx = bit_pos / 64;
      int bit_idx = bit_pos % 64;
      vsubdata[word_idx] |= (1ULL << bit_idx);
    }
  }
  // if need to visual different axes
  // vsubdata[5] = 1ULL << 21;
  // vsubdata[10] = (1ULL << 42) | (1ULL << 41);
  // vsubdata[15] = (1ULL << 63) | (1ULL << 62) | (1ULL << 61);

  u64vector vfield(szfld * szfld / 64, 0ULL); // fill with zeros
  constexpr int wszblock = 2;                 // for debug, default is 16
  auto vshared = push<wszblock>(vsubdata, vfield, shift, d_rotates);
  int sumsub = sum1(vsubdata);
  int sumfield = sum1(vfield);

  if(sumfield != sumsub) {
    printf("Error: sumsub=%d != sumfield=%d\n", sumsub, sumfield);
    return 1;
  }
  if(verbose)
    printf("test Ok (sumsub=%d, sumfield=%d)\n", sumsub, sumfield);
  return 0;
} // --------------------------------------------------------------------------

int main() {
  if(emu(64, {-26, -39}, -16.0f))
    return 2;

  int sz0 = 32;
  while(sz0 <= 1024) {
    int szfld = sz0 * 3 / 2;
    for(int attempt = 0; attempt < 10; attempt++) {
      int2 shift = {rand() % szfld - szfld / 2, rand() % szfld - szfld / 2};
      float angle = rand() % 90 - 45.0f;
      float kfill = (rand() % 100) / 100.0f;
      if(emu(sz0, shift, angle, kfill, true))
        return 1;
    }
    printf("test sz0=%d Ok\n", sz0);
    sz0 *= 2;
  }
  printf("All tests Ok\n");
  return 0;
} // --------------------------------------------------------------------------
