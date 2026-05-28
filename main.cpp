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
constexpr float half_sqrt2 = sqrt2 * 0.5f;

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
  if(val.x < -sz2)
    val.x += sz;
  if(val.x >= sz2)
    val.x -= sz;
  if(val.y < -sz2)
    val.y += sz;
  if(val.y >= sz2)
    val.y -= sz;
  return val;
} // --------------------------------------------------------------------------
int __float2int_rn(float f) {
  return int(roundf(f));
} // --------------------------------------------------------------------------
int2 rotate_device(int x, int y, float2 d_rotates) {
  // Step 1: Horizontal shear
  int idx1 = __float2int_rn(d_rotates.x * y);
  int2 fld0 = {x + idx1, y};

  // Step 2: Vertical shear
  int idy2 = __float2int_rn(d_rotates.y * fld0.x);
  int2 fld1 = {fld0.x, fld0.y + idy2};

  // Step 3: Horizontal shear
  int idx3 = __float2int_rn(d_rotates.x * fld1.y);
  int2 fld2 = {fld1.x + idx3, fld1.y};

  return fld2;
} // --------------------------------------------------------------------------
int2 rotate_device(int2 pos, float2 d_rotates) {
  return rotate_device(pos.x, pos.y, d_rotates);
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
template <int wszblock = 2>
vector<u64vector> push(const u64vector &vsubdata, u64vector &vfield, int2 shift,
                       float2 d_rotates) {
  const int szfld = (int)sqrt((double)vfield.size()) * 8; // sz0 * 3 / 2
  const int hszfld = szfld / 2;                           // szfld / 2
  const int wszfld = szfld / 8;                           // szfld / 8
  const int hwszfld = wszfld / 2;                         // szfld / 16
  const int sz0 = szfld * 2 / 3;                          // 2 ^ N
  const int hsz0 = sz0 / 2;                               // sz0 / 2
  const int wsz0 = sz0 / 8;                               // sz0 / 8
  const int hwsz0 = wsz0 / 2;                             // sz0 / 16
  const int shbound = hwsz0; // hwszfld;// + wszblock / 2;

  constexpr int szblock = wszblock * 8;
  constexpr int wszshared = wszblock * 2;
  constexpr int szshared = wszshared * 8;

  dim3 blockDim(wszblock, wszblock, 1u); // only debug - not fixed
  dim3 gridDim(wszfld / wszblock, wszfld / wszblock, 1u);
  dim3 blockIdx, threadIdx;

  // emulation of shared memory
  vector<u64vector> vshared(gridDim.x * gridDim.y);
  vector<ivector> vsrcshmem(vshared.size()); // sources

  constexpr int szshall = wszshared * wszshared;
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++)
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      auto idx = blockIdx.y * gridDim.x + blockIdx.x;
      vshared[idx].resize(szshall, 0ULL);
      vsrcshmem[idx].resize(szshall, {INT_MAX, INT_MAX});
    }
  ivector vfld_base0(gridDim.x * gridDim.y); // single variable in cuda

  uint64_t *field = vfield.data();

  assert(wsz0 * wsz0 == vsubdata.size());
  const uint64_t *subdata = vsubdata.data(); // as kernel

  // find base of block in field for s_sub[0]
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++) {
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      int idblock = int(blockIdx.y * gridDim.x + blockIdx.x);
      u64vector &s_sub = vshared[idblock]; // shared memory for block

      // center of block (cob) in words
      int2 wcob = {(int)blockIdx.x * wszblock + wszblock / 2,
                   (int)blockIdx.y * wszblock + wszblock / 2};

      int cobx = wcob.x * 8 + shift.x;
      int coby = wcob.y * 8 + shift.y;
      int2 cobc = {cobx - hszfld, coby - hszfld};
      cobc = wrap_toroid0(cobc, szfld); // qwen

      // qwen
      int2 cobrotc = rotate_device(cobc, d_rotates);
      int2 cobrotc_floor = floor8({cobrotc.x + 4, cobrotc.y + 4});
      int2 basec = {cobrotc_floor.x - szblock, cobrotc_floor.y - szblock};

      vfld_base0[idblock] = basec;

      int2 wbasec = {floor8(basec.x) / 8, floor8(basec.y) / 8};

      for(threadIdx.y = 0; threadIdx.y < blockDim.y; threadIdx.y++) {
        for(threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
          // fill shared memory. 4 values per thread
          for(int j = 0; j < 4; j++) {
            int2 wshared = {int(threadIdx.x) * 2 + j % 2,
                            int(threadIdx.y) * 2 + j / 2};
            int2 wsubc_cart = {wbasec.x + wshared.x, wbasec.y + wshared.y};

            // toroidal wrap в координатах subdata относительно центра
            wsubc_cart = wrap_toroid0(wsubc_cart, sz0 / 8);

            int2 wsub_cart = {wsubc_cart.x + hwsz0, wsubc_cart.y + hwsz0};
            auto idw_sub_morton = morton_encode(wsub_cart);
            int idshared = wshared.y * wszshared + wshared.x;
            s_sub[idshared] = subdata[idw_sub_morton];
            vsrcshmem[idblock][idshared] = wsub_cart;
          } // for(int j = 0; j < 4; j++)
        } // threadIdx.x
      } // threadIdx.y
    } // blockIdx.x
  } // blockIdx.y
  dump_base(vfld_base0);
  // dump_shmem(vshared);
  dump_src_shmem(vsrcshmem);

  // copy shared memory to global memory
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++) {
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      const int idblock = int(blockIdx.y * gridDim.x + blockIdx.x);
      int2 &base = vfld_base0[idblock];
      u64vector &s_sub = vshared[idblock];
      for(threadIdx.y = 0; threadIdx.y < blockDim.y; threadIdx.y++) {
        for(threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
          int2 w = {(int)blockIdx.x * wszblock + (int)threadIdx.x,
                    (int)blockIdx.y * wszblock + (int)threadIdx.y};
          int idw = w.y * wszfld + w.x;
          uint64_t tile_field =
              field[idw]; // if don't overwrite bit then old value will be used
          int2 id_bit0 = {w.x * 8, w.y * 8};
          // #pragma unroll
          for(int nbit = 0; nbit < 64; ++nbit) {
            int2 bit = {nbit & 7, nbit >> 3};

            int2 fld = {id_bit0.x + bit.x, id_bit0.y + bit.y};
            int2 fld_shift = {fld.x + shift.x, fld.y + shift.y};
            int2 wcob = {(int)blockIdx.x * wszblock + wszblock / 2,
                         (int)blockIdx.y * wszblock + wszblock / 2};

            int2 cob_shift = {wcob.x * 8 + shift.x, wcob.y * 8 + shift.y};
            int2 cobc = {cob_shift.x - hszfld, cob_shift.y - hszfld};
            cobc = wrap_toroid0(cobc, szfld);
            int2 fldc_raw = {fld_shift.x - hszfld, fld_shift.y - hszfld};

            // wrap относительно центра блока
            int2 df = {fldc_raw.x - cobc.x, fldc_raw.y - cobc.y};
            df = wrap_toroid0(df, szfld);

            int2 fldc = {cobc.x + df.x, cobc.y + df.y};
            int2 rotc = rotate_device(fldc, d_rotates);
            
            bool err_dump = blockIdx.x == 2 && blockIdx.y == 1 && threadIdx.x == 1 &&
               threadIdx.y == 1 && nbit == 7;
            if(err_dump) {
              printf("blockIdx:%d %d threadIdx:%d %d nbit:%d\n", blockIdx.x,
                     blockIdx.y, threadIdx.x, threadIdx.y, nbit);
              printf("fld:%d %d  fld_shift:%d %d  fldc:%d %d\nrotc:%d %d "
                     "[%d,%d)\n",
                     fld.x, fld.y, fld_shift.x, fld_shift.y, fldc.x, fldc.y,
                     rotc.x, rotc.y, -hsz0, hsz0);
            }
            // rotc — проекция на subdata в центрированных координатах.
            // Если вне subdata, ничего не пишем.
            if(rotc.x < -hsz0 || rotc.y < -hsz0 || rotc.x >= hsz0 ||
               rotc.y >= hsz0)
              continue;
            int2 shr = {rotc.x - base.x, rotc.y - base.y};
            int2 wshr = {shr.x / 8, shr.y / 8};

            if(err_dump) 
              printf("shr:%d %d  wshr:%d %d (%d)\n", shr.x, shr.y, wshr.x,
                     wshr.y, wszshared);

            if(wshr.x < 0 || wshr.y < 0 || wshr.x >= wszshared ||
               wshr.y >= wszshared)
              continue;
            int idwshared = wshr.y * wszshared + wshr.x;
            uint64_t val_shared = s_sub[idwshared];
            uint32_t shift_bit = morton_encode(shr) & 63;
            uint32_t val_bit = uint32_t(val_shared >> shift_bit) & 1;
            replace_bit(tile_field, nbit, val_bit);
            if(err_dump) 
              printf("idwshared:%d  val_shared:%zX  shift_bit:%d  val_bit:%d\n",
                     idwshared, val_shared, shift_bit, val_bit);
          }
          field[idw] = tile_field;
        } // threadIdx.x
      } // threadIdx.y
    } // blockIdx.x
  } // blockIdx.y
  dump(vfield);
  return vshared;
} // --------------------------------------------------------------------------

int emu(int sz0, int2 shift, float angle) {
  printf("\nemu: sz0:%d  shift:%d %d  angle:%.2f\n", sz0, shift.x, shift.y,
         angle);
  int ret = 0;
  int szfld = sz0 * 3 / 2;
  float2 d_rotates = get_d_rotates(angle);
  u64vector vsubdata(sz0 * sz0 / 64, ~0ULL); // fill with ones
  // vsubdata[5] = 1ULL << 21;
  // vsubdata[10] = (1ULL << 42) | (1ULL << 41);
  // vsubdata[15] = (1ULL << 63) | (1ULL << 62) | (1ULL << 61);
  u64vector vfield(szfld * szfld / 64, 0ULL); // fill with zeros
  constexpr int wszblock = 2;                 // for debug, default is 16
  auto vshared = push<wszblock>(vsubdata, vfield, shift, d_rotates);
  int sumsub = sum1(vsubdata);
  int sumfield = sum1(vfield);

  if(sumsub != sumfield) {
    printf("Error: sumsub=%d != sumfield=%d\n", sumsub, sumfield);
    if(sz0 <= 32) {
      // dump_shmem(vshared);
      // dump(vfield);
    }
    return 1;
  }
  printf("test Ok\n");
  return 0;
} // --------------------------------------------------------------------------

int main() {
  if(emu(32, {5, 0}, 45.0f))
    return 1;
  printf("All tests Ok\n");
  return 0;
} // --------------------------------------------------------------------------
