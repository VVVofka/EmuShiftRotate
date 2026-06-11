#include "NoShared/rotate.h"
#include <array>
#include <assert.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

using std::array;
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
// for debug
ivector vfld_base0;        // single variable in cuda
vector<ivector> vsrcshmem; // sources

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
  vector<array<int2, 4>> vshared_bounds(gridDim.x *
                                        gridDim.y); // for get min max
  vsrcshmem.resize(vshared.size());                 // sources
  vfld_base0.resize(gridDim.x * gridDim.y);         // single variable in cuda

  constexpr int szshall = wszshared * wszshared;
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++)
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      auto idx = blockIdx.y * gridDim.x + blockIdx.x;
      vshared[idx].resize(szshall, 0ULL);
      vsrcshmem[idx].resize(szshall, {INT_MAX, INT_MAX});
    }

  uint64_t *field = vfield.data();

  assert(wsz0 * wsz0 == vsubdata.size());
  const uint64_t *subdata = vsubdata.data(); // as kernel

  // find base of block in field for s_sub[0]
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++) {
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      int idblock = int(blockIdx.y * gridDim.x + blockIdx.x);
      u64vector &s_sub = vshared[idblock]; // current shared memory for block
      array<int2, 4> &s_bounds = vshared_bounds[idblock];

      // find corners of the block
      for(int corner = 0; corner < 4; corner++) {
        int x = (blockIdx.x * wszblock + (corner & 1) * wszblock) * 8 - 1;
        int y = (blockIdx.y * wszblock + (corner >> 1) * wszblock) * 8 - 1;

        int2 c = {x - hszfld, y - hszfld};
        int2 shc = {c.x - shift.x, c.y - shift.y};
        int2 shcwr = wrap_toroid0(shc, szfld);
        s_bounds[corner] = rotate_device(shcwr, d_rotates);
      }
      // synchronize threads in cuda  __syncthreads(); there

      // get min max of the block
      int bound_min = s_bounds[0], bound_max = s_bounds[0];
      for(int corner = 1; corner < 4; corner++) {
        bound_min = min(bound_min, s_bounds[corner]);
        bound_max = max(bound_max, s_bounds[corner]);
      }
      assert(bound_max.x - bound_min.x < szshared);
      assert(bound_max.y - bound_min.y < szshared);

      int2 wbasemin = floor8(wbound_min) / 8;
      int2 wbasemax = floor8(wbound_max) / 8; 
      for(int wy = wbasemin.y; wy <= wbasemax.y; wy++) {
        for(int wx = wbasemin.x; wx <= wbasemax.x; wx++) {
          // TODO: stop
        }
      }
      int2 base = {basemin.x - szblock, basemin.y - szblock};
      int2 base_wr = wrap_toroid0(base, szfld);

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

      // DEBUG: Check if shared memory covers all possible rotc values
      // For rotc in [-hsz0, hsz0), we need shr in [-szblock, szblock)
      // But wrap_relative can produce values in [-sz0, sz0) which is much
      // larger
      int2 wbasec = {floor8(basec.x) / 8, floor8(basec.y) / 8};

      // For debugging: print block 2 base
      if(idblock == 2) {
        printf("DEBUG block 2: wcob:(%d,%d) cobc:(%d,%d) cobrotc:(%d,%d) "
               "basec:(%d,%d) wbasec:(%d,%d)\n",
               wcob.x, wcob.y, cobc.x, cobc.y, cobrotc.x, cobrotc.y, basec.x,
               basec.y, wbasec.x, wbasec.y);
      }

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
  // dump_base(vfld_base0);
  //  dump_shmem(vshared);
  // dump_src_shmem(vsrcshmem);

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
            // NoShared uses: flds = (szfield + fld + shift) % szfield, then
            // fldc = flds - szfield/2 Main uses: fld_shift = fld + shift, then
            // fldc = fld_shift - hszfld These are equivalent for shift in
            // [-szfield/2, szfield/2)
            int2 fld_shift = {fld.x + shift.x, fld.y + shift.y};
            int2 fldc = {fld_shift.x - hszfld, fld_shift.y - hszfld};
            fldc = wrap_toroid0(fldc, szfld);

            int2 rotc = rotate_device(fldc, d_rotates);

            bool err_dump = 0 && blockIdx.x == 2 && blockIdx.y == 1 &&
                            threadIdx.x == 1 && threadIdx.y == 1 && nbit == 7;
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

            // Calculate shr_raw = rotc - base (without wrap)
            // If shr_raw is in the shared memory range [-szshared, szshared),
            // use it directly If shr_raw wrapped around (via wrap_relative), it
            // means the value is outside our shared memory region
            int2 shr_raw = {rotc.x - base.x,
                            rotc.y - base.y}; // before wrap_relative

            // Check if shr_raw is in the shared memory range BEFORE wrapping
            // shr_raw should be in [-szshared, szshared) for it to be covered
            // by shared memory
            if(shr_raw.x < -szshared || shr_raw.x >= szshared ||
               shr_raw.y < -szshared || shr_raw.y >= szshared) {
              // Value is outside shared memory region - skip
              continue;
            }

            // Now shr_raw is in range, convert to positive coordinates for
            // Morton encoding
            int2 shr = {shr_raw.x + szshared, shr_raw.y + szshared};
            int2 wshr = {shr.x / 8, shr.y / 8};

            bool wshr_check = wshr.x < 0 || wshr.y < 0 || wshr.x >= wszshared ||
                              wshr.y >= wszshared;
            if(wshr_check && (rotc.x >= -hsz0 && rotc.y >= -hsz0 &&
                              rotc.x < hsz0 && rotc.y < hsz0)) {
              // This should NOT happen - rotc is in bounds but wshr is out of
              // shared memory range Print detailed error for block 2
              int blk_id = blockIdx.y * gridDim.x + blockIdx.x;
              if(blk_id == 2) {
                printf("ERROR: blk_id=%d wshr out of bounds! rotc:(%d,%d) "
                       "base:(%d,%d) shr:(%d,%d) wshr:(%d,%d) wszshared:%d\n",
                       blk_id, rotc.x, rotc.y, base.x, base.y, shr.x, shr.y,
                       wshr.x, wshr.y, wszshared);
              }
              continue;
            }
            if(err_dump)
              printf("shr:%d %d  wshr:%d %d (%d)\n", shr.x, shr.y, wshr.x,
                     wshr.y, wszshared);

            if(wshr.x < 0 || wshr.y < 0 || wshr.x >= wszshared ||
               wshr.y >= wszshared)
              continue;

            if(err_dump)
              printf("wshr check Ok\n");

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
  // dump(vfield);
  return vshared;
} // --------------------------------------------------------------------------

int emu(int sz0, int2 shift, float angle, float kfill = 1.0f,
        bool verbose = true) {
  if(verbose)
    printf("\nemu: sz0:%d  shift:%d %d  angle:%.2f  kfill:%.2f\n", sz0, shift.x,
           shift.y, angle, kfill);
  int ret = 0;
  int szfld = sz0 * 3 / 2;
  float2 d_rotates =
      get_d_rotates(-angle); // FIX: negate angle to match NoShared

  // DEBUG: Find the exact error location
  int id_first_err_debug = -1;

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
  int id_first_err = RotateShiftHost::check_raw(angle, shift, vsubdata, vfield);
  if(id_first_err == -1) {
    if(verbose)
      printf("test Ok\n");
    return 0;
  }
  if(sz0 <= 32) {
    dump_base(vfld_base0);
    dump_shmem(vshared);
    dump_src_shmem(vsrcshmem);
    dump(vfield);
  }
  return 1;
} // --------------------------------------------------------------------------
#include <string>
using namespace std;
int main() {
  // return emu(32, {20,30}, 0.0f, 0.5f, true);
  srand(999);
  int sz0 = 32;
  while(sz0 <= 1024) {
    int szfld = sz0 * 3 / 2;
    for(int attempt = 0; attempt < 3; attempt++) {
      int2 shift = {rand() % szfld - szfld / 2, rand() % szfld - szfld / 2};
      float angle = rand() % 90 - 45.0f;
      float kfill = (10 + rand() % 90) / 100.0f;
      if(emu(sz0, shift, angle, kfill, true))
        return 1;
    }
    printf("test sz0=%d Ok\n", sz0);
    sz0 *= 2;
  }
  printf("All tests Ok\n");
  return 0;
} // --------------------------------------------------------------------------
