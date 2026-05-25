#include <assert.h>
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
  printf("sz=%d\n", sz);
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
  const int szfld = (int)sqrt((double)vfield.size()) * 8;
  const int wszfld = szfld / 8;
  const int sz0 = szfld * 2 / 3;
  const int hsz0 = sz0 / 2;
  const int wsz0 = sz0 / 8;

  constexpr int szblock = wszblock * 8;
  constexpr int wszshared = wszblock * 2;
  constexpr int szshared = wszshared * 8;

  dim3 blockDim(wszblock, wszblock, 1u); // only debug - not fixed
  dim3 gridDim(wszfld / wszblock, wszfld / wszblock, 1u);
  dim3 blockIdx, threadIdx;

  // emulation of shared memory
  vector<u64vector> vshared(gridDim.x * gridDim.y);
  constexpr int szshall = wszshared * wszshared;
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++)
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++)
      vshared[blockIdx.y * gridDim.x + blockIdx.x].resize(szshall, 0ULL);

  vector<int2> vfld_base0(gridDim.x * gridDim.y); // single variable in cuda

  uint64_t *field = vfield.data();

  assert(wsz0 * wsz0 == vsubdata.size());
  const uint64_t *subdata = vsubdata.data(); // as kernel

  // find base of block in field for s_sub[0]
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++) {
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      int id_block = int(blockIdx.y * gridDim.x + blockIdx.x);
      u64vector &s_sub = vshared[id_block]; // shared memory for block

      // center of block (cob) in words
      int2 wcob = {(int)blockIdx.x * wszblock + wszblock / 2,
                   (int)blockIdx.y * wszblock + wszblock / 2};
      int cobx = (wcob.x * 8 + shift.x + szfld) % szfld;
      int coby = (wcob.y * 8 + shift.y + szfld) % szfld;
      int2 cobc = {cobx - szfld / 2, coby - szfld / 2};
      int2 rotc = rotate_device(cobc, d_rotates);
      int2 rot = {rotc.x + szfld / 2, rotc.y + szfld / 2};

      vfld_base0[id_block].x = (((rot.x + 4 + szfld) / 8) * 8) - szfld;
      vfld_base0[id_block].y = (((rot.y + 4 + szfld) / 8) * 8) - szfld;

      int2 wsub_cob = {(rotc.x + 4 + szfld / 2) / 8,
                       (rotc.y + 4 + szfld / 2) / 8};
      int2 wsub_down = {wsub_cob.x - wszblock, wsub_cob.y - wszblock};

      for(threadIdx.y = 0; threadIdx.y < blockDim.y; threadIdx.y++) {
        for(threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
          int2 w = {(int)blockIdx.x * wszblock + (int)threadIdx.x,
                    (int)blockIdx.y * wszblock + (int)threadIdx.y};
          int idw = w.y * wszfld + w.x; // id word

          // fill shared memory. 4 values per thread
          for(int j = 0; j < 4; j++) {
            int x = int(threadIdx.x) * 2 + j % 2;
            if(wsub_down.x + x < 0 || wsub_down.x + x >= wsz0)
              continue;

            int y = int(threadIdx.y) * 2 + j / 2;
            if(wsub_down.y + y < 0 || wsub_down.y + y >= wsz0)
              continue;

            auto idw_sub_morton = morton_encode({x, y});
            s_sub[y * wszshared + x] = subdata[idw_sub_morton];
          } // for(int j = 0; j < 4; j++) 
        } // threadIdx.x
      } // threadIdx.y
    } // blockIdx.x
  } // blockIdx.y
  dump_shmem(vshared);

  // copy shared memory to global memory
  for(blockIdx.y = 0; blockIdx.y < gridDim.y; blockIdx.y++) {
    for(blockIdx.x = 0; blockIdx.x < gridDim.x; blockIdx.x++) {
      int idblock = int(blockIdx.y * gridDim.x + blockIdx.x);
      u64vector &s_sub = vshared[idblock];
      for(threadIdx.y = 0; threadIdx.y < blockDim.y; threadIdx.y++) {
        for(threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
          int2 w = {(int)blockIdx.x * wszblock + (int)threadIdx.x,
                    (int)blockIdx.y * wszblock + (int)threadIdx.y};
          int idw = w.y * wszfld + w.x;
          uint64_t tile_field = field[idw];
          int2 id_bit = {w.x << 3, w.y << 3};
          for(int bit = 0; bit < 64; ++bit) {
            int2 local = {bit & 7, bit >> 3};
            int2 fld = {id_bit.x + local.x, id_bit.y + local.y};
            if(fld.x == 8 && fld.y == 8)
              printf("pause\n");

            int2 fld_shift = {fld.x + shift.x, fld.y + shift.y};
            int2 fldc = {fld_shift.x - szfld / 2, fld_shift.y - szfld / 2};
            int2 rotc = rotate_device(fldc, d_rotates);
            if(rotc.x >= -hsz0 && rotc.x < hsz0 && rotc.y >= -hsz0 &&
               rotc.y < hsz0) {
              int2 rot = {rotc.x + szfld / 2, rotc.y + szfld / 2};
              int2 shr = {rot.x - vfld_base0[idblock].x,
                          rot.y - vfld_base0[idblock].y};
              if(shr.x < 0 || shr.y < 0 || shr.x >= szshared ||
                 shr.y >= szshared) {
                // printf("shr:%d %d  b:%u %u t:%u %u  bit:%d(%d %d)  fld:%d
                // %d\n", shr.x,shr.y, bx, by, tx, ty, bit, bit & 7, bit >> 3,
                // fld.x,fld.y); continue;
              }
              int2 wshr = {shr.x / 8, shr.y / 8};
              int idwshared = wshr.y * wszshared + wshr.x;
              uint64_t val_shared = s_sub[idwshared];
              uint32_t nbit = morton_encode(shr) & 63;
              uint32_t val_bit = uint32_t(val_shared >> nbit) & 1;
              replace_bit(tile_field, bit, val_bit);
            }
          }
          field[idw] = tile_field;
        } // threadIdx.x
      } // threadIdx.y
    } // blockIdx.x
  } // blockIdx.y
  return vshared;
} // --------------------------------------------------------------------------

int emu(int sz0, int2 shift, float angle) {
  printf("\nemu: sz0:%d  shift:%d %d  angle:%.2f\n", sz0, shift.x, shift.y,
         angle);
  int ret = 0;
  int szfld = sz0 * 3 / 2;
  float2 d_rotates = get_d_rotates(angle);
  u64vector vsubdata(sz0 * sz0 / 64, ~0ULL);  // fill with ones
  u64vector vfield(szfld * szfld / 64, 0ULL); // fill with zeros
  constexpr int wszblock = 2;                 // for debug, default is 16
  auto sharedmem = push<wszblock>(vsubdata, vfield, shift, d_rotates);
  int sumsub = sum1(vsubdata);
  int sumfield = sum1(vfield);

  if(sumsub != sumfield) {
    printf("Error: sumsub=%d != sumfield=%d\n", sumsub, sumfield);
    if(sz0 <= 32)
      dump(vfield);
    return 1;
  }
  printf("test Ok\n");
  return 0;
} // --------------------------------------------------------------------------

int main() {
  if(emu(32, {0, 0}, 0.0f))
    return 1;
  return 0;
  if(emu(32, {1, 0}, 0.0f))
    return 2;
  if(emu(32, {0, 0}, 41.0f))
    return 3;
  if(emu(32, {1, 0}, 41.0f))
    return 4;
  if(emu(32, {12, -10}, -41.0f))
    return 5;
  printf("All tests Ok\n");
  return 0;
} // --------------------------------------------------------------------------
