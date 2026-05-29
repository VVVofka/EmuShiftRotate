#include "rotate.h"
#include <cassert>
#include <cmath>
#include <cstdio>

int2 RotateShiftHost::rotate(int2 sub, float2 k) {
  int idx1 = (int)roundf(k.x * sub.y);
  int2 fld0 = {sub.x + idx1, sub.y};

  // int idy2 = truncf(k.y * fld0.x);
  int idy2 = (int)roundf(k.y * fld0.x);
  int2 fld1 = {fld0.x, fld0.y + idy2};

  int idx3 = (int)roundf(k.x * fld1.y);
  int2 fld2 = {fld1.x + idx3, fld1.y};

  return fld2;
} // --------------------------------------------------------------------------

void RotateShiftHost::dump_field(const std::vector<int> &vfield, int2 shift) {
  int szfield = int(sqrt(static_cast<double>(vfield.size())));
  int sz0 = int(sqrt(static_cast<double>(vfield.size())) / 1.5);
  assert(std::abs(sqrt(static_cast<double>(sz0 * sz0)) - sz0) < 0.001);
  if(sz0 <= 16) {
    for(int yr = 0; yr < szfield; yr++) {
      int y = (2 * szfield - 1 - yr - shift.y) % szfield;
      for(int xr = 0; xr < szfield; xr++) {
        int x = (xr + shift.x) % szfield;
        int v = vfield[y * szfield + x];
        if(v == -1)
          printf(" ..");
        else
          printf(" %X%X", v & (sz0 - 1), (v / sz0) & (sz0 - 1));
      }
      printf("\n");
    }
  } else {
    for(int yr = 0; yr < szfield; yr++) {
      int y = (2 * szfield - 1 - yr - shift.y) % szfield;
      for(int xr = 0; xr < szfield; xr++) {
        int x = (xr + shift.x) % szfield;
        int v = vfield[y * szfield + x];
        printf("%c", v == -1 ? '.' : 'x');
      }
      printf("\n");
    }
  }
  printf("\n");
} // --------------------------------------------------------------------------

int RotateShiftHost::check(const std::vector<int> &vfield) {
  int sz0 = int(sqrt(static_cast<double>(vfield.size())) / 1.5);
  if(sz0 <= 0 || (sz0 & (sz0 - 1)) != 0) {
    printf("size error sz0=%d\n", sz0);
    return -1; // size error
  }
  int expected_count = sz0 * sz0;
  std::vector<bool> found(expected_count, false);
  int found_count = 0;

  for(int v : vfield) {
    if(v == -1)
      continue;

    if(v < 0 || v >= expected_count) {
      printf("value error - out of range: v=%d(%d)\n", v, expected_count);
      return -2; // value error - out of range
    }
    if(found[v]) {
      printf("value error - duplicate: v=%d\n", v);
      return -3; // value error - duplicate
    }
    found[v] = true;
    found_count++;
  }

  if(found_count != expected_count) {
    printf("expected count error found_count=%d expected_count=%d\n",
           found_count, expected_count);
    return -4; //  expected count error
  }
  return 0; // no errors
} // --------------------------------------------------------------------------

static uint32_t morton_decode(uint32_t code) {
  code &= 0x55555555;
  code = (code ^ (code >> 1)) & 0x33333333;
  code = (code ^ (code >> 2)) & 0x0F0F0F0F;
  code = (code ^ (code >> 4)) & 0x00FF00FF;
  code = (code ^ (code >> 8)) & 0x0000FFFF;
  return code;
} // --------------------------------------------------------------------------
static int2 morton_decode_xy(uint32_t morton) {
  return {(int)morton_decode(morton), (int)morton_decode(morton >> 1)};
} // --------------------------------------------------------------------------
std::vector<int>
RotateShiftHost::convert_raw_morton(const std::vector<uint64_t> &v) {
  assert(v.size());
  int sz0 = 8 * int(sqrt(static_cast<double>(v.size())));
  std::vector<int> vret(sz0 * sz0, -1);
  for(int z = 0; z < vret.size(); z++) {
    int2 xy = morton_decode_xy(z);
    int i = xy.y * sz0 + xy.x;
    if(v[z / 64] & (1ULL << (z % 64)))
      vret[i] = i;
  }
  return vret;
} // --------------------------------------------------------------------------

int RotateShiftHost::check_raw_field(const std::vector<uint64_t> &vrawfield,
                                     const std::vector<int> &vfield) {
  int ret = -1;
  int wsz = int(sqrt(static_cast<double>(vrawfield.size())));
  int szfld = int(sqrt(static_cast<double>(vfield.size())));
  assert(wsz * 8 == szfld);
  int cnterr = 0;
  for(int wy = 0; wy < wsz; wy++) {
    for(int wx = 0; wx < wsz; wx++) {
      uint64_t w = vrawfield[wy * wsz + wx];
      for(int by = 0; by < 8; by++) {
        for(int bx = 0; bx < 8; bx++) {
          int bit = int(w >> (by * 8 + bx)) & 1;
          int y = wy * 8 + by;
          int x = wx * 8 + bx;
          int idx = y * szfld + x;
          if(bit && vfield[idx] < 0 || !bit && vfield[idx] >= 0) {
            if(ret == - 1)
              ret = idx;
            if(++cnterr < 4)
              printf("check_raw_field error x=%d y=%d idx=%d raw=%d tst=%d\n",
                     x, y, idx, bit, vfield[idx]);
          }
        }
      }
    }
  }
  printf("check_raw_field errors=%d\n", cnterr);
  return ret;
} // --------------------------------------------------------------------------

int RotateShiftHost::check_raw(float angle, int2 shift,
                               const std::vector<uint64_t> &vrawsubdata,
                               const std::vector<uint64_t> &vrawfield) {
  std::vector<int> vsubdata = convert_raw_morton(vrawsubdata);
  std::vector<int> vfield = push(angle, shift, vsubdata);
  return check_raw_field(vrawfield, vfield);
} // --------------------------------------------------------------------------

std::vector<int> RotateShiftHost::def_subdata(int sz0) {
  std::vector<int> v(sz0 * sz0);
  for(int y = 0; y < sz0; y++)
    for(int x = 0; x < sz0; x++)
      v[y * sz0 + x] = y * sz0 + x;
  return v;
} // --------------------------------------------------------------------------

void RotateShiftHost::dump_subdata(const std::vector<int> &vsubdata) {
  int sz0 = int(sqrt(static_cast<double>(vsubdata.size())));
  if(sz0 <= 16) {
    for(int yr = 0; yr < sz0; yr++) {
      int y = sz0 - 1 - yr;
      for(int x = 0; x < sz0; x++) {
        int v = vsubdata[y * sz0 + x];
        if(v == -1)
          printf(" ..");
        else
          printf(" %X%X", v & (sz0 - 1), (v / sz0) & (sz0 - 1));
      }
      printf("\n");
    }
  } else {
    for(int yr = 0; yr < sz0; yr++) {
      int y = 2 * sz0 - 1 - yr;
      for(int x = 0; x < sz0; x++) {
        int v = vsubdata[y * sz0 + x];
        printf("%c", v == -1 ? '.' : 'x');
      }
      printf("\n");
    }
  }
  printf("\n");
} // --------------------------------------------------------------------------

std::vector<int> RotateShiftHost::push(float angle, int2 shift,
                                       const std::vector<int> &vsubdata) {
  int sz0 = int(sqrt(static_cast<double>(vsubdata.size())));
  int szfield = sz0 * 3 / 2;
  std::vector<int> vfield(szfield * szfield, -1);
  float rad = -angle * (3.14159265358979323846f / 180.0f);
  float2 d_rotate = {tanf(rad) * (1.f - sqrtf(2)), sinf(rad * 2.f) / sqrtf(2)};

  int2 fld;  // не центрированное поле со сдвигом [0..szfield-1]
  int2 fldc; // центрированное поле [-szfield/2..szfield/2-1]
  int2 subc; // центрированное субполе [-sz0/2..sz0/2-1] (включая точки которые
             // ему не принадлежат)
  int2
      sub; // не центрированное субполе [0..sz0-1] (включая точки которые ему не
  // принадлежат)
  int2 flds; // shifted field
  // printf("\nfld.y;fld.x;fldc.y;fldc.x;sub.y;sub.x;subc.y;subc.x;idfield;idsub;field;"
  //        "subdata\n");
  for(fld.y = 0; fld.y < szfield; fld.y++) {
    flds.y = (szfield + fld.y + shift.y) % szfield;
    fldc.y = flds.y - szfield / 2;
    for(fld.x = 0; fld.x < szfield; fld.x++) {
      flds.x = (szfield + fld.x + shift.x) % szfield;
      fldc.x = flds.x - szfield / 2;
      subc = rotate(fldc, d_rotate);
      sub = {subc.x + sz0 / 2, subc.y + sz0 / 2};
      int idfield = fld.y * szfield + fld.x;
      int idsub = sub.y * sz0 + sub.x;
      if(sub.x < 0 || sub.x >= sz0 || sub.y < 0 || sub.y >= sz0) {
        // printf("%5d;%5d;%6d;%6d;%5d;%5d;%6d;%6d;%7d;%5d;%5d;\n", fld.y,
        // fld.x,
        //        fldc.y, fldc.x, sub.y, sub.x, subc.y, subc.x, idfield, idsub,
        //        vfield[idfield]);
      } else {
        vfield[idfield] = vsubdata[idsub];
        // printf("%5d;%5d;%6d;%6d;%5d;%5d;%6d;%6d;%7d;%5d;%5d;%02X\n", fld.y,
        // fld.x,
        //        fldc.y, fldc.x, sub.y, sub.x, subc.y, subc.x, idfield, idsub,
        //        vfield[idfield], (unsigned)vsubdata[idsub]);
      }
    }
  }
  return vfield;
} // --------------------------------------------------------------------------

std::vector<int> RotateShiftHost::pull(float angle, int2 shift,
                                       const std::vector<int> &vfield) {
  int szfield = int(sqrt(static_cast<double>(vfield.size())));
  int sz0 = szfield * 2 / 3;
  assert(sz0 == sqrt(sz0 * sz0) && sz0 >= 16);
  std::vector<int> vsubdata(sz0 * sz0, -1);
  float rad = angle * (3.14159265358979323846f / 180.0f);
  float2 d_rotate = {tanf(rad) * (1.f - sqrtf(2.0f)),
                     sinf(rad * 2.0f) / sqrtf(2.0f)};
  for(int ysub = 0; ysub < sz0; ysub++) {
    int2 sub = {0, ysub - sz0 / 2};
    for(int xsub = 0; xsub < sz0; xsub++) {
      sub.x = xsub - sz0 / 2;
      int2 fldc2 = rotate(sub, d_rotate);
      int2 fld = {fldc2.x + szfield / 2, fldc2.y + szfield / 2};
      int y = (fld.y + szfield - shift.y) % szfield;
      int x = (fld.x + szfield - shift.x) % szfield;
      assert(x >= 0 && x < szfield && y >= 0 && y < szfield);
      int idsub = ysub * sz0 + xsub;
      int idfield = y * szfield + x;
      vsubdata[idsub] = vfield[idfield];
    }
  }
  return vsubdata;
} // --------------------------------------------------------------------------
