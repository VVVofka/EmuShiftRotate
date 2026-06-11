#include "dbdbg.h"
#include <cassert>

static uint32_t part1by1_32(uint32_t x) {
  x &= 0x0000FFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;
  return x;
} // --------------------------------------------------------------------------
uint32_t morton_encode(int x, int y) {
  return part1by1_32((uint32_t)x) | (part1by1_32((uint32_t)y) << 1);
} // --------------------------------------------------------------------------

void DbDbg::create(const std::vector<uint64_t> &v_subdata) {
  wsz0 = int(sqrt(static_cast<double>(v_subdata.size())));
  wszfield = sz0 * 3 / 2;
  sz0 = sz0 * 8;
  szfield = szfield * 8;
  vrows.resize(sz0 * sz0);

  for(int y = 0; y < sz0; y++) {
    for(int x = 0; x < sz0; x++) {
      int i = y * sz0 + x;
      int z = morton_encode(x, y);
      RowDbDbg &r = vrows[i];
      r.id_subdata_morton = z;
      r.subdata.id = i;
      r.subdata = {x, y};
      r.val = uint32_t(v_subdata[z / 64] >> (z % 64)) & 1;
    }
  }
} // --------------------------------------------------------------------------
void DbDbg::create(const std::vector<int> &v_subdata) {
  sz0 = int(sqrt(static_cast<double>(v_subdata.size())));
  szfield = sz0 * 3 / 2;
  wsz0 = sz0 / 8;
  wszfield = szfield / 8;
  vrows.resize(sz0 * sz0);

  for(int y = 0; y < sz0; y++) {
    for(int x = 0; x < sz0; x++) {
      int2 sub = {x, y};
      int i = y * sz0 + x;
      int z = morton_encode(x, y);
      RowDbDbg &r = vrows[i];
      r.subdata.id = i;
      r.id_subdata_morton = z;
      r.subdata.xy = sub;
      r.val = v_subdata[i];
      r.a.xy = {x - sz0 / 2, y - sz0 / 2};
    }
  }
} // --------------------------------------------------------------------------
RowDbDbg DbDbg::get_row(int id) {
  assert(id >= 0 && id < sz0 * sz0);
  return vrows[id];
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_field(int id) {
  assert(id >= 0 && id < szfield * szfield);
  for(RowDbDbg &r : vrows)
    if(r.field.id == id)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_field(int x, int y) {
  assert(x >= 0 && x < szfield);
  assert(y >= 0 && y < szfield);
  for(RowDbDbg &r : vrows)
    if(r.field.xy.x == x && r.field.xy.y == y)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_b(int id) {
  assert(id >= 0 && id < szfield * szfield);
  for(RowDbDbg &r : vrows)
    if(r.b.id == id)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_b(int x, int y) {
  assert(x >= 0 && x < szfield);
  assert(y >= 0 && y < szfield);
  for(RowDbDbg &r : vrows)
    if(r.b.xy.x == x && r.b.xy.y == y)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_a(int id) {
  assert(id >= 0 && id < szfield * szfield);
  for(RowDbDbg &r : vrows)
    if(r.a.id == id)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_a(int x, int y) {
  assert(x >= 0 && x < szfield);
  assert(y >= 0 && y < szfield);
  for(RowDbDbg &r : vrows)
    if(r.a.xy.x == x && r.a.xy.y == y)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_subdata(int id_morton) {
  assert(id_morton >= 0 && id_morton < sz0 * sz0);
  for(RowDbDbg &r : vrows)
    if(r.id_subdata_morton == id_morton)
      return &r;
  return nullptr;
} // --------------------------------------------------------------------------

RowDbDbg *DbDbg::find_by_subdata(int x, int y) {
  assert(x >= 0 && x < sz0);
  assert(y >= 0 && y < sz0);
  int z = morton_encode(x, y);
  return find_by_subdata(z);
} // --------------------------------------------------------------------------
