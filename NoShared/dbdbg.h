#pragma once
#include <cstdint>
#include <vector>
#include "emu_vector_types.h"

struct IdXY {
  int id = -1;
  int2 xy = {-INT32_MAX, -INT32_MAX};
}; // --------------------------------------------------------------------------

// fill push:
//(subdata) -> center -> (a) -> -rotate -> (b) -> -shift -> (c) -> wrapc ->
//(d) -> center -> (field)
//(d) -> +shift -> (e)
//
// check pull:
//(subdata) <- center <- (a) <- +rotate <- (b) <- wrapc <- (e) <- +shift <-
//(d) <- center <- (field)
struct RowDbDbg {
  int val = -1; // 0 or 1
  int id_subdata_morton = -1;

  // push: source
  // pull: result of unrotate,
  // non centered, [0, sz0)
  IdXY subdata;

  // push: result of center(subdata)
  // pull: result of rotate(b + angle)
  // centered, [-sz0/2, sz0/2)
  IdXY a;

  // push: result of rotate(a - angle)
  // pull: result of of wrapc(e, szfield)
  // centered [-sz0/sqrt(2), sz0/sqrt(2))
  IdXY b;  

  // push: result of b - shift
  // pull: not used
  // centered
  IdXY c;

  // push: result of wrapc(c, szfield)
  // pull: result of center(field)
  // centered
  IdXY d;

  // push(write): result of d + shift
  // pull(read): result of d + shift
  // centered
  IdXY e;

  // push: result decenter(d)
  // pull: source
  // not centered
  IdXY field;
}; // --------------------------------------------------------------------------

class DbDbg {
public:
  int sz0 = -1;      // size of subdata
  int szfield = -1;  // size of field
  int wsz0 = -1;     // size of subdata in words
  int wszfield = -1; // size of field in words
  std::vector<RowDbDbg> vrows;

  /// @brief fill vrows: subdata, id_subdata_morton, val
  /// @param v_subdata 1 бит/ячейка. Код Мортона
  void create(const std::vector<uint64_t> &v_subdata);

  /// @brief fill vrows: subdata, id_subdata_morton, val
  /// @param v_subdata 1 значение/ячейка, код Мортона
  void create(const std::vector<int> &v_subdata);

  RowDbDbg get_row(int id);
  RowDbDbg operator[](int id) { return get_row(id); }

  // find row
  RowDbDbg *find_by_field(int id);
  RowDbDbg *find_by_field(int x, int y);
  RowDbDbg *find_by_b(int id);
  RowDbDbg *find_by_b(int x, int y);
  RowDbDbg *find_by_a(int id);
  RowDbDbg *find_by_a(int x, int y);
  RowDbDbg *find_by_subdata(int id_morton);
  RowDbDbg *find_by_subdata(int x, int y);
}; // --------------------------------------------------------------------------