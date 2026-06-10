#pragma once
#include <cstdint>
#include <vector>
/*
fill push:
(subdata) -> rotate -> (a) -> shift -> (b) -> wrap -> (field)

check pull:
(subdata) <- unrotate <- (a) <- unwrap <- (b) <- unshift <- (field)
*/
struct IdXY {
  int id = -1;
  int x = -INT32_MAX;
  int y = -INT32_MAX;
}; // --------------------------------------------------------------------------

struct RowDbDbg {
  int val = -1; // 0 or 1
  int id_subdata_morton = -1;
  IdXY subdata; // pull: result of unrotate
  IdXY a;       // push: result of rotate, pull: result of unwrap
  IdXY b;       // push: result of shift, pull: result of unshift
  IdXY field;   // push: result wrap
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
  RowDbDbg* find_by_field(int id);
  RowDbDbg* find_by_field(int x, int y);
  RowDbDbg* find_by_b(int id);
  RowDbDbg* find_by_b(int x, int y);
  RowDbDbg* find_by_a(int id);
  RowDbDbg* find_by_a(int x, int y);
  RowDbDbg* find_by_subdata(int id_morton);
  RowDbDbg* find_by_subdata(int x, int y);
}; // --------------------------------------------------------------------------