#pragma once
#include "emu_vector_types.h"
#include <cstdint>
#include <vector>
//#include "dbdbg.h"

// @brief Пространство имен RotateShiftHost для поворота квадратов
// без потери и дублирования данных
// @note Это Host версия - нужна для экспериментов и отладки.
// @note Реализация для CUDA-версии - rotate_shift.h
namespace RotateShiftHost {
//DbDbg dbdbg;

// @brief Основная функция:
// shear поворот по центру против часовой стрелки
// @param pos позиция в входном поле [-szfield/2, szfield/2-1]
// @param d_rotate = {tanf(rad) * (1 - SQRT_2), sinf(rad * 2) / SQRT_2};
int2 rotate(int2 pos, float2 d_rotate);

// @brief Помещает vsubdata в выходное поле vfield с наклоном и сдвигом
// относительно центра
// @param angle угол поворота в градусах [-45, 45]
// @param shift сдвиг в выходном поле [-szfield/2, szfield/2)
// @param vsubdata входное поле размером sz0*sz0, где sz0=2^N
// @return vfield: выходное поле размером szfield*szfield, где szfield=1.5*sz0
std::vector<int> push(float angle, int2 shift,
                      const std::vector<int> &vsubdata);

// @brief Извлекает из vfield квадрат с наклоном и сдвигом относительно центра
// @param angle угол поворота в градусах [-45, 45]
// @param shift сдвиг в выходном поле [-szfield/2, szfield/2)
// @param vfield входное поле размером szfield*szfield, где szfield=1.5*sz0
// @return vsubdata: выходное поле размером sz0*sz0, где sz0=2^N
std::vector<int> pull(float angle, int2 shift, const std::vector<int> &vfield);

// Функции для тестирования

// @brief Создаёт тестовый квадрат для отладки push/pull
// @param sz0 сторона квадрата 2^N
// @note Каждая ячейка содержит текущий индекс [0, sz0*sz0)
std::vector<int> def_subdata(int sz0);

// @brief Выводит субполя в консоль для отладки push/pull
void dump_subdata(const std::vector<int> &vsubdata);

// @brief Вывод поля в консоль для отладки push/pull
// @param shift сдвиг в выходном поле
// @note shift может использоваться для наглядности,
// чтобы сместить сильно сдвинутую фигуру в центр
void dump_field(const std::vector<int> &vfield, int2 shift = {0, 0});

// @brief Проверка поля на корректность
int check(const std::vector<int> &vfield);

// @brief Преобразование Z 1-битовых данных в вектор для push()
std::vector<int> convert_raw_morton(const std::vector<uint64_t> &vrawsubdata);

// @brief Сравнение внешних битовых данных и возвращаемого вектора из push()
// @return true, если ошибок нет
// @param vrawfield внешнее поле размером 8*8 значения / элемент
// @param vfield поле размером 1.5*sz0 возвращаемое из push() 1 значение /
// элемент
// @note Выводит 5 первых ошибок в консоль
int check_raw_field(const std::vector<uint64_t> &vrawfield,
                     const std::vector<int> &vfield);

int check_raw(float angle, int2 shift,
    const std::vector<uint64_t> &vrawsubdata,
               const std::vector<uint64_t> &vrawfield);
} // namespace RotateShiftHost