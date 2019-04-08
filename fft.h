#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>

// Реализуте пропущенные методы
// в качестве Т будет использоваться std::complex<float> // <double> // <long double>
namespace FFT {
// Возвращает корень степени degree из 1
template <typename T>
T GetRoot(size_t degree);

// Выполняет преобразование фурье квадратичным алгоритмом для вектора произвольной длины
template <typename T>
std::vector<T> FourierTransform(const std::vector<T>& data);

// Выполняет обратное преобразование квадратичным алгоритмом для вектора фиксированной длины
template <typename T>
std::vector<T> InverseFourierTransform(const std::vector<T>& data);

// Добивает вектор в конце нулями до длины expected_length,
// выбрасывает std::runtime_error если expected_length < data.size()
template <typename T>
std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length);

// Быстрое преобразование Фурье для вектора длины 2^k
template <typename T>
std::vector<T> FastFourierTransform(const std::vector<T>& data);

// Обратное быстрое преобразование Фурье для вектора длины 2^k
template <typename T>
std::vector<T> FastInverseFourierTransform(const std::vector<T>& data);

//Тесты для namespace FFT
void TestGetRoot();
void TestFourierTransform();
void TestInverseFourierTransform();
void TestAddPadding();
void TestFastFourierTransform();
void TestFastInverseFourierTransform();
} // namespace FFT