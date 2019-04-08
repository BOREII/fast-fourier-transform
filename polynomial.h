#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>


// Операции над многочленами с помощью ффт
template <typename T>
class Polynomial {
 public:
  explicit Polynomial(const std::vector<T>& coefficients)
      : coefficients_(coefficients), degree_(coefficients.size()) {
  }

  // Чтобы можно было написать Polynomial<std::complex<double>> p({1, 2, 1})
  // И получить представление многочлена 1 + 2x + x^2
  Polynomial(const std::initializer_list<T>& coefficients)
      : coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size()) {
  }

  Polynomial(const Polynomial&) = default;
  Polynomial(Polynomial&&) noexcept = default;

  Polynomial& operator=(const Polynomial&) = default;
  Polynomial& operator=(Polynomial&&) noexcept = default;

  size_t GetDegree() const {
      return degree_;
  }

  // Реализуйте следующие методы
  bool operator==(const Polynomial& other) const;

  Polynomial& operator+=(const Polynomial& other);

  Polynomial& operator-=(const Polynomial& other);

  // С Использованием быстрого преобразования фурье
  Polynomial& operator*=(const Polynomial& other);

  // Возведение в степень pow с помощью комбинации FFT и индийского возведения в степень
  Polynomial& operator^=(size_t pow);

  template <typename U>
  friend std::ostream& operator<<(std::ostream& ostr, const Polynomial<U>& polynomial);

  // Используя предыдущее
  const Polynomial operator+(const Polynomial& other) const;

  const Polynomial operator-(const Polynomial& other) const;

  const Polynomial operator*(const Polynomial& other) const;

  const Polynomial operator^(size_t pow) const;

  // И еще один, унарный минус
  template <typename U>
  friend const Polynomial<U> operator- (const Polynomial<U>& other);

  //для решения задачи про подстроки нужны будут коэффициенты
  const std::vector<T> GetCoefficients() const;

 private:
  std::vector<T> coefficients_;
  size_t degree_;
};


//Тесты для Polynomial
namespace PolynomialTests {
void CompareOperator();
void AddAndSubstractOperators();
void OutputStream();
void Multiply();
void Power();
} // namespace PolynomialTests
