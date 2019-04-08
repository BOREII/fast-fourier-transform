#include "substring_matching.h"
#include "fft.h"
#include "polynomial.h"
#include "test_runner.h"

namespace SubstringMatching {
std::vector<size_t> FindSubstrings(const std::string& str,
                                   const std::string& pattern) {
    std::vector<size_t> result;
    result.reserve(str.size() - pattern.size() + 1);

    for (size_t i = 0; i < result.capacity(); ++i) {
        bool substring_found = true;
        for (size_t j = 0; j < pattern.size(); ++j) {
            if (str[i + j] != pattern[j]) {
                substring_found = false;
            }
        }
        if (substring_found) {
            result.push_back(i);
        }
    }

    result.shrink_to_fit();
    return result;
}

std::vector<size_t> FindSubstringsFFT(const std::string& str,
                                      const std::string& pattern) {

    if (str.empty() || pattern.empty() || pattern.size() > str.size()) {
        return {};
    }

    //закодируем строчки числами в ascii кодировке
    //исходная строка - прямая, а подстрока - развернутая
    //кажется, что такой конструктор сделает
    //static_cast<std::complex<float>> к каждому элементу строки
    std::vector<std::complex<long double>> str_(begin(str), end(str));
    std::vector<std::complex<long double>> pattern_(rbegin(pattern), rend(pattern));

    //вычислим сумму квадратов всех элементов для подстроки
    long double square_sum_pattern = std::accumulate(
        begin(pattern_), end(pattern_), 0,
        [](const long double& curr_sum, const auto& elem) {
            return curr_sum + pow(std::real(elem), 2);
        }
    );

    //вычислим произведение многочленов c коэфициентами,
    //равными элементам векторов str_ и pattern_ c помощью ффт
    std::vector<std::complex<long double>> multiply_coefficients =
        (Polynomial<std::complex<long double>>(str_) *
         Polynomial<std::complex<long double>>(pattern_)
         ).GetCoefficients();

    //вычислим сумму квадратов первых pattern_.size() элементов для строки
    long double first_m_square_sum_str = std::accumulate(
        begin(str_), begin(str_) + pattern_.size(), 0,
        [](const long double& curr_sum, const std::complex<long double>& elem) {
          return curr_sum + pow(std::real(elem), 2);
        }
    );

    //теперь вычислим суммы квадратов разностей элементов в строке и в подстроке
    //для всех str_.size() - pattern_.size() + 1 подстрок
    std::vector<long double> sum_square_diff(str_.size() - pattern_.size() + 1);

    sum_square_diff[0] =
        square_sum_pattern -
        2 * std::real(multiply_coefficients[pattern_.size() - 1]) +
        first_m_square_sum_str;

    for (size_t i = 1; i < sum_square_diff.size(); ++i) {
        sum_square_diff[i] = sum_square_diff[i - 1] +
            2 * (std::real(multiply_coefficients[pattern_.size() - 2 + i])) -
            2 * (std::real(multiply_coefficients[pattern_.size() - 1 + i])) -
            std::real(str_[i - 1]) * std::real(str_[i - 1]) +
            std::real(str_[pattern_.size() - 1 + i]) *
            std::real(str_[pattern_.size() - 1 + i]);
    }

    //если квадрат разности равен нулю, то
    //подстрока в строке и исходная подстрока совпали
    //соберем ответ
    std::vector<size_t> result;
    result.reserve(sum_square_diff.size());

    for (size_t i = 0; i < sum_square_diff.size(); ++i) {
        if (std::round(sum_square_diff[i]) == 0) {
            result.push_back(i);
        }
    }

    result.shrink_to_fit();

    return result;
}

std::vector<size_t> FindMatches(const std::string& str,
                                const std::string& pattern) {

    std::vector<size_t> result;
    result.reserve(str.size() - pattern.size() + 1);

    for (size_t i = 0; i < result.capacity(); ++i) {
        bool substring_found = true;
        for (size_t j = 0; j < pattern.size(); ++j) {
            if (str[i + j] != pattern[j] && pattern[j] != '?') {
                substring_found = false;
            }
        }
        if (substring_found) {
            result.push_back(i);
        }
    }

    result.shrink_to_fit();
    return result;
}

std::vector<size_t> FindMatchesFFT(const std::string& str,
                                const std::string& pattern) {
    //алгоритм практически не отличается от предыдущего,
    //тепень нам нужно, чтобы занулилась не просто сумма квадратов разностей, а
    //сумма квадратов разнсотей, домноженных на элементы элемента pattern,
    //'?' будем кодировать нулем, и они будут все занулять, когда строка подходящая


    if (str.empty() || pattern.empty() || pattern.size() > str.size()) {
        return {};
    }

    std::vector<std::complex<long double>> str_(begin(str), end(str));
    std::vector<std::complex<long double>> pattern_(rbegin(pattern), rend(pattern));
    std::vector<std::complex<long double>> str_squared(str_.size());
    std::vector<std::complex<long double>> pattern_squared(pattern_.size());

    std::for_each(begin(pattern_), end(pattern_), [](auto& elem){
         elem = static_cast<char>(std::real(elem)) == '?' ? 0 : elem;
    });

    auto make_squared = [](const auto& elem) {
      return std::pow(elem, 2);
    };

    std::transform(begin(str_), end(str_), begin(str_squared), make_squared);
    std::transform(begin(pattern_), end(pattern_), begin(pattern_squared), make_squared);

    //вычислим сумму кубов всех элементов для подстроки
    long double cube_sum_pattern = std::accumulate(
        begin(pattern_), end(pattern_), 0,
        [](const long double& curr_sum, const auto& elem) {
          return curr_sum + pow(std::real(elem), 3);
        }
    );

    //вычислим c помощью ффт произведение многочленов c коэффициентами,
    //равными квадратам элементам векторов str_ и элементам pattern_
    //а так же с коэффициентами равными элементам str_ и
    //квадратам элементов pattern_
    std::vector<std::complex<long double>> multiply_coefficients_str_squared =
        (Polynomial<std::complex<long double>>(str_squared) *
            Polynomial<std::complex<long double>>(pattern_)
        ).GetCoefficients();
    std::vector<std::complex<long double>> multiply_coefficients_pattern_squared =
        (Polynomial<std::complex<long double>>(str_) *
            Polynomial<std::complex<long double>>(pattern_squared)
        ).GetCoefficients();

    //подсчитаем итоговую сумму
    std::vector<long double> sum_square_diff(str_.size() - pattern_.size() + 1);

    sum_square_diff[0] =
        cube_sum_pattern -
        2 * std::real(multiply_coefficients_pattern_squared[pattern_.size() - 1]) +
        std::real(multiply_coefficients_str_squared[pattern_.size() - 1]);

    for (size_t i = 1; i < sum_square_diff.size(); ++i) {
        sum_square_diff[i] = sum_square_diff[i - 1] +
            2 * (std::real(multiply_coefficients_pattern_squared[pattern_.size() - 2 + i])) -
            2 * (std::real(multiply_coefficients_pattern_squared[pattern_.size() - 1 + i]))  -
            std::real(multiply_coefficients_str_squared[pattern_.size() - 2 + i]) +
            std::real(multiply_coefficients_str_squared[pattern_.size() - 1 + i]);
    }

    //если зануляется элемент суммы,
    //то будет подстроки совпали с учетом джокеров
    std::vector<size_t> result;
    result.reserve(sum_square_diff.size());

    for (size_t i = 0; i < sum_square_diff.size(); ++i) {
        if (std::round(sum_square_diff[i]) == 0) {
            result.push_back(i);
        }
    }

    result.shrink_to_fit();

    return result;

}
}