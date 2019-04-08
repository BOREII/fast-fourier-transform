#include "polynomial.h"
#include "fft.h"

template <typename T>
bool Polynomial<T>::operator==(const Polynomial& other) const {
    //ффт вычисляет приближенные значения, поэтому
    //будем округлять коэффициенты при сравнении

    size_t min_ = std::min(degree_, other.degree_);
    size_t max_ = std::max(degree_, other.degree_);

    //значения в коэффициентах при маленьких степенях должны совпадать
    for (size_t i = 0; i < min_; ++i) {
        if (std::round(std::real(coefficients_[i])) !=
            std::round(std::real(other.coefficients_[i]))
            ) {
            return false;
        }
    }

    //а при больших степенях, у одного из многочленов они нулевые,
    //поэтому у другого они не должны сильно отличаться от нуля
    for (size_t i = min_; i < max_; ++i) {
        if ((degree_ < other.degree_ && std::real(other.coefficients_[i]) != 0) ||
            (degree_ > other.degree_ && std::real(coefficients_[i]) != 0)
            ) {
            return false;
        }
    }

    return true;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& other) {
    size_t min_ = std::min(degree_, other.degree_);
    size_t max_ = std::max(degree_, other.degree_);

    //если текущая степень меньше новой, то
    //необходимо приплюсовать все элементы, в том числе
    //с более высокой степенью
    //иначе, достаточно просто добавить новые элементы, до
    //степени меньше либо равной текущей

    coefficients_.resize(max_, T(0));
    size_t bound = degree_ < other.degree_ ? max_ : min_;

    for (size_t i = 0; i < bound; ++i) {
        coefficients_[i] += other.coefficients_[i];
    }

    degree_ = max_;
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& other) {
    return *this += -other;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial& other) {
    size_t future_degree = (degree_ - 1) + (other.degree_ - 1) + 1;

    size_t new_deg = std::max(degree_, other.degree_);
    while ((new_deg & (new_deg - 1)) != 0) {
        ++new_deg;
    }

    new_deg *= 2;

    std::vector<T> this_values =
        FFT::FastFourierTransform(FFT::AddPadding<T>(coefficients_, new_deg));
    std::vector<T> other_values =
        FFT::FastFourierTransform(FFT::AddPadding<T>(other.coefficients_, new_deg));

    std::vector<T> multiply_values(new_deg);
    for (size_t i = 0; i < new_deg; ++i) {
        multiply_values[i] = this_values[i] * other_values[i];
    }

    *this = Polynomial(FFT::FastInverseFourierTransform(multiply_values));

    degree_ = future_degree;
    coefficients_.resize(degree_);

    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator^=(size_t pow) {
    if (pow == 0) {
        *this = Polynomial({1});
        return *this;
    }
    if (pow == 1) {
        return *this;
    }

    if (pow % 2 == 0) {
        size_t future_degree = (degree_ - 1) * pow + 1;

        *this ^= pow / 2;
        *this *= *this;

        degree_ = future_degree;
        coefficients_.resize(degree_);

        return *this;
    }
    else {
        size_t future_degree = (degree_ - 1) * pow + 1;

        Polynomial curr = *this;
        *this ^= pow - 1;
        *this *= curr;

        degree_ = future_degree;
        coefficients_.resize(degree_);

        return *this;
    }
}


template <typename T>
const Polynomial<T> Polynomial<T>::operator+(const Polynomial& other) const {
    Polynomial<T> result = *this;
    return result += other;
}

template <typename T>
const Polynomial<T> Polynomial<T>::operator-(const Polynomial& other) const {
    Polynomial<T> result = *this;
    return result -= other;
}

template <typename T>
const Polynomial<T> Polynomial<T>::operator*(const Polynomial<T>& other) const {
    Polynomial<T> result = *this;
    return result *= other;
}

template <typename T>
const Polynomial<T> Polynomial<T>::operator^(size_t pow) const {
    Polynomial<T> result = *this;
    return result ^= pow;
}

template <typename T>
const std::vector<T> Polynomial<T>::GetCoefficients() const {
    return coefficients_;
}


template <typename T>
const Polynomial<T> operator-(const Polynomial<T>& other) {
    std::vector<T> result_coefficients(other.degree_);

    std::transform(begin(other.coefficients_), end(other.coefficients_),
                   begin(result_coefficients), [](const T& element) {
          return -element;
        });

    return Polynomial(result_coefficients);
}


template <typename T>
std::ostream& operator<<(std::ostream& ostr, const Polynomial<T>& polynomial) {
    //вообще можно оставить и так, но вроде мы будем умножать многочлены с
    //действительными целыми коэффициентами,
    //поэтому для удобства будем округлять до целых

    const std::vector<T>& coeffs = polynomial.coefficients_;

    if (polynomial.degree_ == 0) {
        return ostr;
    }
    if (polynomial.degree_ == 1) {
        ostr << std::round(std::real(coeffs[0]));
        return ostr;
    }


    size_t first_non_zero_coeff_num =
        std::find_if(begin(coeffs), end(coeffs), [](const T& elem){
          return std::round(std::real(elem)) != 0;
        }) - begin(coeffs);

    int64_t rounded_coef = std::round(std::real(coeffs[first_non_zero_coeff_num]));

    if (first_non_zero_coeff_num != 0) {
        if (rounded_coef != 1) {
            ostr << rounded_coef;
        }
        ostr << "x";
        if (first_non_zero_coeff_num != 1) {
            ostr << "^";
        }
        ostr << first_non_zero_coeff_num;
    } else {
        ostr << rounded_coef;
    }


    if (first_non_zero_coeff_num != polynomial.degree_ - 1) {
        ostr << " + ";
    }

    size_t second_non_zero_coeff_num = std::find_if(
        begin(coeffs) + first_non_zero_coeff_num + 1, end(coeffs),
        [](const T& elem) {
            return std::round(std::real(elem)) != 0;
        }) - begin(coeffs);

    for (size_t coeff_num = second_non_zero_coeff_num;
         coeff_num < polynomial.degree_;
         ++coeff_num) {

        rounded_coef =
            std::round(std::real(coeffs[coeff_num]));

        if (rounded_coef != 0) {
            if (rounded_coef != 1 || rounded_coef != -1) {
                ostr << rounded_coef;
            }
            if (rounded_coef == -1) {
                ostr << "-";
            }
            ostr << "x";
            if (coeff_num != 1) {
                ostr << "^" << coeff_num;
            }
            if (coeff_num != polynomial.degree_ - 1) {
                ostr << " + ";
            }
        }
    }
    return ostr;
}

template class Polynomial<std::complex<float>>;
template class Polynomial<std::complex<double>>;
template class Polynomial<std::complex<long double>>;

template const Polynomial<std::complex<float>>
operator-(const Polynomial<std::complex<float>>& other);
template const Polynomial<std::complex<double>>
operator-(const Polynomial<std::complex<double>>& other);
template const Polynomial<std::complex<long double>>
operator-(const Polynomial<std::complex<long double>>& other);

template std::ostream&
operator<<(std::ostream& ostr, const Polynomial<std::complex<float>>& polynomial);
template std::ostream&
operator<<(std::ostream& ostr, const Polynomial<std::complex<double>>& polynomial);
template std::ostream&
operator<<(std::ostream& ostr, const Polynomial<std::complex<long double>>& polynomial);