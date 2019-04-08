#include "fft.h"

namespace FFT {
    template <typename T>
    T GetRoot(size_t degree) {
        static const typename T::value_type dPI = 2 * acos(-1);
        return {cos(dPI / degree), sin(dPI / degree)};
    }

    template <typename T>
    std::vector<T> FourierTransform(const std::vector<T>& data) {
        size_t degree = data.size();
        T root = GetRoot<T>(degree);

        std::vector<std::vector<T>> Matrix(degree, std::vector<T>(degree));
        for (size_t i = 0; i < degree; ++i) {
            //чтобы не вычислять степень корня каждый раз заново
            //будем умножать каждый следующий элемент в строке
            //на второй слева
            Matrix[i][0] = 1;
            if (degree > 1) {
                Matrix[i][1] = pow(root, i);
            }
            for (size_t j = 2; j < degree; ++j) {
                Matrix[i][j] = Matrix[i][1] * Matrix[i][j-1];
            }
        }

        std::vector<T> result(degree, 0);
        for (size_t i = 0; i < degree; ++i) {
            for (size_t j = 0; j < degree; ++j) {
                result[i] += Matrix[i][j] * data[j];
            }
        }

        return result;
    }

    template <typename T>
    std::vector<T> InverseFourierTransform(const std::vector<T>& data) {
        size_t degree = data.size();
        T inv_root = std::conj(GetRoot<T>(degree)) / std::norm(GetRoot<T>(degree));

        std::vector<std::vector<T>> Matrix(degree, std::vector<T>(degree));
        for (size_t i = 0; i < degree; ++i) {
            //чтобы не вычислять степень корня каждый раз заново
            //будем умножать каждый следующий элемент в строке
            //на второй слева
            Matrix[i][0] = 1;
            if (degree > 1) {
                Matrix[i][1] = pow(inv_root, i);
            }
            for (size_t j = 2; j < degree; ++j) {
                Matrix[i][j] = Matrix[i][1] * Matrix[i][j-1];
            }
        }
        for (size_t i = 0; i < degree; ++i) {
            for (size_t j = 0; j < degree; ++j) {
                Matrix[i][j] /= T(degree);
            }
        }


        std::vector<T> result(degree);
        for (size_t i = 0; i < degree; ++i) {
            for (size_t j = 0; j < degree; ++j) {
                result[i] += Matrix[i][j] * data[j];
            }
        }

        return result;
    }

    template <typename T>
    std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length) {
        if (expected_length < data.size()) {
            std::ostringstream os;
            os << "Exception thrown in AddPadding, expected length is less than current: "
                << expected_length << " < " << data.size() << "\n";
            throw std::runtime_error(os.str());
        }

        std::vector<T> result = data;
        result.resize(expected_length, T(0));
        return result;
    }

    template <typename T>
    std::vector<T> FastFourierTransform(const std::vector<T>& data) {
        size_t degree = data.size();
        if (degree == 1) {
            return data;
        }
        //проверка на степень двойки
        if (degree == 0 || ((degree & (degree - 1)) != 0)) {
            std::ostringstream os;
            os << "Exception thrown in FastFourierTransform,"
                  " expected degree is not the power of two: "
                << degree << "\n";
            throw std::runtime_error(os.str());
        }

        T root = GetRoot<T>(degree);
        T curr_root = 1;

        std::vector<T> even_coef(degree / 2), odd_coef(degree / 2);
        for (size_t i = 0; i < degree / 2; ++i) {
            even_coef[i] += data[2 * i];
            odd_coef[i] += data[2 * i + 1];
        }

        std::vector<T> even_data = FastFourierTransform(even_coef);
        std::vector<T> odd_data = FastFourierTransform(odd_coef);

        std::vector<T> result(degree);

        for (size_t i = 0; i < degree / 2; ++i) {
            result[i] += even_data[i] + curr_root * odd_data[i];
            result[i + degree / 2] +=  even_data[i] - curr_root * odd_data[i];
            curr_root *= root;
        }

        return result;
    }

    template <typename T>
    std::vector<T> FastInverseFourierTransform(const std::vector<T>& data) {
        size_t degree = data.size();
        if (degree == 1) {
            return data;
        }
        //проверка на степень двойки
        if (degree == 0 || ((degree & (degree - 1)) != 0)) {
            std::ostringstream os;
            os << "Exception thrown in FastInverseFourierTransform,"
                  " expected degree is not the power of two: "
                  << degree << "\n";
            throw std::runtime_error(os.str());
        }

        T inv_root = std::conj(GetRoot<T>(degree)) / std::norm(GetRoot<T>(degree));
        T curr_inv_root = 1;

        std::vector<T> even_coef(degree / 2), odd_coef(degree / 2);
        for (size_t i = 0; i < degree / 2; ++i) {
            even_coef[i] += data[2 * i];
            odd_coef[i] += data[2 * i + 1];
        }

        std::vector<T> even_data = FastInverseFourierTransform(even_coef);
        std::vector<T> odd_data = FastInverseFourierTransform(odd_coef);

        std::vector<T> result(degree);

        for (size_t i = 0; i < degree / 2; ++i) {
            result[i] += even_data[i] + curr_inv_root * odd_data[i];
            result[i + degree / 2] +=  even_data[i] - curr_inv_root * odd_data[i];
            //делим на два на каждом уровне рекурсии,
            //тогда в конце произойдет деление на 2^(log_2 n) = n
            result[i] /= T(2);
            result[i + degree / 2] /= T(2);
            curr_inv_root *= inv_root;
        }

        return result;
    }


    template std::complex<float> GetRoot<std::complex<float>>(size_t degree);
    template std::complex<double> GetRoot<std::complex<double>>(size_t degree);
    template std::complex<long double> GetRoot<std::complex<long double>>(size_t degree);

    template std::vector<std::complex<float>>
    FourierTransform(const std::vector<std::complex<float>>& data);
    template std::vector<std::complex<double>>
    FourierTransform(const std::vector<std::complex<double>>& data);
    template std::vector<std::complex<long double>>
    FourierTransform(const std::vector<std::complex<long double>>& data);

    template std::vector<std::complex<float>>
    InverseFourierTransform(const std::vector<std::complex<float>>& data);
    template std::vector<std::complex<double>>
    InverseFourierTransform(const std::vector<std::complex<double>>& data);
    template std::vector<std::complex<long double>>
    InverseFourierTransform(const std::vector<std::complex<long double>>& data);

    template std::vector<std::complex<float>>
    AddPadding(const std::vector<std::complex<float>>& data, size_t expected_length);
    template std::vector<std::complex<double>>
    AddPadding(const std::vector<std::complex<double>>& data, size_t expected_length);
    template std::vector<std::complex<long double>>
    AddPadding(const std::vector<std::complex<long double>>& data, size_t expected_length);

    template std::vector<std::complex<float>>
    FastFourierTransform(const std::vector<std::complex<float>>& data);
    template std::vector<std::complex<double>>
    FastFourierTransform(const std::vector<std::complex<double>>& data);
    template std::vector<std::complex<long double>>
    FastFourierTransform(const std::vector<std::complex<long double>>& data);

    template std::vector<std::complex<float>>
    FastInverseFourierTransform(const std::vector<std::complex<float>>& data);
    template std::vector<std::complex<double>>
    FastInverseFourierTransform(const std::vector<std::complex<double>>& data);
    template std::vector<std::complex<long double>>
    FastInverseFourierTransform(const std::vector<std::complex<long double>>& data);
}