#include "fft.h"
#include "test_runner.h"

namespace FFT {
void TestGetRoot() {
    using std::complex;

    const long double error = 1.0e-06;

    ASSERT_ERROR(GetRoot<complex<float>>(1), complex<float>(1, 0), error);
    ASSERT_ERROR(GetRoot<complex<double>>(2), complex<double>(-1, 0), error);
    ASSERT_ERROR(
        GetRoot<complex<long double>>(3), complex<long double>(-0.5, sqrt(3) / 2), error
    );
    ASSERT_ERROR(
        GetRoot<complex<long double>>(4), complex<long double>(0, 1), error
    );
    ASSERT_ERROR(
        GetRoot<complex<long double>>(8),
        complex<long double>(sqrt(2) / 2, sqrt(2) / 2), error
    );
}

void TestFourierTransform() {
    using std::complex;
    using std::vector;

    const long double error = 1.0e-15;

    vector<complex<float>> v1 = {4};
    ASSERT_VECTOR(FourierTransform<complex<float>>({4}), v1, error);
    vector<complex<double>> v2 = {9, -1};
    ASSERT_VECTOR(FourierTransform<complex<double>>({4, 5}), v2, error);
    vector<complex<long double>> v3 = {3, 0, 0};
    ASSERT_VECTOR(FourierTransform<complex<long double >>({1, 1, 1}), v3, error);
    vector<complex<long double>> v4 = {
        {6, 0}, {-1.5, -sqrt(3) / 2}, {-1.5, sqrt(3) / 2}
    };
    ASSERT_VECTOR(FourierTransform<complex<long double>>({1, 2, 3}), v4, error);
}

void TestInverseFourierTransform() {
    using std::complex;
    using std::vector;

    const long double error = 1.0e-15;

    vector<complex<float>> v1 = {4};
    ASSERT_VECTOR(InverseFourierTransform<complex<float>>({4}), v1, error);
    vector<complex<double>> v2 = {4, 5};
    ASSERT_VECTOR(InverseFourierTransform<complex<double>>({9, -1}), v2, error);
    vector<complex<long double>> v3 = {1, 1, 1};
    ASSERT_VECTOR(InverseFourierTransform<complex<long double >>({3, 0, 0}), v3, error);
    vector<complex<long double>> v4 = {
        {6, 0}, {-1.5, -sqrt(3) / 2}, {-1.5, sqrt(3) / 2}
    };
    vector<complex<long double>> v5 = {1, 2, 3};
    ASSERT_VECTOR(InverseFourierTransform<complex<long double>>(v4), v5, error);

    //более сложная проверка, обратное преобразование от прямого
    const long double error_ = 1.0e-12;
    vector<complex<long double>> v6 = {71, 125, 161, 123, 356, 124, 5512, 123, 3, 0, 123};
    ASSERT_VECTOR(
        InverseFourierTransform<complex<long double>>(FourierTransform<complex<long double>>(v6)),
        v6, error_
    );
}

void TestAddPadding() {
    using std::complex;
    using std::vector;

    const long double error = 1.0e-15;

    vector<complex<float>> v1 = {4, 0};
    ASSERT_VECTOR(AddPadding<complex<float>>({4, 0}, 2), v1, error);
    vector<complex<double>> v2 = {4, 5, 0, 0, 0};
    ASSERT_VECTOR(AddPadding<complex<double>>({4, 5}, 5), v2, error);

    //проверка длины
    //try {
    //    AddPadding<complex<long double>>({1, 1, 1}, 2);
    //} catch (std::exception& e) {
    //    std::cerr << e.what();
    //}

}

void TestFastFourierTransform() {
    using std::complex;
    using std::vector;

    const long double error = 1.0e-5;

    //проверка на степень двойки
    //try {
    //    FastFourierTransform<complex<long double>>({1, 1, 1});
    //} catch (std::exception& e) {
    //    std::cerr << e.what();
    //}

    vector<complex<float>> v1 = {0, 1, 2, 3};
    ASSERT_VECTOR(FourierTransform<complex<float>>(v1),
                  FastFourierTransform<complex<float>>(v1), error);

    vector<complex<double>> v2 = {5, 3, 2, 4};
    ASSERT_VECTOR(FourierTransform<complex<double>>(v2),
                  FastFourierTransform<complex<double>>(v2), error);

    vector<complex<long double>> v3 = {1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_VECTOR(FourierTransform<complex<long double>>(v3),
                  FastFourierTransform<complex<long double>>(v3), error);
    vector<complex<long double>> v4 = {8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_VECTOR(FourierTransform<complex<long double>>(v4),
                  FastFourierTransform<complex<long double>>(v4), error);
}

void TestFastInverseFourierTransform() {
    using std::complex;
    using std::vector;

    const long double error = 1.0e-5;

    //проверка на степень двойки
    //try {
    //    FastInverseFourierTransform<complex<long double>>({1, 1, 1});
    //} catch (std::exception& e) {
    //    std::cerr << e.what();
    //}

    vector<complex<float>> v1 = {9, -1, 5, 3};
    ASSERT_VECTOR(InverseFourierTransform<complex<float>>(v1),
                  FastInverseFourierTransform<complex<float>>(v1), error);

    vector<complex<double>> v2 = {3, 4, 5, -6};
    ASSERT_VECTOR(InverseFourierTransform<complex<double>>(v2),
                  FastInverseFourierTransform<complex<double>>(v2), error);

    vector<complex<long double>> v3 = {15, 22, 31, 46, 51, 36, 72, 84};
    ASSERT_VECTOR(InverseFourierTransform<complex<long double>>(v3),
                  FastInverseFourierTransform<complex<long double>>(v3), error);

    vector<complex<long double>> v4 = {8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8,
                                       8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_VECTOR(InverseFourierTransform<complex<long double>>(v4),
                  FastInverseFourierTransform<complex<long double>>(v4), error);

    vector<complex<long double>> v5 = {8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8,
                                       8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8,
                                       8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8,
                                       8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_VECTOR(InverseFourierTransform<complex<long double>>(v5),
                  FastInverseFourierTransform<complex<long double>>(v5), error);

    //более сложная проверка, обратное преобразование от прямого
    ASSERT_VECTOR(
        FastInverseFourierTransform<complex<long double>>(
            FastFourierTransform<complex<long double>>(v5)
        ), v5, error
    );
}
}