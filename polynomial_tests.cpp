#include "fft.h"
#include "test_runner.h"
#include "polynomial.h"

namespace PolynomialTests {
void CompareOperator() {
    using std::complex;
    using std::vector;

    ASSERT(
        Polynomial<complex<float>>({1, 2, 3, 4, 0, 0, 0}) ==
            Polynomial<complex<float>>({1, 2, 3, 4, 0, 0})
    );

    ASSERT(!(Polynomial<complex<double>>({-1}) == Polynomial<complex<double>>({1})));

    vector<complex<long double>> v5 = {5, 3, 2, 4};
    ASSERT(
        Polynomial(FFT::FourierTransform<complex<long double>>(v5)) ==
            Polynomial(FFT::FastFourierTransform<complex<long double>>(v5))
    );
}

void AddAndSubstractOperators() {
    using std::complex;
    using std::vector;
    Polynomial<complex<long double>> p1({1, 2, 3, 4, 2, 3, 16}),
        p2({1, 2, 3, 4, -1, 7}),
        p3({2, 4, 6, 8, 1, 10, 16});
    ASSERT_EQUAL(p1 + p2, p3);
    ASSERT_EQUAL(p1, p3 - p2);
    ASSERT_EQUAL(p2, p3 - p1);

    p1 += p2;
    ASSERT_EQUAL(p1, p3);
    p1 -= p2;
    p3 -= p1;
    ASSERT_EQUAL(p2, p3);

    Polynomial<complex<long double>> p4({1, 2, 3}), p5({-1, -2, -3});

    ASSERT_EQUAL(p4, -p5);
    ASSERT_EQUAL(p4, -(-p4));
}

void OutputStream() {
    Polynomial<std::complex<long double>> p({1, 1});
    p ^= 4;
    std::ostringstream os;
    os << p;
    ASSERT_EQUAL(
        os.str(),
        "1 + 4x + 6x^2 + 4x^3 + x^4"
    );

}

void Multiply() {
    using std::complex;

    Polynomial<complex<long double>> p1({1, 2, 3}), p2({3, 2, 1}), p3({3, 8, 14, 8, 3});
    ASSERT_EQUAL(p1 * p2, p3);

    Polynomial<complex<long double>>
        p4({7, 2, 13, 6}),
        p5({3, 1, 0,  0, 5, 8, 9}),
        p6({21, 13, 41, 31, 41, 66, 144, 152, 165, 54});

    ASSERT_EQUAL(p4 * p5, p6);
    ASSERT_EQUAL(p4 * p5 * p5, p6 * p5);
    ASSERT_EQUAL(p4 * p5 * p5 * p6, p6 * p5 * p6);
}

void Power() {
    using std::complex;
    Polynomial<complex<long double>>
        p1({1, 2}), p2({1, 4, 4}), p3({1, 6, 12, 8}), p4({1, 8, 24, 32, 16});
    ASSERT_EQUAL(p1^2, p2);
    ASSERT_EQUAL(p1^3, p3);
    ASSERT_EQUAL(p1^4, p4);
    ASSERT_EQUAL(p1 * p1 * p2 * p2, (p1 ^ 2) * (p2^2));

    Polynomial<complex<long double>> p5({1, 2, 3, 4}), p6({
                                                              1, 6, 21, 56, 111, 174, 219, 204, 144, 64
                                                          });
    ASSERT_EQUAL(p5^3, p6);
    ASSERT_EQUAL((p5^2) * (p6^3), p6 * p6 * p6 * p5 * p5);

    p5 ^= 6;
    p6 ^= 2;
    ASSERT_EQUAL(p5, p6);
}
}