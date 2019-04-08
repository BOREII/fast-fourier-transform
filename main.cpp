#include "test_runner.h"
#include "fft.h"
#include "polynomial.h"
#include "substring_matching.h"
#include "profile.h"

//Запуск всех тестов
void TestAll() {
    TestRunner tr;
    RUN_TEST(tr, FFT::TestGetRoot);
    RUN_TEST(tr, FFT::TestFourierTransform);
    RUN_TEST(tr, FFT::TestInverseFourierTransform);
    RUN_TEST(tr, FFT::TestAddPadding);
    RUN_TEST(tr, FFT::TestFastFourierTransform);
    RUN_TEST(tr, FFT::TestFastInverseFourierTransform);
    RUN_TEST(tr, PolynomialTests::CompareOperator);
    RUN_TEST(tr, PolynomialTests::AddAndSubstractOperators);
    RUN_TEST(tr, PolynomialTests::OutputStream);
    RUN_TEST(tr, PolynomialTests::Multiply);
    RUN_TEST(tr, PolynomialTests::Power);
    RUN_TEST(tr, SubstringMatching::TestSubstrings);
    RUN_TEST(tr, SubstringMatching::TestSubstringsHeyJude);
    RUN_TEST(tr, SubstringMatching::TestMatches);
    RUN_TEST(tr, SubstringMatching::TestMatchesHeyJude);
}

int main() {
    //TestAll();

    Polynomial<std::complex<long double>> p1({1, 2, 3}), p2({-4, 3, 0, 6});

    std::cout << (p1 * p2);

    return 0;
}