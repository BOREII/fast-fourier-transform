#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>
#include <numeric>

// Задачи, решаемые с помощью умножения многочленов
// Если вы напишете решение, работающее для произольных строк над ascii кодировкой - укажете это и
// возможно, получите небольшой бонусный балл
namespace SubstringMatching {
//Квадратичный алгоритм
std::vector<size_t> FindSubstrings(const std::string& str, const std::string& pattern);

//С помощью ффт
std::vector<size_t> FindSubstringsFFT(const std::string& str, const std::string& pattern);

//Квадратичный алгоритм
std::vector<size_t> FindMatches(const std::string& str, const std::string& pattern);
//С помощью ффт
std::vector<size_t> FindMatchesFFT(const std::string& str, const std::string& pattern);

//Тесты
void TestSubstrings();
void TestSubstringsHeyJude();
void TestMatches();
void TestMatchesHeyJude();
} // namespace SubstringMatching