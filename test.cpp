#include "sort-cpp/insertion.hpp"
#include "sort-cpp/impractical.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace sortcpp;



int main() {
	std::random_device rd;
    std::mt19937 gen(rd());
	std::vector<int> a{2, 0, 134, 1, 8, 7, 6, 1145, 1919, 801, 3, 2, 77, 45, 12, 34, 66, 111, -114514, 4326, 3215, 2151324, 3152, 21154136, 426547543, 33};
	// std::sort(a.begin(), a.end());
	// auto key = [&](int i) { return i + 114514; };
	insertion::simplified_library_sort(a.begin(), a.end());
	// impractical::hyper_stooge_sort(a.begin(), a.end());
	for (auto &i: a) std::cout << i << ' ';
	std::cout << std::endl;
	// std::cout << a.size();
	return 0;
}