#pragma once

#include <random>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include "templates.hpp"

namespace sortcpp {
    namespace misc {
        // Brunt Pancake Sort
        template<class RandomAccessIterator, class Compare>
        void brunt_pancake_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::pancake_sorting::brunt_pancake_sort;
            brunt_pancake_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void brunt_pancake_sort(RandomAccessIterator first, RandomAccessIterator last) {
            brunt_pancake_sort(first, last, std::less<void>());
        }


        // Pancake Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void pancake_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::pancake_sorting::pancake_insertion_sort;
            pancake_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pancake_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pancake_insertion_sort(first, last, std::less<void>());
        }

        template<class RandomAccessIterator, class Compare>
        void pancake_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::pancake_sorting::pancake_sort;
            pancake_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pancake_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pancake_sort(first, last, std::less<void>());
        }
    }
}