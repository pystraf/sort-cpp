#pragma once

#include <random>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <stack>
#include <queue>
#include <utility>
#include "templates.hpp"

namespace sortcpp {
    namespace insertion {
        // Double Insertion Sort With Binary Search
        template<class RandomAccessIterator, class Compare>
        void binary_double_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::binary_double_insertion_sort;
            binary_double_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void binary_double_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            binary_double_insertion_sort(first, last, std::less<void>());
        }


        // Insertion Sort With Binary Search
        template<class RandomAccessIterator, class Compare>
        void binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::binary_insertion_sort;
            binary_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            binary_insertion_sort(first, last, std::less<void>());
        }

        // TODO: Block Insertion Sort (Will write after implementing grail sort)


        // Double Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void double_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::double_insertion_sort;
            double_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void double_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            double_insertion_sort(first, last, std::less<void>());
        }


        // Inserion Sort
        template<class RandomAccessIterator, class Compare>
        void insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::insertion_sort;
            insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            insertion_sort(first, last, std::less<void>());
        }


        // Library Sort
        template<class RandomAccessIterator, class Compare>
        void library_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::library_sorting::library_sort;
            library_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void library_sort(RandomAccessIterator first, RandomAccessIterator last) {
            library_sort(first, last, std::less<void>());
        }


        // Patience Sort
        template<class RandomAccessIterator, class Compare>
        void patience_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::patience_sorting::patience_sort;
            patience_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void patience_sort(RandomAccessIterator first, RandomAccessIterator last) {
            patience_sort(first, last, std::less<void>());
        }


        // Recursive Shell Sort
        template<class RandomAccessIterator, class Compare>
        void recursive_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::shell_sorting::recursive_shell_sort;
            recursive_shell_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void recursive_shell_sort(RandomAccessIterator first, RandomAccessIterator last) {
            recursive_shell_sort(first, last, std::less<void>());
        }


        // Shell Sort
        template<class RandomAccessIterator, class Compare>
        void shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::shell_sorting::shell_sort;
            shell_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void shell_sort(RandomAccessIterator first, RandomAccessIterator last) {
            shell_sort(first, last, std::less<void>());
        }


        // Simplified Library Sort
        template<class RandomAccessIterator, class Compare>
        void simplified_library_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::simplified_library_sorting::library_sort;
            library_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void simplified_library_sort(RandomAccessIterator first, RandomAccessIterator last) {
            simplified_library_sort(first, last, std::less<void>());
        }
    }
}