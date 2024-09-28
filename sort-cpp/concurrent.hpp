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
    namespace concurrent {
        // Bitonic Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void bitonic_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::bitonic_sort_iterative;
            bitonic_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bitonic_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            bitonic_sort_iterative(first, last, std::less<void>());
        }


        // Bitonic Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void bitonic_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::bitonic_sort_recursive;
            bitonic_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bitonic_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            bitonic_sort_recursive(first, last, std::less<void>());
        }


        // Bose Nelson Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void bose_nelson_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::bose_nelson_sort_iterative;
            bose_nelson_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bose_nelson_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            bose_nelson_sort_iterative(first, last, std::less<void>());
        }


        // Bose Nelson Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void bose_nelson_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::bose_nelson_sort_recursive;
            bose_nelson_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bose_nelson_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            bose_nelson_sort_recursive(first, last, std::less<void>());
        }


        // Crease Sort
        template<class RandomAccessIterator, class Compare>
        void crease_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::crease_sort;
            crease_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void crease_sort(RandomAccessIterator first, RandomAccessIterator last) {
            crease_sort(first, last, std::less<void>());
        }


        // Diamond Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void diamond_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::diamond_sort_iterative;
            diamond_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void diamond_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            diamond_sort_iterative(first, last, std::less<void>());
        }


        // Diamond Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void diamond_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::diamond_sort_recursive;
            diamond_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void diamond_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            diamond_sort_recursive(first, last, std::less<void>());
        }


        // Fold Sort
        template<class RandomAccessIterator, class Compare>
        void fold_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::fold_sort;
            fold_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void fold_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fold_sort(first, last, std::less<void>());
        }


        // Matrix Sort
        template<class RandomAccessIterator, class Compare>
        void matrix_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::matrix_sort;
            matrix_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void matrix_sort(RandomAccessIterator first, RandomAccessIterator last) {
            matrix_sort(first, last, std::less<void>());
        }


        // Odd-Even Merge Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void odd_even_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::odd_even_merge_sort_iterative;
            odd_even_merge_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void odd_even_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_merge_sort_iterative(first, last, std::less<void>());
        }


        // Odd-Even Merge Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void odd_even_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::odd_even_merge_sort_recursive;
            odd_even_merge_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void odd_even_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_merge_sort_recursive(first, last, std::less<void>());
        }


        // Pairwise Merge Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void pairwise_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::pairwise_merge_sort_iterative;
            pairwise_merge_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pairwise_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_merge_sort_iterative(first, last, std::less<void>());
        }


        // Pairwise Merge Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void pairwise_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::pairwise_merge_sort_recursive;
            pairwise_merge_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pairwise_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_merge_sort_recursive(first, last, std::less<void>());
        }


        // Pairwise Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void pairwise_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::pairwise_sort_iterative;
            pairwise_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pairwise_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_sort_iterative(first, last, std::less<void>());
        }


        // Pairwise Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void pairwise_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::pairwise_sort_recursive;
            pairwise_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pairwise_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_sort_recursive(first, last, std::less<void>());
        }

        
        // Weave Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void weave_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::weave_sort_iterative;
            weave_sort_iterative(first, last, comp);
        }

        template<class RandomAccessIterator>
        void weave_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            weave_sort_iterative(first, last, std::less<void>());
        }


        // Weave Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void weave_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::weave_sort_recursive;
            weave_sort_recursive(first, last, comp);
        }

        template<class RandomAccessIterator>
        void weave_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            weave_sort_recursive(first, last, std::less<void>());
        }


        // Apollyon Sort
        template<class RandomAccessIterator, class Compare>
        void apollyon_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::apollyon_sort;
            apollyon_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void apollyon_sort(RandomAccessIterator first, RandomAccessIterator last) {
            apollyon_sort(first, last, std::less<void>());
        }


        // Optimized Odd-Even Merge Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_odd_even_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::concurrent_sorting::optimized_oddeven_merge_sort;
            optimized_oddeven_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_odd_even_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_odd_even_merge_sort(first, last, std::less<void>());
        }
    }
}
