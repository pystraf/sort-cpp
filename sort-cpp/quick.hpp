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
#include "templates.hpp"

namespace sortcpp {
    namespace quick {
        // Ternary Quick Sort (LL ptrs)
        template<class RandomAccessIterator, class Compare>
        void ternary_quick_sort_ll(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::ternary_quick_sorting::ternary_quick_sort_ll;
            ternary_quick_sort_ll(first, last, comp);
        }

        template<class RandomAccessIterator>
        void ternary_quick_sort_ll(RandomAccessIterator first, RandomAccessIterator last) {
            ternary_quick_sort_ll(first, last, std::less<void>());
        }


        // Ternary Quick Sort (LR ptrs)
        template<class RandomAccessIterator, class Compare>
        void ternary_quick_sort_lr(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::ternary_quick_sorting::ternary_quick_sort_lr;
            ternary_quick_sort_lr(first, last, comp);
        }

        template<class RandomAccessIterator>
        void ternary_quick_sort_lr(RandomAccessIterator first, RandomAccessIterator last) {
            ternary_quick_sort_lr(first, last, std::less<void>());
        }


        // Cube Root Quick Sort
        template<class RandomAccessIterator, class Compare>
        void cube_root_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::quick_sorting::cube_root_quick_sort;
            cube_root_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void cube_root_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cube_root_quick_sort(first, last, std::less<void>());
        }


        // Shell Unstable Singularity Quick Sort
        template<class RandomAccessIterator, class Compare>
        void shell_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::singularity_quick_sorting::shell_unstable_singularity_quick_sort;
            shell_unstable_singularity_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void shell_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            shell_unstable_singularity_quick_sort(first, last, std::less<void>());
        }


        // Singularity Quick Sort
        template<class RandomAccessIterator, class Compare>
        void singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::singularity_quick_sorting::singularity_quick_sort;
            singularity_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            singularity_quick_sort(first, last, std::less<void>());
        }


        // Unbounded Singularity Quick Sort
        template<class RandomAccessIterator, class Compare>
        void unbounded_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::singularity_quick_sorting::unbounded_singularity_quick_sort;
            unbounded_singularity_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void unbounded_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unbounded_singularity_quick_sort(first, last, std::less<void>());
        }


        // Unstable Singularity Quick Sort
        template<class RandomAccessIterator, class Compare>
        void unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::singularity_quick_sorting::unstable_singularity_quick_sort;
            unstable_singularity_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unstable_singularity_quick_sort(first, last, std::less<void>());
        }


        // Unbounded Unstable Singularity Quick Sort
        template<class RandomAccessIterator, class Compare>
        void unbounded_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::singularity_quick_sorting::unbounded_unstable_singularity_quick_sort;
            unbounded_unstable_singularity_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void unbounded_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unbounded_unstable_singularity_quick_sort(first, last, std::less<void>());
        }
    }
}