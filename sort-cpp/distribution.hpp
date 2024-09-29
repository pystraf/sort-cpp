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
    namespace distribution {

        // American Flag Sort
        template<class RandomAccessIterator, class Key>
        void american_flag_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int buckets = 128) {
            using sortcpp::_internal::american_flag_sorting::american_flag_sort;
            american_flag_sort(first, last, buckets, key);
        }


        // Binary Quick Sort (Iterative)
        template<class RandomAccessIterator, class Key>
        void binary_quick_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::binary_quick_sorting::binary_quick_sort;
            using sortcpp::_internal::binary_quick_sorting::get_max_bit;

            int max_bit = get_max_bit(first, last, key);
            binary_quick_sort(first, last - 1, max_bit - 1, key);
        }
        

        // Binary Quick Sort (Recursive)
        template<class RandomAccessIterator, class Key>
        void binary_quick_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::binary_quick_sorting::binary_quick_sort_recursive;
            using sortcpp::_internal::binary_quick_sorting::get_max_bit;

            int max_bit = get_max_bit(first, last, key);
            binary_quick_sort_recursive(first, last - 1, max_bit - 1, key);
        }
        

        // Counting Sort
        template<class RandomAccessIterator, class Key>
        void counting_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::counting_sorting::counting_sort;
            counting_sort(first, last, key);
        }


        // Flash Sort
        template<class RandomAccessIterator, class Key>
        void flash_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::flash_sorting::flash_sort;
            flash_sort(first, last, key);
        }


        // In-Place LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void inplace_lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::radix_sorting::inplace_lsd_radix_sort;
            inplace_lsd_radix_sort(first, last, base, key);
        }


        // LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::radix_sorting::lsd_radix_sort;
            lsd_radix_sort(first, last, base, key);
        }

        // MSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void msd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::radix_sorting::msd_radix_sort;
            msd_radix_sort(first, last, base, key);
        }


        // Rotate LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void rotate_lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::radix_sorting::rotate_lsd_radix_sort;
            rotate_lsd_radix_sort(first, last, base, key);
        }


        // Rotate MSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void rotate_msd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::radix_sorting::rotate_msd_radix_sort;
            rotate_msd_radix_sort(first, last, base, key);
        }

        
        // Stackless American Flag Sort
        template<class RandomAccessIterator, class Key>
        void stackless_american_flag_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int buckets = 128) {
            using sortcpp::_internal::stackless_american_flag_sorting::american_flag_sort;
            american_flag_sort(first, last, buckets, key);
        }


        // Stackless Binary Quick Sort
        template<class RandomAccessIterator, class Key>
        void stackless_binary_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::stackless_binary_quick_sorting::binary_quick_sort;
            binary_quick_sort(first, last, key);
        }


        // Static Sort
        template<class RandomAccessIterator, class Key>
        void static_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::static_sorting::static_sort;
            static_sort(first, last, key);
        }
    }
}
