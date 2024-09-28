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
    namespace merge {
        // Andrey Sort
        template<class RandomAccessIterator, class Compare>
        void andrey_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::andrey_sorting::andrey_sort;
            andrey_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void andrey_sort(RandomAccessIterator first, RandomAccessIterator last) {
            andrey_sort(first, last, std::less<void>());
        }


        // Block Swap Merge Sort
        template<class RandomAccessIterator, class Compare>
        void block_swap_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::multi_swap_merge_sort;
            multi_swap_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void block_swap_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            block_swap_merge_sort(first, last, std::less<void>());
        }


        // Bottom-Up Merge Sort
        template<class RandomAccessIterator, class Compare>
        void bottom_up_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::bottom_up_merge_sort;
            bottom_up_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bottom_up_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bottom_up_merge_sort(first, last, std::less<void>());
        }


        // Buffered Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void buffered_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::buffered_stooge_sort;
            buffered_stooge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void buffered_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            buffered_stooge_sort(first, last, std::less<void>());
        }


        // Improved In-Place Merge Sort
        template<class RandomAccessIterator, class Compare>
        void improved_inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::inplace_merge_sorting::improved_inplace_merge_sort;
            improved_inplace_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void improved_inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            improved_inplace_merge_sort(first, last, std::less<void>());
        }


        // In-Place Merge Sort
        template<class RandomAccessIterator, class Compare>
        void inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::inplace_merge_sorting::inplace_merge_sort;
            inplace_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            inplace_merge_sort(first, last, std::less<void>());
        }


        // Iterative Top-Down Merge Sort
        template<class RandomAccessIterator, class Compare>
        void iterative_top_down_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::iterative_top_down_merge_sorting::merge_sort;
            merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void iterative_top_down_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            iterative_top_down_merge_sort(first, last, std::less<void>());
        }


        // TODO: Lazy Stable Sort


        // Merge Sort
        template<class RandomAccessIterator, class Compare>
        void merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::merge_sort;
            merge_sort(first, last, false, comp);
        }

        template<class RandomAccessIterator>
        void merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            merge_sort(first, last, std::less<void>());
        }



        // New Shuffle Merge Sort
        template<class RandomAccessIterator, class Compare>
        void new_shuffle_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::new_shuffle_merge_sorting::merge_sort;
            merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void new_shuffle_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            new_shuffle_merge_sort(first, last, std::less<void>());
        }


        // Pattern-Defeating Merge Sort (Simplified Tim Sort)
        template<class RandomAccessIterator, class Compare>
        void pd_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::pd_merge_sorting::pd_merge_sort;
            pd_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pd_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pd_merge_sort(first, last, std::less<void>());
        }

        
        // TODO: Quad Sort


        // Rotate Merge Sort
        template<class RandomAccessIterator, class Compare>
        void rotate_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::rotate_merge_sort;
            rotate_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void rotate_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            rotate_merge_sort(first, last, std::less<void>());
        }

        
        // Strand Sort
        template<class RandomAccessIterator, class Compare>
        void strand_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::strand_sort;
            strand_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void strand_sort(RandomAccessIterator first, RandomAccessIterator last) {
            strand_sort(first, last, std::less<void>());
        }


        // Twin Sort
        template<class RandomAccessIterator, class Compare>
        void twin_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::twin_sorting::twin_sort;
            twin_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void twin_sort(RandomAccessIterator first, RandomAccessIterator last) {
            twin_sort(first, last, std::less<void>());
        }


        // Weaved Merge Sort
        template<class RandomAccessIterator, class Compare>
        void weaved_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::weaved_merge_sort;
            weaved_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void weaved_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            weaved_merge_sort(first, last, std::less<void>());
        }


        // Extra Sorts

        // Index Merge Sort
        template<class RandomAccessIterator, class Compare>
        void index_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::index_merge_sorting::index_merge_sort;
            index_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void index_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            index_merge_sort(first, last, std::less<void>());
        }


        // Mob Merge Sort
        template<class RandomAccessIterator, class Compare>
        void mob_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::mob_merge_sort;
            mob_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void mob_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            mob_merge_sort(first, last, std::less<void>());
        }


        // Natural Merge Sort
        template<class RandomAccessIterator, class Compare>
        void natural_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::natural_merge_sort;
            natural_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void natural_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            natural_merge_sort(first, last, std::less<void>());
        }


        // Optimized Pancake Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_pancake_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::pancake_sorting::pancake_merge_sort;
            pancake_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_pancake_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_pancake_sort(first, last, std::less<void>());
        }


        // Partial Merge Sort
        template<class RandomAccessIterator, class Compare>
        void partial_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::merge_sorting::partial_merge_sort;
            partial_merge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void partial_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            partial_merge_sort(first, last, std::less<void>());
        }
    }
}