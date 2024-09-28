#pragma once

#include <random>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <deque>
#include <queue>
#include "templates.hpp"

namespace sortcpp {
    namespace selection {
        // Asynchronous Sort
        template<class RandomAccessIterator, class Compare>
        void asynchronous_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::insertion_sort;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

            int length = last - first;
            auto min = *first, max = *first;
            std::vector<T> ext(length);

            for (auto i = first; i < last; i++) {
                ext[i - first] = *i;
                if (comp(*i, min)) min = *i;
                if (comp(max, *i)) max = *i;
            }
            max++;

            auto cur = min;
            auto i = first;

            while (i < last) {
                for (int j = 0; j < length; j++)
                    if (!comp(cur, ext[j])) {
                        *i = ext[j];
                        ext[j] = max;
                        i++;
                    }
                cur++;
            }

            insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void asynchronous_sort(RandomAccessIterator first, RandomAccessIterator last) {
            asynchronous_sort(first, last, std::less<void>());
        }



        // Base-N Max Heap Sort
        template<class RandomAccessIterator, class Compare>
        void base_n_max_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, int base = 4) {
            using sortcpp::_internal::base_n_heap_sorting::base_n_heap_sort;
            base_n_heap_sort(first, last, base, comp);
        }

        template<class RandomAccessIterator>
        void base_n_max_heap_sort(RandomAccessIterator first, RandomAccessIterator last, int base = 4) {
            base_n_max_heap_sort(first, last, std::less<void>());
        }


        // Bingo Sort
        template<class RandomAccessIterator, class Compare>
        void bingo_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::equals;

            auto max = last - 1;
            auto nxt = *max;

            for (auto i = max - 1; i >= first; i--)
                if (*i > nxt) nxt = *i;

            while (max > first && *max == nxt) max--;

            while (max > first) {
                auto val = nxt;
                nxt = *max;

                for (auto j = max - 1; j >= first; j--) {
                    if (equals(*j, val, comp)) {
                        std::iter_swap(j, max);
                        max--;
                    }
                    else if (comp(nxt, *j)) nxt = *j;
                }

                while (max > first && *max == nxt) max--;
            }
        }

        template<class RandomAccessIterator>
        void bingo_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bingo_sort(first, last, std::less<void>());
        }


        // Binomial Heap Sort
        template<class RandomAccessIterator, class Compare>
        void binomial_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binomial_heap_sorting::binomial_heap_sort;
            binomial_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void binomial_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            binomial_heap_sort(first, last, std::less<void>());
        }


        // Binomial Smooth Sort
        template<class RandomAccessIterator, class Compare>
        void binomial_smooth_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binomial_smooth_sorting::binomial_smooth_sort;
            binomial_smooth_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void binomial_smooth_sort(RandomAccessIterator first, RandomAccessIterator last) {
            binomial_smooth_sort(first, last, std::less<void>());
        }


        // Bottom-up Heap Sort
        template<class RandomAccessIterator, class Compare>
        void bottom_up_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::bottom_up_heap_sorting::bottom_up_heap_sort;
            bottom_up_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bottom_up_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bottom_up_heap_sort(first, last, std::less<void>());
        }


        // Classic Tournament Sort
        template<class RandomAccessIterator, class Compare>
        void classic_tournament_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            using sortcpp::_internal::classic_tournament_sorting::tournament_sort;
            tournament_sort<T>(first, last, comp);
        }

        template<class RandomAccessIterator>
        void classic_tournament_sort(RandomAccessIterator first, RandomAccessIterator last) {
            classic_tournament_sort(first, last, std::less<void>());
        }


        // Cycle Sort
        template<class RandomAccessIterator, class Compare>
        void cycle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::cycle_sorting::cycle_sort;
            cycle_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void cycle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cycle_sort(first, last, std::less<void>());
        }


        // Double Selection Sort
        template<class RandomAccessIterator, class Compare>
        void double_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto left = first, right = last - 1;
            auto smallest = first, biggest = first;

            while (left <= right) {
                for (auto i = left; i <= right; i++) {
                    if (comp(*biggest, *i)) biggest = i;
                    if (comp(*i, *smallest)) smallest = i;
                }

                if (biggest == left) biggest = smallest;

                std::iter_swap(left, smallest);
                std::iter_swap(right, biggest);

                left++;
                right--;

                smallest = left;
                biggest = right;
            }
        }

        template<class RandomAccessIterator>
        void double_selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            double_selection_sort(first, last, std::less<void>());
        }


        // Flipped Min Heap Sort
        template<class RandomAccessIterator, class Compare>
        void filpped_min_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::filpped_min_heap_sorting::filpped_heap_sort;
            filpped_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void filpped_min_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            filpped_min_heap_sort(first, last, std::less<void>());
        }


        // Lazy Heap Sort
        template<class RandomAccessIterator, class Compare>
        void lazy_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::lazy_heap_sorting::lazyheap_sort;
            lazyheap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void lazy_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            lazy_heap_sort(first, last, std::less<void>());
        }


        // Max Heap Sort
        template<class RandomAccessIterator, class Compare>
        void max_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::heap_sorting::heap_sort;
            heap_sort(first, last, true, comp);
        }

        template<class RandomAccessIterator>
        void max_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            max_heap_sort(first, last, std::less<void>());
        }


        // Min Heap Sort
        template<class RandomAccessIterator, class Compare>
        void min_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::heap_sorting::heap_sort;
            heap_sort(first, last, false, comp);
        }

        template<class RandomAccessIterator>
        void min_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            min_heap_sort(first, last, std::less<void>());
        }


        // Max-Min Heap Sort
        template<class RandomAccessIterator, class Compare>
        void min_max_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::min_max_heap_sorting::min_max_heap_sort;
            min_max_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void min_max_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            min_max_heap_sort(first, last, std::less<void>());
        }


        // Poplar Heap Sort
        template<class RandomAccessIterator, class Compare>
        void poplar_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::poplar_heap_sorting::make_poplar_heap;
            using sortcpp::_internal::poplar_heap_sorting::sort_poplar_heap;

            make_poplar_heap(first, last, comp);
            sort_poplar_heap(first, last, comp);
        }

        template<class RandomAccessIterator>
        void poplar_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            poplar_heap_sort(first, last, std::less<void>());
        }


        // Selection Sort
        template<class RandomAccessIterator, class Compare>
        void selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i < last - 1; i++) {
                auto lowest_index = i;

                for (auto j = i + 1; j < last; j++)
                    if (comp(*j, *lowest_index)) lowest_index = j;

                std::iter_swap(i, lowest_index);
            }
        }

        template<class RandomAccessIterator>
        void selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            selection_sort(first, last, std::less<void>());
        }


        // Smooth Sort
        template<class RandomAccessIterator, class Compare>
        void smooth_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::smooth_sorting::smooth_sort;
            smooth_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void smooth_sort(RandomAccessIterator first, RandomAccessIterator last) {
            smooth_sort(first, last, std::less<void>());
        }


        // Stable Cycle Sort
        template<class RandomAccessIterator, class Compare>
        void stable_cycle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::cycle_sorting::stable_cycle_sort;
            stable_cycle_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void stable_cycle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_cycle_sort(first, last, std::less<void>());
        }


        // Stable Selection Sort
        template<class RandomAccessIterator, class Compare>
        void stable_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i < last - 1; i++) {
                auto min = i;
                for (auto j = i + 1; j < last; j++)
                    if (comp(*j, *min)) min = j;

                auto tmp = *min;
                auto pos = min;
                while (pos > i) {
                    *pos = *(pos - 1);
                    pos--;
                }
                *pos = tmp;
            }
        }

        template<class RandomAccessIterator>
        void stable_selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_selection_sort(first, last, std::less<void>());
        }


        // Ternary Heap Sort
        template<class RandomAccessIterator, class Compare>
        void ternary_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::ternary_heap_sorting::ternary_heap_sort;
            ternary_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void ternary_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            ternary_heap_sort(first, last, std::less<void>());
        }


        // Tournament Sort
        template<class RandomAccessIterator, class Compare>
        void tournament_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            using sortcpp::_internal::tournament_sorting::tournament_sort;
            tournament_sort<T>(first, last, comp);
        }

        template<class RandomAccessIterator>
        void tournament_sort(RandomAccessIterator first, RandomAccessIterator last) {
            tournament_sort(first, last, std::less<void>());
        }


        // Triangular Heap Sort
        template<class RandomAccessIterator, class Compare>
        void triangular_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::triangular_heap_sorting::triangular_heap_sort;
            triangular_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void triangular_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            triangular_heap_sort(first, last, std::less<void>());
        }


        // Weak Heap Sort
        template<class RandomAccessIterator, class Compare>
        void weak_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::weakheap_sorting::weakheap_sort;
            weakheap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void weak_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            weak_heap_sort(first, last, std::less<void>());
        }


        // Extra Sorts

        
        // Cocktail Peel Sort
        template<class RandomAccessIterator, class Compare>
        void cocktail_peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::peel_sorting::cocktail_peel_sort;
            cocktail_peel_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void cocktail_peel_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cocktail_peel_sort(first, last, std::less<void>());
        }


        // Dequeue Sort
        template<class RandomAccessIterator, class Compare>
        void dequeue_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::queue_sorting::deque_sort;
            deque_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void dequeue_sort(RandomAccessIterator first, RandomAccessIterator last) {
            dequeue_sort(first, last, std::less<void>());
        }


        // Ecolo Sort
        template<class RandomAccessIterator, class Compare>
        void ecolo_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto left = first + 1, right = last, i = first + 1;
            int way = 1;

            while (left < right) {
                if (way == 1) i = left;
                else i = right;

                while ((way == 1 && i < right) || (way == -1 && i > left)) {
                    if (comp(*(i - 1), *(left - 1))) std::iter_swap(left - 1, i - 1);
                    if (comp(*(right - 1), *(i - 1))) std::iter_swap(right - 1, i - 1);
                    i += way;
                }

                left++;
                right--;
                way *= -1;
            }
        }

        template<class RandomAccessIterator>
        void ecolo_sort(RandomAccessIterator first, RandomAccessIterator last) {
            ecolo_sort(first, last, std::less<void>());
        }


        // Fall Sort
        template<class RandomAccessIterator, class Compare>
        void fall_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fall_sorting::fall_sort;
            fall_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void fall_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fall_sort(first, last, std::less<void>());
        }


        // Forced Stable Heap Sort
        template<class RandomAccessIterator, class Compare>
        void forced_stable_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::forced_stable_heap_sorting::stable_heap_sort;
            stable_heap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void forced_stable_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            forced_stable_heap_sort(first, last, std::less<void>());
        }


        // Improved Multi Selection Sort
        template<class RandomAccessIterator, class Compare>
        void improved_multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::multi_selection_sorting::improved_multi_selection_sort;
            improved_multi_selection_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void improved_multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            improved_multi_selection_sort(first, last, std::less<void>());
        }

        // Multi Selection Sort
        template<class RandomAccessIterator, class Compare>
        void multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::multi_selection_sorting::multi_selection_sort;
            multi_selection_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            multi_selection_sort(first, last, std::less<void>());
        }


        // Optimized Lazy Heap Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_lazy_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::lazy_heap_sorting::optimized_lazyheap_sort;
            optimized_lazyheap_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_lazy_heap_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_lazy_heap_sort(first, last, std::less<void>());
        }


        // Peel Sort
        template<class RandomAccessIterator, class Compare>
        void peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::peel_sorting::peel_sort;
            peel_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void peel_sort(RandomAccessIterator first, RandomAccessIterator last) {
            peel_sort(first, last, std::less<void>());
        }


        // Queue Sort
        template<class RandomAccessIterator, class Compare>
        void queue_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::queue_sorting::queue_sort;
            queue_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void queue_sort(RandomAccessIterator first, RandomAccessIterator last) {
            queue_sort(first, last, std::less<void>());
        }


        // Reverse Peel Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::peel_sorting::reverse_peel_sort;
            reverse_peel_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_peel_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_peel_sort(first, last, std::less<void>());
        }


        // Reverse Sandpaper Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::sandpaper_sorting::reverse_sandpaper_sort;
            reverse_sandpaper_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_sandpaper_sort(first, last, std::less<void>());
        }


        // Reverse Selection Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = last - 1; i >= first; i--) {
                auto highest_index = first;

                for (auto j = first + 1; j < i + 1; j++) {
                    if (comp(*highest_index, *j)) highest_index = j;
                }
                std::iter_swap(i, highest_index);
            }
        }

        template<class RandomAccessIterator>
        void reverse_selection_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_selection_sort(first, last, std::less<void>());
        }


        // Sandpaper Sort
        template<class RandomAccessIterator, class Compare>
        void sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::sandpaper_sorting::sandpaper_sort;
            sandpaper_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last) {
            sandpaper_sort(first, last, std::less<void>());
        }


        // Stable Fall Sort
        template<class RandomAccessIterator, class Compare>
        void stable_fall_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fall_sorting::stable_fall_sort;
            stable_fall_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void stable_fall_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_fall_sort(first, last, std::less<void>());
        }
    }
}