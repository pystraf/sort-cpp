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
    namespace exchange {
        // Gnome Sort with Binary Search
        template<class RandomAccessIterator, class Compare>
        void binary_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first + 1; i < last; i++) {
                auto lo = first, hi = i;
                while (lo < hi) {
                    auto mid = lo + (hi - lo) / 2;
                    if (comp(*i, *mid)) hi = mid;
                    else lo = mid + 1;
                }
                
                auto j = i;
                while (j > lo) {
                    std::iter_swap(j, j - 1);
                    j--;
                }
            }
        }
        
        template<class RandomAccessIterator>
        void binary_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            binary_gnome_sort(first, last, std::less<void>());
        }
        
        
        // Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = last - 1; i > first; i--) {
                bool sorted = true;
                for (auto j = first; j < i; j++)
                    if (comp(*(j + 1), *j)) {
                        std::iter_swap(j, j + 1);
                        sorted = false;
                    }
                if (sorted) break;
            }
        }
        
        template<class RandomAccessIterator>
        void bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bubble_sort(first, last, std::less<void>());
        }
        
        
        // Circle Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void iterative_circle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::iterative_circle_sorting::circle_sort_routine;
            int length = last - first;
            int size = 1;
            for (; size < length; size *= 2);
            
            bool swapped = false;
            do {
                swapped = circle_sort_routine(first, last, size, comp);
            } while (swapped);
        }
        
        template<class RandomAccessIterator>
        void iterative_circle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            iterative_circle_sort(first, last, std::less<void>());
        }
        
        
        // Circle Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void recursive_circle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::circle_sorting::circle_sort_routine;
            int length = last - first;
            int size = 1;
            for (; size < length; size *= 2);
            
            bool swapped = false;
            do {
                swapped = circle_sort_routine(first, first + size - 1, last, false, comp);
            } while (swapped);
        }
        
        template<class RandomAccessIterator>
        void recursive_circle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            recursive_circle_sort(first, last, std::less<void>());
        }
        
        
        // Circloid Sort
        template<class RandomAccessIterator, class Compare>
        void circloid_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::circloid_sorting::circle_pass;
            while (circle_pass(first, last - 1, comp));
        }
        
        template<class RandomAccessIterator>
        void circloid_sort(RandomAccessIterator first, RandomAccessIterator last) {
            circloid_sort(first, last, std::less<void>());
        }
        
        
        // Classic 3-Smooth Comb Sort
        template<class RandomAccessIterator, class Compare>
        void classic_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::three_smooth_comb_sorting::classic_3smooth_comb_sort;
            classic_3smooth_comb_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void classic_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last) {
            classic_3smooth_comb_sort(first, last, std::less<void>());
        }
        
        
        // Cocktail Shaker Sort
        template<class RandomAccessIterator, class Compare>
        void cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto i = first;
            auto mid = first + (last - first) / 2;
            while (i < mid) {
                bool sorted = true;
                
                auto pos = last - i - 1 + first;
                for (auto j = i; j < pos; j++)
                    if (comp(*(j + 1), *j)) {
                        std::iter_swap(j, j + 1);
                        sorted = false;
                    }
                
                for (auto j = pos; j > i; j--)
                    if (comp(*j, *(j - 1))) {
                        std::iter_swap(j, j - 1);
                        sorted = false;
                    }
                
                if (sorted) break;
                else i++;
            }
        }
        
        template<class RandomAccessIterator>
        void cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cocktail_shaker_sort(first, last, std::less<void>());
        }
        
        
        // Comb Sort
        template<class RandomAccessIterator, class Compare>
        void comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, double shrink = 1.3) {
            using sortcpp::_internal::comb_sorting::comb_sort;
            comb_sort(first, last, comp, shrink, false);
        }
        
        template<class RandomAccessIterator>
        void comb_sort(RandomAccessIterator first, RandomAccessIterator last, double shrink = 1.3) {
            comb_sort(first, last, std::less<void>(), shrink);
        }
    
        
        // Complete Graph Sort
        template<class RandomAccessIterator, class Compare>
        void complete_graph_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::complete_graph_sorting::complete_graph_sort;
            complete_graph_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void complete_graph_sort(RandomAccessIterator first, RandomAccessIterator last) {
            complete_graph_sort(first, last, std::less<void>());
        }
    
        
        // Dual-Pivot Quick Sort
        template<class RandomAccessIterator, class Compare>
        void dual_pivot_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::quick_sorting::dual_pivot_loop;
            dual_pivot_loop(first, last - 1, 3, comp);
        }
        
        template<class RandomAccessIterator>
        void dual_pivot_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            dual_pivot_quick_sort(first, last, std::less<void>());
        }
        
        
        // Forced Stable Quick Sort
        template<class RandomAccessIterator, class Compare>
        void forced_stable_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::forced_stable_quick_sorting::stable_quick_sort;
            stable_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void forced_stable_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            forced_stable_quick_sort(first, last, std::less<void>());
        }
    
    
        // Fun Sort
        template<class RandomAccessIterator, class Compare>
        void fun_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fun_sorting::fun_sort;
            fun_sort(first, last, comp);
        }
    
        template<class RandomAccessIterator>
        void fun_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fun_sort(first, last, std::less<void>());
        }
    
        
        // Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first + 1; i < last; ) {
                if (!comp(*i, *(i - 1))) i++;
                else {
                    std::iter_swap(i, i - 1);
                    if (i > first + 1) i--;
                }
            }
        }
        
        template<class RandomAccessIterator>
        void gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            gnome_sort(first, last, std::less<void>());
        }
    
        // Quick Sort (LL ptrs)
        template<class RandomAccessIterator, class Compare>
        void quick_sort_ll(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::quick_sorting::quick_ll_loop;
            quick_ll_loop(first, last - 1, comp);
        }
        
        template<class RandomAccessIterator>
        void quick_sort_ll(RandomAccessIterator first, RandomAccessIterator last) {
            quick_sort_ll(first, last, std::less<void>());
        }
    
        
        // Quick Sort (LR ptrs)
        template<class RandomAccessIterator, class Compare>
        void quick_sort_lr(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::quick_sorting::quick_lr_loop;
            quick_lr_loop(first, last - 1, comp);
        }
        
        template<class RandomAccessIterator>
        void quick_sort_lr(RandomAccessIterator first, RandomAccessIterator last) {
            quick_sort_lr(first, last, std::less<void>());
        }
    
        
        // Odd-Even Sort (a.k.a Brick Sort)
        template<class RandomAccessIterator, class Compare>
        void odd_even_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            bool sorted = false;
            while (!sorted) {
                sorted = true;
                for (auto i = first + 1; i < last - 1; i += 2)
                    if (comp(*(i + 1), *i)) {
                        std::iter_swap(i, i + 1);
                        sorted = false;
                    }
                for (auto i = first; i < last - 1; i += 2)
                    if (comp(*(i + 1), *i)) {
                        std::iter_swap(i, i + 1);
                        sorted = false;
                    }
            }
        }
    
        template<class RandomAccessIterator>
        void odd_even_sort(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_sort(first, last, std::less<void>());
        }
    
    
        // Optimized Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            int consec_sorted;
            for (auto i = last - 1; i > first; i -= consec_sorted) {
                consec_sorted = 1;
                for (auto j = first; j < i; j++) {
                    if (comp(*(j + 1), *j)) {
                        std::iter_swap(j, j + 1);
                        consec_sorted = 1;
                    } else consec_sorted++;
                }
            }
        }
    
        template<class RandomAccessIterator>
        void optimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_bubble_sort(first, last, std::less<void>());
        }
        
        
        // Optimized Cocktail Shaker Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto start = first, end = last - 1; start < end; ) {
                int consec_sorted = 1;
                for (auto i = start; i < end; i++) {
                    if (comp(*(i + 1), *i)) {
                        std::iter_swap(i, i + 1);
                        consec_sorted = 1;
                    } else consec_sorted++;
                }
                end -= consec_sorted;
                
                consec_sorted = 1;
                for (auto i = end; i > start; i--) {
                    if (comp(*i, *(i - 1))) {
                        std::iter_swap(i - 1, i);
                        consec_sorted = 1;
                    } else consec_sorted++;
                }
                start += consec_sorted;
            }
        }
        
        template<class RandomAccessIterator>
        void optimized_cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_cocktail_shaker_sort(first, last, std::less<void>());
        }
        
        
        // Optimized Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_sorting::smart_gnome_sort;
            smart_gnome_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void optimized_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_gnome_sort(first, last, std::less<void>());
        }
    
        
        // Optimized Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::optimized_stooge_sort;
            optimized_stooge_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void optimized_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_stooge_sort(first, last, std::less<void>());
        }
    
    
        // Optistooge Sort
        template<class RandomAccessIterator, class Compare>
        void optistooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::optimized_stooge_sort_studio;
            optimized_stooge_sort_studio(first, last, comp);
        }
    
        template<class RandomAccessIterator>
        void optistooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optistooge_sort(first, last, std::less<void>());
        }
    
        
        // Slope Sort
        template<class RandomAccessIterator, class Compare>
        void slope_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first + 1, j = first + 1; i < last; i++, j++) {
                for (auto k = i - 1; k >= first; k--, i--)
                    if (comp(*i, *k)) std::iter_swap(i, k);
                i = j;
            }
        }
    
        template<class RandomAccessIterator>
        void slope_sort(RandomAccessIterator first, RandomAccessIterator last) {
            slope_sort(first, last, std::less<void>());
        }
        
        
        
        // Stable Quick Sort
        template<class RandomAccessIterator, class Compare>
        void stable_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            using sortcpp::_internal::stable_quick_sorting::stable_quick_loop;
            stable_quick_loop<T>(first, last - 1, false, comp);
        }
        
        template<class RandomAccessIterator>
        void stable_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_quick_sort(first, last, std::less<void>());
        }
        
        
        // Swapless Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void swapless_bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            RandomAccessIterator pre;
            for (auto i = last; i > first; i = pre) {
                pre = first;
                auto pos = first;
                auto tmp = *first;
                for (auto j = first + 1; j < i; j++) {
                    if (comp(*j, tmp)) {
                        *(j - 1) = *j;
                        pre = j;
                    } else {
                        if (pos + 1 < j) *(j - 1) = tmp;
                        pos = j;
                        tmp = *j;
                    }
                }
                *(i - 1) = tmp;
            }
        }
        
        template<class RandomAccessIterator>
        void swapless_bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            swapless_bubble_sort(first, last, std::less<void>());
        }
    
    
        // Table Sort
        template<class RandomAccessIterator, class Compare>
        void table_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::table_sorting::table_sort;
            table_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void table_sort(RandomAccessIterator first, RandomAccessIterator last) {
            table_sort(first, last, std::less<void>());
        }
        
        
        // 3-Smooth Comb Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void iterative_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::three_smooth_comb_sorting::iterative_3smooth_comb_sort;
            iterative_3smooth_comb_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void iterative_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last) {
            iterative_3smooth_comb_sort(first, last, std::less<void>());
        }
        
        
        
        // 3-Smooth Comb Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void recursive_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::three_smooth_comb_sorting::recursive_3smooth_comb_sort;
            recursive_3smooth_comb_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void recursive_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last) {
            recursive_3smooth_comb_sort(first, last, std::less<void>());
        }
        
        
        // Unoptimized Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void unoptimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            bool sorted = false;
            
            while (!sorted) {
                sorted = true;
                for (auto i = first; i < last - 1; i++) 
                    if (comp(*(i + 1), *i)) {
                        std::iter_swap(i, i + 1);
                        sorted = false;
                    }	
            }
        }
        
        template<class RandomAccessIterator>
        void unoptimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unoptimized_bubble_sort(first, last, std::less<void>());
        }
        
        
        // Unoptimized Cocktail Shaker Sort
        template<class RandomAccessIterator, class Compare>
        void unoptimized_cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto i = first;
            auto mid = first + (last - first) / 2;
            while (i < mid) {
                for (auto j = i; j < last - i - 1 + first; j++) 
                    if (comp(*(j + 1), *j)) std::iter_swap(j, j + 1);
                
                for (auto j = last - i - 1 + first; j > i; j--)
                    if (comp(*j, *(j - 1))) std::iter_swap(j, j - 1);
                
                i++;
            }
        }
        
        template<class RandomAccessIterator>
        void unoptimized_cocktail_shaker_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unoptimized_cocktail_shaker_sort(first, last, std::less<void>());
        }


        // Extra Sorts

        // Asteraceae Sort
        template<class RandomAccessIterator, class Compare>
        void asteraceae_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            int length = last - first;
            
            auto i = first + 1;
            auto first_swap = first + 2;
            bool any_swaps = true, last_swap = false;

            while (any_swaps) {
                if (first_swap - 1 == first) i = first + 1;
                else i = first_swap - 1;
                any_swaps = last_swap = false;

                while (i < last) {
                    if (comp(*i, *(i - 1))) {
                        std::iter_swap(i - 1, i);
                        i++;
                        if (!any_swaps) first_swap = i - 1;
                        any_swaps = last_swap = true;
                    }
                    else {
                        if (last_swap) i += std::floor(std::sqrt(length));
                        else i++;
                        last_swap = false;
                    }
                }
            }
        }

        template<class RandomAccessIterator>
        void asteraceae_sort(RandomAccessIterator first, RandomAccessIterator last) {
            asteraceae_sort(first, last, std::less<void>());
        }


        // Chinotto Sort
        template<class RandomAccessIterator, class Compare>
        void chinotto_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            bool done = false;
            int gap = 1;

            while (!done) {
                auto i = first;
                done = true;

                while (i + gap < last) {
                    if (comp(*(i + gap), *i)) {
                        done = false;
                        for (auto j = i; j < i + gap; j++) std::iter_swap(j, j + 1);
                        gap++;
                    }
                    else if (gap >= 2) gap--;
                    i++;
                }

                while (i - gap > first) {
                    if (comp(*i, *(i - gap))) {
                        done = false;
                        for (auto j = i; j > i - gap; j--) std::iter_swap(j, j - 1);
                        gap++;
                    }
                    else if (gap >= 2) gap--;
                    i--;
                }
            }
        }

        template<class RandomAccessIterator>
        void chinotto_sort(RandomAccessIterator first, RandomAccessIterator last) {
            chinotto_sort(first, last, std::less<void>());
        }


        // Clamber Sort
        template<class RandomAccessIterator, class Compare>
        void clamber_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first + 1; i < last; i++)
                for (auto j = first; j < i; j++)
                    if (comp(*i, *j)) std::iter_swap(i, j);
        }

        template<class RandomAccessIterator>
        void clamber_sort(RandomAccessIterator first, RandomAccessIterator last) {
            clamber_sort(first, last, std::less<void>());
        }


        


        // Cocktail Push Sort
        template<class RandomAccessIterator, class Compare>
        void cocktail_push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            bool any_swaps = true;
            auto i = first + 1;
            int gap = 1;

            while (any_swaps) {
                any_swaps = false;
                i = first + 1;
                gap = 1;

                while (i + gap <= last) {
                    if (comp(*(i - 1 + gap), *(i - 1))) {
                        for (int j = 1; j <= gap; j++) std::iter_swap(i - 1, i - 1 + j);
                        any_swaps = true;
                        gap++;
                    }
                    else i++;
                }

                i = last;
                gap = 1;
                while (i - gap > first) {
                    if (comp(*(i - 1), *(i - 1 - gap))) {
                        for (int j = 1; j <= gap; j++) std::iter_swap(i - 1, i - 1 - j);
                        any_swaps = true;
                        gap++;
                    }
                    else i--;
                }
            }
        }

        template<class RandomAccessIterator>
        void cocktail_push_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cocktail_push_sort(first, last, std::less<void>());
        }


        // Dandelion Sort
        template<class RandomAccessIterator, class Compare>
        void dandelion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto b = first; b < last; ) {
                auto pointer = b;
                bool any_swap = false;

                while (pointer < last - 1 && comp(*(pointer + 1), *pointer)) {
                    std::iter_swap(pointer, pointer + 1);
                    any_swap = true;
                    pointer++;
                }

                if (any_swap) {
                    if (b > first) b--;
                    continue;
                }
                b++;
            }
        }

        template<class RandomAccessIterator>
        void dandelion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            dandelion_sort(first, last, std::less<void>());
        }


        // Fibonacci Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void fibonacci_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::_fib_search;

            for (auto i = first + 1; i < last; i++) {
                auto tmp = *i;
                auto pos = _fib_search(first, i - 1, tmp, comp);
                auto j = i;
                while (j > pos) {
                    std::iter_swap(j, j - 1);
                    j--;
                }
            }
        }

        template<class RandomAccessIterator>
        void fibonacci_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fibonacci_gnome_sort(first, last, std::less<void>());
        }


        // Float Sort
        template<class RandomAccessIterator, class Compare>
        void float_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            bool sorted = false;
            while (!sorted) {
                auto h = first;
                sorted = true;

                for (auto g = last - 1; g > first; g--) {
                    auto i = h;
                    auto j = h + 1;

                    while (i >= first && comp(*j, *i)) {
                        std::iter_swap(i, j);
                        sorted = false;
                        i--;
                        j--;
                    }

                    if (i >= first) {
                        i++;
                        j++;

                        while (j < last && comp(*j, *i)) {
                            std::iter_swap(i, j);
                            sorted = false;
                            i++;
                            j++;
                        }
                    }

                    h++;
                }
            }
        }

        template<class RandomAccessIterator>
        void float_sort(RandomAccessIterator first, RandomAccessIterator last) {
            float_sort(first, last, std::less<void>());
        }


        // Gambit Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void gambit_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_sorting::gambit_gnome_sort;
            gambit_gnome_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void gambit_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            gambit_gnome_sort(first, last, std::less<void>());
        }


        // Gnome Weave Sort (High Prime)
        template<class RandomAccessIterator, class Compare>
        void gnome_weave_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_weave_sorting::gnome_weave_sort;
            gnome_weave_sort(first, last, true, comp);
        }

        template<class RandomAccessIterator>
        void gnome_weave_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last) {
            gnome_weave_sort_high_prime(first, last, std::less<void>());
        }


        // Gnome Weave Sort (Low Prime)
        template<class RandomAccessIterator, class Compare>
        void gnome_weave_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_weave_sorting::gnome_weave_sort;
            gnome_weave_sort(first, last, false, comp);
        }

        template<class RandomAccessIterator>
        void gnome_weave_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last) {
            gnome_weave_sort_low_prime(first, last, std::less<void>());
        }


        // Head Pull Sort
        template<class RandomAccessIterator, class Compare>
        void head_pull_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto i = first + 1, pull = first + 1;
            while (i + 1 <= last) {
                if (comp(*i, *(i - 1))) {
                    pull = i;
                    while (pull > first) {
                        std::iter_swap(pull, pull - 1);
                        pull--;
                    }
                    i = first + 1;
                }
                else i++;
            }
        }

        template<class RandomAccessIterator>
        void head_pull_sort(RandomAccessIterator first, RandomAccessIterator last) {
            head_pull_sort(first, last, std::less<void>());
        }


        // TODO: Lazy Stable Quick Sort


        // Iterative Quick Sort
        template<class RandomAccessIterator, class Compare>
        void iterative_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::iterative_quick_sorting::iterative_quick_sort;
            iterative_quick_sort(first, last - 1, comp);
        }

        template<class RandomAccessIterator>
        void iterative_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            iterative_quick_sort(first, last, std::less<void>());
        }


        // Linked Iterative Quick Sort
        template<class RandomAccessIterator, class Compare>
        void linked_iterative_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::iterative_quick_sorting::linked_iterative_quick_sort;
            linked_iterative_quick_sort(first, last - 1, comp);
        }

        template<class RandomAccessIterator>
        void linked_iterative_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            linked_iterative_quick_sort(first, last, std::less<void>());
        }


        // Quick Sort (LL ptrs with Middle Pivot)
        template<class RandomAccessIterator, class Compare>
        void quick_sort_ll_middle(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::quick_sorting::median_quick_ll_loop;
            median_quick_ll_loop(first, last, comp);
        }

        template<class RandomAccessIterator>
        void quick_sort_ll_middle(RandomAccessIterator first, RandomAccessIterator last) {
            quick_sort_ll_middle(first, last, std::less<void>());
        }


        // More Optimized Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void more_optimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            int c = 1;
            RandomAccessIterator s, f = first + 1;
            bool a = false;

            for (auto j = last - 1; j > first; j -= c) {
                if (f - 1 < first) s = first;
                else s = f - 1;
                
                a = false;
                c = 1;

                for (auto i = s; i < j; i++) {
                    if (comp(*(i + 1), *i)) {
                        std::iter_swap(i, i + 1);
                        if (!a) f = i;
                        a = true;
                        c = 1;
                    }
                    else c++;
                }
            }
        }

        template<class RandomAccessIterator>
        void more_optimized_bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            more_optimized_bubble_sort(first, last, std::less<void>());
        }


        // Odd Even Weave Sort (High Prime)
        template<class RandomAccessIterator, class Compare>
        void odd_even_weave_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::odd_even_weave_sorting::odd_even_weave_sort;
            odd_even_weave_sort(first, last, true, comp);
        }

        template<class RandomAccessIterator>
        void odd_even_weave_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_weave_sort_high_prime(first, last, std::less<void>());
        }


        // Odd Even Weave Sort (Low Prime)
        template<class RandomAccessIterator, class Compare>
        void odd_even_weave_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::odd_even_weave_sorting::odd_even_weave_sort;
            odd_even_weave_sort(first, last, false, comp);
        }

        template<class RandomAccessIterator>
        void odd_even_weave_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_weave_sort_low_prime(first, last, std::less<void>());
        }


        // Optimized Zipper Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_zipper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::zipper_sorting::optimized_zipper_sort;
            optimized_zipper_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_zipper_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_zipper_sort(first, last, std::less<void>());
        }


        // Pattern-Defeating Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void pattern_defeating_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_sorting::pd_gnome_sort;
            pd_gnome_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pattern_defeating_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pattern_defeating_gnome_sort(first, last, std::less<void>());
        }


        // Push Sort
        template<class RandomAccessIterator, class Compare>
        void push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::push_sorting::push_sort;
            push_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void push_sort(RandomAccessIterator first, RandomAccessIterator last) {
            push_sort(first, last, std::less<void>());
        }


        // Reverse Bubble Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_bubble_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i < last - 1; i++) {
                bool sorted = true;
                for (auto j = last - 1; j > i; j--) 
                    if (comp(*j, *(j - 1))) {
                        std::iter_swap(j - 1, j);
                        sorted = false;
                    }
                if (sorted) break;
            }
        }

        template<class RandomAccessIterator>
        void reverse_bubble_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_bubble_sort(first, last, std::less<void>());
        }


        // Reverse Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_sorting::reverse_gnome_sort;
            reverse_gnome_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_gnome_sort(first, last, std::less<void>());
        }


        // Reverse Push Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::push_sorting::reverse_push_sort;
            reverse_push_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_push_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_push_sort(first, last, std::less<void>());
        }


        // Split Center Sort
        template<class RandomAccessIterator, class Compare>
        void split_center_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            int length = last - first;
            int way = 1;
            auto i = first + 1;

            for (auto r = first + 1; r < last; r++) {
                i = first + length / 2;
                while (i < last && i > first) {
                    if (comp(*i, *(i - 1))) std::iter_swap(i, i - 1);
                    i += way;
                }
                way *= -1;
            }
        }

        template<class RandomAccessIterator>
        void split_center_sort(RandomAccessIterator first, RandomAccessIterator last) {
            split_center_sort(first, last, std::less<void>());
        }


        // Stable Quick Sort (Middle Pivot)
        template<class RandomAccessIterator, class Compare>
        void stable_quick_sort_middle(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            using sortcpp::_internal::stable_quick_sorting::stable_quick_loop;
            stable_quick_loop<T>(first, last - 1, true, comp);
        }

        template<class RandomAccessIterator>
        void stable_quick_sort_middle(RandomAccessIterator first, RandomAccessIterator last) {
            stable_quick_sort_middle(first, last, std::less<void>());
        }


        // Stackless Quick Sort
        template<class RandomAccessIterator, class Compare>
        void stackless_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::iterative_quick_sorting::stackless_quick_sort;
            stackless_quick_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void stackless_quick_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stackless_quick_sort(first, last, std::less<void>());
        }


        // Strange Push Sort
        template<class RandomAccessIterator, class Compare>
        void strange_push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, int base = 2) {
            using sortcpp::_internal::push_sorting::strange_push_sort;
            strange_push_sort(first, last, base, comp);
        }

        template<class RandomAccessIterator>
        void strange_push_sort(RandomAccessIterator first, RandomAccessIterator last, int base = 2) {
            strange_push_sort(first, last, std::less<void>(), base);
        }


        // Swap Map Sort
        template<class RandomAccessIterator, class Compare>
        void swap_map_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            std::vector<RandomAccessIterator> map;
            while (true) {
                for (auto i = first; i < last - 1; i++)
                    if (comp(*(i + 1), *i)) map.push_back(i);
                
                if (map.empty()) break;

                for (int i = 0; i < map.size(); i++)
                    std::iter_swap(map[i], map[i] + 1);

                map.clear();
            }
        }

        template<class RandomAccessIterator>
        void swap_map_sort(RandomAccessIterator first, RandomAccessIterator last) {
            swap_map_sort(first, last, std::less<void>());
        }


        // Tri Search Gnome Sort
        template<class RandomAccessIterator, class Compare>
        void tri_search_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::gnome_sorting::tri_gnome_sort;
            tri_gnome_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void tri_search_gnome_sort(RandomAccessIterator first, RandomAccessIterator last) {
            tri_search_gnome_sort(first, last, std::less<void>());
        }


        // Unbelievable Sort
        template<class RandomAccessIterator, class Compare>
        void unbelievable_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i < last; i++)
                for (auto j = first; j < last; j++)
                    if (comp(*i, *j)) std::iter_swap(i, j);
        }

        template<class RandomAccessIterator>
        void unbelievable_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unbelievable_sort(first, last, std::less<void>());
        }


        // Wiggle Sort
        template<class RandomAccessIterator, class Compare>
        void wiggle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::wiggle_sorting::wiggle_loop;
            wiggle_loop(first, last, comp);
        }

        template<class RandomAccessIterator>
        void wiggle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            wiggle_sort(first, last, std::less<void>());
        }


        // Zipper Sort
        template<class RandomAccessIterator, class Compare>
        void zipper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::zipper_sorting::zipper_sort;
            zipper_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void zipper_sort(RandomAccessIterator first, RandomAccessIterator last) {
            zipper_sort(first, last, std::less<void>());
        }
    }
}