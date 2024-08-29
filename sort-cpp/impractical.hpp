#pragma once

#include <random>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <stack>
#include "templates.hpp"

namespace sortcpp {
    namespace impractical {
        // Bubble Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void bubble_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_range_sorted;
            using sortcpp::_internal::bogo_sorting::randiter;
            
            while (!is_range_sorted(first, last, comp)) {
                auto i = randiter(first, last - 1, rng);
                if (comp(*(i + 1), *i)) std::iter_swap(i, i + 1);
            }
        }
        
        template<class RandomAccessIterator, class RandomNumberGenerator>
        void bubble_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng) {
            bubble_bogo_sort(first, last, rng, std::less<void>());
        }
    
        
        // Exchange Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void exchange_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_range_sorted;
            using sortcpp::_internal::bogo_sorting::randiter;
            
            while (!is_range_sorted(first, last, comp)) {
                auto i = randiter(first, last, rng);
                auto j = randiter(first, last, rng);
                
                if (i < j) {
                    if (comp(*j, *i)) std::iter_swap(i, j);
                } else {
                    if (comp(*i, *j)) std::iter_swap(i, j);
                }
            }
        }
        
        template<class RandomAccessIterator, class RandomNumberGenerator>
        void exchange_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng) {
            exchange_bogo_sort(first, last, rng, std::less<void>());
        }
    
        
        // Quad Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void quad_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::quad_stooge_sort;
            quad_stooge_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void quad_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            quad_stooge_sort(first, last, std::less<void>());
        }
    
        
        // Shove Sort
        template<class RandomAccessIterator, class Compare>
        void shove_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto i = first;
            while (i < last - 1) {
                if (comp(*(i + 1), *i)) {
                    for (auto f = i; f < last - 1; f++) std::iter_swap(f, f + 1);
                    if (i > first) i--;
                    continue;
                }
                i++;
            }
        }
    
        template<class RandomAccessIterator>
        void shove_sort(RandomAccessIterator first, RandomAccessIterator last) {
            shove_sort(first, last, std::less<void>());
        }
    
        
        // Silly Sort
        template<class RandomAccessIterator, class Compare>
        void silly_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::slow_sorting::silly_sort;
            silly_sort(first, last, comp);
        }
    
        template<class RandomAccessIterator>
        void silly_sort(RandomAccessIterator first, RandomAccessIterator last) {
            silly_sort(first, last, std::less<void>());
        }
    
    
        // Slow Sort
        template<class RandomAccessIterator, class Compare>
        void slow_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::slow_sorting::slow_sort;
            slow_sort(first, last, comp);
        }
    
        template<class RandomAccessIterator>
        void slow_sort(RandomAccessIterator first, RandomAccessIterator last) {
            slow_sort(first, last, std::less<void>());
        }
    
        
        // Snuffle Sort
        template<class RandomAccessIterator, class Compare>
        void snuffle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::slow_sorting::snuffle_sort;
            snuffle_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void snuffle_sort(RandomAccessIterator first, RandomAccessIterator last) {
            snuffle_sort(first, last, std::less<void>());
        }
        
        
        // Stable Permutation Sort
        template<class RandomAccessIterator, class Compare>
        void stable_permutation_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::permutation_sorting::stable_permutation_sort;
            stable_permutation_sort(first, last, comp);
        }
    
        template<class RandomAccessIterator>
        void stable_permutation_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_permutation_sort(first, last, std::less<void>());
        }
    
    
        // Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::stooge_sort;
            stooge_sort(first, last, comp);
        }
        
        template<class RandomAccessIterator>
        void stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stooge_sort(first, last, std::less<void>());
        }		
        
        
        // Bad Sort
        template<class RandomAccessIterator, class Compare>
        void bad_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            for (auto i = first; i < last; i++) {
                auto shortest = i;
                
                for (auto j = i; j < last; j++) {
                    bool is_shortest = true;
                    for (auto k = j + 1; k < last; k++) 
                        if (comp(*k, *j)) {
                            is_shortest = false;
                            break;
                        }
                    
                    if (is_shortest) {
                        shortest = j;
                        break;
                    }
                }
                
                std::iter_swap(i, shortest);
            }
        }
        
        template<class RandomAccessIterator>
        void bad_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bad_sort(first, last, std::less<void>());
        }


        // Hanoi Sort
        template<class RandomAccessIterator, class Compare>
        void hanoi_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::hanoi_sorting::hanoi_sort;
            hanoi_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void hanoi_sort(RandomAccessIterator first, RandomAccessIterator last) {
            hanoi_sort(first, last, std::less<void>());
        }

        // Bogo Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void bogo_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::bogobogo_sort_loop;
            bogobogo_sort_loop(first, last, rng, comp);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void bogo_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            bogo_bogo_sort(first, last, rng, std::less<void>());
        }

        // Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_range_sorted;
            using sortcpp::_internal::bogo_sorting::bogo_swap;

            while (!is_range_sorted(first, last, comp)) bogo_swap(first, last, rng);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            bogo_sort(first, last, rng, std::less<void>());
        }


        // Bozo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void bozo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_range_sorted;
            using sortcpp::_internal::bogo_sorting::randiter;

            while (!is_range_sorted(first, last, comp)) {
                auto i = randiter(first, last, rng);
                auto j = randiter(first, last, rng);
                std::iter_swap(i, j);
            }
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void bozo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            bozo_sort(first, last, rng, std::less<void>());
        }


        // Cocktail Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void cocktail_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_min_sorted;
            using sortcpp::_internal::bogo_sorting::is_max_sorted;
            using sortcpp::_internal::bogo_sorting::bogo_swap;
            
            auto min = first, max = last;
            while (min < max - 1) {
                if (is_min_sorted(min, max, comp)) min++;
                if (is_max_sorted(min, max, comp)) max--;
                bogo_swap(min, max, rng);
            }
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void cocktail_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            cocktail_bogo_sort(first, last, rng, std::less<void>());
        }


        // Deterministic Bogo Sort (Permutation Sort)
        template<class RandomAccessIterator, class Compare>
        void deterministic_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::permutation_sorting::permutation_sort;
            permutation_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void deterministic_bogo_sort(RandomAccessIterator first, RandomAccessIterator last) {
            deterministic_bogo_sort(first, last, std::less<void>());
        }


        // Guess Sort
        template<class RandomAccessIterator, class Compare>
        void guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::guess_sorting::guess_sort;
            guess_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void guess_sort(RandomAccessIterator first, RandomAccessIterator last) {
            guess_sort(first, last, std::less<void>());
        }


        // Less Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void less_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::is_min_sorted;
            using sortcpp::_internal::bogo_sorting::bogo_swap;

            for (auto i = first; i < last; i++)
                while (!is_min_sorted(i, last, comp)) bogo_swap(i, last, rng);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void less_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            less_bogo_sort(first, last, rng, std::less<void>());
        }


        // Median Quick Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void median_quick_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::median_quick_bogo_loop;
            median_quick_bogo_loop(first, last, rng, comp);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void median_quick_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            median_quick_bogo_sort(first, last, rng, std::less<void>());
        }


        // Merge Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void merge_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::merge_bogo_loop;
            merge_bogo_loop(first, last, rng, comp);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void merge_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            merge_bogo_sort(first, last, rng, std::less<void>());
        }


        // Optimized Guess Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::guess_sorting::optimized_guess_sort;
            optimized_guess_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_guess_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_guess_sort(first, last, std::less<void>());
        }


        // Quick Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void quick_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::quick_bogo_loop;
            quick_bogo_loop(first, last, rng, comp);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void quick_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            quick_bogo_sort(first, last, rng, std::less<void>());
        }

        
        // Random Guess Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void random_guess_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::guess_sorting::random_guess_sort;
            random_guess_sort(first, last, rng, comp);
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void random_guess_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            random_guess_sort(first, last, rng, std::less<void>());
        }

        
        // Selection Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void selection_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::randiter;
            using sortcpp::_internal::bogo_sorting::is_min_sorted;

            for (auto i = first; i < last; i++)
                while (!is_min_sorted(i, last, comp)) std::iter_swap(i, randiter(i, last, rng));
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void selection_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            selection_bogo_sort(first, last, rng, std::less<void>());
        }



        // Smart Guess Sort
        template<class RandomAccessIterator, class Compare>
        void smart_guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::guess_sorting::smart_guess_sort;
            smart_guess_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void smart_guess_sort(RandomAccessIterator first, RandomAccessIterator last) {
            smart_guess_sort(first, last, std::less<void>());
        }


        // Extra Sorts

        // Awkward Sort
        template<class RandomAccessIterator, class Compare>
        void awkward_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::awkward_sort;
            awkward_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void awkward_sort(RandomAccessIterator first, RandomAccessIterator last) {
            awkward_sort(first, last, std::less<void>());
        }


        // Fire Sort
        template<class RandomAccessIterator, class Compare>
        void fire_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fire_sorting::fire_sort;
            fire_sort(first, last, last - first, comp);
        }

        template<class RandomAccessIterator>
        void fire_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fire_sort(first, last, std::less<void>());
        }


        // In-Order Shove Sort
        template<class RandomAccessIterator, class Compare>
        void in_order_shove_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto left = first, right = first + 1, pull = first;
            while (left < last) {
                right = left + 1;
                while (right < last) {
                    if (comp(*right, *left)) {
                        pull = left;
                        while (pull + 1 < last) {
                            std::iter_swap(pull, pull + 1);
                            pull++;
                        }
                        right = left + 1;
                    }
                    else right++;
                }
                left++;
            }
        }

        template<class RandomAccessIterator>
        void in_order_shove_sort(RandomAccessIterator first, RandomAccessIterator last) {
            in_order_shove_sort(first, last, std::less<void>());
        }


        // Markov Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void markov_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng, Compare comp) {
            using sortcpp::_internal::gnome_sorting::markov_sort;
            markov_sort(first, last, rng, comp);
        }
        
        template<class RandomAccessIterator, class RandomNumberGenerator>
        void markov_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            markov_sort(first, last, rng, std::less<void>());
        }


        // Merry-Go-Round Sort
        template<class RandomAccessIterator, class Compare>
        void merry_go_round_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fire_sorting::merry_go_round_sort;
            merry_go_round_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void merry_go_round_sort(RandomAccessIterator first, RandomAccessIterator last) {
            merry_go_round_sort(first, last, std::less<void>());
        }


        // Napoleon Sort
        template<class RandomAccessIterator, class Compare>
        void napoleon_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::napoleon_sorting::tilsit;
            using sortcpp::_internal::napoleon_sorting::napoleon;

            tilsit(first, last - 1, comp);
            napoleon(first, last - 1, comp);
        }

        template<class RandomAccessIterator>
        void napoleon_sort(RandomAccessIterator first, RandomAccessIterator last) {
            napoleon_sort(first, last, std::less<void>());
        }


        // Noisy Sort
        template<class RandomAccessIterator, class Compare>
        void noisy_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, int base = 16) {
            auto left = first + 1, right = first + 1, verifyi = first + 1;
            bool verfiy_pass = false;

            while (!verfiy_pass) {
                right = verifyi + 1;
                while (right <= last) {
                    left = verifyi;
                    while (left <= right && right <= last) {
                        if (comp(*(right - 1), *(left - 1))) {
                            std::iter_swap(left - 1, right - 1);
                            if (right - 1 > verifyi) right--;
                            left = verifyi;
                        }
                        else left++;
                    }
                    right += base;
                }

                if (verifyi - 1 > first) verifyi--;
                verfiy_pass = true;
                while (verifyi < last && verfiy_pass) {
                    if (!comp(*verifyi, *(verifyi - 1))) verifyi++;
                    else verfiy_pass = false;
                }
            }
        }

        template<class RandomAccessIterator>
        void noisy_sort(RandomAccessIterator first, RandomAccessIterator last, int base = 16) {
            noisy_sort(first, last, std::less<void>(), base);
        }


        // Optimized Bubble Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void optimized_bubble_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::randiter;
            
            auto start = first, end = last - 1;
            while (!comp(*(start + 1), *start) && start <= end) start++;
            while (!comp(*end, *(end - 1)) && start <= end) end--;

            while (start <= end) {
                auto index = randiter(start, end, rng);
                if (comp(*(index + 1), *index)) {
                    std::iter_swap(index, index + 1);
                    if (index == start) {
                        if (start > first) start--;
                        while (!comp(*(start + 1), *start) && start <= end) start++;
                    }

                    if (index == end - 1) {
                        if (end < last - 1) end++;
                        while (!comp(*end, *(end - 1)) && start <= end) end--;
                    }
                }
            }
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void optimized_bubble_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng) {
            optimized_bubble_bogo_sort(first, last, rng, std::less<void>());
        }


        // Playground Sort
        template<class RandomAccessIterator, class Compare>
        void playground_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::playground_sorting::playground_sort;
            playground_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void playground_sort(RandomAccessIterator first, RandomAccessIterator last) {
            playground_sort(first, last, std::less<void>());
        }


        // Cocktail Grate Sort
        template<class RandomAccessIterator, class Compare>
        void cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::cocktail_grate_sort;
            cocktail_grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cocktail_grate_sort(first, last, std::less<void>());
        }


        // Grate Sort
        template<class RandomAccessIterator, class Compare>
        void grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::grate_sort;
            grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            grate_sort(first, last, std::less<void>());
        }

        // Optimized Cocktail Grate Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::optimized_cocktail_grate_sort;
            optimized_cocktail_grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_cocktail_grate_sort(first, last, std::less<void>());
        }


        // Optimized Grate Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::optimized_grate_sort;
            optimized_grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_grate_sort(first, last, std::less<void>());
        }


        // Optimized Reverse Grate Sort
        template<class RandomAccessIterator, class Compare>
        void optimized_reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::optimized_reverse_grate_sort;
            optimized_reverse_grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void optimized_reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            optimized_reverse_grate_sort(first, last, std::less<void>());
        }


        // Reverse Grate Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::grate_sorting::reverse_grate_sort;
            reverse_grate_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_grate_sort(first, last, std::less<void>());
        }


        // Safe Bogo Sort
        template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
        void safe_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
            using sortcpp::_internal::bogo_sorting::randiter;
            using sortcpp::_internal::bogo_sorting::find_last_sorted;

            auto p = find_last_sorted(first, last, comp);
            while (p < last - 1) {
                std::iter_swap(p, randiter(p, last, rng));
                p = find_last_sorted(first, last, comp);
            }
        }

        template<class RandomAccessIterator, class RandomNumberGenerator>
        void safe_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng) {
            safe_bogo_sort(first, last, rng, std::less<void>());
        }


        // Stable Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void stable_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::stable_stooge_sort;
            stable_stooge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void stable_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stable_stooge_sort(first, last, std::less<void>());
        }


        // Stupid Fire Sort
        template<class RandomAccessIterator, class Compare>
        void stupid_fire_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fire_sorting::fire_sort;
            fire_sort(first, last, 1, comp);
        }

        template<class RandomAccessIterator>
        void stupid_fire_sort(RandomAccessIterator first, RandomAccessIterator last) {
            stupid_fire_sort(first, last, std::less<void>());
        }


        // Ternary Slow Sort
        template<class RandomAccessIterator, class Compare>
        void ternary_slow_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::slow_sorting::ternary_slow_sort;
            ternary_slow_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void ternary_slow_sort(RandomAccessIterator first, RandomAccessIterator last) {
            ternary_slow_sort(first, last, std::less<void>());
        }


        // X-Pattern Sort
        template<class RandomAccessIterator, class Compare>
        void x_pattern_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::fire_sorting::x_pattern_sort;
            x_pattern_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void x_pattern_sort(RandomAccessIterator first, RandomAccessIterator last) {
            x_pattern_sort(first, last, std::less<void>());
        }


        // Hyper Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void hyper_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::stooge_sorting::hyper_stooge_sort;
            hyper_stooge_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void hyper_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            hyper_stooge_sort(first, last, std::less<void>());
        }

    }
}