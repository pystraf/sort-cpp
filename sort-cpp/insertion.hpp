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


        // Extra Sorts


        // Adaptive Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void adaptive_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::adaptive_insertion_sorting::adaptive_insertion_sort;
            adaptive_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void adaptive_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            adaptive_insertion_sort(first, last, std::less<void>());
        }


        // Adaptive Binary Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void adaptive_binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::adaptive_insertion_sorting::adaptive_binary_insertion_sort;
            adaptive_binary_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void adaptive_binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            adaptive_binary_insertion_sort(first, last, std::less<void>());
        }


        // Bidirectional Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void bidirectional_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::bidirectional_insertion_sort;
            bidirectional_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void bidirectional_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bidirectional_insertion_sort(first, last, std::less<void>());
        }


        // Cocktail Shell Sort
        template<class RandomAccessIterator, class Compare>
        void cocktail_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::shell_sorting::cocktail_shell_sort;
            cocktail_shell_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void cocktail_shell_sort(RandomAccessIterator first, RandomAccessIterator last) {
            cocktail_shell_sort(first, last, std::less<void>());
        }


        // Exponential Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void exponential_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::exp_insertion_sort;
            exp_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void exponential_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            exponential_insertion_sort(first, last, std::less<void>());
        }


        // Fibonacci Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void fibonacci_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::fib_insertion_sort;
            fib_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void fibonacci_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fibonacci_insertion_sort(first, last, std::less<void>());
        }


        // Gambit Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void gambit_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::gambit_insertion_sort;
            gambit_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void gambit_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            gambit_insertion_sort(first, last, std::less<void>());
        }


        // Pattern-Defeating Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void pd_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::adaptive_insertion_sorting::pd_insertion_sort;
            pd_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pd_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pd_insertion_sort(first, last, std::less<void>());
        }


        // Par Shell Sort
        template<class RandomAccessIterator, class Compare>
        void par_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, int constant = 2) {
            using sortcpp::_internal::par_shell_sorting::par_shell_sort;
            par_shell_sort(first, last, constant, comp);
        }

        template<class RandomAccessIterator>
        void par_shell_sort(RandomAccessIterator first, RandomAccessIterator last, int constant = 2) {
            par_shell_sort(first, last, std::less<void>(), constant);
        }


        // Pattern-Defeating Binary Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void pd_binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::adaptive_insertion_sorting::pd_binary_insertion_sort;
            pd_binary_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pd_binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pd_binary_insertion_sort(first, last, std::less<void>());
        }


        // Pattern-Defeating Exponential Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void pd_exponential_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::adaptive_insertion_sorting::pd_exp_insertion_sort;
            pd_exp_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void pd_exponential_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pd_exponential_insertion_sort(first, last, std::less<void>());
        }


        // Rendezvous Sort
        template<class RandomAccessIterator, class Compare>
        void rendezvous_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::rendezvous_sort;
            rendezvous_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void rendezvous_sort(RandomAccessIterator first, RandomAccessIterator last) {
            rendezvous_sort(first, last, std::less<void>());
        }


        // Reverse Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void reverse_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::reverse_insertion_sort;
            reverse_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void reverse_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            reverse_insertion_sort(first, last, std::less<void>());
        }


        // Room Sort
        template<class RandomAccessIterator, class Compare>
        void room_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::room_sort;
            room_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void room_sort(RandomAccessIterator first, RandomAccessIterator last) {
            room_sort(first, last, std::less<void>());
        }


        // Shell Sort (High Prime)
        template<class RandomAccessIterator, class Compare>
        void shell_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::shell_sorting::shell_sort_primes;
            shell_sort_primes(first, last, true, comp);
        }

        template<class RandomAccessIterator>
        void shell_sort_high_prime(RandomAccessIterator first, RandomAccessIterator last) {
            shell_sort_high_prime(first, last, std::less<void>());
        }


        // Shell Sort (Low Prime)
        template<class RandomAccessIterator, class Compare>
        void shell_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::shell_sorting::shell_sort_primes;
            shell_sort_primes(first, last, false, comp);
        }

        template<class RandomAccessIterator>
        void shell_sort_low_prime(RandomAccessIterator first, RandomAccessIterator last) {
            shell_sort_low_prime(first, last, std::less<void>());
        }


        // Tri Search Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void tri_search_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::binary_insertion_sorting::tri_insertion_sort;
            tri_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void tri_search_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            tri_search_insertion_sort(first, last, std::less<void>());
        }


        // Unstable Insertion Sort
        template<class RandomAccessIterator, class Compare>
        void unstable_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::insertion_sorting::unstable_insertion_sort;
            unstable_insertion_sort(first, last, comp);
        }

        template<class RandomAccessIterator>
        void unstable_insertion_sort(RandomAccessIterator first, RandomAccessIterator last) {
            unstable_insertion_sort(first, last, std::less<void>());
        }
    }
}