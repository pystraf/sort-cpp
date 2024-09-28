#pragma once

#include <random>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <queue>
#include <deque>
#include <stack>

namespace sortcpp {
    namespace _internal {
        // Determines if a equals b in the special rule of the sort.
        // a - The first element
        // b - The second element
        // comp - The comparison function
        template<class T, class Compare>
        bool equals(const T& a, const T& b, Compare comp) {
            return (!comp(a, b)) && (!comp(b, a));            // a == b
        }

        // Reverses the range [first, last)
        // first - The beginning of the range
        // last - The end of the range
        template<class RandomAccessIterator>
        void reversal(RandomAccessIterator first, RandomAccessIterator last) {
            if (first == last) return;                          // Empty range.
            last--;                                             // Convert to [first, last].
            while (first < last) {
                std::iter_swap(first, last);                    // Swap the elements.
                first++;
                last--;
            }
        }

        // Return the length of the binary representation of `x` within the sign bit.
        // We consider 0 as a zero digit number.
        // x - The number.
        template<class T>
        int bit_length(const T& x) {
            return x ? (int)(std::log2(x)) + 1 : 0;                 // = ceil(log2(x))
        }

        // Return the number of trailing zeros in the binary representation of `x`.
        // We consider 0 as a zero digit number. So ctz(0) = sizeof(T) * 8.
        // x - The number.
        template<class T>
        int ctz(const T& x) {
            if (x == 0) return sizeof(T) * 8;                       // For 0.
            return bit_length(~x & (x - 1));                        // = The length of the binary representation of `~x & (x - 1)`.
        }

        // Returns the smallest number of powers of 2 that are greater than or equal to `x`.
        // x - The number.
        template<class T>
        T ceil_pow2(const T& x) {
            return 1 << (int)(std::ceil(std::log2(x)));              // = 2 ^ ceil(log2(x))
        }



        // Get the kth bit from the right of `x`.
        // x - The number.
        // k - The bit position.
        template<class T>
        bool get_bit(const T& x, int k) {
            return (x >> k) & 1;
        }

        // Get the kth digit from the right of `x` in the given radix.
        // x - The number.
        // power - The digit position.
        // radix - The base of the number system.
        template<class T>
        T get_digit(const T& a, int power, int radix) {
            return (T)(a / std::pow(radix, power)) % radix;
        }

        // Find `value` in the range [first, last) using binary search.
        // If `value` is not found, return the position where it should be inserted.
        // first - The beginning of the range.
        // last - The end of the range.
        // value - The value to find.
        // comp - The comparison function.
        template<class RandomAccessIterator, class T, class Compare>
        RandomAccessIterator binary_find(RandomAccessIterator first, RandomAccessIterator last, const T& value, Compare comp) {
            if (first == last) return first;                       // Empty range.
            auto lo = first, hi = last;
            while (lo < hi) {
                auto mid = lo + (hi - lo) / 2;
                if (comp(value, *mid)) hi = mid;                   // If the value is less than the middle element, narrow the search to the lower half.
                else lo = mid + 1;                                 // Otherwise, narrow the search to the upper half.
            }
            // If the value is not found, `lo` is the position where it should be inserted.
            return lo;
        }


        template<class RandomAccessIterator, class Register>
        void transcribe(RandomAccessIterator start, Register& registers) {
            auto total = start;
            for (int index = 0; index < registers.size(); index++) {
                for (int i = 0; i < registers[index].size(); i++) {
                    *total = registers[index][i];
                    total++;
                }
                registers[index].clear();
            }
        }

        template<class RandomAccessIterator, class Register>
        void transcribe_msd(RandomAccessIterator start, Register& registers) {
            auto total = start;
            int temp = 0;

            for (auto reg: registers) total += reg.size();

            for (int index = registers.size() - 1; index >= 0; index--) {
                for (int i = registers[index].size() - 1; i >= 0; i--)
                    *(total - (temp++) - 1) = registers[index][i];
            }
        }

        template<class T, class RandomAccessIterator, class Register>
        void fancy_transcribe(RandomAccessIterator start, RandomAccessIterator end, Register& registers) {
            int length = end - start;
            std::vector<T> tmp(length);
            std::vector<bool> tmp_write(length, false);

            int radix = registers.size();
            transcribe(tmp.begin(), registers);

            for (int i = 0; i < length; i++) {
                int reg = i % radix;
                int pos = (reg * (length / radix)) + (i / radix);

                if (!tmp_write[pos]) {
                    *(start + pos) = tmp[pos];
                    tmp_write[pos] = true;
                }
            }

            for (int i = 0; i < length; i++) 
                if (!tmp_write[i]) *(start + i) = tmp[i];
        }

        template<class RandomAccessIterator, class T, class Compare>
        RandomAccessIterator _fib_search(RandomAccessIterator start, RandomAccessIterator end, const T& item, Compare comp) {
            int fib_m2 = 0, fib_m1 = 0, fib_m = 1;
            while (fib_m <= end - start) {
                fib_m2 = fib_m1;
                fib_m1 = fib_m;
                fib_m = fib_m1 + fib_m2;
            }

            auto offset = start - 1;
            while (fib_m > 1) {
                auto i = std::min(offset + fib_m2, end);

                if (comp(item, *i)) {
                    fib_m = fib_m2;
                    fib_m1 -= fib_m2;
                    fib_m2 = fib_m - fib_m1;
                }
                else {
                    fib_m = fib_m1;
                    fib_m1 = fib_m2;
                    fib_m2 = fib_m - fib_m1;
                    offset = i;
                }
            }

            auto pos = ++offset;
            if (!comp(item, *pos)) pos++;
            return pos;
        }

        template<class RandomAccessIterator, class T, class Compare>
        RandomAccessIterator _tri_search(RandomAccessIterator l, RandomAccessIterator h, const T& val, Compare comp) {
            auto mid = l + (h - l) / 2;
            if (comp(val, *l)) return l;
            else {
                if (comp(val, *h)) {
                    if (comp(val, *mid)) return _tri_search(l + 1, mid - 1, val, comp);
                    else return _tri_search(mid + 1, h - 1, val, comp);
                }
                else return h + 1;
            }
        }

        namespace iterative_circle_sorting {
            // Sort the range [first, last) once use the circle sort algorithm with iterative.
            // first - The beginning of the range.
            // last - The real end of the range.
            // size - Size of current iteration.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool circle_sort_routine(RandomAccessIterator first, RandomAccessIterator last, int size, Compare comp) {
                bool swapped = false;                              // If the range is swapped.
                for (int gap = size / 2; gap > 0; gap /= 2)        // Let gap = size / 2, size / 4, size / 8, ..., 1.
                    for (auto start = first; start + gap < last; start += 2 * gap) {
                        // Sort the range [start, start + 2 * gap - 1] once.
                        auto high = start + 2 * gap - 1;
                        auto low = start;

                        while (low < high) {
                            if (high < last && comp(*high, *low)) {// Swap the elements if the high element is less than the low element.
                                std::iter_swap(low, high);
                                swapped = true;
                            }

                            low++;
                            high--;
                        }
                    }
                return swapped;
            }
        }

        namespace circle_sorting {
            // Sort the range [lo, hi] once use the circle sort algorithm with recursion.
            // lo - The beginning of the range.
            // hi - The end of the range currently (inclusive).
            // last - The real end of the range.
            // swapped - If the range is swapped.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool circle_sort_routine(RandomAccessIterator lo, RandomAccessIterator hi, RandomAccessIterator last, bool swapped, Compare comp) {
                if (lo == hi) return swapped;                      // Empty range.

                auto high = hi, low = lo;
                auto mid = lo + (hi - lo) / 2;

                while (lo < hi) {
                    if (hi < last && comp(*hi, *lo)) {             // Swap the elements if the high element is less than the low element.
                        std::iter_swap(lo, hi);
                        swapped = true;
                    }

                    lo++;
                    hi--;
                }

                // Sort the range [lo, mid] and [mid + 1, hi] recursively.
                swapped = circle_sort_routine(low, mid, last, swapped, comp);
                if (mid + 1 < last) swapped = circle_sort_routine(mid + 1, high, last, swapped, comp);
                return swapped;
            }
        }

        namespace circloid_sorting {
            template<class RandomAccessIterator, class Compare>
            bool circle(RandomAccessIterator left, RandomAccessIterator right, Compare comp) {
                auto a = left, b = right;
                bool swapped = false;

                while (a < b) {
                    if (comp(*b, *a)) {
                        std::iter_swap(a, b);
                        swapped = true;
                    }

                    a++;
                    b--;
                    if (a == b) b++;
                }
                return swapped;
            }

            template<class RandomAccessIterator, class Compare>
            bool circle_pass(RandomAccessIterator left, RandomAccessIterator right, Compare comp) {
                if (left >= right) return false;
                auto mid = left + (right - left) / 2;
                bool l = circle_pass(left, mid, comp);
                bool r = circle_pass(mid + 1, right, comp);
                return circle(left, right, comp) || l || r;
            }
        }

        namespace insertion_sorting {
            using sortcpp::_internal::reversal;

            // Sort the range [first, last) use the insertion sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last; i++) {              // Sort the range [first, i).
                    auto tmp = *i;                                 // Store the current element.
                    auto pos = i - 1;
                    while (pos >= first && comp(tmp, *pos)) {      // Move the elements that are greater than the current element.
                        *(pos + 1) = *pos;
                        pos--;
                    }

                    *(pos + 1) = tmp;                              // Insert the current element to the correct position.
                }
            }


            template<class RandomAccessIterator, class Compare>
            void double_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto left = first + (last - first) / 2 - 1, right = left + 1;
                if (comp(*right, *left)) std::iter_swap(left, right);

                left--;
                right++;

                while (left >= first && right < last) {
                    if (comp(*right, *left)) {
                        auto left_item = *right, right_item = *left;

                        auto pos = left + 1;
                        while (pos <= right && !comp(left_item, *pos)) {
                            *(pos - 1) = *pos;
                            pos++;
                        }
                        *(pos - 1) = left_item;

                        pos = right - 1;
                        while (pos >= left && !comp(*pos, right_item)) {
                            *(pos + 1) = *pos;
                            pos--;
                        }
                        *(pos + 1) = right_item;
                    }
                    else {
                        auto left_item = *left, right_item = *right;

                        auto pos = left + 1;
                        while (comp(*pos, left_item)) {
                            *(pos - 1) = *pos;
                            pos++;
                        }
                        *(pos - 1) = left_item;

                        pos = right - 1;
                        while (comp(right_item, *pos)) {
                            *(pos + 1) = *pos;
                            pos--;
                        }
                        *(pos + 1) = right_item;
                    }

                    left--;
                    right++;
                }

                if (right < last) {
                    auto pos = right - 1;
                    auto current = *right;
                    while (comp(current, *pos)) {
                        *(pos + 1) = *pos;
                        pos--;
                    }
                    *(pos + 1) = current;
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            void insert_fw(RandomAccessIterator i, const T& current, Compare comp) {
                auto pos = i - 1;
                while (comp(current, *pos)) {
                    *(pos + 1) = *pos;
                    pos--;
                }
                *(pos + 1) = current;
            }

            template<class RandomAccessIterator, class T, class Compare>
            void insert_bw(RandomAccessIterator i, const T& current, Compare comp) {
                auto pos = i - 1;
                while (!comp(current, *pos)) {
                    *(pos + 1) = *pos;
                    pos--;
                }
                *(pos + 1) = current;
            }

            template<class RandomAccessIterator, class Compare>
            void bidirectional_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                bool dir = true;
                for (auto i = a + 1; i < b; i++) {
                    auto current = *i;
                    if (dir) {
                        if (comp(current, *a)) {
                            reversal(a, i);
                            dir = !dir;
                        }
                        else {
                            insert_fw(i, current, comp);
                        }
                    }
                    else {
                        if (!comp(current, *a)) {
                            reversal(a, i);
                            dir = !dir;
                        }
                        else {
                            insert_bw(i, current, comp);
                        }
                    }
                }
                if (!dir) {
                    reversal(a, b);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _rendezvous_sort(RandomAccessIterator first, int length, Compare comp) {
                int len = length;
                RandomAccessIterator i, j, ticker;
                int f, gap;
                while (len >= 1) {
                    i = first;
                    j = first + len;

                    while (length >= j - first) {
                        ticker = i;
                        f = 0;

                        while (!(i < first + 1 || !comp(*j, *i))) {
                            std::iter_swap(i, j);
                            f = 1;
                            gap = 1 + (j - i);
                            i -= gap;
                            j -= gap;
                        }

                        i++;
                        j++;

                        if (f == 1) {
                            len /= 0.25;
                            i = first;
                            j = first + len;
                        }
                    }

                    len /= 8;
                }
                i = first;
                j = first + 1;

                while (length >= j - first) {
                    ticker = i;

                    while (!(i < first || !comp(*j, *i))) {
                        std::iter_swap(i, j);
                        i--;
                        j--;
                    }

                    i = ticker + 1;
                    j = ticker + 2;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void rendezvous_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _rendezvous_sort(first, last - first - 1, comp);
            }


            template<class RandomAccessIterator, class Compare>
            void reverse_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = last - 1; i >= first; i--) {
                    auto current = *i;
                    auto pos = i + 1;

                    while (pos <= last - 1 && comp(*pos, current)) {
                        *(pos - 1) = *pos;
                        pos++;
                    }
                    *(pos - 1) = current;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void insert_to(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto val = *a;
                a--;
                while (a >= b && comp(val, *a)) {
                    *(a + 1) = *a;
                    a--;
                }
                *(a + 1) = val;
            }

            template<class RandomAccessIterator, class Compare>
            void room_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int room_len = (int)std::sqrt(length) + 1;

                RandomAccessIterator end, i;
                for (end = last; end > first + room_len; end -= room_len) {
                    insertion_sort(first, first + room_len, comp);
                    for (auto i = first + room_len; i < end; i++) insert_to(i, i - room_len, comp);
                }
                insertion_sort(first, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void unstable_insertion_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                for (auto i = start + 1; i < end; i++) {
                    if (comp(*i, *start)) std::iter_swap(i, start);
                    auto tmp = *i;
                    auto j = i - 1;
                    for (; comp(tmp, *j); j--) *(j + 1) = *j;
                    *(j + 1) = tmp;
                }
            }
        }

        

        namespace complete_graph_sorting {
            template<class RandomAccessIterator, class Compare>
            void comp_swap(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            template<class RandomAccessIterator, class Compare>
            void split(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, Compare comp) {
                if (b - a < 2) return;

                auto c = first;
                int len1 = (b - a) / 2;
                bool odd = (b - a) & 1;

                if (odd) {
                    if (m - a > b - m) c = a++;
                    else c = --b;
                }

                for (int s = 0; s < len1; s++) {
                    auto i = a;
                    for (int j = s; j < len1; j++) comp_swap(i++, m + j, comp);
                    for (int j = 0; j < s; j++) comp_swap(i++, m + j, comp);
                }

                if (odd) {
                    if (c < m) for (int j = 0; j < len1; j++) comp_swap(c, m + j, comp);
                    else for (int j = 0; j < len1; j++) comp_swap(a + j, c, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void complete_graph_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int n = last - first;
                int d = 2, end = 1 << ((int)std::log2(n - 1) + 1);

                while (d <= end) {
                    auto i = first;
                    int dec = 0;

                    while (i < last) {
                        auto j = i;
                        dec += n;

                        while (dec >= d) {
                            dec -= d;
                            j++;
                        }

                        auto k = j;
                        dec += n;

                        while (dec >= d) {
                            dec -= d;
                            k++;
                        }
                        split(first, i, j, k, comp);
                        i = k;
                    }

                    d *= 2;
                }
            }
        }

        namespace comb_sorting {
            // Sort the range [first, last) use the comb sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range.
            // comp - The comparison function.
            // shrink - The shrink factor.
            // hybird - If use the insertion sort algorithm when the gap is small.
            template<class RandomAccessIterator, class Compare>
            void comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp, double shrink, bool hybird) {
                using sortcpp::_internal::insertion_sorting::insertion_sort;

                bool swapped = false;                                     // If the range is swapped.
                int length = last - first;                                // The length of the range.
                int gap = length;                                         // The gap.
                int threshold = std::min(8, (int)(length * 0.03125));     // The threshold for the insertion sort algorithm = min(8, 1/32 * length).

                while ((gap > 1) || swapped) {
                    if (gap > 1) gap = (int)(gap / shrink);               // Decrease the gap.
                    swapped = false;                                      // There are currently no swaps.

                    for (auto i = first; (i + gap) < last; i++) {
                        if (hybird && gap <= threshold) {                 // If the gap is too small, use the insertion sort algorithm and exit.
                            gap = 0;                                      // Let the gap be 0 to exit the loop.
                            insertion_sort(first, last, comp);
                            break;
                        }

                        if (comp(*(i + gap), *i)) {                       // If the elements are not in order.
                            std::iter_swap(i, i + gap);                   // Swap the elements.
                            swapped = true;
                        }
                    }
                }
            }
        }

        namespace three_smooth_comb_sorting {
            // Determine whether `n` is a 3-smooth number.
            // n - The number to be determined.
            bool is_3smooth(int n) {
                while (n % 6 == 0) n /= 6;
                while (n % 3 == 0) n /= 3;
                while (n % 2 == 0) n /= 2;
                return n == 1;
            }


            // Sort the range [first, last) using the 3-smooth comb sort algorithm.
            template<class RandomAccessIterator, class Compare>
            void classic_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                for (int g = length - 1; g > 0; g--)
                    if (is_3smooth(g))
                        for (auto i = first + g; i < last; i++)
                            if (comp(*i, *(i - g))) std::iter_swap(i - g, i);
            }

            template<class RandomAccessIterator, class Compare>
            void comp_swap(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            template<class RandomAccessIterator, class Compare>
            void iterative_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int pow2 = (int)(std::log(length - 1) / std::log(2));

                for (int k = pow2; k >= 0; k--) {
                    int pow3 = (int)((std::log(length) - k * std::log(2)) / std::log(3));
                    for (int j = pow3; j >= 0; j--) {
                        int gap = (int)(std::pow(2, k) * std::pow(3, j));
                        for (auto i = first; i + gap < last; i++) comp_swap(i, i + gap, comp);
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _power_of_three(RandomAccessIterator pos, int gap, RandomAccessIterator end, Compare comp) {
                if (pos + gap > end) return;

                _power_of_three(pos, gap * 3, end, comp);
                _power_of_three(pos + gap, gap * 3, end, comp);
                _power_of_three(pos + 2 * gap, gap * 3, end, comp);

                for (auto i = pos; i + gap < end; i += gap) comp_swap(i, i + gap, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _recursive_comb_sort_loop(RandomAccessIterator pos, int gap, RandomAccessIterator end, Compare comp) {
                if (pos + gap > end) return;

                _recursive_comb_sort_loop(pos, gap * 2, end, comp);
                _recursive_comb_sort_loop(pos + gap, gap * 2, end, comp);
                _power_of_three(pos, gap, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void recursive_3smooth_comb_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _recursive_comb_sort_loop(first, 1, last, comp);
            }
        }

        namespace shove_sorting {
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
        }

        namespace bogo_sorting {
            // A seeder for the random number generators.
            std::random_device rd;

            // Return a random integer in the range [start, end) use the specified random number generator.
            // start - The start of the range.
            // end - The end of the range (exclusive).
            // rng - The random number generator.
            template<class RandomNumberGenerator>
            int randint(int start, int end, RandomNumberGenerator &rng) {
                // rng.seed(rd());   // Reseed the random number generator.
                std::uniform_int_distribution<> distr(start, end - 1);

                return distr(rng);
            }

            // Return a random boolean use the specified random number generator.
            // rng - The random number generator.
            template<class RandomNumberGenerator>
            bool randbool(RandomNumberGenerator &rng) {
                return randint(0, 1, rng);
            }

            // Return a random iterator in the range [first, last) use the specified random number generator.
            // And reseed the random number generator.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // rng - The random number generator.
            template<class RandomAccessIterator, class RandomNumberGenerator>
            RandomAccessIterator randiter(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng) {
                return first + randint(0, last - first, rng);
            }

            // Shuffle the range [first, last) use the specified random number generator.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // rng - The random number generator.
            template<class RandomAccessIterator, class RandomNumberGenerator>
            void bogo_swap(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng) {
                for (auto i = first; i < last; i++) std::iter_swap(i, randiter(i, last, rng));
            }

            // Generate a random sequence in the range [first, last) with size ones use the specified random number generator.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // size - The number of ones.
            // rng - The random number generator.
            template<class RandomAccessIterator, class RandomNumberGenerator>
            void bogo_combo(RandomAccessIterator first, RandomAccessIterator last, int size, RandomNumberGenerator &rng) {
                for (auto i = first; i < last; i++) (*i) = 0;
                for (auto i = last - size; i < last; i++) {
                    auto j = randiter(first, i + 1, rng);
                    auto pos = (*j == 0) ? j : i;
                    *pos = 1;
                }
            }


            // Determine whether the range [first, last) is sorted with the specified comparison function.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool is_range_sorted(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last - 1; i++)
                    if (comp(*(i + 1), *i)) return false;
                return true;
            }

            // Determine whether the range [first, last) is partitioned with the specified pivot and comparison function.
            // first - The beginning of the range.
            // pivot - The pivot.
            // last - The end of the range (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool is_range_partitioned(RandomAccessIterator first, RandomAccessIterator pivot, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < pivot; i++)
                    if (comp(*pivot, *i)) return false;
                for (auto i = pivot; i < last; i++)
                    if (comp(*i, *pivot)) return false;
                return true;
            }

            // Determine whether the minimum of the range [first, last) is on the left side with the specified comparison function.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool is_min_sorted(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                return is_range_partitioned(first, first, last, comp);
            }

            // Determine whether the maximum of the range [first, last) is on the right side with the specified comparison function.
            // first - The beginning of the range.
            // last - The end of the range (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool is_max_sorted(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                return is_range_partitioned(first, last - 1, last, comp);
            }

            // Determine whether the range [first, last) is split with the specified middle and comparison function.
            // first - The beginning of the range.
            // mid - The middle.
            // last - The end of the range (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool is_range_split(RandomAccessIterator first, RandomAccessIterator mid, RandomAccessIterator last, Compare comp) {
                auto low_max = first;
                for (auto i = first + 1; i < mid; i++)
                    if (comp(*low_max, *i)) low_max = i;

                for (auto i = mid; i < last; i++)
                    if (comp(*i, *low_max)) return false;

                return true;
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void bubble_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
                using sortcpp::_internal::bogo_sorting::is_range_sorted;
                using sortcpp::_internal::bogo_sorting::randiter;

                while (!is_range_sorted(first, last, comp)) {
                    auto i = randiter(first, last - 1, rng);
                    if (comp(*(i + 1), *i)) std::iter_swap(i, i + 1);
                }
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void exchange_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
                using sortcpp::_internal::bogo_sorting::is_range_sorted;
                using sortcpp::_internal::bogo_sorting::randiter;

                while (!is_range_sorted(first, last, comp)) {
                    auto i = randiter(first, last, rng);
                    auto j = randiter(first, last, rng);

                    if (i < j) {
                        if (comp(*j, *i)) std::iter_swap(i, j);
                    }
                    else {
                        if (comp(*i, *j)) std::iter_swap(i, j);
                    }
                }
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void bogobogo_sort_loop(RandomAccessIterator start, RandomAccessIterator end, RandomNumberGenerator& rng, Compare comp) {
                if (end - start <= 1) return;

                bogobogo_sort_loop(start, end - 1, rng, comp);
                while (comp(*(end - 1), *(end - 2))) {
                    bogo_swap(start, end, rng);
                    bogobogo_sort_loop(start, end - 1, rng, comp);
                }
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void median_quick_bogo_loop(RandomAccessIterator start, RandomAccessIterator end, RandomNumberGenerator &rng, Compare comp) {
                if (start >= end - 1) return;

                auto mid = start + (end - start) / 2;
                while (!is_range_split(start, mid, end, comp)) bogo_swap(start, end, rng);

                median_quick_bogo_loop(start, mid, rng, comp);
                median_quick_bogo_loop(mid, end, rng, comp);
            }

            template<class T, class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void bogo_weave(RandomAccessIterator first, RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, 
                            std::vector<T> &tmp, RandomNumberGenerator& rng, Compare comp) {
                int length = end - start;
                std::vector<int> combo(length);
                bogo_combo(combo.begin(), combo.end(), end - mid, rng);

                int low = start - first, high = mid - first;
                for (int i = 0; i < length; i++) {
                    if (combo[i] == 0) *(start + i) = tmp[low++];
                    else *(start + i) = tmp[high++];
                }
            }

            template<class T, class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void _merge_bogo_loop(RandomAccessIterator first, RandomAccessIterator start, RandomAccessIterator end, 
                                 std::vector<T> &tmp, RandomNumberGenerator& rng, Compare comp) {
                if (start >= end - 1) return;

                auto mid = start + (end - start) / 2;

                _merge_bogo_loop(first, start, mid, tmp, rng, comp);
                _merge_bogo_loop(first, mid, end, tmp, rng, comp);

                for (auto i = start; i < end; i++) tmp[i - first] = *i;
                while (!is_range_sorted(start, end, comp)) bogo_weave<T>(first, start, mid, end, tmp, rng, comp);
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void merge_bogo_loop(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator &rng, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                int length = last - first;
                std::vector<T> tmp(length);

                _merge_bogo_loop<T>(first, first, last, tmp, rng, comp);
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            RandomAccessIterator quick_bogo_swap(RandomAccessIterator start, RandomAccessIterator pivot, RandomAccessIterator end, 
                                                 RandomNumberGenerator &rng, Compare comp) {
                for (auto i = start; i < end; i++) {
                    auto j = randiter(i, end, rng);
                    if (pivot == i) pivot = j;
                    else if (pivot == j) pivot = i;
                    std::iter_swap(i, j);
                }
                return pivot;
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void quick_bogo_loop(RandomAccessIterator start, RandomAccessIterator end, RandomNumberGenerator &rng, Compare comp) {
                if (start >= end - 1) return;

                auto pivot = start;
                while (!is_range_partitioned(start, pivot, end, comp))
                    pivot = quick_bogo_swap(start, pivot, end, rng, comp);

                quick_bogo_loop(start, pivot, rng, comp);
                quick_bogo_loop(pivot + 1, end, rng, comp);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator find_last_sorted(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto i = first + 1;
                for (; i < last && !comp(*i, *(i - 1)); i++);
                return i - 1;
            }
        }

        namespace grate_sorting {
            template<class RandomAccessIterator, class Compare>
            void cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bool sorted = false;
                while (!sorted) {
                    sorted = true;
                    RandomAccessIterator i;
                    for (i = first; i < last - 1; i++)
                        for (auto j = last - 1; j > i; j--)
                            if (comp(*j, *i)) {
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }

                    if (sorted) break;

                    for (i = first; i < last - 1; i++)
                        for (auto j = i + 1; j < last; j++)
                            if (comp(*j, *i)) {
                                std::iter_swap(i, j);
                                break;
                            }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bool sorted = false;
                while (!sorted) {
                    sorted = true;
                    for (auto i = first; i < last - 1; i++)
                        for (auto j = last - 1; j > i; j--)
                            if (comp(*j, *i)) {
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_cocktail_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto a_bound = last - 1;
                auto left = first;
                auto a_right = last - 1;
                auto a_first_swap = first;
                auto a_last_swap = first;
                auto a_testi = last - 1;
                auto b_bound = last - 1;
                auto b_first_swap = first;
                auto b_last_swap = first;
                bool sorted = false, a_higher_val = false;

                while (!sorted) {
                    if (!sorted) {
                        while (!a_higher_val) {
                            if (a_testi < left) {
                                a_bound--;
                                a_testi = a_right;
                                if (a_bound < a_right) a_higher_val = true;
                            }
                            else {
                                if (comp(*a_bound, *a_testi)) a_higher_val = true;
                                else a_testi--;
                            }
                        }
                    }

                    sorted = true;
                    for (auto i = left; i < a_right; i++) {
                        for (auto j = a_bound; j > i; j--) {
                            if (comp(*j, *i)) {
                                if (sorted) a_first_swap = i;
                                a_last_swap = i;
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }
                        }
                    }

                    b_bound--;
                    a_testi = a_right;
                    a_higher_val = false;
                    left = a_first_swap;
                    a_right = a_last_swap + 1;

                    if (sorted) break;
                    sorted = true;
                    for (auto i = left; i < b_bound; i++) {
                        for (auto j = i + 1; j <= b_bound; j++) {
                            if (comp(*j, *i)) {
                                if (sorted) b_first_swap = i;
                                b_last_swap = i;
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }
                        }
                    }

                    b_bound = b_last_swap;
                    left = b_first_swap;
                    a_bound = b_bound;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto bound = last - 1;
                auto left = first;
                auto right = last - 1;
                auto first_swap = first;
                auto last_swap = first;
                auto testi = last - 1;
                bool sorted = false;
                bool higher_val = false;

                while (!sorted) {
                    if (!sorted) {
                        while (!higher_val) {
                            if (testi < left) {
                                bound--;
                                testi = right;
                                if (bound < right) higher_val = true;
                            }
                            else {
                                if (comp(*bound, *testi)) higher_val = true;
                                else testi--;
                            }
                        }
                    }

                    sorted = true;
                    for (auto i = left; i < right; i++) {
                        for (auto j = bound; j > i; j--) {
                            if (comp(*j, *i)) {
                                if (sorted) first_swap = i;
                                last_swap = i;
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }
                        }
                    }

                    bound--;
                    testi = right;
                    higher_val = false;
                    left = first_swap;
                    right = last_swap + 1;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto bound = last - 1;
                auto left = first;
                auto first_swap = first;
                auto last_swap = first;
                bool sorted = false;

                while (!sorted) {
                    sorted = true;
                    for (auto i = left; i < bound; i++) {
                        for (auto j = i + 1; j <= bound; j++) {
                            if (comp(*j, *i)) {
                                if (sorted) first_swap = i;
                                last_swap = i;
                                sorted = false;
                                std::iter_swap(i, j);
                                break;
                            }
                        }
                    }

                    bound = last_swap;
                    left = first_swap;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void reverse_grate_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bool sorted = false;
                while (!sorted) {
                    sorted = true;
                    for (auto i = first; i < last - 1; i++) {
                        for (auto j = i + 1; j < last; j++) {
                            if (comp(*j, *i)) {
                                sorted = false;
                                std::iter_swap(i, j);
                            }
                        }
                    }
                }
            }
        }

        namespace heap_sorting {
            // Determine whether a is less than b with is_max = false, or a is greater than b with is_max = true.
            // a - The first element.
            // b - The second element.
            // is_max - Indicates whether a large root heap is used.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool heap_compare(RandomAccessIterator a, RandomAccessIterator b, bool is_max, Compare comp) {
                if (is_max) return comp(*a, *b);
                else return comp(*b, *a);
            }


            // Sift down the root of the heap.
            // array - The real start position of the array.
            // root - The position of the root.
            // dist - The length of the heap.
            // start - The beginning position of the heap.
            // is_max - Indicates whether a large root heap is used.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void sift_down(RandomAccessIterator array, int root, int dist, int start, bool is_max, Compare comp) {
                while (root <= dist / 2) {
                    int leaf = 2 * root;                                                                                     // The left child of the root.
                    if (leaf < dist && heap_compare(array + start + leaf - 1, array + start + leaf, is_max, comp)) leaf++;   // If the right child is larger, choose the right child.
                    if (heap_compare(array + start + root - 1, array + start + leaf - 1, is_max, comp)) {                    // If the root is larger than the child, swap them.
                        std::iter_swap(array + start + root - 1, array + start + leaf - 1);
                        root = leaf;                                                                                         // Perform similar operations on the subtree.
                    }
                    else break;
                }
            }

            // Make a heap in range [low, high) of the array.
            // If is_max = true, it is a large root heap, otherwise it is a small root heap.
            // array - The real start position of the array.
            // low - The beginning position of the heap.
            // high - The end position of the heap.
            // is_max - Indicates whether a large root heap is used.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void heapify(RandomAccessIterator array, int low, int high, bool is_max, Compare comp) {
                int length = high - low;
                for (int i = length / 2; i >= 1; i--) sift_down(array, i, length, low, is_max, comp);
            }

            // Run heap sort in range [start, end) of the array.
            // If is_max = true, it will be a large heap sort, otherwise it will be a small heap sort.
            // array - The real start position of the array.
            // start - The beginning position of the heap sort.
            // end - The end position of the heap sort (exclusive).
            // is_max - Indicates whether a large root heap is used.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void heap_sort(RandomAccessIterator array, int start, int end, bool is_max, Compare comp) {
                using sortcpp::_internal::reversal;

                heapify(array, start, end, is_max, comp);                          // Make a heap between [start, end).
                for (int i = end - start; i > 1; i--) {                            // Pop elements from the heap one by one.
                    std::iter_swap(array + start, array + start + i - 1);          // Pop the top element of the heap to the last of the array.
                    sift_down(array, 1, i - 1, start, is_max, comp);               // Adjust the current heap.
                }
                if (!is_max) reversal(array + start, array + end);                 // If it's a min heap sort, the array is reverse sorted. So reverse it.
            }

            // Run heap sort in range [first, last) of the array.
            // If is_max = true, it will be a large heap sort, otherwise it will be a small heap sort.
            // first - The beginning position of the array.
            // last - The end position of the array (exclusive).
            // is_max - Indicates whether a large root heap is used.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void heap_sort(RandomAccessIterator first, RandomAccessIterator last, bool is_max, Compare comp) {
                heap_sort(first, 0, last - first, is_max, comp);
            }

            
        }

        namespace binary_insertion_sorting {
            using sortcpp::_internal::_fib_search;
            using sortcpp::_internal::_tri_search;
            using sortcpp::_internal::insertion_sorting::insertion_sort;

            // Run binary insertion sort in range [first, last) of the array.
            // first - The beginning position of the array.
            // last - The end position of the array (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void binary_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last; i++) {         // Sort the range [first, i).
                    // Find the insertion position of the current element.
                    auto num = *i;
                    auto lo = first, hi = i;

                    while (lo < hi) {
                        auto mid = lo + (hi - lo) / 2;
                        if (comp(num, *mid)) hi = mid;
                        else lo = mid + 1;
                    }

                    // Move elements of arr[0..i-1] that are greater than the current element to one position ahead of their current position.
                    auto j = i - 1;
                    while (j >= lo) {
                        *(j + 1) = *j;
                        j--;
                    }
                    *lo = num;                                // Insert the current element at the found position.
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator left_binary_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(*m, val)) a = m + 1;
                    else b = m;
                }
                return a;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator right_binary_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(val, *m)) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class RandomAccessIterator, class T>
            void insert_to_left(RandomAccessIterator a, RandomAccessIterator b, const T& temp) {
                while (a > b) {
                    *a = *(a - 1);
                    a--;
                }
                *b = temp;
            }

            template<class RandomAccessIterator, class T>
            void insert_to_right(RandomAccessIterator a, RandomAccessIterator b, const T& temp) {
                while (a < b) {
                    *a = *(a + 1);
                    a++;
                }
                *a = temp;
            }

            template<class RandomAccessIterator, class Compare>
            void binary_double_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (b - a < 2) return;

                auto j = a + (b - a - 2) / 2 + 1;
                auto i = a + (b - a - 1) / 2;

                if (j > i && comp(*j, *i)) std::iter_swap(i, j);

                i--;
                j++;

                while (j < b) {
                    if (comp(*j, *i)) {
                        auto l = *j, r = *i;
                        auto m = right_binary_search(i + 1, j, l, comp);
                        insert_to_right(i, m - 1, l);
                        insert_to_left(j, left_binary_search(m, j, r, comp), r);
                    }
                    else {
                        auto l = *i, r = *j;
                        auto m = left_binary_search(i + 1, j, l, comp);
                        insert_to_right(i, m - 1, l);
                        insert_to_left(j, right_binary_search(m, j, r, comp), r);
                    }
                    i--;
                    j++;
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator right_exp_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                int i = 1;
                while (b - i >= a && comp(val, *(b - i))) i *= 2;
                return right_binary_search(std::max(a, b - i + 1), b - i / 2, val, comp);
            }

            template<class RandomAccessIterator>
            void insert_to(RandomAccessIterator a, RandomAccessIterator b) {
                if (a > b) {
                    auto temp = *a;
                    do {
                        *a = *(a - 1);
                        a--;
                    } while (a > b);
                    *b = temp;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void exp_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                for (auto i = a + 1; i < b; i++)
                    insert_to(i, right_exp_search(a, i, *i, comp));
            }

            template<class RandomAccessIterator, class Compare>
            void fib_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                for (auto i = a + 1; i < b; i++) {
                    auto tmp = *i;
                    auto position = _fib_search(a, i - 1, tmp, comp);
                    auto j = i - 1;
                    while (j >= position) {
                        *(j + 1) = *j;
                        j--;
                    }
                    *(j + 1) = tmp;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void tri_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                for (auto i = a + 1; i < b; i++) {
                    auto tmp = *i;
                    auto position = _tri_search(a, i - 1, tmp, comp);
                    auto j = i - 1;
                    while (j >= position) {
                        *(j + 1) = *j;
                        j--;
                    }
                    *(j + 1) = tmp;
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator bin_search(RandomAccessIterator begin, RandomAccessIterator end, const T& target, Compare comp) {
                while (true) {
                    int delta = end - begin;
                    if (delta <= 0) break;

                    auto p = begin + delta / 2;
                    if (equals(*p, target, comp)) return p;
                    if (comp(target, *p)) {
                        end = p;
                        continue;
                    }
                    begin = p + 1;
                }
                return end;
            }

            template<class RandomAccessIterator, class Compare>
            void bin_insert(RandomAccessIterator start, RandomAccessIterator end, int len, Compare comp) {
                int offset = 1;
                for (; offset * offset < len; offset *= 2);
                for (auto b_start = start, b_end = end, i = start + offset; i < end; i++) {
                    auto target = bin_search(b_start, b_end, *i, comp);
                    auto tmp = *i;
                    auto j = i;
                    while (j > target && !comp(*j, tmp)) {
                        std::iter_swap(j, j - 1);
                        j--;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void gambit_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                bin_insert(a, b, b - a, comp);
                insertion_sort(a, b, comp);
            }

            
        }

        namespace adaptive_insertion_sorting {
            using sortcpp::_internal::reversal;
            using sortcpp::_internal::equals;
            using sortcpp::_internal::binary_insertion_sorting::right_binary_search;
            using sortcpp::_internal::binary_insertion_sorting::right_exp_search;
            using sortcpp::_internal::binary_insertion_sorting::insert_to;

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator find_run(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = a + 1;
                
                if (comp(*i, *(i - 1))) {
                    do i++;
                    while (i < b && comp(*i, *(i - 1)));
                    reversal(a, i);
                }
                else {
                    do i++;
                    while (i < b && !comp(*i, *(i - 1)));
                }
                
                return i;
            }

            template<class RandomAccessIterator, class Compare>
            void move_front(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, Compare comp) {
                auto mini = m;
                auto minv = *mini;

                for (auto i = m + 1; i < b; i++) {
                    if (comp(*i, minv)) {
                        mini = i;
                        minv = *i;
                    }
                }

                m--;
                while (mini > m) {
                    *mini = *(mini - 1);
                    mini--;
                }

                a--;
                while (mini > a && comp(minv, *mini)) {
                    *(mini + 1) = *mini;
                    mini--;
                }
                *(mini + 1) = minv;
            }

            template<class RandomAccessIterator, class Compare>
            void adaptive_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = find_run(a, b, comp);
                if (i < b) {
                    move_front(a, i++, b, comp);
                    while (i < b) {
                        auto current = *i;
                        auto pos = i - 1;
                        
                        while (comp(current, *pos)) {
                            *(pos + 1) = *pos;
                            pos--;
                        }
                        *(pos + 1) = current;
                        i++;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void adaptive_binary_insertion_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                int count = 0;
                for (auto i = start + 1; i < end; i++) {
                    auto num = *i;
                    int v = (2 * count / (i - start)) + 1;
                    auto lo = std::max(i - v, start), hi = i;

                    while (lo >= start && comp(num, *lo)) {
                        lo -= v;
                        hi -= v;
                    }
                    lo++;
                    if (lo < start) lo = start;

                    while (lo < hi) {
                        auto mid = lo + (hi - lo) / 2;

                        if (comp(num, *mid)) hi = mid;
                        else lo = mid + 1;
                    }

                    count += (i - lo);

                    auto j = i - 1;

                    if (j >= lo) {
                        while (j >= lo) {
                            *(j + 1) = *j;
                            j--;
                        }
                        *lo = num;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pd_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = find_run(a, b, comp);

                while (i < b) {
                    auto current = *i;
                    auto pos = i - 1;
                    while (pos >= a && comp(current, *pos)) {
                        *(pos + 1) = *pos;
                        pos--;
                    }
                    *(pos + 1) = current;
                    i++;
                }
            }

            
            template<class RandomAccessIterator, class Compare>
            void stable_segment_reversal(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                if (end - start < 3) std::iter_swap(start, end);
                else reversal(start, end + 1);
                RandomAccessIterator i = start, left, right;
                while (i < end) {
                    left = i;
                    while (equals(*i, *(i + 1), comp) && i < end) i++;
                    right = i;
                    if (left != right) {
                        if (right - left < 3) std::iter_swap(left, right);
                        else reversal(left, right + 1);
                    }
                    i++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator pd(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto forward = start;
                int cmp = comp(*(forward + 1), *forward) - comp(*forward, *(forward + 1));
                bool less_unique = false;
                while (cmp <= 0 && forward + 1 < end) {
                    forward++;
                    if (cmp == 0) less_unique = true;
                    if (forward + 1 < end) cmp = comp(*(forward + 1), *forward) - comp(*forward, *(forward + 1));
                }

                auto reverse = start;
                if (forward == start) {
                    bool different = false;
                    cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                    while (cmp >= 0 && reverse + 1 < end) {
                        if (cmp == 0) less_unique = true;
                        else different = true;
                        reverse++;
                        if (reverse + 1 < end) cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                    }
                    if (reverse > start && different) {
                        if (less_unique) stable_segment_reversal(start, reverse, comp);
                        else if (reverse < start + 3) std::iter_swap(start, reverse);
                        else reversal(start, reverse + 1);
                    }
                }
                return std::max(forward, reverse);
            }

            template<class RandomAccessIterator, class Compare>
            void pd_binary_insertion_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto pattern = pd(start, end, comp);
                for (auto i = pattern + 1; i < end; i++) {
                    auto item = *i;
                    auto left = right_binary_search(start, i, item, comp);

                    bool w = false;
                    for (auto right = i; right > left; right--) {
                        *right = *(right - 1);
                        w = true;
                    }

                    if (w) *left = item;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pd_exp_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = find_run(a, b, comp);
                while (i < b) {
                    insert_to(i, right_exp_search(a, i, *i, comp));
                    i++;
                }
            }
        }

        namespace shell_sorting {
            // Extra Ciura Gaps
            constexpr int incs[] = { 1, 4, 10, 23, 57, 132, 301, 701, 1750, 3938, 8861 };
            constexpr int incs_count = sizeof(incs) / sizeof(int);

            // Return the n-th Extra Ciura Gap.
            // n - The index of the gap.
            int ciura(int n) {
                if (n < incs_count) return incs[n];
                return (int)std::pow(2.25, n);
            }

            template<class RandomAccessIterator, class Compare>
            void gapped_insertion_sort(RandomAccessIterator a, RandomAccessIterator b, int gap, Compare comp) {
                for (auto i = a + gap; i < b; i++) {
                    auto key = *i;
                    auto j = i - gap;

                    while (j >= a && comp(key, *j)) {
                        *(j + gap) = *j;
                        j -= gap;
                    }
                    *(j + gap) = key;
                }
            }


            // Run shell sort in range [first, last) of the array.
            // first - The beginning position of the array.
            // last - The end position of the array (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first, gap = 1, k;
                for (k = 0; gap < length; k++) gap = ciura(k);

                while (--k >= 0) {
                    gap = ciura(k);
                    gapped_insertion_sort(first, last, gap, comp);
                }
            }

            // Run quick shell sort in range [first, last) of the array (for small arrays).
            // first - The beginning position of the array.
            // last - The end position of the array (exclusive).
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void quick_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                static const int q_incs[] = { 48, 21, 7, 3, 1 };
                static const int q_incs_count = sizeof(q_incs) / sizeof(int);
                for (int k = 0; k < q_incs_count; k++)
                    if (q_incs[k] < length) gapped_insertion_sort(first, last, q_incs[k], comp);
            }

            template<class RandomAccessIterator, class Compare>
            void cocktail_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first, gap = 1, k;
                for (k = 0; gap < length; k++) gap = ciura(k);

                bool dir = true;
                while (--k >= 0) {
                    gap = ciura(k);

                    if (dir) {
                        for (auto i = first + gap; i < last; i++) {
                            auto tmp = *i;
                            auto j = i;
                            
                            while (j >= first + gap && comp(tmp, *(j - gap))) {
                                *j = *(j - gap);
                                j -= gap;
                            }

                            *j = tmp;
                        }
                    }
                    else {
                        for (auto i = last - gap; i >= first; i--) {
                            auto tmp = *i;
                            auto j = i;

                            while (j < last - gap && comp(*(j + gap), tmp)) {
                                *j = *(j + gap);
                                j += gap;
                            }

                            *j = tmp;
                        }
                    }

                    dir = !dir;
                }
            }

            
            template<class RandomAccessIterator, class Compare>
            void _recusive_shell_sort(RandomAccessIterator start, RandomAccessIterator end, int g, Compare comp) {
                if (start + g <= end) {
                    _recusive_shell_sort(start, end, 3 * g, comp);
                    _recusive_shell_sort(start + g, end, 3 * g, comp);
                    _recusive_shell_sort(start + 2 * g, end, 3 * g, comp);
                    gapped_insertion_sort(start, end, g, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void recursive_shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _recusive_shell_sort(first, last, 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void shell_sort_primes(RandomAccessIterator first, RandomAccessIterator last, bool high, Compare comp) {
                int length = last - first;
                int prime_test_running = 1, prime_testi = 2, gap = length;
                bool prime_test = false;

                while (gap != 1) {
                    prime_test_running = gap;
                    while (high? prime_test_running != 1: prime_test_running == gap) {
                        prime_test = false;
                        prime_testi = 2;
                        while (!prime_test) {
                            if (prime_test_running % prime_testi == 0) {
                                prime_test_running /= prime_testi;
                                prime_test = true;
                            }
                            else prime_testi++;
                        }
                    }

                    gap /= prime_testi;
                    gapped_insertion_sort(first, last, gap, comp);
                }
            }
        }

        namespace par_shell_sorting {
            using sortcpp::_internal::shell_sorting::gapped_insertion_sort;

            constexpr double mult = 1.5;
            double threshold(int x) { return std::sqrt(x); }

            template<class RandomAccessIterator, class Compare>
            int get_par(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int len = last - first;
                std::vector<bool> is_max(len, false);

                auto maximum = *first;
                for (auto i = first + 1; i < last; i++) 
                    if (comp(maximum, *i)) {
                        maximum = *i;
                        is_max[i - first] = true;
                    }

                auto i = last - 1, j = last - 1;
                int p = 1;

                while (j >= first && i >= first + p) {
                    while (!is_max[j - first] && j > first) j--;
                    maximum = *j;
                    while (!comp(*i, maximum) && i >= first + p) i--;
                    if (comp(*i, *j) && p < i - j) p = i - j;
                    j--;
                }
                return p;
            }

            template<class RandomAccessIterator, class Compare>
            int shell_pass(RandomAccessIterator first, RandomAccessIterator last, int gap, int last_gap, Compare comp) {
                if (gap >= last_gap) return last_gap;
                if (gap == last_gap - 1 && gap != 1) return last_gap;

                gapped_insertion_sort(first, last, gap, comp);
                return gap;
            }

            template<class RandomAccessIterator, class Compare>
            void par_shell_sort(RandomAccessIterator first, RandomAccessIterator last, int constant_div, Compare comp) {
                double true_div = constant_div;
                int len = last - first;
                int last_par = len, last_gap = len;

                while (true) {
                    int par = get_par(first, last, comp);
                    if (par >= last_par) par = last_par - (int)true_div;
                    if (par / (int)true_div <= 1) {
                        last_gap = shell_pass(first, last, 1, last_gap, comp);
                        break;
                    }
                    last_gap = shell_pass(first, last, (int)((par / (int)true_div) + par % (int)true_div), last_gap, comp);
                    if (last_par - par <= threshold(last_par)) true_div *= mult;
                    last_par = par;
                }
            }
        }

        namespace merge_sorting {
            // Do merge sort in range [first, last) of the array with a temporary buffer.
            // start - The beginning position of the array.
            // mid - The middle position of the array.
            // end - The end position of the array (exclusive).
            // tmp - The beginning position of the temporary buffer.
            // binary - Whether to use binary insertion sort for small arrays.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Buffer, class Compare>
            void internal_merge_sort(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, 
                                     Buffer tmp, bool binary, Compare comp) {
                using sortcpp::_internal::binary_insertion_sorting::binary_insertion_sort;

                if (end - start < 32 && binary) return;
                if (end - start < 64 && binary) {
                    binary_insertion_sort(start, end, comp);
                    return;
                }

                internal_merge_sort(start, start + (mid - start) / 2, mid, tmp, binary, comp);
                internal_merge_sort(mid, mid + (end - mid) / 2, end, tmp, binary, comp);

                auto low = start, high = mid;
                for (int nxt = 0; nxt < end - start; nxt++) {
                    if (low >= mid && high >= end) break;

                    if (low < mid && high >= end) {
                        *(tmp + nxt) = *low;
                        low++;
                    }
                    else if (low >= mid && high < end) {
                        *(tmp + nxt) = *high;
                        high++;
                    }
                    else if (comp(*high, *low)) {
                        *(tmp + nxt) = *high;
                        high++;
                    }
                    else {
                        *(tmp + nxt) = *low;
                        low++;
                    }
                }

                for (int i = 0; i < end - start; i++) *(start + i) = *(tmp + i);
            }

            // Do merge sort in range [first, last) of the array.
            // first - The beginning position of the array.
            // last - The end position of the array (exclusive).
            // binary - Whether to use binary insertion sort for small arrays.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void merge_sort(RandomAccessIterator first, RandomAccessIterator last, bool binary, Compare comp) {
                using sortcpp::_internal::binary_insertion_sorting::binary_insertion_sort;
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                int length = last - first;
                if (length < 32) {
                    binary_insertion_sort(first, last, comp);
                    return;
                }

                std::vector<T> tmp(length);
                auto start = first, end = last, mid = start + length / 2;
                internal_merge_sort(start, mid, end, tmp.begin(), binary, comp);
            }

            template<class RandomAccessIterator>
            void multi_swap(RandomAccessIterator a, RandomAccessIterator b, int len) {
                for (int i = 0; i < len; i++) std::iter_swap(a + i, b + i);
            }

            template<class RandomAccessIterator, class Compare>
            int binary_search_mid(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, Compare comp) {
                int a = 0, b = std::min(mid - start, end - mid), m = a + (b - a) / 2;

                while (b > a) {
                    if (comp(*(mid + m), *(mid - m - 1))) a = m + 1;
                    else b = m;
                    m = a + (b - a) / 2;
                }
                return m;
            }

            template<class RandomAccessIterator, class Compare>
            void multi_swap_merge(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, Compare comp) {
                int m = binary_search_mid(start, mid, end, comp);
                while (m > 0) {
                    multi_swap(mid - m, mid, m);
                    multi_swap_merge(mid, mid + m, end, comp);

                    end = mid;
                    mid -= m;

                    m = binary_search_mid(start, mid, end, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void multi_swap_merge_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                int len = b - a;
                RandomAccessIterator i;

                for (int j = 1; j < len; j *= 2) {
                    for (i = a; i + 2 * j <= b; i += 2 * j) multi_swap_merge(i, i + j, i + 2 * j, comp);
                    if (i + j < b) multi_swap_merge(i, i + j, b, comp);
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void bottom_up_merge(RandomAccessIterator first, RandomAccessIterator last, RandomAccessIterator index, 
                                 int merge_size, std::vector<T> &scratch, int &copy_length, Compare comp) {
                auto left = index, mid = left + merge_size / 2, right = mid, end = std::min(last, index + merge_size);
                int scratch_index = left - first;

                if (right < end) {
                    while (left < mid && right < end) {
                        if (comp(*right, *left)) {
                            scratch[scratch_index++] = *right;
                            right++;
                        }
                        else {
                            scratch[scratch_index++] = *left;
                            left++;
                        }
                    }

                    if (left < mid) {
                        while (left < mid) {
                            scratch[scratch_index++] = *left;
                            left++;
                        }
                    }

                    if (right < end) {
                        while (right < end) {
                            scratch[scratch_index++] = *right;
                            right++;
                        }
                    }
                }
                else copy_length = left - first;
            }

            template<class RandomAccessIterator, class Compare>
            void bottom_up_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                int length = last - first;
                std::vector<T> scratch(length);
                int copy_length, merge_size = 2;

                while (merge_size <= length) {
                    copy_length = length;
                    for (auto i = first; i < last; i += merge_size) bottom_up_merge(first, last, i, merge_size, scratch, copy_length, comp);
                    for (int i = 0; i < copy_length; i++) *(first + i) = scratch[i];
                    merge_size *= 2;
                }

                if ((merge_size / 2) != length) {
                    bottom_up_merge(first, last, first, merge_size, scratch, copy_length, comp);
                    for (int i = 0; i < length; i++) *(first + i) = scratch[i];
                }
            }

            template<class RandomAccessIterator>
            void _rotate(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) {
                int l = m - a, r = b - m;

                while (l > 0 && r > 0) {
                    if (r < l) {
                        multi_swap(m - r, m, r);
                        b -= r;
                        m -= r;
                        l -= r;
                    }
                    else {
                        multi_swap(a, m, l);
                        a += l;
                        m += l;
                        r -= l;
                    }
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator _binary_search(RandomAccessIterator a, RandomAccessIterator b, const T& value, bool left, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    bool cmp = left ? (!comp(*m, value)) : comp(value, *m);
                    if (cmp) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            void rotate_merge(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, Compare comp) {
                RandomAccessIterator m1, m2, m3;

                if (m - a >= b - m) {
                    m1 = a + (m - a) / 2;
                    m2 = _binary_search(m, b, *m1, true, comp);
                    m3 = m1 + (m2 - m);
                }
                else {
                    m2 = m + (b - m) / 2;
                    m1 = _binary_search(a, m, *m2, false, comp);
                    m3 = (m2++) - (m - m1);
                }
                rotate(m1, m, m2);

                if (m2 - (m3 + 1) > 0 && b - m2 > 0) rotate_merge(m3 + 1, m2, b, comp);
                if (m1 - a > 0 && m3 - m1 > 0) rotate_merge(a, m1, m3, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void rotate_merge_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                int len = b - a;
                RandomAccessIterator i;

                for (int j = 1; j < len; j *= 2) {
                    for (i = a; i + 2 * j <= b; i += 2 * j) rotate_merge(i, i + j, i + 2 * j, comp);
                    if (i + j < b) rotate_merge(i, i + j, b, comp);
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void merge_to(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, std::vector<T>& sublist, Compare comp) {
                int i = 0, s = m - a;

                while (i < s && m < b) {
                    if (comp(sublist[i], *m)) *(a++) = sublist[i++];
                    else *(a++) = *(m++);
                }

                while (i < s) *(a++) = sublist[i++];
            }

            template<class RandomAccessIterator, class Compare>
            void strand_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                
                int length = last - first;
                std::vector<T> sublist(length);

                auto j = last, k = j;
                while (j > first) {
                    sublist[0] = *first;
                    k--;

                    int i = 0;
                    for (auto p = first, m = first + 1; m < j; m++) {
                        if (comp(*m, sublist[i])) *(p++) = *m;
                        else {
                            sublist[++i] = *m;
                            k--;
                        }
                    }

                    merge_to(k, j, last, sublist, comp);
                    j = k;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void _weaved_merge_sort(RandomAccessIterator residue, RandomAccessIterator end, int modulus, Compare comp) {
                if (residue + modulus >= end) return;

                auto low = residue, high = residue + modulus;
                int dmodulus = modulus << 1;

                _weaved_merge_sort<T>(low, end, dmodulus, comp);
                _weaved_merge_sort<T>(high, end, dmodulus, comp);

                std::vector<T> tmp(end - residue);
                int nxt = 0;

                for (; low < end && high < end; nxt += modulus) {
                    if (comp(*high, *low) || equals(*low, *high, comp) && low > high) {
                        tmp[nxt] = *high;
                        high += dmodulus;
                    }
                    else {
                        tmp[nxt] = *low;
                        low += dmodulus;
                    }
                }

                if (low >= end) {
                    while (high < end) {
                        tmp[nxt] = *high;
                        nxt += modulus;
                        high += dmodulus;
                    }
                }
                else {
                    while (low < end) {
                        tmp[nxt] = *low;
                        nxt += modulus;
                        low += dmodulus;
                    }
                }

                for (auto i = residue; i < end; i += modulus) *i = tmp[i - residue];
            }

            template<class RandomAccessIterator, class Compare>
            void weaved_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                _weaved_merge_sort<T>(first, last, 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void bubble(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                int c = 1;
                bool a = false;
                RandomAccessIterator s, f = start + (end - start) / 2;
                
                for (auto j = end - 1; j > start; j -= c) {
                    if (f - 1 < start) s = start;
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

            template<class RandomAccessIterator, class Compare>
            void mob_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int n = last - first, len = 2;
                RandomAccessIterator index = first;

                while (len <= n) {
                    index = first;
                    while (index + len <= last) {
                        if (len == 2) {
                            if (comp(*(index + 1), *index)) std::iter_swap(index, index + 1);
                        }
                        else bubble(index, index + len, comp);
                        index += len;
                    }
                    if (index != last) bubble(index, last, comp);
                    len *= 2;
                }
                bubble(first, last, comp);
            }

            template<class T, class RandomAccessIterator, class Compare>
            void natural_merge(RandomAccessIterator left, RandomAccessIterator right, RandomAccessIterator stop, 
                               std::vector<T>& merged, Compare comp) {
                bool first = true;
                int index = 0, start = 0;
                auto orig_right = right, orig_left = left;

                while (left < orig_right && right < stop)
                    if (comp(*right, *left)) {
                        first = false;
                        merged[index++] = *(right++);
                    }
                    else
                        if (first) {
                            start++;
                            index++;
                            left++;
                        }
                        else merged[index++] = *(left++);
                
                while (left < orig_right) merged[index++] = *(left++);
                for (int i = start; i < index; i++) *(i + orig_left) = merged[i];
            }

            template<class RandomAccessIterator, class Compare>
            void natural_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                int length = last - first;
                std::vector<T> merged(length);

                bool done = false;
                auto start = first, stop = last - 1;

                while (!done) {
                    auto prev = first, left = last;
                    done = true;
                    for (auto i = start; i < stop; i++)
                        if (comp(*(i + 1), *i)) {
                            if (left == last) {
                                left = prev;
                                prev = i + 1;
                            }
                            else {
                                natural_merge(left, prev, i + 1, merged, comp);
                                if (done) start = i;
                                prev = i + 1;
                                left = last;
                                done = false;
                            }
                        }
                    if (left != last) {
                        natural_merge(left, prev, last, merged, comp);
                        done = false;
                        stop = left;
                    }
                }
            }
            
            template<class T, class RandomAccessIterator, class Compare>
            void partial_merge(RandomAccessIterator left_start, RandomAccessIterator right_start, RandomAccessIterator end, 
                               std::vector<T>& copied, Compare comp) {
                for (int i = 0; i < right_start - left_start; i++) copied[i] = *(i + left_start);

                auto left = left_start;
                auto right = right_start;
                for (int nxt = 0; nxt < end - left_start; nxt++) {
                    if (left >= right_start && right >= end) break;

                    if (left < right_start && right >= end) {
                        *(nxt + left_start) = copied[(left++) - left_start];
                    }
                    else if (left >= right_start && right < end) {
                        *(nxt + left_start) = *(right++);
                    }
                    else if (!comp(*right, copied[left - left_start])) {
                        *(nxt + left_start) = copied[(left++) - left_start];
                    }
                    else {
                        *(nxt + left_start) = *(right++);
                    }
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void merge_run(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end,
                           std::vector<T>& copied, Compare comp) {
                if (start == mid) return;

                merge_run<T>(start, start + (mid - start) / 2, mid, copied, comp);
                merge_run<T>(mid, mid + (end - mid) / 2, end, copied, comp);
                partial_merge<T>(start, mid, end, copied, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void partial_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                int length = last - first;
                std::vector<T> copied(length / 2);

                merge_run<T>(first, first + length / 2, last, copied, comp);
            }
        }
        

        namespace iterative_top_down_merge_sorting {
            template<class RandomAccessIterator, class Compare>
            void internal_merge(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                std::vector<T> tmp(end - start);

                auto low = start, high = mid;
                int nxt = 0;

                for (; low < mid && high < end; nxt++) {
                    if (comp(*high, *low)) {
                        tmp[nxt] = *high;
                        high++;
                    }
                    else {
                        tmp[nxt] = *low;
                        low++;
                    }
                }

                if (low >= mid) {
                    while (high < end) {
                        tmp[nxt] = *high;
                        high++;
                        nxt++;
                    }
                }
                else {
                    while (low < mid) {
                        tmp[nxt] = *low;
                        low++;
                        nxt++;
                    }
                }

                for (auto i = start; i < end; i++) *i = tmp[i - start];
            }

            template<class RandomAccessIterator, class Merger, class Compare>
            void run_sort_small(RandomAccessIterator first, RandomAccessIterator last, Merger merger, Compare comp) {
                int length = last - first;
                for (int subarray_count = ceil_pow2(length); subarray_count > 1; subarray_count >>= 1)
                    for (int i = 0; i < subarray_count; i += 2)
                        merger(
                            first + length * i / subarray_count, 
                            first + length * (i + 1) / subarray_count, 
                            first + length * (i + 2) / subarray_count,
                            comp
                        );
            }

            template<class RandomAccessIterator, class Merger, class Compare>
            void run_sort_large(RandomAccessIterator first, RandomAccessIterator last, Merger merger, Compare comp) {
                int length = last - first;
                for (int subarray_count = ceil_pow2(length), whole_i = length / subarray_count, frac_i = length % subarray_count; subarray_count > 1; ) {
                    int frac = 0;
                    for (auto whole = first; whole < last; ) {
                        auto start = whole;
                        whole += whole_i;
                        frac += frac_i;
                        if (frac >= subarray_count) {
                            whole++;
                            frac -= subarray_count;
                        }

                        auto mid = whole;
                        whole += whole_i;
                        frac += frac_i;
                        if (frac >= subarray_count) {
                            whole++;
                            frac -= subarray_count;
                        }

                        merger(start, mid, whole, comp);
                    }

                    subarray_count >>= 1;
                    whole_i <<= 1;
                    if (frac_i >= subarray_count) {
                        whole_i++;
                        frac_i -= subarray_count;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                if (length < 1 << 15) run_sort_small(first, last, internal_merge<RandomAccessIterator, Compare>, comp);
                else run_sort_large(first, last, internal_merge<RandomAccessIterator, Compare>, comp);
            }
        }

        namespace new_shuffle_merge_sorting {
            using sortcpp::_internal::ceil_pow2;
            using sortcpp::_internal::equals;
            using sortcpp::_internal::iterative_top_down_merge_sorting::run_sort_small;
            using sortcpp::_internal::iterative_top_down_merge_sorting::run_sort_large;

            template<class RandomAccessIterator>
            void _rotate_equal(RandomAccessIterator a, RandomAccessIterator b, int size) {
                for (int i = 0; i < size; i++) std::iter_swap(a + i, b + i);
            }

            template<class RandomAccessIterator>
            void _rotate(RandomAccessIterator mid, int a, int b) {
                while (a > 0 && b > 0) {
                    if (a > b) {
                        _rotate_equal(mid - b, mid, b);
                        mid -= b;
                        a -= b;
                    }
                    else {
                        _rotate_equal(mid - a, mid, a);
                        mid += a;
                        b -= a;
                    }
                }
            }

            template<class RandomAccessIterator>
            void _shuffle_easy(RandomAccessIterator start, int size) {
                for (int i = 1; i < size; i *= 3) {
                    auto val = *(start + i - 1);
                    for (int j = i * 2 % size; j != i; j = j * 2 % size) {
                        auto nval = *(start + j - 1);
                        *(start + j - 1) = val;
                        val = nval;
                    }
                    *(start + i - 1) = val;
                }
            }

            template<class RandomAccessIterator>
            void _shuffle(RandomAccessIterator start, RandomAccessIterator end) {
                while (end - start > 1) {
                    int n = (end - start) / 2;
                    int l = 1;
                    while (l * 3 - 1 <= 2 * n) l *= 3;
                    int m = (l - 1) / 2;

                    _rotate(start + n, n - m, m);
                    _shuffle_easy(start, l);
                    start += l - 1;
                }
            }

            template<class RandomAccessIterator>
            void _rotate_shuffled_equal(RandomAccessIterator a, RandomAccessIterator b, int size) {
                for (int i = 0; i < size; i += 2) std::iter_swap(a + i, b + i);
            }

            template<class RandomAccessIterator>
            void _rotate_shuffled(RandomAccessIterator mid, int a, int b) {
                while (a > 0 && b > 0) {
                    if (a > b) {
                        _rotate_shuffled_equal(mid - b, mid, b);
                        mid -= b;
                        a -= b;
                    }
                    else {
                        _rotate_shuffled_equal(mid - a, mid, a);
                        mid += a;
                        b -= a;
                    }
                }
            }

            template<class RandomAccessIterator>
            void _rotate_shuffled_outer(RandomAccessIterator mid, int a, int b) {
                if (a > b) {
                    _rotate_shuffled_equal(mid - b, mid + 1, b);
                    mid -= b;
                    a -= b;
                    _rotate_shuffled(mid, a, b);
                }
                else {
                    _rotate_shuffled_equal(mid - a, mid + 1, a);
                    mid += a + 1;
                    b -= a;
                    _rotate_shuffled(mid, a, b);
                }
            }

            template<class RandomAccessIterator>
            void _unshuffle_easy(RandomAccessIterator start, int size) {
                for (int i = 1; i < size; i *= 3) {
                    int prev = i;
                    auto val = *(start + i - 1);
                    for (int j = i * 2 % size; j != i; j = j * 2 % size) {
                        *(start + prev - 1) = *(start + j - 1);
                        prev = j;
                    }
                    *(start + prev - 1) = val;
                }
            }

            template<class RandomAccessIterator>
            void _unshuffle(RandomAccessIterator start, RandomAccessIterator end) {
                while (end - start > 1) {
                    int n = (end - start) / 2;
                    int l = 1;
                    while (l * 3 - 1 <= 2 * n) l *= 3;
                    int m = (l - 1) / 2;

                    _rotate_shuffled_outer(start + 2 * m, 2 * m, 2 * n - 2 * m);
                    _unshuffle_easy(start, l);
                    start += l - 1;
                }
            }
            
            template<class RandomAccessIterator, class Compare>
            void _merge_up(RandomAccessIterator start, RandomAccessIterator end, bool type, Compare comp) {
                auto i = start, j = i + 1;
                while (j < end) {
                    if (comp(*i, *j) || !type && equals(*i, *j, comp)) {
                        i++;
                        if (i == j) {
                            j++;
                            type = !type;
                        }
                    }
                    else if (end - j == 1) {
                        _rotate(j, j - i, 1);
                        break;
                    }
                    else {
                        int r = 0;
                        if (type) while (j + 2 * r < end && !comp(*i, *(j + 2 * r))) r++;
                        else while (j + 2 * r < end && comp(*(j + 2 * r), *i)) r++;
                        j--;
                        _unshuffle(j, j + 2 * r);
                        _rotate(j, j - i, r);
                        i += r + 1;
                        j += 2 * r + 1;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void new_shuffle_merge(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, Compare comp) {
                if (mid - start <= end - mid) {
                    _shuffle(start, end);
                    _merge_up(start, end, true, comp);
                }
                else {
                    _shuffle(start + 1, end);
                    _merge_up(start, end, false, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                if (length < 1 << 15) run_sort_small(first, last, new_shuffle_merge<RandomAccessIterator, Compare>, comp);
                else run_sort_large(first, last, new_shuffle_merge<RandomAccessIterator, Compare>, comp);
            }
        }

        namespace pd_merge_sorting {
            using sortcpp::_internal::reversal;

            template<class T, class RandomAccessIterator, class Compare>
            void pd_merge_up(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, 
                             std::vector<T>& copied, Compare comp) {
                for (int i = 0; i < mid - start; i++) copied[i] = *(start + i);

                int buffer_ptr = 0;
                auto left = start, right = mid;

                while (left < right && right < end) {
                    if (comp(*right, copied[buffer_ptr])) {
                        *left = *right;
                        left++;
                        right++;
                    }
                    else {
                        *left = copied[buffer_ptr];
                        left++;
                        buffer_ptr++;
                    }
                }

                while (left < right) {
                    *left = copied[buffer_ptr];
                    left++;
                    buffer_ptr++;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void pd_merge_down(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, 
                               std::vector<T>& copied, Compare comp) {
                for (int i = 0; i < end - mid; i++) copied[i] = *(mid + i);

                int buffer_ptr = end - mid - 1;
                auto left = mid - 1, right = end - 1;

                while (right > left && left >= start) {
                    if (comp(copied[buffer_ptr], *left)) {
                        *right = *left;
                        right--;
                        left--;
                    }
                    else {
                        *right = copied[buffer_ptr];
                        right--;
                        buffer_ptr--;
                    }
                }

                while (right > left) {
                    *right = copied[buffer_ptr];
                    right--;
                    buffer_ptr--;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void pd_merge(RandomAccessIterator left_start, RandomAccessIterator right_start, RandomAccessIterator end, 
                          std::vector<T>& copied, Compare comp) {
                if (end - right_start < right_start - left_start)
                    pd_merge_down<T>(left_start, right_start, end, copied, comp);
                else
                    pd_merge_up<T>(left_start, right_start, end, copied, comp);
            }

            template<class RandomAccessIterator, class Compare>
            bool compare(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                return !comp(*b, *a);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator identify_run(RandomAccessIterator index, RandomAccessIterator max_index, Compare comp) {
                auto start_index = index;
                if (index >= max_index) return max_index + 1;

                bool cmp = compare(index, index + 1, comp);
                index++;

                while (index < max_index) {
                    bool check_cmp = compare(index, index + 1, comp);
                    if (check_cmp != cmp) break;
                    index++;
                }

                if (!cmp) reversal(start_index, index + 1);
                if (index >= max_index) return max_index + 1;
                return index + 1;
            }

            template<class RandomAccessIterator, class Compare>
            std::vector<RandomAccessIterator> find_runs(RandomAccessIterator start, RandomAccessIterator max_index, int &run_count, Compare comp) {
                std::vector<RandomAccessIterator> runs((max_index - start) / 2 + 2);
                run_count = 0;

                auto last_run = start;
                while (last_run != max_index + 1) {
                    runs[run_count++] = last_run;
                    auto new_run = identify_run(last_run, max_index, comp);
                    last_run = new_run;
                }
                return runs;
            }

            template<class RandomAccessIterator, class Compare>
            void pd_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                
                int length = last - first, run_count;
                std::vector<T> copied(length);

                auto runs = find_runs(first, last, run_count, comp);
                while (run_count > 1) {
                    for (int i = 0; i < run_count - 1; i += 2) {
                        auto end = i + 2 >= run_count ? last : runs[i + 2];
                        pd_merge(runs[i], runs[i + 1], end, copied, comp);
                    }
                    for (int i = 1, j = 2; i < run_count; i++, j += 2, run_count--)
                        runs[i] = runs[j];
                }
            }
        }

        

        namespace inplace_merge_sorting {
            template<class RandomAccessIterator>
            void improved_push(RandomAccessIterator p, RandomAccessIterator a, RandomAccessIterator b) {
                if (a == b) return;

                auto temp = *p;
                *p = *a;

                for (auto i = a + 1; i < b; i++) *(i - 1) = *i;
                *(b - 1) = temp;
            }

            template<class RandomAccessIterator, class Compare>
            void _improved_inplace_merge(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, Compare comp) {
                auto i = a, j = m;

                while (i < m && j < b) {
                    if (comp(*j, *i)) j++;
                    else improved_push(i++, m, j);
                }

                while (i < m) improved_push(i++, m, b);
            }

            template<class RandomAccessIterator, class Compare>
            void improved_inplace_merge_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto m = a + (b - a) / 2;

                if (b - a > 2) {
                    if (b - a > 3) improved_inplace_merge_sort(a, m, comp);
                    improved_inplace_merge_sort(m, b, comp);
                }
                _improved_inplace_merge(a, m, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void push(RandomAccessIterator low, RandomAccessIterator high, Compare comp) {
                for (auto i = low; i < high; i++)
                    if (comp(*(i + 1), *i)) std::iter_swap(i, i + 1);
            }

            template<class RandomAccessIterator, class Compare>
            void _inplace_merge(RandomAccessIterator min, RandomAccessIterator max, RandomAccessIterator mid, Compare comp) {
                auto i = min;
                while (i <= mid) {
                    if (comp(*(mid + 1), *i)) {
                        std::iter_swap(i, mid + 1);
                        push(mid + 1, max, comp);
                    }
                    i++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _inplace_merge_sort(RandomAccessIterator min, RandomAccessIterator max, Compare comp) {
                if (min == max) return;
                if (max - min == 1) {
                    if (comp(*max, *min)) std::iter_swap(min, max);
                    return;
                }

                auto mid = min + (max - min) / 2;
                _inplace_merge_sort(min, mid, comp);
                _inplace_merge_sort(mid + 1, max, comp);
                _inplace_merge(min, max, mid, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _inplace_merge_sort(first, last - 1, comp);
            }
        }

        namespace andrey_sorting {
            template<class RandomAccessIterator, class Compare>
            void insert_sort(RandomAccessIterator a, int b, Compare comp) {
                while (b > 1) {
                    int k = 0;
                    for (int i = 1; i < b; i++)
                        if (comp(*(a + i), *(a + k))) k = i;
                    std::iter_swap(a, a + k);
                    a++;
                    b--;
                }
            }

            template<class RandomAccessIterator>
            void aswap(RandomAccessIterator arr1, RandomAccessIterator arr2, int l) {
                while (l-- > 0) std::iter_swap(arr1++, arr2++);
            }

            template<class RandomAccessIterator, class Compare>
            int backward_merge(RandomAccessIterator arr1, int l1, RandomAccessIterator arr2, int l2, Compare comp) {
                auto arr0 = arr2 + l1;
                for (;;) {
                    if (comp(*arr2, *arr1)) {
                        std::iter_swap(arr1--, arr0--);
                        if (--l1 == 0) return 0;
                    }
                    else {
                        std::iter_swap(arr2--, arr0--);
                        if (--l2 == 0) break;
                    }
                }
                int res = l1;
                do std::iter_swap(arr1--, arr0--);
                while (--l1 != 0);
                return res;
            }

            template<class RandomAccessIterator, class Compare>
            void rotate_merge(RandomAccessIterator a, int l, int r, Compare comp) {
                for (int i = 0; i < l; i += r) {
                    int q = i;
                    for (int j = i + r; j < l; j += r)
                        if (comp(*(a + j), *(a + q))) q = j;
                    
                    if (q != i) aswap(a + i, a + q, r);

                    if (i != 0) {
                        aswap(a + l, a + i, r);
                        backward_merge(a + (l + r - 1), r, a + (i - 1), r, comp);
                    }
                }
            }

            int rbnd(int len) {
                len = len / 2;
                int k = 0;
                for (int i = 1; i < len; i *= 2) {
                    k++;
                }
                len /= k;
                for (k = 1; k <= len; k *= 2);
                return k;
            }

            template<class RandomAccessIterator, class Compare>
            void sort_rec(RandomAccessIterator a, int len, Compare comp) {
                if (len < 12) {
                    insert_sort(a, len, comp);
                    return;
                }

                int r = rbnd(len);
                int lr = (len / r - 1) * r;

                for (int p = 2; p <= lr; p += 2) {
                    if (comp(*(a + (p - 1)), *(a + (p - 2)))) std::iter_swap(a + (p - 2), a + (p - 1));
                    if ((p & 2) != 0) continue;

                    aswap(a + (p - 2), a + p, 2);

                    int m = len - p;
                    int q = 2;

                    for (;;) {
                        int q0 = 2 * q;
                        if (q0 > m || (p & q0) != 0) break;
                        backward_merge(a + (p - q - 1), q, a + (p + q - 1), q, comp);
                        q = q0;
                    }

                    backward_merge(a + (p + q - 1), q, a + (p - q - 1), q, comp);
                    int q1 = q;
                    q *= 2;

                    while ((q & p) == 0) {
                        q *= 2;
                        rotate_merge(a + (p - q), q, q1, comp);
                    }
                }

                int q1 = 0;
                for (int q = r; q < lr; q *= 2) {
                    if ((lr & q) != 0) {
                        q1 += q;
                        if (q1 != q) {
                            rotate_merge(a + (lr - q1), q1, r, comp);
                        }
                    }
                }

                int s = len - lr;
                sort_rec(a + lr, s, comp);
                aswap(a, a + lr, s);
                s += backward_merge(a + (s - 1), s, a + (lr - 1), lr - s, comp);
                sort_rec(a, s, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void andrey_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                sort_rec(first, last - first, comp);
            }
        }

        namespace rotations {
            template<class RandomAccessIterator>
            void swap_blocks_backwards(RandomAccessIterator a, RandomAccessIterator b, int len) {
                for (int i = 0; i < len; i++) std::iter_swap(a + len - i - 1, b + len - i - 1);
            }

            template<class RandomAccessIterator>
            void block_swap(RandomAccessIterator a, RandomAccessIterator b, int len) {
                for (int i = 0; i < len; i++) std::iter_swap(a + i, b + i);
            }

            template<class RandomAccessIterator>
            void shift_forwards(RandomAccessIterator start, int length) {
                auto temp = *start;
                for (int i = 0; i < length; i++) *(start + i) = *(start + i + 1);
                *(start + length) = temp;
            }

            template<class RandomAccessIterator>
            void shift_backwards(RandomAccessIterator start, int length) {
                auto temp = *(start + length);
                for (int i = length; i > 0; i--) *(start + i) = *(start + i - 1);
                *start = temp;
            }

            int map_index(int index, int n, int length) {
                return (index - n + length) % length;
            }

            template<class RandomAccessIterator, class T>
            const T& internal_swap(RandomAccessIterator a, const T& v) {
                auto old = *a;
                *a = v;
                return old;
            }

            template<class RandomAccessIterator>
            void gries_mills(RandomAccessIterator pos, int len_a, int len_b) {
                while (len_a != 0 && len_b != 0) {
                    if (len_a <= len_b) {
                        block_swap(pos, pos + len_a, len_a);
                        pos += len_a;
                        len_b -= len_a;
                    }
                    else {
                        block_swap(pos + (len_a - len_b), pos + len_a, len_b);
                        len_a -= len_b;
                    }
                }
            }

            template<class RandomAccessIterator>
            void three_reversal(RandomAccessIterator pos, int len_a, int len_b) {
                using sortcpp::_internal::reversal;
                reversal(pos, pos + len_a);
                reversal(pos + len_a, pos + len_a + len_b);
                reversal(pos, pos + len_a + len_b);
            }

            template<class RandomAccessIterator>
            void holy_gries_mills(RandomAccessIterator pos, int len_a, int len_b) {
                while (len_a > 1 && len_b > 1) {
                    while (len_a <= len_b) {
                        block_swap(pos, pos + len_a, len_a);
                        pos += len_a;
                        len_b -= len_a;
                    }
                    if (len_a <= 1 || len_b <= 1) break;

                    while (len_b <= len_a) {
                        swap_blocks_backwards(pos + len_a - len_b, pos + len_a, len_b);
                        len_a -= len_b;
                    }
                }

                if (len_a == 1) shift_forwards(pos, len_b);
                else if (len_b == 1) shift_backwards(pos, len_a);
            }

            template<class RandomAccessIterator>
            void helium(RandomAccessIterator pos, int len_a, int len_b) {
                while (len_b > 1 && len_a > 1) {
                    if (len_b < len_a) {
                        block_swap(pos, pos + len_a, len_b);
                        pos += len_b;
                        len_a -= len_b;
                    }
                    else {
                        swap_blocks_backwards(pos, pos + len_b, len_a);
                        len_b -= len_a;
                    }
                }

                if (len_b == 1) shift_backwards(pos, len_a);
                else if (len_a == 1) shift_forwards(pos, len_b);
            }

            template<class RandomAccessIterator>
            void cycle_reverse(RandomAccessIterator pos, int len_a, int len_b) {
                using sortcpp::_internal::reversal;
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                if (len_a < 1 || len_b < 1) return;

                auto a = pos;
                auto b = pos + len_a - 1;
                auto c = pos + len_a;
                auto d = pos + len_a + len_b - 1;
                T swap;

                while (a < b && c < d) {
                    swap = *b;
                    *(b--) = *a;
                    *(a++) = *c;
                    *(c++) = *d;
                    *(d--) = swap;
                }

                while (a < b) {
                    swap = *b;
                    *(b--) = *a;
                    *(a++) = *d;
                    *(d--) = swap;
                }

                while (c < d) {
                    swap = *c;
                    *(c++) = *d;
                    *(d--) = *a;
                    *(a++) = swap;
                }

                if (a < d) reversal(a, d + 1);
            }

            template<class RandomAccessIterator>
            void juggling(RandomAccessIterator pos, int len_a, int len_b) {
                int length = len_a + len_b;
                len_a %= length;

                if (len_a == 0) return;

                int cnt = 0, index = 0, start_index = 0;
                auto value = *(pos + index);

                for (; cnt < length; cnt++) {
                    int next_index = map_index(index, len_a, length);

                    value = internal_swap(pos + next_index, value);

                    if (next_index == start_index) {
                        start_index = map_index(index, 1, length);
                        value = *(pos + index);
                    }
                    else index = next_index;
                }
            }

        }

        namespace binary_quick_sorting {
            // Represents a task for binary quick sorting.
            template<class RandomAccessIterator>
            struct Task {
                RandomAccessIterator p, r;
                int bit;
                Task() {}

                // p - the beginning of the range.
                // r - the end of the range (inclusive).
                // bit - the bit to partition by.
                Task(RandomAccessIterator _p, RandomAccessIterator _r, int _bit) : p(_p), r(_r), bit(_bit) {}
            };

            // Partition the range [p, r] by the bit-th bit of the key.
            // p - the beginning of the range.
            // r - the end of the range (inclusive).
            // bit - the bit to partition by.
            // key - the function that returns the key of the element.
            template<class RandomAccessIterator, class Key>
            RandomAccessIterator binary_partition(RandomAccessIterator p, RandomAccessIterator r, int bit, Key key) {
                using sortcpp::_internal::get_bit;

                auto i = p - 1, j = r + 1;

                while (true) {
                    i++;
                    while (i <= r && !get_bit(key(*i), bit)) i++;

                    j--;
                    while (j >= p && get_bit(key(*j), bit)) j--;

                    if (i < j) std::iter_swap(i, j);
                    else return j;
                }
                return j;
            }

            
            // Sort the range [p, r] use the binary quick sort algorithm (recursive).
            // p - the beginning of the range.
            // r - the end of the range (inclusive).
            // bit - the highest bit of the key.
            // key - the function that returns the key of the element.
            template<class RandomAccessIterator, class Key>
            void binary_quick_sort_recursive(RandomAccessIterator p, RandomAccessIterator r, int bit, Key key) {
                if (p < r && bit >= 0) {
                    auto q = binary_partition(p, r, bit, key);
                    binary_quick_sort_recursive(p, q, bit - 1, key);
                    binary_quick_sort_recursive(q + 1, r, bit - 1, key);
                }
            }


            // Sort the range [p, r] use the binary quick sort algorithm (iterative with a queue).
            // p - the beginning of the range.
            // r - the end of the range (inclusive).
            // bit - the highest bit of the key.
            // key - the function that returns the key of the element.
            template<class RandomAccessIterator, class Key>
            void binary_quick_sort(RandomAccessIterator p, RandomAccessIterator r, int bit, Key key) {
                std::queue<Task<RandomAccessIterator>> tasks;
                tasks.emplace(p, r, bit);

                while (!tasks.empty()) {
                    auto task = tasks.front();
                    tasks.pop();

                    if (task.p < task.r && task.bit >= 0) {
                        auto q = binary_partition(task.p, task.r, task.bit, key);
                        tasks.emplace(task.p, q, task.bit - 1);
                        tasks.emplace(q + 1, task.r, task.bit - 1);
                    }
                }
            }

            // Get the length of the binary representation of the key of the largest element in the range [first, last).
            // first - the beginning of the range.
            // last - the end of the range (exclusive).
            // key - the function that returns the key of the element.
            template<class RandomAccessIterator, class Key>
            int get_max_bit(RandomAccessIterator first, RandomAccessIterator last, Key key) {
                using sortcpp::_internal::bit_length;
                int max = 0;
                for (auto i = first; i < last; i++) {
                    auto d = key(*i);
                    int bits = bit_length(d);
                    if (bits > max) max = bits;
                }
                return max;
            }
        }

        namespace forced_stable_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            bool stable_comp(RandomAccessIterator array, std::vector<int> &key, RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                return comp(*b, *a) || equals(*a, *b, comp) && key[a - array] > key[b - array];
            }

            template<class RandomAccessIterator>
            void stable_swap(RandomAccessIterator array, std::vector<int>& key, RandomAccessIterator a, RandomAccessIterator b) {
                std::iter_swap(a, b);
                std::swap(key[a - array], key[b - array]);
            }

            template<class RandomAccessIterator, class Compare>
            void stable_sift_down(RandomAccessIterator array, std::vector<int>& key, int root, int dist, int start, Compare comp) {
                while (root <= dist / 2) {
                    int leaf = 2 * root;
                    if (leaf < dist && stable_comp(array, key, array + start + leaf, array + start + leaf - 1, comp)) leaf++;
                
                    if (stable_comp(array, key, array + start + leaf - 1, array + start + root - 1, comp)) {
                        stable_swap(array, key, array + start + leaf - 1, array + start + root - 1);
                        root = leaf;
                    }
                    else break;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void stable_heapify(RandomAccessIterator array, std::vector<int> &key, int low, int high, Compare comp) {
                int length = high - low;
                for (int i = length / 2; i >= 1; i--) stable_sift_down(array, key, i, length, low, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _stable_heap_sort(RandomAccessIterator array, std::vector<int>& key, int start, int end, Compare comp) {
                stable_heapify(array, key, start, end, comp);
                for (int i = end - start; i > 1; i--) {
                    stable_swap(array, key, array + start, array + start + i - 1);
                    stable_sift_down(array, key, 1, i - 1, start, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void stable_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> key(length);
                std::iota(key.begin(), key.end(), 0);

                _stable_heap_sort(first, key, 0, length, comp);
            }
        }

        namespace forced_stable_quick_sorting {
            using sortcpp::_internal::forced_stable_heap_sorting::stable_comp;
            using sortcpp::_internal::forced_stable_heap_sorting::stable_swap;

            template<class RandomAccessIterator, class Compare>
            void stable_median_of_three(RandomAccessIterator array, std::vector<int> &key, RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto m = a + (b - 1 - a) / 2;

                if (stable_comp(array, key, a, m, comp)) stable_swap(array, key, a, m);

                if (stable_comp(array, key, m, b - 1, comp)) {
                    stable_swap(array, key, m, b - 1);
                    if (stable_comp(array, key, a, m, comp)) return;
                }
                stable_swap(array, key, a, m);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _stable_partition(RandomAccessIterator array, std::vector<int>& key, RandomAccessIterator a, RandomAccessIterator b, RandomAccessIterator p, Compare comp) {
                auto i = a - 1, j = b;

                while (true) {
                    do i++;
                    while (i < j && !stable_comp(array, key, i, p, comp));

                    do j--;
                    while (j >= i && stable_comp(array, key, j, p, comp));


                    if (i < j) stable_swap(array, key, i, j);
                    else return j;
                }
                return i;
            }

            template<class RandomAccessIterator, class Compare>
            void _stable_quick_sort(RandomAccessIterator array, std::vector<int>& key, RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (b - a < 3) {
                    if (b - a == 2 && stable_comp(array, key, a, a + 1, comp)) stable_swap(array, key, a, a + 1);
                    return;
                }

                stable_median_of_three(array, key, a, b, comp);
                auto p = _stable_partition(array, key, a + 1, b, a, comp);
                stable_swap(array, key, a, p);

                _stable_quick_sort(array, key, a, p, comp);
                _stable_quick_sort(array, key, p + 1, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void stable_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> key(length);
                std::iota(key.begin(), key.end(), 0);

                _stable_quick_sort(first, key, first, last, comp);
            }
        }
        
        namespace fun_sorting {
            using sortcpp::_internal::equals;

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator bin_search(RandomAccessIterator start, RandomAccessIterator end, const T& value, Compare comp) {
                while (start < end) {
                    auto mid = start + (end - start) / 2;
                    if (comp(*mid, value)) start = mid + 1;
                    else end = mid;
                }
                return start;
            }

            template<class RandomAccessIterator, class Compare>
            void fun_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first + 1; i < last; i++) {
                    bool done = false;
                    do {
                        done = true;
                        auto pos = bin_search(first, last, *i, comp);
                        if (!equals(*pos, *i, comp)) {
                            if (i < pos - 1) std::iter_swap(i, pos - 1);
                            else if (i > pos) std::iter_swap(i, pos);
                            done = false;
                        }
                    } while (!done);
                }
            }
        }

        namespace base_n_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            void base_n_sift_down(RandomAccessIterator array, int base, int node, int stop, Compare comp) {
                int left = node * base + 1;
                if (left < stop) {
                    int max_index = left;
                    for (int i = left + 1; i < left + base; i++) {
                        if (i >= stop) break;
                        if (comp(*(array + max_index), *(array + i))) max_index = i;
                    }

                    if (comp(*(array + node), *(array + max_index))) {
                        std::iter_swap(array + node, array + max_index);
                        base_n_sift_down(array, base, max_index, stop, comp);
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _base_n_heap_sort(RandomAccessIterator first, int length, int base, Compare comp) {
                for (int i = length - 1; i >= 0; i--) base_n_sift_down(first, base, i, length, comp);
                for (int i = length - 1; i > 0; i--) {
                    std::iter_swap(first, first + i);
                    base_n_sift_down(first, base, 0, i, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void base_n_heap_sort(RandomAccessIterator first, RandomAccessIterator last, int base, Compare comp) {
                _base_n_heap_sort(first, last - first, base, comp);
            }
        }

        namespace binomial_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            void _binomial_heap_sort(RandomAccessIterator array, int length, Compare comp) {
                int max_node, focus, index, depth;
                for (index = 2; index <= length; index++) {
                    max_node = index;
                    do {
                        focus = max_node;
                        for (depth = 1; (focus & depth) == 0; depth *= 2)
                            if (comp(*(array + max_node - 1), *(array + focus - depth - 1))) max_node = (focus - depth);

                        if (focus != max_node) std::iter_swap(array + focus - 1, array + max_node - 1);

                    } while (focus != max_node);
                }

                for (index = length; index > 2; index--) {
                    max_node = index;
                    focus = index;

                    for (depth = 1; focus != 0; depth *= 2)
                        if ((focus & depth) != 0) {
                            if (comp(*(array + max_node - 1), *(array + focus - 1))) max_node = focus;
                            focus -= depth;
                        }

                    if (max_node != index) {
                        focus = index;
                        do {
                            std::iter_swap(array + focus - 1, array + max_node - 1);
                            focus = max_node;
                            for (depth = 1; (focus & depth) == 0; depth *= 2)
                                if (comp(*(array + max_node - 1), *(array + focus - depth - 1))) max_node = (focus - depth);
                        } while (focus != max_node);
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void binomial_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _binomial_heap_sort(first, last - first, comp);
            }
        }

        namespace binomial_smooth_sorting {
            int height(int node) {
                int count = 0;
                while ((node >> count) % 2 == 1) count++;
                return count;
            }

            template<class RandomAccessIterator, class Compare>
            void thrift(RandomAccessIterator array, int node, bool parent, bool root, Compare comp) {
                root = root && (node >= (1 << height(node)));
                if (!root && !parent) return;
                int choice = height(node) - (root ? 0 : 1);
                if (parent)
                    for (int child = choice - 1; child >= 0; child--)
                        if (!comp(*(array + node - (1 << child)), *(array + node - (1 << choice)))) choice = child;
                if (!comp(*(array + node), *(array + node - (1 << choice)))) return;
                std::iter_swap(array + node, array + node - (1 << choice));
                thrift(array, node - (1 << choice), (node - (1 << choice)) % 2 == 1, choice == height(node), comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _binomial_smooth_sort(RandomAccessIterator first, int length, Compare comp) {
                int node;
                for (node = 1; node < length; node++)
                    thrift(first, node, node & 1, (node + (1 << height(node)) >= length), comp);
                for (node -= (node - 1) % 2; node > 2; node -= 2)
                    for (int child = height(node) - 1; child >= 0; child--)
                        thrift(first, node - (1 << child), false, true, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void binomial_smooth_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _binomial_smooth_sort(first, last - first, comp);
            }
        }

        namespace bottom_up_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            void bottom_up_sift_down(RandomAccessIterator array, int i, int b, Compare comp) {
                int j = i;
                for (; 2 * j + 1 < b; j = 2 * j + 2 < b ? (comp(*(array + 2 * j + 1), *(array + 2 * j + 2)) ? 2 * j + 2 : 2 * j + 1) : 2 * j + 1);
                for (; comp(*(array + j), *(array + i)); j = (j - 1) / 2);
                for (; j > i; j = (j - 1) / 2) std::iter_swap(array + i, array + j);
            }

            template<class RandomAccessIterator, class Compare>
            void _bottom_up_heap_sort(RandomAccessIterator first, int length, Compare comp) {
                for (int i = (length - 1) / 2; i >= 0; i--) bottom_up_sift_down(first, i, length, comp);
                for (int i = length - 1; i > 0; i--) {
                    std::iter_swap(first, first + i);
                    bottom_up_sift_down(first, 0, i, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void bottom_up_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _bottom_up_heap_sort(first, last - first, comp);
            }
        }

        namespace classic_tournament_sorting {
            int ceil_pow2(int n) {
                int r = 1;
                while (r < n) r *= 2;
                return r;
            }

            template<class RandomAccessIterator, class Compare>
            bool tree_compare(RandomAccessIterator array, std::vector<int>& tree, int a, int b, Compare comp) {
                return !comp(*(array + tree[b]), *(array + tree[a]));
            }

            template<class RandomAccessIterator, class Compare>
            void build_tree(RandomAccessIterator array, std::vector<int>& tree, int &size, int n, Compare comp) {
                size = ceil_pow2(n) - 1;
                int mod = n & 1;
                int tree_size = n + size + mod;

                tree.assign(tree_size, -1);
                for (int i = size; i < tree_size - mod; i++) tree[i] = i - size;

                for (int i, j = size, k = tree_size - mod; j > 0; j /= 2, k /= 2) {
                    for (i = j; i + 1 < k; i += 2) {
                        int val = tree_compare(array, tree, i, i + 1, comp) ? tree[i] : tree[i + 1];
                        tree[i / 2] = val;
                    }
                    if (i < k) tree[i / 2] = tree[i];
                }
            }

            template<class RandomAccessIterator>
            auto peek(RandomAccessIterator array, std::vector<int>& tree) {
                return *(array + tree[0]);
            }

            template<class RandomAccessIterator, class Compare>
            auto find_next(RandomAccessIterator array, std::vector<int>& tree, int size, Compare comp) {
                int root = tree[0] + size;

                for (int i = root; i > 0; i = (i - 1) / 2) tree[i] = -1;

                for (int i = root; i > 0; ) {
                    int j = i + ((i & 1) << 1) - 1;

                    int c1 = tree[i] >> 31;
                    int c2 = tree[j] >> 31;

                    int n_val = (c1 & ((c2 & -1) + (~c2 & tree[j]))) + (~c1 & ((c2 & tree[i]) + (~c2 & -2)));

                    if (n_val == -2) {
                        if (i < j) n_val = tree_compare(array, tree, i, j, comp) ? tree[i] : tree[j];
                        else n_val = tree_compare(array, tree, j, i, comp) ? tree[j] : tree[i];
                    }

                    i = (i - 1) / 2;
                    if (n_val != -1) tree[i] = n_val;
                }
                return peek(array, tree);
            }

            template<class T, class RandomAccessIterator, class Compare>
            void _tournament_sort(RandomAccessIterator first, int length, Compare comp) {
                std::vector<int> tree;
                std::vector<T> tmp(length);

                int size;
                build_tree(first, tree, size, length, comp);

                tmp[0] = peek(first, tree);

                for (int i = 1; i < length; i++) {
                    auto val = find_next(first, tree, size, comp);
                    tmp[i] = val;
                }

                for (auto p : tmp) {
                    *first = p;
                    first++;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void tournament_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _tournament_sort<T>(first, last - first, comp);
            }
        }

        namespace cycle_sorting {
            using sortcpp::_internal::equals;

            // Unstable

            template<class T, class RandomAccessIterator, class Compare>
            RandomAccessIterator count_lesser(RandomAccessIterator a, RandomAccessIterator b, const T& t, Compare comp) {
                auto r = a;
                for (auto i = a + 1; i < b; i++) r += comp(*i, t);
                return r;
            }

            template<class RandomAccessIterator, class Compare>
            void cycle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last - 1; i++) {
                    auto t = *i;
                    auto r = count_lesser(i, last, t, comp);

                    if (r != i) {
                        do {
                            while (equals(*r, t, comp)) r++;

                            auto t1 = *r;
                            *r = t;
                            t = t1;

                            r = count_lesser(i, last, t, comp);
                        } while (r != i);

                        *i = t;
                    }
                }
            }


            // Stable

            const int WLEN = 3;

            bool get_bit(std::vector<int> &bits, int idx) {
                int b = (bits[idx >> WLEN]) >> (idx & ((1 << WLEN) - 1)) & 1;
                return b == 1;
            }

            void flag(std::vector<int>& bits, int idx) {
                bits[idx >> WLEN] |= 1 << (idx & ((1 << WLEN) - 1));
            }

            template<class RandomAccessIterator, class Compare>
            int destination(RandomAccessIterator array, std::vector<int>& bits, int a, int b1, int b, Compare comp) {
                int d = a, e = 0;

                for (int i = a + 1; i < b; i++) {
                    if (comp(*(array + i), *(array + a))) d++;
                    else if (i < b1 && !get_bit(bits, i) && equals(*(array + i), *(array + a), comp)) e++;
                }

                while (get_bit(bits, d) || e-- > 0) d++;
                return d;
            }

            template<class RandomAccessIterator, class Compare>
            void _stable_cycle_sort(RandomAccessIterator first, std::vector<int>& bits, int length, Compare comp) {
                for (int i = 0; i < length - 1; i++)
                    if (!get_bit(bits, i)) {
                        int j = i;
                        do {
                            int k = destination(first, bits, i, j, length, comp);
                            std::iter_swap(first + i, first + k);
                            flag(bits, k);
                            j = k;
                        } while (j != i);
                    }
            }

            template<class RandomAccessIterator, class Compare>
            void stable_cycle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> bits(((length - 1) >> WLEN) + 1);

                _stable_cycle_sort(first, bits, length, comp);
            }
        }

        namespace peel_sorting {
            template<class RandomAccessIterator, class Compare>
            void cocktail_peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int stacked = 0;
                for (auto left = first; left < last; left++) {
                    stacked = 0;
                    for (auto right = last - 1; right > left; right--) {
                        if (comp(*(right + stacked), *left)) {
                            auto item = *(right + stacked);
                            for (auto pull = right + stacked; pull > left; pull--) *pull = *(pull - 1);
                            *left = item;
                            stacked++;
                        }
                    }

                    left++;
                    for (auto right = left + 1; right < last; right++) {
                        if (comp(*right, *left)) {
                            auto item = *right;
                            for (auto pull = right; pull > left; pull--) *pull = *(pull - 1);
                            *left = item;
                        }
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto left = first; left < last; left++) {
                    int stacked = 0;
                    for (auto right = last - 1; right > left; right--) {
                        if (comp(*(right + stacked), *left)) {
                            auto item = *(right + stacked);
                            for (auto pull = right + stacked; pull > left; pull--) *pull = *(pull - 1);
                            *left = item;
                            stacked++;
                        }
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void reverse_peel_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int stacked = 0;
                for (auto left = first; left < last; left++) {
                    for (auto right = left + stacked + 1; right < last; right++) {
                        if (right == left + stacked + 1) stacked = 0;
                        if (comp(*right, *left)) {
                            auto item = *right;
                            for (auto pull = right; pull > left; pull--) *pull = *(pull - 1);
                            *left = item;
                            stacked++;
                        }
                    }
                }
            }
        }

        namespace sandpaper_sorting {
            template<class RandomAccessIterator, class Compare>
            void reverse_sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last - 1; i++) {
                    for (auto j = last - 1; j > i; j--) {
                        if (comp(*j, *i)) std::iter_swap(i, j);
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void sandpaper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first; i < last - 1; i++) {
                    for (auto j = i + 1; j < last; j++) {
                        if (comp(*j, *i)) std::iter_swap(i, j);
                    }
                }
            }


        }

        namespace fall_sorting {
            template<class RandomAccessIterator, class Compare>
            void fall_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto left = first + 1, right = first + 2, highest_low = first;

                while (left <= last) {
                    right = left + 1;
                    highest_low = first;

                    while (right <= last) {
                        if (comp(*(right - 1), *(left - 1))) {
                            if (highest_low == first) highest_low = right;
                            else if (comp(*(highest_low - 1), *(right - 1))) highest_low = right;
                        }
                        right++;
                    }

                    if (highest_low == first) left++;
                    else std::iter_swap(left - 1, highest_low - 1);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void stable_fall_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto left = first + 1, right = first + 2, highest_low = first;
                int stacked = 0;

                while (left <= last) {
                    right = left + 1 + stacked;
                    highest_low = first;
                    while (right <= last) {
                        if (comp(*(right - 1), *(left - 1))) {
                            if (highest_low == first) highest_low = right;
                            else if (comp(*(highest_low - 1), *(right - 1))) highest_low = right;
                        }
                        right++;
                    }
                    if (highest_low == first) {
                        left++;
                        stacked = 0;
                    }
                    else {
                        auto item = *(highest_low - 1);
                        auto pull = highest_low;
                        while (pull > left) {
                            *(pull - 1) = *(pull - 2);
                            pull--;
                        }
                        *(left - 1) = item;
                        stacked++;
                    }
                }
            }
        }

        namespace filpped_min_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            void filpped_sift_down(RandomAccessIterator array, int length, int root, int dist, Compare comp) {
                while (root <= dist / 2) {
                    int leaf = 2 * root;
                    if (leaf < dist && comp(*(array + length - leaf - 1), *(array + length - leaf))) leaf++;
                    if (comp(*(array + length - leaf), *(array + length - root))) {
                        std::iter_swap(array + length - root, array + length - leaf);
                        root = leaf;
                    }
                    else break;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _filpped_heap_sort(RandomAccessIterator first, int length, Compare comp) {
                for (int i = length / 2; i >= 1; i--) filpped_sift_down(first, length, i, length, comp);
                for (int i = length; i > 1; i--) {
                    std::iter_swap(first + length - 1, first + length - i);
                    filpped_sift_down(first, length, 1, i - 1, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void filpped_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _filpped_heap_sort(first, last - first, comp);
            }
        }

        namespace lazy_heap_sorting {
            template<class RandomAccessIterator, class Compare>
            void max_to_front(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto max = a;
                for (auto i = a + 1; i < b; i++)
                    if (comp(*max, *i)) max = i;
                std::iter_swap(max, a);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator find_min(RandomAccessIterator p, RandomAccessIterator a, RandomAccessIterator b, int s, Compare comp) {
                auto min = p;
                for (auto i = a; i < b; i += s)
                    if (comp(*i, *min)) min = i;
                return min;
            }

            template<class RandomAccessIterator, class Compare>
            void lazyheap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int s = std::sqrt(length - 1) + 1;

                for (auto i = first; i < last; i += s) max_to_front(i, std::min(i + s, last), comp);

                for (auto j = last; j > first; ) {
                    auto max = first;

                    for (auto i = max + s; i < j; i += s)
                        if (!comp(*i, *max)) max = i;

                    std::iter_swap(max, --j);
                    max_to_front(max, std::min(max + s, j), comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_lazyheap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int s = std::sqrt(length - 1) + 1;

                auto f = first + (length - 1) % s + 1;
                auto fmin = find_min(first, first + 1, f, 1, comp);

                for (auto j = f; j < last; j += s) {
                    auto min = find_min(j, j + 1, std::min(j + s, last), 1, comp);
                    if (j != min) std::iter_swap(min, j);
                }

                for (auto j = first; j < last; ) {
                    auto min = find_min(fmin, f, last, s, comp);

                    if (min == fmin) {
                        if (j != min) std::iter_swap(min, j);
                        if (++j == f) f += s;
                        fmin = find_min(j, j + 1, f, 1, comp);
                    }
                    else {
                        if (j == fmin) fmin = find_min(j + 1, j + 2, f, 1, comp);
                        auto nmin = find_min(j, min + 1, std::min(min + s, last), 1, comp);

                        if (nmin == j) std::iter_swap(min, j);
                        else {
                            auto t = *j;
                            *j = *min;
                            *min = *nmin;
                            *nmin = t;
                        }

                        if (++j == f) f += s;
                    }
                }
            }
        }

        namespace min_max_heap_sorting {
            using sortcpp::_internal::bit_length;

            template<class RandomAccessIterator, class Compare>
            bool compare(RandomAccessIterator a, RandomAccessIterator b, bool is_gt, Compare comp) {
                return (is_gt ? comp(*b, *a) : comp(*a, *b));
            }

            bool is_min_level(int start, int index) {
                index = index - start + 1;
                return bit_length(index) & 1;
            }

            template<class RandomAccessIterator, class Compare>
            void down_heap(RandomAccessIterator a, int start, int end, int i, Compare comp) {
                bool cf = !is_min_level(start, i);
                int left = 2 * i + 1;
                while (left < end) {
                    int right = left + 1, nexti = left;
                    for (int c : {right, 2 * left + 1, 2 * left + 2, 2 * right + 1, 2 * right + 2}) {
                        if (c >= end) break;
                        if (compare(a + c, a + nexti, cf, comp)) nexti = c;
                    }

                    if (nexti <= right) {
                        if (compare(a + nexti, a + i, cf, comp)) std::iter_swap(a + nexti, a + i);
                        return;
                    }
                    else {
                        if (compare(a + nexti, a + i, cf, comp)) {
                            std::iter_swap(a + nexti, a + i);
                            int parent = (nexti - 1) / 2;
                            if (compare(a + parent, a + nexti, cf, comp)) std::iter_swap(a + nexti, a + parent);
                        }
                        else return;
                    }

                    i = nexti;
                    left = 2 * i + 1;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void min_max_heapify(RandomAccessIterator a, int start, int end, Compare comp) {
                for (int i = (end - 1) / 2; i >= start; i--) down_heap(a, start, end, i, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void store_max(RandomAccessIterator a, int start, int& end, Compare comp) {
                if (end <= start + 1) return;
                int imax = start + 1;
                if (end > imax + 1 && compare(a + imax, a + imax + 1, false, comp)) imax++;
                end--;
                std::iter_swap(a + imax, a + end);
                if (imax < end) down_heap(a, start, end, imax, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _min_max_heap_sort(RandomAccessIterator a, int start, int end, Compare comp) {
                min_max_heapify(a, start, end, comp);
                for (int i = end - 1; i > start; i--) store_max(a, start, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void min_max_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _min_max_heap_sort(first, 0, last - first, comp);
            }
        }

        namespace poplar_heap_sorting {
            int hyperfloor(int n) {
                return (int)std::pow(2, std::floor(std::log2(n)));
            }

            template<class RandomAccessIterator, class Compare>
            void unchecked_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto cur = first + 1; cur < last; cur++) {
                    auto sift = cur;
                    auto sift_1 = cur - 1;

                    if (comp(*sift, *sift_1)) {
                        auto tmp = *sift;
                        do {
                            *sift = *sift_1;
                        } while (--sift != first && comp(tmp, *--sift_1));
                        *sift = tmp;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                if (first == last) return;
                unchecked_insertion_sort(first, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void sift(RandomAccessIterator first, int size, Compare comp) {
                if (size < 2) return;

                auto root = first + (size - 1);
                auto child_root1 = root - 1;
                auto child_root2 = first + (size / 2 - 1);

                while (true) {
                    auto max_root = root;
                    if (comp(*max_root, *child_root1)) max_root = child_root1;
                    if (comp(*max_root, *child_root2)) max_root = child_root2;
                    if (max_root == root) return;

                    std::iter_swap(root, max_root);

                    size /= 2;
                    if (size < 2) return;

                    root = max_root;
                    child_root1 = root - 1;
                    child_root2 = max_root - (size - size / 2);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pop_heap_with_size(RandomAccessIterator first, RandomAccessIterator last, int size, Compare comp) {
                int poplar_size = hyperfloor(size + 1) - 1;
                auto last_root = last - 1;
                auto bigger = last_root;
                int bigger_size = poplar_size;

                auto it = first;
                while (true) {
                    auto root = it + poplar_size - 1;
                    if (root == last_root) break;
                    if (comp(*bigger, *root)) {
                        bigger = root;
                        bigger_size = poplar_size;
                    }
                    it = root + 1;
                    size -= poplar_size;
                    poplar_size = hyperfloor(size + 1) - 1;
                }

                if (bigger != last_root) {
                    std::iter_swap(bigger, last_root);
                    sift(bigger - (bigger_size - 1), bigger_size, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void make_poplar_heap(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int size = last - first;
                if (size < 2) return;

                int small_poplar_size = 15;
                if (size <= small_poplar_size) {
                    unchecked_insertion_sort(first, last, comp);
                    return;
                }

                int poplar_level = 1;

                auto it = first;
                auto next = it + small_poplar_size;
                while (true) {
                    unchecked_insertion_sort(it, next, comp);
                    int poplar_size = small_poplar_size;

                    for (int i = (poplar_level & -poplar_level) >> 1; i > 0; i >>= 1) {
                        it -= poplar_size;
                        poplar_size = 2 * poplar_size + 1;
                        sift(it, poplar_size, comp);
                        next++;
                    }

                    if ((last - next) <= small_poplar_size) {
                        unchecked_insertion_sort(next, last, comp);
                        return;
                    }

                    it = next;
                    next += small_poplar_size;
                    poplar_level++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void sort_poplar_heap(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int size = last - first;
                if (size < 2) return;

                do {
                    pop_heap_with_size(first, last, size, comp);
                    last--;
                    size--;
                } while (size > 1);
            }
        }

        namespace smooth_sorting {
            using sortcpp::_internal::ctz;

            const int LP[] = { 1, 1, 3, 5, 9, 15, 25, 41, 67, 109, 177, 287,
                465, 753, 1219, 1973, 3193, 5167, 8361, 13529, 21891, 35421, 57313,
                92735, 150049, 242785, 392835, 635621, 1028457, 1664079, 2692537,
                4356617, 7049155, 11405773, 18454929, 29860703, 48315633, 78176337,
                126491971, 204668309, 331160281, 535828591, 866988873
            };

            template<class RandomAccessIterator, class Compare>
            void sift(RandomAccessIterator array, int pshift, int head, Compare comp) {
                auto val = *(array + head);

                while (pshift > 1) {
                    int rt = head - 1;
                    int lf = head - 1 - LP[pshift - 2];

                    if ((!comp(val, *(array + lf))) && (!comp(val, *(array + rt)))) break;
                    if (!comp(*(array + lf), *(array + rt))) {
                        *(array + head) = *(array + lf);
                        head = lf;
                        pshift -= 1;
                    }
                    else {
                        *(array + head) = *(array + rt);
                        head = rt;
                        pshift -= 2;
                    }
                }
                *(array + head) = val;
            }


            template<class RandomAccessIterator, class Compare>
            void trinkle(RandomAccessIterator array, int p, int pshift, int head, bool is_trusty, Compare comp) {
                auto val = *(array + head);

                while (p != 1) {
                    int stepson = head - LP[pshift];
                    if (!comp(val, *(array + stepson))) break;

                    if (!is_trusty && pshift > 1) {
                        int rt = head - 1;
                        int lf = head - 1 - LP[pshift - 2];

                        if ((!comp(*(array + rt), *(array + stepson))) || (!comp(*(array + lf), *(array + stepson))))
                            break;
                    }

                    *(array + head) = *(array + stepson);
                    head = stepson;

                    int trail = ctz(p & ~1);
                    p >>= trail;
                    pshift += trail;
                    is_trusty = false;
                }

                if (!is_trusty) {
                    *(array + head) = val;
                    sift(array, pshift, head, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _smooth_sort(RandomAccessIterator array, int lo, int hi, bool full, Compare comp) {
                int head = lo;
                int p = 1;
                int pshift = 1;

                while (head < hi) {
                    if ((p & 3) == 3) {
                        sift(array, pshift, head, comp);
                        p >>= 2;
                        pshift += 2;
                    }
                    else {
                        if (LP[pshift - 1] >= hi - head) trinkle(array, p, pshift, head, false, comp);
                        else sift(array, pshift, head, comp);

                        if (pshift == 1) {
                            p <<= 1;
                            pshift--;
                        }
                        else {
                            p <<= (pshift - 1);
                            pshift = 1;
                        }
                    }
                    p |= 1;
                    head++;
                }

                if (full) {
                    trinkle(array, p, pshift, head, false, comp);
                    while (pshift != 1 || p != 1) {
                        if (pshift <= 1) {
                            int trail = ctz(p & ~1);
                            p >>= trail;
                            pshift += trail;
                        }
                        else {
                            p <<= 2;
                            p ^= 7;
                            pshift -= 2;

                            trinkle(array, p >> 1, pshift + 1, head - LP[pshift] - 1, true, comp);
                            trinkle(array, p, pshift, head - 1, true, comp);
                        }
                        head--;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void smooth_heapify(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _smooth_sort(first, 0, last - first - 1, false, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void smooth_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _smooth_sort(first, 0, last - first - 1, true, comp);
            }
        }

        namespace ternary_heap_sorting {
            int left_branch(int i) { return 3 * i + 1; };
            int middle_branch(int i) { return 3 * i + 2; };
            int right_branch(int i) { return 3 * i + 3; };

            template<class RandomAccessIterator, class Compare>
            void ternary_max_heapify(RandomAccessIterator array, int i, int heap_size, Compare comp) {
                int left_child = left_branch(i);
                int middle_child = middle_branch(i);
                int right_child = right_branch(i);
                int largest;

                largest = left_child <= heap_size && comp(*(array + i), *(array + left_child)) ? left_child : i;
                if (middle_child <= heap_size && comp(*(array + largest), *(array + middle_child))) largest = middle_child;
                if (right_child <= heap_size && comp(*(array + largest), *(array + right_child))) largest = right_child;

                if (largest != i) {
                    std::iter_swap(array + i, array + largest);
                    ternary_max_heapify(array, largest, heap_size, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            int build_max_ternary_heap(RandomAccessIterator array, int length,  Compare comp) {
                int heap_size = length - 1;
                for (int i = (length - 1) / 3; i >= 0; i--) ternary_max_heapify(array, i, heap_size, comp);
                return heap_size;
            }

            template<class RandomAccessIterator, class Compare>
            void _ternary_heap_sort(RandomAccessIterator first, int length, Compare comp) {
                int heap_size = build_max_ternary_heap(first, length, comp);
                for (int i = length - 1; i >= 0; i--) {
                    std::iter_swap(first, first + i);
                    heap_size--;
                    ternary_max_heapify(first, 0, heap_size, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void ternary_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _ternary_heap_sort(first, last - first, comp);
            }
        }

        namespace tournament_sorting {
            template<class RandomAccessIterator, class Compare>
            bool compare(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                return comp(*a, *b);
            }
            
            bool is_player(int i) { return i <= 0; }
            void set_winner(std::vector<int>& matches, int root, int winner) { matches[root] = winner; }
            void set_winners(std::vector<int>& matches, int root, int winners) { matches[root + 1] = winners; }
            void set_losers(std::vector<int>& matches, int root, int losers) { matches[root + 2] = losers; }
            int get_winner(std::vector<int>& matches, int root) { return matches[root]; }
            int get_winners(std::vector<int>& matches, int root) { return matches[root + 1]; }
            int get_losers(std::vector<int>& matches, int root) { return matches[root + 2]; }
            void set_match(std::vector<int>& matches, int root, int winner, int winners, int losers) {
                set_winner(matches, root, winner);
                set_winners(matches, root, winners);
                set_losers(matches, root, losers);
            }

            int get_player(std::vector<int>& matches, int i) { return i <= 0 ? std::abs(i) : get_winner(matches, i); }
            int make_player(int i) { return -i; }

            template<class RandomAccessIterator, class Compare>
            int make_match(RandomAccessIterator array, std::vector<int> &matches, int top, int bot, int root, Compare comp) {
                int top_w = get_player(matches, top);
                int bot_w = get_player(matches, bot);

                if (!compare(array + bot_w, array + top_w, comp)) set_match(matches, root, top_w, top, bot);
                else set_match(matches, root, bot_w, bot, top);
                return root;
            }

            template<class RandomAccessIterator, class Compare>
            int knockout(RandomAccessIterator array, std::vector<int>& matches, int i, int k, int root, Compare comp) {
                if (i == k) return make_player(i);

                int j = (i + k) / 2;
                return make_match(
                    array, matches, 
                    knockout(array, matches, i, j, 2 * root, comp),
                    knockout(array, matches, j + 1, k, 2 * root + 3, comp),
                    root, comp
                );
            }

            template<class RandomAccessIterator, class Compare>
            int rebuild(RandomAccessIterator array, std::vector<int>& matches, int root, Compare comp) {
                if (is_player(get_winners(matches, root))) return get_losers(matches, root);
                set_winners(matches, root, rebuild(array, matches, get_winners(matches, root), comp));
                if (compare(
                    array + get_player(matches, get_losers(matches, root)),
                    array + get_player(matches, get_winners(matches, root)), comp)) {

                    set_winner(matches, root, get_player(matches, get_losers(matches, root)));
                    int tmp = get_losers(matches, root);
                    set_losers(matches, root, get_winners(matches, root));
                    set_winners(matches, root, tmp);
                }
                else set_winner(matches, root, get_player(matches, get_winners(matches, root)));
                return root;
            }

            template<class RandomAccessIterator, class Compare>
            auto tournament_pop(RandomAccessIterator array, std::vector<int>& matches, int &tourney, Compare comp) {
                auto result = *(array + get_player(matches, tourney));
                tourney = is_player(tourney) ? 0 : rebuild(array, matches, tourney, comp);
                return result;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void _tournament_sort(RandomAccessIterator array, std::vector<int> &matches, int length, int &tourney, Compare comp) {
                std::vector<T> copy(length);
                for (int i = 0; i < length; i++) copy[i] = tournament_pop(array, matches, tourney, comp);
                for (auto& x : copy) {
                    *array = x;
                    array++;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void tournament_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> matches(6 * length);
                int tourney = knockout(first, matches, 0, length - 1, 3, comp);
                _tournament_sort<T>(first, matches, length, tourney, comp);
            }
        }

        namespace triangular_heap_sorting {
            int triangular_root(int val) {
                return ((int)std::sqrt((8 * val + 1)) - 1) / 2;
            }

            template<class RandomAccessIterator, class Compare>
            void triangular_sift_down(RandomAccessIterator first, RandomAccessIterator last, int root, Compare comp) {
                auto temp = *(first + root);
                int len = triangular_root(root);
                auto left = first + root + len + 1;
                auto right = left + 1;

                while (left < last) {
                    if (right >= last) {
                        if (comp(temp, *left)) *(first + root) = *left;
                        break;
                    }

                    auto max = (!comp(*left, *right)) ? left : right;

                    if (comp(temp, *max)) {
                        *(first + root) = *max;
                        root = max - first;
                        len = triangular_root(root);
                        left = first + root + len + 1;
                        right = left + 1;
                        continue;
                    }
                    break;
                }
                *(first + root) = temp;
            }

            template<class RandomAccessIterator, class Compare>
            void triangular_heapify(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                for (int i = length - 1; i >= 0; i--) triangular_sift_down(first, last, i, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void triangular_heap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                triangular_heapify(first, last, comp);

                for (int i = 1; i < length - 1; i++) {
                    std::iter_swap(first, last - i);
                    triangular_sift_down(first, last - i, 0, comp);
                }
                if (comp(*(first + 1), *first)) std::iter_swap(first, first + 1);
            }
        }

        namespace weakheap_sorting {
            int get_bitwise_flag(std::vector<int>& bits, int x) {
                return ((bits[x >> 3] >> (x & 7)) & 1);
            }

            void toggle_bitwise_flag(std::vector<int> &bits, int x) {
                int flag = bits[x >> 3];
                flag ^= (1 << (x & 7));
                bits[x >> 3] = flag;
            }

            template<class RandomAccessIterator, class Compare>
            void weakheap_merge(RandomAccessIterator first, RandomAccessIterator i, RandomAccessIterator j, std::vector<int> &bits, Compare comp) {
                if (comp(*i, *j)) {
                    toggle_bitwise_flag(bits, j - first);
                    std::iter_swap(i, j);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _weakheap_sort(RandomAccessIterator first, int n, std::vector<int>& bits, Compare comp) {
                int i, j, x, y, g_parent;

                for (i = n - 1; i > 0; i--) {
                    j = i;
                    while ((j & 1) == get_bitwise_flag(bits, j >> 1)) j >>= 1;
                    g_parent = j >> 1;
                    weakheap_merge(first, first + g_parent, first + i, bits, comp);
                }

                for (i = n - 1; i >= 2; i--) {
                    std::iter_swap(first, first + i);

                    x = 1;
                    while ((y = 2 * x + get_bitwise_flag(bits, x)) < i) x = y;
                    while (x > 0) {
                        weakheap_merge(first, first, first + x, bits, comp);
                        x >>= 1;
                    }
                }
                std::iter_swap(first, first + 1);
            }

            template<class RandomAccessIterator, class Compare>
            void weakheap_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int n = last - first;
                int bits_length = (n + 7) / 8;
                std::vector<int> bits(bits_length);

                _weakheap_sort(first, n, bits, comp);
            }
        }

        namespace pancake_sorting {
            using sortcpp::_internal::reversal;
            using sortcpp::_internal::bogo_sorting::is_range_sorted;
            using sortcpp::_internal::merge_sorting::_binary_search;

            template<class RandomAccessIterator, class Compare>
            void brunt_pancake_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = last - 1; i > first; i--) {
                    auto max = first;

                    for (auto j = max + 1; j <= i; j++)
                        if (!comp(*j, *max)) max = j;

                    if (max != i) {
                        reversal(first, max + 1);
                        reversal(first, i + 1);
                        reversal(first, i);
                        reversal(first, max);
                    }
                }
            }

            template<class T, class Compare>
            bool le(const T& a, const T& b, Compare comp) { return !comp(b, a); }

            template<class T, class Compare>
            bool gt(const T& a, const T& b, Compare comp) { return comp(b, a); }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator monobound_fw(RandomAccessIterator start, RandomAccessIterator end, const T& value, Compare comp) {
                int top, mid;
                top = end - start;

                while (top > 1) {
                    mid = top / 2;
                    if (le(value, *(end - mid), comp)) end -= mid;
                    top -= mid;
                }

                if (le(value, *(end - 1), comp)) return end - 1;
                return end;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator monobound_bw(RandomAccessIterator start, RandomAccessIterator end, const T& value, Compare comp) {
                int top, mid;
                top = end - start;

                while (top > 1) {
                    mid = top / 2;
                    if (gt(*(start + mid), value, comp)) start += mid;
                    top -= mid;
                }

                if (gt(*start, value, comp)) return start + 1;
                return start;
            }

            template<class RandomAccessIterator, class Compare>
            bool front(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                int length = end - start;
                bool dir = true;

                if (gt(*start, *(start + 1), comp)) reversal(start, start + 2);

                if (length > 2) {
                    if (gt(*(start + 1), *(start + 2), comp)) {
                        if (gt(*start, *(start + 2), comp)) {
                            reversal(start, start + 2);
                            return false;
                        }
                        else {
                            reversal(start, start + 3);
                            reversal(start, start + 2);
                        }
                        return false;
                    }
                    else return true;
                }
                return dir;
            }

            template<class RandomAccessIterator, class Compare>
            void pancake_insertion_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bool dir = front(first, last, comp);
                for (auto i = first + 3; i < last; i++) {
                    if (dir) {
                        if (le(*(i - 1), *i, comp)) continue;
                        else if (gt(*first, *i, comp)) {
                            reversal(first, i);
                            dir = !dir;
                        }
                        else {
                            auto idx = monobound_fw(first, i, *i, comp);
                            reversal(first, i + 1);
                            auto end = i - idx + first;
                            reversal(first, end + 1);
                            reversal(first, end);
                            dir = !dir;
                        }
                    }
                    else {
                        if (gt(*(i - 1), *i, comp)) continue;
                        else if (le(*first, *i, comp)) {
                            reversal(first, i);
                            dir = !dir;
                        }
                        else {
                            auto idx = monobound_bw(first, i, *i, comp);
                            reversal(first, i + 1);
                            auto end = i - idx + first;
                            reversal(first, end + 1);
                            reversal(first, end);
                            dir = !dir;
                        }
                    }
                }

                if (!dir) reversal(first, last);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator find_max(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto index = start;
                auto max = *start;

                for (auto i = start + 1; i < end; i++)
                    if (comp(max, *i)) {
                        max = *i;
                        index = i;
                    }
                return index;
            }

            template<class RandomAccessIterator, class Compare>
            void pancake_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = last - 1; i >= first; i--) {
                    if (!is_range_sorted(first, i + 1, comp)) {
                        auto index = find_max(first, i + 1, comp);
                        if (index == first) reversal(first, i + 1);
                        else if (index != i) {
                            reversal(first, index + 1);
                            reversal(first, i + 1);
                        }
                    }
                    else break;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void cupcake_wrapper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto right = last;
                for (auto limit = last - 1; limit >= first; limit--) {
                    for (right = limit; right > first; right--)
                        if (comp(*right, *first)) reversal(first, right + 1);
                    reversal(first, last);
                    reversal(first, last - 1);
                }
            }


            using sortcpp::_internal::bogo_sorting::is_max_sorted;
            using sortcpp::_internal::bogo_sorting::randiter;

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void pancake_bogo_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& gen, Compare comp) {
                for (auto i = last - 1; i >= first; i--) {
                    while (!is_max_sorted(first, i + 1, comp)) reversal(first, randiter(first, i + 1, gen) + 1);
                }
            }


            template<class RandomAccessIterator>
            void cursed_rotate(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) {
                reversal(first, a);
                reversal(first, m);
                reversal(first, b);
                reversal(first, b - m + a);
            }

            template<class RandomAccessIterator, class Compare>
            void pancake_merge(RandomAccessIterator first, RandomAccessIterator m, RandomAccessIterator b, Compare comp) {
                RandomAccessIterator m1, m2, m3;

                if (m >= b - m + first) {
                    m1 = first + (m - first) / 2;
                    m2 = _binary_search(m, b, *m1, false, comp);
                    m3 = m1 + (m2 - m);
                }
                else {
                    m2 = m + (b - m) / 2;
                    m1 = _binary_search(first, m, *m2, true, comp);
                    m3 = (m2++) - (m - m1);
                }
                cursed_rotate(first, m1, m, m2);

                if (m1 > first && m3 > m1) pancake_merge(first, m1, m3, comp);

                m3++;
                if (m2 > m3 && b > m2) {
                    cursed_rotate(first, first, m3, b);
                    pancake_merge(first, m2 - m3 + first, b - m3 + first, comp);
                    cursed_rotate(first, first, b - m3 + first, b);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pancake_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int n = last - first;
                if (n > 1) {
                    if (!comp(*(first + 1), *first)) {
                        auto i = first + 2;
                        for (; i < last && !comp(*i, *(i - 1)); i++);
                        if (i == last) return;
                    }
                    else {
                        auto i = first + 2;
                        for (; i < last && comp(*i, *(i - 1)); i++);

                        if (i == last) {
                            reversal(first, last);
                            return;
                        }
                    }

                    int m = n / 2;

                    pancake_merge_sort(first, first + m, comp);
                    cursed_rotate(first, first, first + m, first + n);
                    m = n - m;

                    pancake_merge_sort(first, first + m, comp);
                    pancake_merge(first, first + m, first + n, comp);
                }
            }
        }

        namespace multi_selection_sorting {
            using sortcpp::_internal::pancake_sorting::monobound_fw;
            using sortcpp::_internal::pancake_sorting::monobound_bw;

            template<class RandomAccessIterator>
            void move_front(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) {
                auto stop = std::max(mid, start + (end - start) / 2);
                while (end > stop) {
                    std::iter_swap(start, end);
                    start++;
                    end--;
                }
            }

            template<class RandomAccessIterator>
            void move_back(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) {
                while (mid > start) {
                    std::iter_swap(mid, end);
                    mid--;
                    end--;
                }
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator select_smallest(RandomAccessIterator left, RandomAccessIterator run_end, RandomAccessIterator right, Compare comp) {
                auto run_start = left;
                auto i = run_end + 1;

                while (i < right) {
                    if (!comp(*run_end, *i)) {
                        if (run_end < i - 1) run_start = i;
                        run_end = i;
                    }
                    else if (run_end != run_start && comp(*i, *run_start)) run_start = monobound_bw(run_start, run_end + 1, *i, comp);
                    i++;
                }

                move_front(left, run_start - 1, run_end);
                return left + (run_end - run_start + 1);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator select_largest(RandomAccessIterator left, RandomAccessIterator run_end, RandomAccessIterator right, Compare comp) {
                auto run_start = left;
                auto i = run_end + 1;

                while (i < right) {
                    if (!comp(*i, *run_end)) {
                        if (run_end < i - 1) run_start = i;
                        run_end = i;
                    }
                    else if (run_end != run_start && comp(*run_start, *i)) run_start = monobound_fw(run_start, run_end + 1, *i, comp);
                    i++;
                }

                if (run_end != right - 1) move_back(run_start - 1, run_end, right - 1);
                return right - (run_end - run_start + 1);
            }

            template<class RandomAccessIterator, class Compare>
            void multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto left = first, right = last;
                while (left < right - 1) {
                    if (comp(*(left + 1), *left)) left = select_smallest(left, left + 1, right, comp);
                    else right = select_largest(left, left + 1, right, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void improved_multi_selection_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                auto left = first, right = last;
                bool dir = comp(*(left + 1), *left);
                if (comp(*(left + 1), *left)) left = select_smallest(left, left + 1, right, comp);
                else right = select_largest(left, left + 1, right, comp);

                bool same_dir = dir;
                int same_dir_cnt = 1, same_dir_max = 1;
                while (same_dir_max * same_dir_max < length) same_dir_max *= 2;

                while (left < right - 1) {
                    auto run_end = left + 1;
                    dir = comp(*(left + 1), *left);
                    if (dir == same_dir) {
                        same_dir_cnt++;
                        if (same_dir_cnt == same_dir_max) {
                            same_dir_cnt = 0;
                            dir = !dir;
                            run_end--;
                        }
                    }
                    else same_dir_cnt = 0;

                    if (dir) left = select_smallest(left, run_end, right, comp);
                    else right = select_largest(left, run_end, right, comp);
                }
            }
        }

        namespace quick_sorting {
            using sortcpp::_internal::insertion_sorting::insertion_sort;

            template<class RandomAccessIterator, class Compare>
            void _dual_pivot_qsort_loop(RandomAccessIterator left, RandomAccessIterator right, int divisor, Compare comp) {
                int length = right - left;

                if (length < 4) {
                    insertion_sort(left, right + 1, comp);
                    return;
                }

                int third = length / divisor;
                auto med1 = left + third;
                auto med2 = right - third;

                if (med1 <= left) med1 = left + 1;
                if (med2 >= right) med2 = right - 1;

                if (comp(*med1, *med2)) {
                    std::iter_swap(med1, left);
                    std::iter_swap(med2, right);
                }
                else {
                    std::iter_swap(med1, right);
                    std::iter_swap(med2, left);
                }

                auto pivot1 = *left, pivot2 = *right;
                auto less = left + 1, great = right - 1;

                for (auto k = less; k <= great; k++) {
                    if (comp(*k, pivot1)) std::iter_swap(k, less++);
                    else if (comp(pivot2, *k)) {
                        while (k < great && comp(pivot2, *great)) great--;
                        std::iter_swap(k, great--);
                        if (comp(*k, pivot1)) std::iter_swap(k, less++);
                    }
                }

                int dist = great - less;
                if (dist < 13) divisor++;

                std::iter_swap(less - 1, left);
                std::iter_swap(great + 1, right);

                _dual_pivot_qsort_loop(left, less - 2, divisor, comp);
                if (pivot1 < pivot2) _dual_pivot_qsort_loop(less, great, divisor, comp);
                _dual_pivot_qsort_loop(great + 2, right, divisor, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void dual_pivot_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _dual_pivot_qsort_loop(first, last - 1, 3, comp);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _partition(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                auto pivot = *hi;
                auto i = lo;

                for (auto j = lo; j < hi; j++) {
                    if (comp(*j, pivot)) {
                        std::iter_swap(i, j);
                        i++;
                    }
                }

                std::iter_swap(i, hi);
                return i;
            }

            template<class RandomAccessIterator, class Compare>
            void quick_ll_loop(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                if (lo < hi) {
                    auto p = _partition(lo, hi, comp);
                    quick_ll_loop(lo, p - 1, comp);
                    quick_ll_loop(p + 1, hi, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void quick_lr_loop(RandomAccessIterator p, RandomAccessIterator r, Compare comp) {
                auto pivot = p + (r - p + 1) / 2;
                auto x = *pivot;
                auto i = p, j = r;

                while (i <= j) {
                    while (comp(*i, x)) i++;
                    while (comp(x, *j)) j--;
                    if (i <= j) {
                        std::iter_swap(i, j);
                        i++;
                        j--;
                    }
                }

                if (p < j) quick_lr_loop(p, j, comp);
                if (i < r) quick_lr_loop(i, r, comp);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _partition_middle(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = a, j = i, m = a + (b - a) / 2;

                while (j < m) {
                    if (!comp(*m, *j)) std::iter_swap(i++, j);
                    j++;
                }

                std::iter_swap(i, m);
                j = m + 1;
                m = i++;

                while (j < b) {
                    if (comp(*j, *m)) std::iter_swap(i++, j);
                    j++;
                }

                std::iter_swap(--i, m);
                return i;
            }

            template<class RandomAccessIterator, class Compare>
            void median_quick_ll_loop(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (b - a > 1) {
                    auto p = _partition_middle(a, b, comp);
                    median_quick_ll_loop(a, p, comp);
                    median_quick_ll_loop(p + 1, b, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void cube_root_quick_sort(RandomAccessIterator start, RandomAccessIterator stop, Compare comp) {
                int len = stop - start;
                if (len >= 2) {
                    int root = (int)std::cbrt(len);
                    auto new_start = start + root;
                    cube_root_quick_sort(start, new_start, comp);

                    std::vector<RandomAccessIterator> pivots(root);
                    for (int i = 0; i < root; i++) pivots[i] = i + start;

                    for (auto i = new_start; i < stop; i++) {
                        int left = 0, right = root;
                        while (left < right) {
                            int mid = (right - left) / 2 + left;
                            if (comp(*i, *pivots[mid])) right = mid;
                            else left = mid + 1;
                        }
                        auto pos = i;
                        for (int j = root - 1; j >= left; j--) {
                            std::iter_swap(pivots[j] + 1, pos);
                            std::iter_swap(pos = pivots[j], pivots[j] + 1);
                            pivots[j]++;
                        }
                    }

                    cube_root_quick_sort(start, pivots[0], comp);
                    for (int i = 1; i < root; i++)
                        cube_root_quick_sort(pivots[i - 1] + 1, pivots[i], comp);
                    cube_root_quick_sort(pivots[root - 1] + 1, stop, comp);
                }
            }

            
        }

        namespace ternary_quick_sorting {
            using sortcpp::_internal::equals;

            template<class RandomAccessIterator>
            struct pivot_pair {
                RandomAccessIterator first, second;
                pivot_pair() {}
                pivot_pair(RandomAccessIterator f, RandomAccessIterator s) : first(f), second(s) {}
            };

            template<class RandomAccessIterator, class Compare>
            bool compare(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                return comp(*lo, *hi);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator select_pivot(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                auto mid = lo + (hi - lo) / 2;

                if (equals(*lo, *mid, comp)) return lo;
                if (equals(*lo, *(hi - 1), comp) || equals(*mid, *(hi - 1), comp)) return hi - 1;

                return compare(lo, mid, comp)
                    ? (compare(mid, hi - 1, comp) ? mid : (compare(lo, hi - 1, comp) ? hi - 1 : lo))
                    : (compare(hi - 1, mid, comp) ? mid : (compare(lo, hi - 1, comp) ? lo : hi - 1));
            }

            template<class RandomAccessIterator, class Compare>
            pivot_pair<RandomAccessIterator> _partition_ll(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                auto p = select_pivot(lo, hi, comp);
                auto pivot = *p;

                std::iter_swap(p, hi - 1);

                auto i = lo, k = hi - 1;

                for (auto j = lo; j < k; j++) {
                    if (equals(*j, pivot, comp)) {
                        std::iter_swap(--k, j);
                        --j;
                    }
                    else if (comp(*j, pivot)) {
                        std::iter_swap(i++, j);
                    }
                }

                auto j = i + (hi - k);

                for (int s = 0; s < hi - k; ++s) {
                    std::iter_swap(i + s, hi - 1 - s);
                }

                return pivot_pair<RandomAccessIterator>(i, j);
            }

            template<class RandomAccessIterator, class Compare>
            void ternary_quick_sort_ll(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                if (lo + 1 < hi) {
                    pivot_pair<RandomAccessIterator> mid = _partition_ll(lo, hi, comp);
                    ternary_quick_sort_ll(lo, mid.first, comp);
                    ternary_quick_sort_ll(mid.second, hi, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _ternary_qsort(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                if (hi <= lo) return;

                auto piv = select_pivot(lo, hi + 1, comp);
                std::iter_swap(piv, hi);

                auto pivot = *hi;

                auto i = lo, j = hi - 1;
                auto p = lo, q = hi - 1;

                for (;;) {
                    while (i <= j && !comp(pivot, *i)) {
                        if (equals(*i, pivot, comp)) {
                            std::iter_swap(i, p++);
                        }
                        ++i;
                    }

                    while (i <= j && !comp(*j, pivot)) {
                        if (equals(*j, pivot, comp)) {
                            std::iter_swap(j, q--);
                        }
                        --j;
                    }

                    if (i > j) break;

                    std::iter_swap(i++, j--);
                }

                std::iter_swap(i, hi);

                auto num_less = i - p;
                auto num_greater = q - j;

                j = i - 1; i = i + 1;

                auto pe = lo + std::min(p - lo, num_less);
                for (auto k = lo; k < pe; k++, j--) {
                    std::iter_swap(k, j);
                }

                auto qe = hi - 1 - std::min(hi - 1 - q, num_greater - 1);
                for (auto k = hi - 1; k > qe; k--, i++) {
                    std::iter_swap(i, k);
                }

                _ternary_qsort(lo, lo + num_less - 1, comp);
                _ternary_qsort(hi - num_greater + 1, hi, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void ternary_quick_sort_lr(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                _ternary_qsort(lo, hi - 1, comp);
            }
        }

        namespace singularity_quick_sorting {
            using sortcpp::_internal::reversal;
            using sortcpp::_internal::shell_sorting::shell_sort;
            using sortcpp::_internal::adaptive_insertion_sorting::pd_binary_insertion_sort;
            using sortcpp::_internal::adaptive_insertion_sorting::stable_segment_reversal;

            int _log2(int x) {
                int n = 1;
                while (1 << n < x) n++;
                if (1 << n > x) n--;
                return n;
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator unstable_pd(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto reverse = start;
                bool different = false;
                int cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                while (cmp >= 0 && reverse + 1 < end) {
                    if (cmp == 1) different = true;
                    reverse++;
                    cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                }

                if (reverse > start && different) {
                    if (reverse < start + 4) std::iter_swap(start, reverse);
                    else reversal(start, reverse + 1);
                }
                return reverse;
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator pd(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto reverse = start;
                bool less_unique = false;
                bool different = false;

                int cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                while (cmp >= 0 && reverse + 1 < end) {
                    if (cmp == 0) less_unique = true;
                    else different = true;
                    reverse++;
                    cmp = comp(*(reverse + 1), *reverse) - comp(*reverse, *(reverse + 1));
                }
                if (reverse > start && different) {
                    if (less_unique) stable_segment_reversal(start, reverse, comp);
                    else if (reverse < start + 3) std::iter_swap(start, reverse);
                    else reversal(start, reverse + 1);
                }
                return reverse;
            }

            template<class RandomAccessIterator, class Compare>
            void _shell_unstable_sqsort(RandomAccessIterator start, RandomAccessIterator offset, RandomAccessIterator end, 
                                        int depth, int rep, int depth_limit, int insert_limit, int rep_limit, Compare comp) {
                if (depth == depth_limit || rep == 4) {
                    shell_sort(start - 1, end, comp);
                    return;
                }
                if (end - start > insert_limit) {
                    auto left = offset;
                    while (!comp(*left, *(left - 1)) && left < end) left++;
                    if (left < end) {
                        auto pivot = *(left - 1);
                        auto pivot_pos = left - 1;
                        auto original_pos = left - 1;
                        auto right = left + 1;
                        auto item = *left;
                        bool broke_loop = false;
                        bool broken_cond = false;
                        bool found_item = false;

                        while (right <= end) {
                            if (comp(*(right - 1), pivot)) {
                                if (right - left == 1) {
                                    if (!found_item) item = *(left - 1);
                                    found_item = true;
                                    *(left - 1) = *left;
                                }
                                else broke_loop = true;
                                if (broke_loop && !broken_cond) {
                                    *(left - 1) = item;
                                    broken_cond = true;
                                }
                                if (right - left > 1) std::iter_swap(left - 1, right - 1);
                                if (pivot_pos == left - 1) pivot_pos = right - 1;
                                left++;
                            }
                            right++;
                        }

                        if (right > end && !broke_loop) *(left - 1) = item;
                        if (pivot_pos != left - 1) std::iter_swap(pivot_pos, left - 1);

                        bool lsmall = left - start < end - (left + 1);
                        if (lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _shell_unstable_sqsort(
                                    start, 
                                    original_pos - 1 > start ? original_pos - 1 : start, 
                                    left - 1, depth + 1, rep + 1,
                                    depth_limit, insert_limit, rep_limit, comp
                                );

                            else _shell_unstable_sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, 0, 
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (end - (left + 1) > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _shell_unstable_sqsort(
                                    left + 1, left + 1, end, depth + 1, rep + 1, 
                                    depth_limit, insert_limit, rep_limit, comp
                                );

                            else _shell_unstable_sqsort(
                                left + 1, left + 1, end, depth + 1, 0, 
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (!lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _shell_unstable_sqsort(
                                    start, 
                                    original_pos - 1 > start ? original_pos - 1 : start, 
                                    left - 1, depth + 1, rep + 1,
                                    depth_limit, insert_limit, rep_limit, comp
                                );
                            else _shell_unstable_sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, 0,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                    }
                }
                else shell_sort(start - 1, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void shell_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int depth_limit = std::min((int)std::sqrt(length), 2 * _log2(length));
                int insert_limit = std::max((depth_limit / 2) - 1, 15);
                int rep_limit = std::max((depth_limit / 4), 2);

                auto real_start = unstable_pd(first, last, comp);
                if (real_start + 1 < last) 
                    _shell_unstable_sqsort(
                        first + 1, real_start + 1, last, 0, 0, 
                        depth_limit, insert_limit, rep_limit, comp
                    );
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator _binary_search(RandomAccessIterator a, RandomAccessIterator b, const T& value, Compare comp) {
                while (a < b) {
                    auto m = a + ((b - a) / 2);
                    if (comp(value, *m)) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            void binsert(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                pd_binary_insertion_sort(start - 1, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _sqsort(RandomAccessIterator start, RandomAccessIterator offset, RandomAccessIterator end,
                         int depth, int rep, int depth_limit, int insert_limit, int rep_limit, Compare comp) {

                if (end - start > insert_limit && depth < depth_limit && rep < 4) {
                    auto left = offset;
                    while (!comp(*left, *(left - 1)) && left < end) left++;
                    if (left < end) {
                        auto right = left + 1;
                        auto pull = end;
                        auto pivot = *(left - 1);
                        auto original_pos = left - 1;
                        bool broke_loop = false;
                        bool broken_cond = false;

                        while (right <= end) {
                            if (comp(*(right - 1), pivot)) {
                                if (right - left == 1) *(left - 1) = *left;
                                else broke_loop = true;
                                if (broke_loop && !broken_cond) {
                                    *(left - 1) = pivot;
                                    broken_cond = true;
                                }
                                if (right - left > 1) {
                                    pull = right - 1;
                                    auto item = *pull;

                                    while (pull >= left) {
                                        *pull = *(pull - 1);
                                        pull--;
                                    }
                                    *pull = item;
                                }
                                left++;
                            }
                            right++;
                        }
                        if (right > end && !broke_loop) *(left - 1) = pivot;

                        bool lsmall = left - start < end - (left + 1);
                        if (lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) _sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start,
                                left - 1, depth + 1, rep + 1,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                            else _sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, 0,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (end - (left + 1) > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) _sqsort(
                                left + 1, left + 1, end, depth + 1, rep + 1,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                            else _sqsort(
                                left + 1, left + 1, end, depth + 1, 0,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (!lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) _sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, rep + 1,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                            else _sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, 0,
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                    }
                }
                else binsert(start, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int depth_limit = std::min((int)std::sqrt(length), 2 * _log2(length));
                int insert_limit = std::max((depth_limit / 2) - 1, 15);
                int rep_limit = std::max((depth_limit / 4), 2);

                auto real_start = pd(first, last, comp);
                if (real_start + 1 < last) 
                    _sqsort(
                        first + 1, real_start + 1, last, 0, 0, 
                        depth_limit, insert_limit, rep_limit, comp
                    );
            }

            template<class RandomAccessIterator, class Compare>
            bool _unbounded_sqsort(RandomAccessIterator beg, RandomAccessIterator start, RandomAccessIterator offset,
                                   RandomAccessIterator end, RandomAccessIterator full, int depth, bool sorted, Compare comp) {

                auto left = offset;
                while (!comp(*left, *(left - 1)) && left < end) left++;

                if (left < end && depth < 2047 && !sorted) {
                    auto right = left + 1;
                    auto pull = end;
                    auto pivot = *(left - 1);
                    auto original_pos = left - 1;
                    bool broke_loop = false;
                    bool broken_cond = false;

                    while (right <= end) {
                        if (comp(*(right - 1), pivot)) {
                            if (right - left == 1) *(left - 1) = *left;
                            else broke_loop = true;
                            if (broke_loop && !broken_cond) {
                                *(left - 1) = pivot;
                                broken_cond = true;
                            }
                            if (right - left > 1) {
                                pull = right - 1;
                                auto item = *pull;

                                while (pull >= left) {
                                    *pull = *(pull - 1);
                                    pull--;
                                }
                                *pull = item;
                            }
                            left++;
                        }
                        right++;
                    }

                    if (right > end && !broke_loop) *(left - 1) = pivot;
                    bool lsmall = left - start < end - (left + 1);
                    if (lsmall && (left - 1) - start > 0) 
                        sorted = _unbounded_sqsort(
                            beg, start, 
                            original_pos - 1 > start ? original_pos - 1 : start, 
                            left - 1, full, depth + 1, sorted, comp
                        );

                    if (end - (left + 1) > 0) sorted = _unbounded_sqsort(beg, left + 1, left + 1, end, full, depth + 1, sorted, comp);
                    if (!lsmall && (left - 1) - start > 0) 
                        sorted = _unbounded_sqsort(
                            beg, start, original_pos - 1 > start ? original_pos - 1 : start, 
                            left - 1, full, depth + 1, sorted, comp
                        );
                }
                else if (!sorted) {
                    singularity_quick_sort(beg, full, comp);
                    sorted = true;
                }
                return sorted;
            }

            template<class RandomAccessIterator, class Compare>
            void unbounded_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _unbounded_sqsort(first, first + 1, first + 1, last, last, 0, false, comp);
            }

            template<class RandomAccessIterator>
            void _shuffle(RandomAccessIterator start, RandomAccessIterator end) {
                using sortcpp::_internal::bogo_sorting::rd;
                using sortcpp::_internal::bogo_sorting::randiter;

                std::mt19937_64 rng(rd());
                for (auto i = start; i < end; i++) {
                    auto index = randiter(i, end, rng);
                    if (index != i) std::iter_swap(i, index);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void _unstable_sqsort(RandomAccessIterator start, RandomAccessIterator offset, RandomAccessIterator end,
                                  int depth, int real_depth, int rep, int depth_limit, int insert_limit, int rep_limit, Compare comp) {
                if (end - start > insert_limit && (depth == depth_limit || rep == 4)) {
                    _shuffle(start - 1, end - 1);
                    offset = start;
                    depth = 0;
                    rep = 0;
                }

                if (end - start > insert_limit && depth < depth_limit) {
                    auto left = offset;
                    while (!comp(*left, *(left - 1)) && left < end) left++;

                    if (left < end) {
                        auto pivot = *(left - 1);
                        auto pivot_pos = left - 1;
                        auto original_pos = left - 1;
                        auto right = left + 1;
                        auto item = *left;
                        bool broke_loop = false;
                        bool broken_cond = false;
                        bool found_item = false;

                        while (right <= end) {
                            if (comp(*(right - 1), pivot)) {
                                if (right - left == 1) {
                                    if (!found_item) item = *(left - 1);
                                    found_item = true;
                                    *(left - 1) = *left;
                                }
                                else broke_loop = true;
                                if (broke_loop && !broken_cond) {
                                    *(left - 1) = item;
                                    broken_cond = true;
                                }

                                if (right - left > 1) std::iter_swap(left - 1, right - 1);
                                if (pivot_pos == left - 1) pivot_pos = right - 1;
                                left++;
                            }
                            right++;
                        }

                        if (right > end && !broke_loop) *(left - 1) = item;
                        if (pivot_pos != left - 1) std::iter_swap(pivot_pos, left - 1);

                        bool lsmall = left - start < end - (left + 1);
                        if (lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _unstable_sqsort(
                                    start, 
                                    original_pos - 1 > start ? original_pos - 1 : start, 
                                    left - 1, depth + 1, real_depth + 1, rep + 1,
                                    depth_limit, insert_limit, rep_limit, comp
                                );
                            else _unstable_sqsort(
                                start, 
                                original_pos - 1 > start ? original_pos - 1 : start, 
                                left - 1, depth + 1, real_depth + 1, 0, 
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (end - (left + 1) > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _unstable_sqsort(
                                    left + 1, left + 1, end, depth + 1, real_depth + 1, rep + 1,
                                    depth_limit, insert_limit, rep_limit, comp
                                );
                            else _unstable_sqsort(
                                left + 1, left + 1, end, depth + 1, real_depth + 1, 0, 
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                        if (!lsmall && (left - 1) - start > 0) {
                            if (end - rep_limit <= left || left <= start + rep_limit) 
                                _unstable_sqsort(
                                    start, 
                                    original_pos - 1 > start ? original_pos - 1 : start, 
                                    left - 1, depth + 1, real_depth + 1, rep + 1,
                                    depth_limit, insert_limit, rep_limit, comp
                                );
                            else _unstable_sqsort(
                                start,
                                original_pos - 1 > start ? original_pos - 1 : 
                                start, left - 1, depth + 1, real_depth + 1, 0, 
                                depth_limit, insert_limit, rep_limit, comp
                            );
                        }
                    }
                }
                else binsert(start, end, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int depth_limit = std::min((int)std::sqrt(length), 2 * _log2(length));
                int insert_limit = std::max((depth_limit / 2) - 1, 15);
                int rep_limit = std::max((depth_limit / 4), 2);

                auto real_start = unstable_pd(first, last, comp);
                if (real_start + 1 < last)
                    _unstable_sqsort(
                        first + 1, real_start + 1, last, 0, 0, 0,
                        depth_limit, insert_limit, rep_limit, comp
                    );
            }


            template<class RandomAccessIterator, class Compare>
            bool _unbounded_unstable_sqsort(RandomAccessIterator beg, RandomAccessIterator start, RandomAccessIterator offset,
                                            RandomAccessIterator end, RandomAccessIterator full, int depth, bool sorted, Compare comp) {
                auto left = offset;
                while (!comp(*left, *(left - 1)) && left < end) left++;

                if (left < end && depth < 2048 && !sorted) {
                    auto pivot = *(left - 1);
                    auto pivot_pos = left - 1;
                    auto original_pos = left - 1;
                    auto right = left + 1;
                    auto item = *left;
                    bool broke_loop = false;
                    bool broken_cond = false;
                    bool found_item = false;

                    while (right <= end) {
                        if (comp(*(right - 1), pivot)) {
                            if (right - left == 1) {
                                if (!found_item) item = *(left - 1);
                                found_item = true;
                                *(left - 1) = *left;
                            }
                            else broke_loop = true;
                            if (broke_loop && !broken_cond) {
                                *(left - 1) = item;
                                broken_cond = true;
                            }
                            if (right - left > 1) std::iter_swap(left - 1, right - 1);
                            if (pivot_pos == left - 1) pivot_pos = right - 1;
                            left++;
                        }
                        right++;
                    }
                    if (right > end && !broke_loop) *(left - 1) = item;
                    if (pivot_pos != left - 1) std::iter_swap(pivot_pos, left - 1);
                    bool lsmall = left - start < end - (left + 1);

                    if (lsmall && (left - 1) - start > 0) 
                        sorted = _unbounded_unstable_sqsort(
                            beg, start, 
                            original_pos - 1 > start ? original_pos - 1 : start, 
                            left - 1, full, depth + 1, sorted, comp
                        );
                    if (end - (left + 1) > 0) sorted = _unbounded_unstable_sqsort(beg, left + 1, left + 1, end, full, depth + 1, sorted, comp);
                    if (!lsmall && (left - 1) - start > 0) 
                        sorted = _unbounded_unstable_sqsort(
                            beg, start, 
                            original_pos - 1 > start ? original_pos - 1 : start, 
                            left - 1, full, depth + 1, sorted, comp
                        );
                }
                else if (depth == 2048 && !sorted) {
                    unstable_singularity_quick_sort(beg, full, comp);
                    sorted = true;
                }
                return sorted;
            }

            template<class RandomAccessIterator, class Compare>
            void unbounded_unstable_singularity_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _unbounded_unstable_sqsort(first, first + 1, first + 1, last, last, 0, false, comp);
            }
        }


        namespace iterative_quick_sorting {
            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _partition(RandomAccessIterator lo, RandomAccessIterator hi, Compare comp) {
                auto pivot = *hi;
                auto i = lo - 1;
                for (auto j = lo; j <= hi - 1; j++)
                    if (!comp(pivot, *j)) {
                        i++;
                        std::iter_swap(i, j);
                    }
                std::iter_swap(i + 1, hi);
                return i + 1;
            }

            template<class RandomAccessIterator, class Compare>
            void iterative_quick_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                std::stack<RandomAccessIterator> stack;
                stack.push(start);
                stack.push(end);

                while (!stack.empty()) {
                    end = stack.top(); stack.pop();
                    start = stack.top(); stack.pop();

                    auto p = _partition(start, end, comp);
                    if (p - 1 > start) {
                        stack.push(start);
                        stack.push(p - 1);
                    }
                    if (p + 1 < end) {
                        stack.push(p + 1);
                        stack.push(end);
                    }
                }
            }

            template<class RandomAccessIterator>
            struct Task {
                RandomAccessIterator p, r;
                Task() {}
                Task(RandomAccessIterator _p, RandomAccessIterator _r) : p(_p), r(_r) {}
            };

            template<class RandomAccessIterator, class Compare>
            void linked_iterative_quick_sort(RandomAccessIterator p, RandomAccessIterator r, Compare comp) {
                std::queue<Task<RandomAccessIterator>> tasks;
                tasks.emplace(p, r);

                while (!tasks.empty()) {
                    auto task = tasks.front();
                    tasks.pop();

                    auto i = task.p, j = task.r;
                    auto pivot = task.p + (task.r - task.p) / 2;
                    auto x = *pivot;

                    while (i <= j) {
                        while (comp(*i, x)) i++;
                        while (comp(x, *j)) j--;

                        if (i <= j) {
                            std::iter_swap(i, j);
                            i++;
                            j--;
                        }
                    }

                    if (task.p < j) tasks.emplace(task.p, j);
                    if (i < task.r) tasks.emplace(i, task.r);
                }
            }

            using sortcpp::_internal::equals;

            template<class RandomAccessIterator, class Compare>
            void _median_of_three(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto m = a + (b - 1 - a) / 2;
                if (comp(*m, *a)) std::iter_swap(a, m);
                if (comp(*(b - 1), *m)) {
                    std::iter_swap(m, b - 1);
                    if (comp(*m, *a)) return;
                }
                std::iter_swap(a, m);
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _stackless_partition(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = a, j = b;
                _median_of_three(a, b, comp);

                do {
                    do i++;
                    while (i < j && comp(*i, *a));

                    do j--;
                    while (j >= i && !comp(*j, *a));

                    if (i < j) std::iter_swap(i, j);
                    else {
                        std::iter_swap(a, j);
                        return j;
                    }

                } while (true);
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator _left_bin_search(RandomAccessIterator a, RandomAccessIterator b, RandomAccessIterator p, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(*m, *p)) a = m + 1;
                    else b = m;
                }
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            void stackless_quick_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto max = *a;
                for (auto i = a + 1; i < b; i++)
                    if (comp(max, *i)) max = *i;

                for (auto i = b - 1; i >= a; i--)
                    if (equals(*i, max, comp)) std::iter_swap(i, --b);

                auto b1 = b;
                do {
                    while (b1 - a > 2) {
                        auto p = _stackless_partition(a, b1, comp);
                        std::iter_swap(p, b);
                        b1 = p;
                    }

                    if (b1 - a == 2 && comp(*(a + 1), *a)) std::iter_swap(a, a + 1);

                    a = b1 + 1;
                    if (a >= b) {
                        if (a - 1 < b) std::iter_swap(a - 1, b);
                        return;
                    }

                    b1 = _left_bin_search(a, b, a - 1, comp);
                    std::iter_swap(a - 1, b);

                    while (a < b1 && equals(*(a - 1), *a, comp)) a++;
                } while (true);
            }
        }

        namespace stable_quick_sorting {
            template<class T, class RandomAccessIterator>
            void _copy(const std::vector<T>& nums, RandomAccessIterator start) {
                for (auto num : nums) {
                    *start = num;
                    start++;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            RandomAccessIterator _stable_partition(RandomAccessIterator start, RandomAccessIterator end, RandomAccessIterator pivot, Compare comp) {
                auto pivot_value = *pivot;

                std::vector<T> left, right;

                for (auto i = start; i <= end; i++) {
                    if (i == pivot) continue;
                    if (comp(*i, pivot_value)) left.push_back(*i);
                    else right.push_back(*i);
                }

                _copy(left, start);

                auto new_pivot_index = start + left.size();
                *new_pivot_index = pivot_value;

                _copy(right, new_pivot_index + 1);

                return new_pivot_index;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void stable_quick_loop(RandomAccessIterator start, RandomAccessIterator end, bool middle_pivot, Compare comp) {
                if (start < end) {
                    auto pivot = start;
                    if (middle_pivot) pivot = start + (end - start - 1) / 2;

                    auto pivot_index = _stable_partition<T>(start, end, pivot, comp);
                    stable_quick_loop<T>(start, pivot_index - 1, middle_pivot, comp);
                    stable_quick_loop<T>(pivot_index + 1, end, middle_pivot, comp);
                }
            }
        }

        namespace table_sorting {
            using sortcpp::_internal::equals;

            template<class RandomAccessIterator, class Compare>
            bool stable_comp(RandomAccessIterator array, std::vector<int> &table, int a, int b, Compare comp) {
                auto x = *(array + table[a]);
                auto y = *(array + table[b]);
                return comp(y, x) || (equals(x, y, comp) && table[a] > table[b]);
            }

            template<class RandomAccessIterator, class Compare>
            void stable_median_of_three(RandomAccessIterator array, std::vector<int>& table, int a, int b, Compare comp) {
                int m = a + (b - 1 - a) / 2;

                if (stable_comp(array, table, a, m, comp)) std::swap(table[a], table[m]);

                if (stable_comp(array, table, m, b - 1, comp)) {
                    std::swap(table[m], table[b - 1]);
                    if (stable_comp(array, table, a, m, comp)) return;
                }

                std::swap(table[a], table[m]);
            }

            template<class RandomAccessIterator, class Compare>
            int _stable_partition(RandomAccessIterator array, std::vector<int>& table, int a, int b, int p, Compare comp) {
                int i = a - 1, j = b;
                while (true) {
                    do i++;
                    while (i < j && !stable_comp(array, table, i, p, comp));

                    do j--;
                    while (j >= i && stable_comp(array, table, j, p, comp));

                    if (i < j) std::swap(table[i], table[j]);
                    else return j;
                }
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            void _quick_sort(RandomAccessIterator array, std::vector<int>& table, int a, int b, Compare comp) {
                if (b - a < 3) {
                    if (b - a == 2 && stable_comp(array, table, a, a + 1, comp)) std::swap(table[a], table[a + 1]);
                    return;
                }

                stable_median_of_three(array, table, a, b, comp);
                int p = _stable_partition(array, table, a + 1, b, a, comp);
                std::swap(table[a], table[p]);

                _quick_sort(array, table, a, p, comp);
                _quick_sort(array, table, p + 1, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void _table_sort(RandomAccessIterator first, int length, std::vector<int>& table, Compare comp) {
                _quick_sort(first, table, 0, length, comp);
                for (int i = 0; i < table.size(); i++)
                    if (i != table[i]) {
                        auto t = *(first + i);
                        int j = i, nxt = table[i];

                        do {
                            *(first + j) = *(first + nxt);
                            table[j] = j;

                            j = nxt;
                            nxt = table[nxt];
                        } while (nxt != i);

                        *(first + j) = t;
                        table[j] = j;
                    }
            }

            template<class RandomAccessIterator, class Compare>
            void table_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> table(length);
                std::iota(table.begin(), table.end(), 0);
                _table_sort(first, length, table, comp);
            }
        }

        namespace gnome_sorting {
            using sortcpp::_internal::equals;
            using sortcpp::_internal::reversal;
            using sortcpp::_internal::bogo_sorting::randbool;
            using sortcpp::_internal::binary_insertion_sorting::bin_insert;

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

            template<class RandomAccessIterator, class Compare>
            void smart_gnome_loop(RandomAccessIterator lower, RandomAccessIterator upper, Compare comp) {
                auto pos = upper;
                while (pos > lower && comp(*pos, *(pos - 1))) {
                    std::iter_swap(pos - 1, pos);
                    pos--;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void smart_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = first + 1; i < last; i++) smart_gnome_loop(first, i, comp);
            }


            template<class RandomAccessIterator, class Compare>
            void gambit_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bin_insert(first, last, last - first, comp);
                smart_gnome_sort(first, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void pd_gnome_sort(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto i = a + 1;
                if (comp(*i, *(i - 1))) {
                    do i++;
                    while (i < b && comp(*i, *(i - 1)));
                    reversal(a, i);
                }
                else {
                    do i++;
                    while (i < b && !comp(*i, *(i - 1)));
                }

                while (i < b) {
                    auto pos = i;
                    while (pos > a && comp(*pos, *(pos - 1))) {
                        std::iter_swap(pos, pos - 1);
                        pos--;
                    }
                    i++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void reverse_gnome_loop(RandomAccessIterator lower, RandomAccessIterator upper, Compare comp) {
                auto pos = lower;
                while (pos < upper && comp(*(pos + 1), *pos)) {
                    std::iter_swap(pos, pos + 1);
                    pos++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void reverse_gnome_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                for (auto i = last - 1; i >= first; i--) reverse_gnome_loop(i, last - 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void tri_gnome_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                using sortcpp::_internal::_tri_search;
                for (auto i = start + 1; i < end; i++) {
                    auto num = *i;
                    auto lo = start;

                    lo = _tri_search(start, i - 1, num, comp);
                    auto j = i;
                    while (j > lo) {
                        std::iter_swap(j, j - 1);
                        j--;
                    }
                }
            }

            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void markov_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
                for (auto i = first; i < last - 1; i++) {
                    auto walk = i + 1;
                    while ((walk == first ? false : comp(*walk, *(walk - 1)))
                        || (walk > i ? false : comp(*(walk + 1), *walk))) {
                        int c = (walk == first || walk <= i && randbool(rng)) ? 1 : -1;
                        std::iter_swap(walk, walk + c);
                        walk += c;
                    }
                }
            }
        }

        namespace push_sorting {
            template<class RandomAccessIterator, class Compare>
            void push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
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
                }
            }

            template<class RandomAccessIterator, class Compare>
            void reverse_push_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bool any_swaps = true;
                auto i = last;
                int gap = 1;
                while (any_swaps) {
                    any_swaps = false;
                    i = last;
                    gap = 1;
                    while (i - gap > first) {
                        if (comp(*(i - 1), *(i - 1 - gap))) {
                            for (int j = 1; j <= gap; j++) std::iter_swap(i - 1 - j, i - 1);
                            any_swaps = true;
                            gap++;
                        }
                        else i--;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void strange_push_sort(RandomAccessIterator first, RandomAccessIterator last, int base, Compare comp) {
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
                            gap *= base;
                        }
                        else i++;
                    }
                }
            }
        }

        namespace wiggle_sorting {
            template<class RandomAccessIterator, class Compare>
            void wiggle_sort(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                if (end - start < 2) return;
                auto left = start, right = end;
                auto mid = left + (right - left) / 2;

                bool start_left = true;
                auto j = mid;

                for (auto i = left; i < mid; i++) {
                    for (auto k = mid; k < right; k++) {
                        if (!comp(*i, *j)) std::iter_swap(i, j);
                        if (start_left) j++;
                        else j--;
                    }
                    if (start_left) {
                        j--;
                        start_left = false;
                    }
                    else {
                        j++;
                        start_left = true;
                    }
                }

                wiggle_sort(start, mid, comp);
                wiggle_sort(mid, end, comp);
            }
        }

        namespace stooge_sorting {
            // Swap the elements in [left, right], [left, right - 1], ..., [left, left + 1]
            // first - The real beginning of the range.
            // left - The beginning of the current range.
            // right - The end of the current range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void forward(RandomAccessIterator left, RandomAccessIterator right, Compare comp) {
                auto first = left;
                while (left < right) {
                    auto index = right;

                    while (left < index) {
                        if (comp(*index, *left)) std::iter_swap(left, index);
                        left++;
                        index--;
                    }

                    left = first;
                    right--;
                }
            }


            // Swap the elements in [left, right], [left + 1, right], ..., [right - 1, right]
            // left - The beginning of the current range.
            // right - The end of the current range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void backward(RandomAccessIterator left, RandomAccessIterator right, Compare comp) {
                auto tmp = right;

                while (left < right) {
                    auto index = left;

                    while (index < right) {
                        if (comp(*right, *index)) std::iter_swap(index, right);
                        index++;
                        right--;
                    }

                    left++;
                    right = tmp;
                }
            }


            // Sort the range [first, last) using the optimized stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void optimized_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto left = first, right = last - 1;
                while (left < right) {
                    if (comp(*right, *left)) std::iter_swap(left, right);
                    left++;
                    right--;
                }

                forward(first, last - 2, comp);
                backward(first + 1, last - 1, comp);
            }
            
            // If `a` is greater than `b`, swap them.
            // a - The first element.
            // b - The second element.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool comp_swap(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (comp(*b, *a)) {
                    std::iter_swap(a, b);
                    return true;
                }
                return false;
            }


            // Sort the range [first, last) using the studio version of the optimized stooge sort algorithm.
            // first - The real beginning of the range.
            // a - The beginning of the current range.
            // m - The middle position of the current range.
            // b - The end of the current range. (exclusive)
            // merge - Whether to merge the sorted ranges.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool optistooge_loop(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, bool merge, Compare comp) {
                if (a >= m) return false;
                if (b - a == 2) return comp_swap(a, m, comp);

                bool lchange = false, rchange = false;

                auto a2 = first + ((a - first) * 2 + (b - first)) / 3;
                auto b2 = first + ((a - first) + 2 * (b - first) + 2) / 3;

                if (m < b2) {
                    lchange = optistooge_loop(first, a, m, b2, merge, comp);

                    if (merge) {
                        rchange = optistooge_loop(first, std::max(a + (b2 - first) - (m - first), a2), b2, b, true, comp);
                        if (rchange) optistooge_loop(first, a + (b2 - first) - (m - first), a2, first + 2 * (a2 - first) - (a - first), true, comp);
                    }
                    else {
                        rchange = optistooge_loop(first, a2, b2, b, false, comp);
                        if (rchange) optistooge_loop(first, a, a2, first + 2 * (a2 - first) - (a - first), true, comp);
                    }
                }
                else {
                    rchange = optistooge_loop(first, a2, m, b, merge, comp);
                    if (rchange) optistooge_loop(first, a, a2, a2 + (b - first) - (m - first), true, comp);
                }
                return lchange || rchange;
            }


            // Sort the range [first, last) using the studio version of the optimized stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void optimized_stooge_sort_studio(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                optistooge_loop(first, first, first + 1, last, false, comp);
            }

            // Sort the range [pos, pos + len) using the quad stooge sort algorithm.
            // pos - The beginning of current range.
            // len - The length of current range.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void quad_stooge_loop(RandomAccessIterator pos, int len, Compare comp) {
                if (len >= 2 && comp(*(pos + len - 1), *pos)) std::iter_swap(pos, pos + len - 1);
                if (len <= 2) return;

                int len1 = len / 2;
                int len2 = (len + 1) / 2;
                int len3 = (len1 + 1) / 2 + (len2 + 1) / 2;

                quad_stooge_loop(pos, len1, comp);
                quad_stooge_loop(pos + len1, len2, comp);
                quad_stooge_loop(pos + len1 / 2, len3, comp);
                quad_stooge_loop(pos + len1, len2, comp);
                quad_stooge_loop(pos, len1, comp);

                if (len > 3) quad_stooge_loop(pos + len1 / 2, len3, comp);
            }


            // Sort the range [first, last) using the quad stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void quad_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                quad_stooge_loop(first, last - first, comp);
            }


            // Sort the range [i, j] use the stooge sorting algorithm.
            // i - The beginning of the range.
            // j - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _stooge_loop(RandomAccessIterator i, RandomAccessIterator j, Compare comp) {
                if (comp(*j, *i)) std::iter_swap(i, j);
                if (j - i + 1 >= 3) {
                    int t = (j - i + 1) / 3;

                    _stooge_loop(i, j - t, comp);
                    _stooge_loop(i + t, j, comp);
                    _stooge_loop(i, j - t, comp);
                }
            }

            // Sort the range [first, last) using the stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _stooge_loop(first, last - 1, comp);
            }


            // Sort the range [pos, pos + len) use the awkward sorting algorithm.
            // pos - The beginning of the range.
            // len - The length of the range.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _awkward_loop(RandomAccessIterator pos, int len, Compare comp) {
                if (len == 1) return;

                _awkward_loop(pos, len / 2, comp);
                _awkward_loop(pos + len / 2, len / 2 + len % 2, comp);

                for (int i = 0; i < len / 2; i++) {
                    auto a = pos + i;
                    auto b = pos + len / 2 + len % 2 + i;
                    if (comp(*b, *a)) std::iter_swap(a, b);
                }

                _awkward_loop(pos + len / 4, len / 2 + len % 2, comp);
                _awkward_loop(pos, len / 2, comp);
                _awkward_loop(pos + len / 2, len / 2 + len % 2, comp);
            }


            // Sort the range [first, last) using the awkward sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void awkward_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _awkward_loop(first, last - first, comp);
            }


            // Sort the range [start, end] use the stable stooge sorting algorithm.
            // start - The beginning of the range.
            // end - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _stable_stooge_loop(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                if (end - start + 1 == 2) {
                    if (comp(*end, *start)) std::iter_swap(start, end);
                }
                else if (end - start + 1 > 2) {
                    auto third = (end - start + 1) / 3;
                    _stable_stooge_loop(start, end - third, comp);
                    _stable_stooge_loop(start + third, end, comp);
                    _stable_stooge_loop(start, end - third, comp);
                }
            }

            // Sort the range [first, last) using the stable stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void stable_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _stable_stooge_loop(first, last - 1, comp);
            }

            // Sort the range [i, j] use the hyper stooge sorting algorithm.
            // i - The beginning of the range.
            // j - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _hyper_stooge_loop(RandomAccessIterator i, RandomAccessIterator j, Compare comp) {
                if (comp(*j, *i)) std::iter_swap(i, j);
                if (j - i > 1) {
                    _hyper_stooge_loop(i, j - 1, comp);
                    _hyper_stooge_loop(i + 1, j, comp);
                    _hyper_stooge_loop(i, j - 1, comp);
                }
            }

            // Sort the range [first, last) using the hyper stooge sort algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void hyper_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _hyper_stooge_loop(first, last - 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void buffered_stooge_sort(RandomAccessIterator start, RandomAccessIterator stop, Compare comp) {
                if (stop - start > 1) {
                    if (stop - start == 2 && comp(*(stop - 1), *start)) std::iter_swap(start, stop - 1);
                    if (stop - start > 2) {
                        auto third = start + std::ceil((stop - start) / 3.0);
                        auto two_third = start + std::ceil((stop - start) / 3.0 * 2.0);

                        if (two_third - third < third - start) two_third--;
                        if ((stop - start - 2) % 3 == 0) two_third--;

                        buffered_stooge_sort(third, two_third, comp);
                        buffered_stooge_sort(two_third, stop, comp);

                        auto left = third, right = two_third, buffer_start = start;
                        while (left < two_third && right < stop) {
                            if (comp(*right, *left)) {
                                std::iter_swap(buffer_start, right);
                                right++;
                            }
                            else {
                                std::iter_swap(buffer_start, left);
                                left++;
                            }
                            buffer_start++;
                        }

                        while (right < stop) {
                            std::iter_swap(buffer_start, right);
                            right++;
                            buffer_start++;
                        }

                        buffered_stooge_sort(two_third, stop, comp);
                        left = two_third - 1;
                        right = stop - 1;
                        while (right > left && left >= start) {
                            if (comp(*right, *left)) {
                                for (auto i = left; i < right; i++) std::iter_swap(i, i + 1);
                                left--;
                            }
                            right--;
                        }
                    }
                }
            };
        }

        namespace slow_sorting {
            // Sort the range [i, j] use the silly sorting algorithm.
            // i - The beginning of the range.
            // j - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _silly_loop(RandomAccessIterator i, RandomAccessIterator j, Compare comp) {
                if (i < j) {
                    auto m = i + (j - i) / 2;

                    _silly_loop(i, m, comp);
                    _silly_loop(m + 1, j, comp);

                    if (comp(*(m + 1), *i)) std::iter_swap(i, m + 1);

                    _silly_loop(i + 1, j, comp);
                }
            }

            // Sort the range [first, last) use the silly sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void silly_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _silly_loop(first, last - 1, comp);
            }

            // Sort the range [i, j] use the slow sorting algorithm.
            // i - The beginning of the range.
            // j - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _slow_loop(RandomAccessIterator i, RandomAccessIterator j, Compare comp) {
                if (i >= j) return;

                auto m = i + (j - i) / 2;
                _slow_loop(i, m, comp);
                _slow_loop(m + 1, j, comp);

                if (comp(*j, *m)) std::iter_swap(m, j);

                _slow_loop(i, j - 1, comp);
            }

            // Sort the range [first, last) use the slow sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void slow_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _slow_loop(first, last - 1, comp);
            }

            // Sort the range [start, stop] use the snuffle sorting algorithm.
            // start - The beginning of the range.
            // stop - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _snuffle_loop(RandomAccessIterator start, RandomAccessIterator stop, Compare comp) {
                if (stop - start + 1 >= 2) {
                    if (comp(*stop, *start)) std::iter_swap(start, stop);
                    if (stop - start + 1 >= 3) {
                        auto mid = start + (stop - start) / 2;
                        for (int i = 0; i < std::ceil((stop - start + 1) / 2.0); i++) {
                            _snuffle_loop(start, mid, comp);
                            _snuffle_loop(mid, stop, comp);
                        }
                    }
                }
            }

            // Sort the range [first, last) use the snuffle sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void snuffle_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _snuffle_loop(first, last - 1, comp);
            }


            // If `a` is greater than `b`, swap them.
            // a - The first element.
            // b - The second element.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void comp_swap(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            // Sort the range [start, end] use the ternary slow sorting algorithm.
            // start - The beginning of the range.
            // end - The end of the range. (inclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void _ternary_slow_loop(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                if (end - start + 1 == 2) comp_swap(start, end, comp);
                else if (end - start + 1 > 2) {
                    int third = (end - start + 1) / 3;
                    auto mid1 = start + third;
                    auto mid2 = start + 2 * third;

                    _ternary_slow_loop(start, mid1 - 1, comp);
                    _ternary_slow_loop(mid1, mid2 - 1, comp);
                    _ternary_slow_loop(mid2 - 1, end, comp);

                    comp_swap(mid1 - 1, end, comp);
                    comp_swap(mid1 - 1, mid2 - 1, comp);
                    comp_swap(mid2 - 1, end, comp);

                    _ternary_slow_loop(start, end - 1, comp);
                }
            }

            // Sort the range [first, last) use the ternary slow sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void ternary_slow_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                _ternary_slow_loop(first, last - 1, comp);
            }
        }

        namespace permutation_sorting {
            using sortcpp::_internal::bogo_sorting::is_range_sorted;

            // Stably find a possible permutation of the range [first, last) that is sorted.
            // first - The beginning of the range.
            // idx - The index table that stabilizes the sort.
            // len - The length of the current range.
            // last - The real end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool stable_permute(RandomAccessIterator first, std::vector<int> &idx, int len, RandomAccessIterator last, Compare comp) {
                if (len < 2) return is_range_sorted(first, last, comp);

                for (int i = len - 2; i >= 0; i--) {
                    if (stable_permute(first, idx, len - 1, last, comp)) return true;

                    std::iter_swap(first + idx[i], first + idx[len - 1]);
                    std::swap(idx[i], idx[len - 1]);
                }
                if (stable_permute(first, idx, len - 1, last, comp)) return true;

                int t = idx[len - 1];
                for (int i = len - 1; i > 0; i--) idx[i] = idx[i - 1];
                idx[0] = t;

                auto t2 = *(first + idx[0]);
                for (int i = 1; i < len; i++) *(first + idx[i - 1]) = *(first + idx[i]);
                *(first + idx[len - 1]) = t2;

                return false;
            }

            // Stably sort the range [first, last) using permutation sort.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void stable_permutation_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> idx(length);
                std::iota(idx.begin(), idx.end(), 0);
                stable_permute(first, idx, length, last, comp);
            }


            // Find a possible permutation of the range [first, last) that is sorted.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // depth - The current depth of the recursion.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool permute(RandomAccessIterator first, RandomAccessIterator last, int depth, Compare comp) {
                int length = last - first;
                if (depth >= length - 1) return is_range_sorted(first, last, comp);

                for (int i = length - 1; i > depth; i--) {
                    if (permute(first, last, depth + 1, comp)) return true;

                    if ((length - depth) % 2 == 0) std::iter_swap(first + depth, first + i);
                    else std::iter_swap(first + depth, last - 1);
                }

                return permute(first, last, depth + 1, comp);
            }


            // Sort the range [first, last) using permutation sort.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void permutation_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                permute(first, last, 0, comp);
            }
        }

        namespace guess_sorting {
            using sortcpp::_internal::equals;
            using sortcpp::_internal::bogo_sorting::randint;


            // If `is_gt` is true, Determine if `a` is greater than `b`.
            // If `is_gt` is false, Determine if `a` is less than `b`.
            // a - The first element to compare.
            // b - The second element to compare.
            // is_gt - Represents whether to compare greater than or less than.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            bool compare(RandomAccessIterator a, RandomAccessIterator b, bool is_gt, Compare comp) {
                return is_gt ? comp(*b, *a) : comp(*a, *b);
            }


            // Reorder the range [first, first + length) according to the given index table.
            // first - The beginning of the range.
            // length - The length of the range.
            // indexes - The index table.
            template<class RandomAccessIterator>
            void rearrange(RandomAccessIterator first, int length, std::vector<int>& indexes) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                std::vector<T> aux(length);
                for (int i = 0; i < length; i++) aux[i] = *(first + indexes[i]);
                for (int i = 0; i < length; i++) *(first + i) = aux[i];
            }


            // Guess an possible index table that sorts the range [first, first + length).
            // first - The beginning of the range.
            // length - The length of the range.
            // loops - The index table used in loops.
            // indexes - The result index table.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void guess(RandomAccessIterator first, int length, std::vector<int>& loops, std::vector<int>& indexes, Compare comp) {
                while (true) {
                    int total = 0;
                    for (int i = 0; i < length; i++)
                        for (int j = 0; j < length; j++)
                            if (loops[i] == loops[j]) total++;

                    for (int i = 0; i < length; i++)
                        for (int j = 0; j < length; j++)
                            if (compare(first + loops[i], first + loops[j], i < j, comp)) total++;


                    if (total == length)
                        for (int i = 0; i < length; i++) indexes[i] = loops[i];

                    int pos = 0;
                    while (pos < length) {
                        if (loops[pos] < length - 1) {
                            loops[pos]++;
                            break;
                        }
                        else {
                            loops[pos] = 0;
                            pos++;
                        }
                    }
                    if (pos == length) break;
                }
            }


            // Sort the range [first, last) use the guess sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> loops(length), indexes(length);
                guess(first, length, loops, indexes, comp);
                rearrange(first, length, indexes);
            }


            // Guess an possible index table that sorts the range [first, first + length) use the optimized guessing algorithm.
            // first - The beginning of the range.
            // length - The length of the range.
            // loops - The index table.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void optimized_guess(RandomAccessIterator first, int length, std::vector<int>& loops, Compare comp) {
                while (true) {
                    bool sorted = true;
                    for (int i = 0; i < length - 1; ++i) {
                        auto left = first + loops[i], right = first + loops[i + 1];
                        if (comp(*left, *right) || equals(*left, *right, comp) && loops[i] < loops[i + 1]) continue;
                        sorted = false;
                        break;
                    }
                    if (sorted) break;

                    for (int pos = 0; pos < length; pos++) {
                        if (loops[pos] < length - 1) {
                            loops[pos]++;
                            break;
                        }
                        else loops[pos] = 0;
                    }
                }
            }
            

            // Sort the range [first, last) use the optimized guess sorting algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void optimized_guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> loops(length);
                optimized_guess(first, length, loops, comp);
                rearrange(first, length, loops);
            }


            // Guess an possible index table that sorts the range [first, first + length) use the random guessing algorithm.
            // first - The beginning of the range.
            // length - The length of the range.
            // loops - The index table.
            // rng - The random number generator.
            // comp - The comparison function.
            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void random_guess(RandomAccessIterator first, int length, std::vector<int>& loops, RandomNumberGenerator& rng, Compare comp) {
                while (true) {
                    bool sorted = true;
                    for (int i = 0; i < length - 1; i++) {
                        auto left = first + loops[i], right = first + loops[i + 1];
                        if (comp(*left, *right) || equals(*left, *right, comp) && loops[i] < loops[i + 1]) continue;
                        sorted = false;
                        break;
                    }
                    if (sorted) break;

                    for (int i = 0; i < length; i++) loops[i] = randint(0, length, rng);
                }
            }

            // Sort the range [first, last) use the random guessing algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // rng - The random number generator.
            // comp - The comparison function.
            template<class RandomAccessIterator, class RandomNumberGenerator, class Compare>
            void random_guess_sort(RandomAccessIterator first, RandomAccessIterator last, RandomNumberGenerator& rng, Compare comp) {
                int length = last - first;
                std::vector<int> loops(length);
                random_guess(first, length, loops, rng, comp);
                rearrange(first, length, loops);
            }


            // Guess an possible index table that sorts the range [first, first + length) use the smart guessing algorithm.
            // first - The beginning of the range.
            // length - The length of the range.
            // loops - The index table.
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void smart_guess(RandomAccessIterator first, int length, std::vector<int>& loops, Compare comp) {
                while (true) {
                    bool sorted = true;
                    int i = length - 2;
                    for (; i >= 0; i--) {
                        auto left = first + loops[i], right = first + loops[i + 1];
                        if (comp(*left, *right) || equals(*left, *right, comp) && loops[i] < loops[i + 1]) continue;
                        sorted = false;
                        break;
                    }
                    if (sorted) break;

                    for (int pos = 0; pos < length; pos++) {
                        if (pos >= i && loops[pos] < length - 1) {
                            loops[pos]++;
                            break;
                        }
                        else loops[pos] = 0;
                    }
                }
            }


            // Sort the range [first, last) use the smart guessing algorithm.
            // first - The beginning of the range.
            // last - The end of the range. (exclusive)
            // comp - The comparison function.
            template<class RandomAccessIterator, class Compare>
            void smart_guess_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> loops(length);
                smart_guess(first, length, loops, comp);
                rearrange(first, length, loops);
            }
        }

        namespace hanoi_sorting {
            using sortcpp::_internal::equals;

            bool vaild_number_moves(int moves) {
                if (moves == 0) return true;
                if (moves % 2 == 0) return false;
                return vaild_number_moves(moves / 2);
            }

            int get_height(int moves_plus1) {
                if (moves_plus1 == 1) return 0;
                return get_height(moves_plus1 / 2) + 1;
            }

            template<class T, class Compare>
            bool end_con_met(std::stack<T> &stack2, int end_con, int moves, int target_moves, const T &target, Compare comp) {
                if (!vaild_number_moves(moves)) return false;
                switch (end_con) {
                case 1: return (stack2.empty() || !comp(stack2.top(), target));
                case 2: return moves == target_moves;
                case 3: return stack2.empty();
                default: return false;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            int move_from_main(std::stack<T>& stack, bool check_unsorted, RandomAccessIterator &sp, 
                               RandomAccessIterator &unsorted, RandomAccessIterator last, Compare comp) {
                int duplicates = 1;
                stack.push(*sp);
                sp++;

                bool end_on_length = (sp >= last) || (check_unsorted && sp >= unsorted);
                while (!end_on_length && equals(*sp, stack.top(), comp)) {
                    duplicates++;
                    stack.push(*sp);
                    sp++;
                    end_on_length = (sp >= last) || (check_unsorted && sp >= unsorted);
                }
                return duplicates;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void move_to_main(std::stack<T>& stack, RandomAccessIterator &sp, Compare comp) {
                sp--;
                *sp = stack.top();
                stack.pop();
                while (!stack.empty() && equals(*sp, stack.top(), comp)) {
                    sp--;
                    *sp = stack.top();
                    stack.pop();
                }
            }

            template<class T, class Compare>
            void move_between_stacks(std::stack<T>& from, std::stack<T>& to, Compare comp) {
                to.push(from.top());
                from.pop();

                while (!from.empty() && equals(from.top(), to.top(), comp)) {
                    to.push(from.top());
                    from.pop();
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            int hanoi(std::stack<T> &stack2, std::stack<T> &stack3, 
                      int start_stack, bool go_right, int end_con, int &target_moves, const T &target,
                      RandomAccessIterator& sp, RandomAccessIterator& unsorted, RandomAccessIterator last, Compare comp) {
                int moves = 0;
                int min_pole_loc = start_stack;

                if (!end_con_met(stack2, end_con, moves, target_moves, target, comp)) {
                    moves++;
                    switch (min_pole_loc) {
                    case 1:
                        if (go_right) {
                            move_from_main<T>(stack2, true, sp, unsorted, last, comp);
                            min_pole_loc = 2;
                        }
                        else {
                            move_from_main<T>(stack3, true, sp, unsorted, last, comp);
                            min_pole_loc = 3;
                        }
                        break;
                    case 2:
                        if (go_right) {
                            move_between_stacks<T>(stack2, stack3, comp);
                            min_pole_loc = 3;
                        }
                        else {
                            move_to_main<T>(stack2, sp, comp);
                            min_pole_loc = 1;
                        }
                        break;
                    case 3:
                        if (go_right) {
                            move_to_main<T>(stack3, sp, comp);
                            min_pole_loc = 1;
                        }
                        else {
                            move_between_stacks<T>(stack3, stack2, comp);
                            min_pole_loc = 2;
                        }
                        break;
                    }
                }

                while (!end_con_met(stack2, end_con, moves, target_moves, target, comp)) {
                    moves += 2;
                    switch (min_pole_loc) {
                    case 1:
                        if (!stack2.empty() && (stack3.empty() || comp(stack2.top(), stack3.top())))
                            move_between_stacks<T>(stack2, stack3, comp);
                        else
                            move_between_stacks<T>(stack3, stack2, comp);
                        if (go_right) {
                            move_from_main<T>(stack2, true, sp, unsorted, last, comp);
                            min_pole_loc = 2;
                        }
                        else {
                            move_from_main<T>(stack3, true, sp, unsorted, last, comp);
                            min_pole_loc = 3;
                        }
                        break;
                    case 2:
                        if (stack3.empty() || (sp < unsorted && comp(*sp, stack3.top())))
                            move_from_main<T>(stack3, true, sp, unsorted, last, comp);
                        else move_to_main<T>(stack3, sp, comp);
                        if (go_right) {
                            move_between_stacks<T>(stack2, stack3, comp);
                            min_pole_loc = 3;
                        }
                        else {
                            move_to_main<T>(stack2, sp, comp);
                            min_pole_loc = 1;
                        }
                        break;
                    case 3:
                        if (stack2.empty() || (sp < unsorted && comp(*sp, stack2.top())))
                            move_from_main<T>(stack2, true, sp, unsorted, last, comp);
                        else move_to_main<T>(stack2, sp, comp);
                        if (go_right) {
                            move_to_main<T>(stack3, sp, comp);
                            min_pole_loc = 1;
                        }
                        else {
                            move_between_stacks<T>(stack3, stack2, comp);
                            min_pole_loc = 2;
                        }
                        break;
                    }
                }

                return moves;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void remove_from_main_stack(std::stack<T>& stack2, std::stack<T>& stack3, int& target_moves, T& target,
                                        RandomAccessIterator& sp, RandomAccessIterator& unsorted, RandomAccessIterator last, bool stable, Compare comp) {
                target = *sp;

                int moves = hanoi(stack2, stack3, 2, true, 1, target_moves, target, sp, unsorted, last, comp);
                int height = get_height(moves + 1);
                target_moves = moves;

                if (!stable) {
                    bool even_height = height % 2 == 0;
                    if (even_height) hanoi(stack2, stack3, 1, true, 2, target_moves, target, sp, unsorted, last, comp);
                    unsorted += move_from_main(stack2, false, sp, unsorted, last, comp);
                    hanoi(stack2, stack3, 3, even_height, 2, target_moves, target, sp, unsorted, last, comp);
                }
                else {
                    if (height % 2 == 1) {
                        unsorted += move_from_main<T>(stack2, false, sp, unsorted, last, comp);
                        hanoi(stack2, stack3, 3, false, 2, target_moves, target, sp, unsorted, last, comp);
                    }
                    else {
                        hanoi(stack2, stack3, 1, false, 2, target_moves, target, sp, unsorted, last, comp);
                        hanoi(stack2, stack3, 2, false, 2, target_moves, target, sp, unsorted, last, comp);
                        unsorted += move_from_main<T>(stack2, false, sp, unsorted, last, comp);
                        hanoi(stack2, stack3, 3, true, 2, target_moves, target, sp, unsorted, last, comp);
                    }
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void return_to_main_stack(std::stack<T>& stack2, std::stack<T>& stack3, int& target_moves, T& target,
                                      RandomAccessIterator& sp, RandomAccessIterator& unsorted, RandomAccessIterator last, bool stable, Compare comp) {
                int moves = hanoi(stack2, stack3, 2, true, 3, target_moves, target, sp, unsorted, last, comp);
                int height = get_height(moves + 1);

                if (!stable) {
                    if (height % 2 == 1) {
                        target_moves = moves;
                        hanoi(stack2, stack3, 3, true, 2, target_moves, target, sp, unsorted, last, comp);
                    }
                }
                else {
                    bool odd_height = height % 2 == 1;
                    target_moves = moves;

                    if (odd_height) hanoi(stack2, stack3, 3, false, 2, target_moves, target, sp, unsorted, last, comp);
                    else {
                        target_moves /= 2;
                        hanoi(stack2, stack3, 1, false, 2, target_moves, target, sp, unsorted, last, comp);
                        hanoi(stack2, stack3, 3, false, 2, target_moves, target, sp, unsorted, last, comp);
                        odd_height = true;
                    }

                    while (!stack2.empty()) {
                        target_moves /= 2;
                        hanoi(stack2, stack3, 2, !odd_height, 2, target_moves, target, sp, unsorted, last, comp);
                        move_to_main<T>(stack2, sp, comp);
                        hanoi(stack2, stack3, 3, odd_height, 2, target_moves, target, sp, unsorted, last, comp);
                        odd_height = !odd_height;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void hanoi_sort(RandomAccessIterator first, RandomAccessIterator last, bool stable, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                std::stack<T> stack2, stack3;

                auto target = T{};
                int target_moves = 0;
                auto sp = first, unsorted = first;

                while (unsorted < last) remove_from_main_stack<T>(stack2, stack3, target_moves, target, sp, unsorted, last, stable, comp);
                return_to_main_stack<T>(stack2, stack3, target_moves, target, sp, unsorted, last, stable, comp);
            }
        }

        namespace napoleon_sorting {
            using sortcpp::_internal::equals;

            template<class RandomAccessIterator, class Compare>
            void tilsit(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                int length = end - start + 1;
                auto mid = start + length / 2;
                for (auto i = start, j = end; i < mid; i++, j--)
                    if (comp(*j, *i)) std::iter_swap(i, j);
            }

            template<class T, class RandomAccessIterator, class Compare>
            RandomAccessIterator look_east(const T& prev, RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                for (auto i = start; i <= end; i++)
                    if (comp(*i, prev)) return i;
                return start - 1;
            }

            template<class T, class RandomAccessIterator, class Compare>
            RandomAccessIterator look_west(const T& prev, RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                for (auto i = start; i >= end; i--)
                    if (comp(prev, *i)) return i;
                return start + 1;
            }

            template<class T, class RandomAccessIterator, class Compare>
            RandomAccessIterator recruit(const T& identical_to, RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                if (start < end) {
                    for (auto i = start; i <= end; i++)
                        if (equals(identical_to, *i, comp)) return i - 1;
                    return end;
                }

                for (auto i = start; i >= end; i--)
                    if (equals(identical_to, *i, comp)) return i + 1;
                return end;
            }

            template<class RandomAccessIterator>
            void attack(RandomAccessIterator start_a, RandomAccessIterator start_b, int swaps_left) {
                while (swaps_left != 0) {
                    std::iter_swap(start_a, start_b);
                    start_a++;
                    start_b++;
                    swaps_left--;
                }
            }

            template<class RandomAccessIterator>
            void march(RandomAccessIterator index, int len1, int len2) {
                while (len1 != 0 && len2 != 0) {
                    if (len1 <= len2) {
                        attack(index, index + len1, len1);
                        index += len1;
                        len2 -= len1;
                    }
                    else {
                        attack(index + len1, index + len1 - len2, len2);
                        len1 -= len2;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void conquer(RandomAccessIterator index, RandomAccessIterator target, Compare comp) {
                int block_size = 1;
                while (index + block_size - 1 < target) {
                    auto march_to = recruit(*index, index + block_size, target, comp);
                    int space_between = march_to - index - block_size + 1;
                    march(index, block_size, space_between);
                    index += space_between;
                    block_size++;
                }

                while (index - block_size + 1 > target) {
                    auto march_to = recruit(*index, index - block_size, target, comp);
                    int space_between = index - block_size - march_to + 1;
                    march(march_to, space_between, block_size);
                    index -= space_between;
                    block_size++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void napoleon(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto lo = start, hi = end;
                auto prev = *start;
                bool go_smaller = true;
                bool none_found = false;

                while (hi > lo) {
                    if (go_smaller) {
                        auto next = look_east(prev, lo + 1, hi, comp);
                        if (next == lo) {
                            if (none_found) {
                                lo++;
                                hi--;
                            }
                            prev = (none_found ? *lo : *hi);
                            go_smaller = none_found;
                            none_found = !none_found;
                        }
                        else if (next == hi) {
                            prev = *hi;
                            conquer(next, hi, comp);
                            none_found = false;
                        }
                        else {
                            prev = *next;
                            conquer(next, hi, comp);
                            none_found = false;
                            go_smaller = false;
                        }
                    }
                    else {
                        auto next = look_west(prev, hi - 1, lo, comp);
                        if (next == hi) {
                            if (none_found) {
                                lo++;
                                hi--;
                            }
                            prev = *lo;
                            go_smaller = true;
                            none_found = !none_found;
                        }
                        else if (next == lo) {
                            prev = *lo;
                            conquer(lo, hi, comp);
                            none_found = false;
                        }
                        else {
                            prev = *next;
                            conquer(next, lo, comp);
                            none_found = false;
                            go_smaller = true;
                        }
                    }
                }
            }
        }

        namespace fire_sorting {
            template<class T, class Compare>
            bool twist_comp(const T& a, const T& b, int twist, Compare comp) {
                if (twist == -1) return comp(a, b);
                return comp(b, a);
            }

            template<class RandomAccessIterator, class Compare>
            void fire_sort(RandomAccessIterator first, RandomAccessIterator last, int increment, Compare comp) {
                auto i = first + 1, testi = first + 1;
                int twist_check = 0, twist_wait = 0, twist = -1;
                bool test_pass = false, test_reverse = false, any_swaps = false;

                while (!test_pass) {
                    if (twist_wait < 1) {
                        twist_check += increment;
                        twist_wait = twist_check;
                        twist *= -1;
                    }
                    else twist_wait--;

                    any_swaps = false;
                    while (i + 1 <= last && i >= first + 1 && !any_swaps) {
                        if (twist_comp(*(i - 1), *i, twist, comp)) {
                            std::iter_swap(i - 1, i);
                            i -= twist;
                            any_swaps = true;
                        }
                        else i += twist;
                    }

                    if (i < first + 1) {
                        i = last - 1;
                        testi = first + 1;
                        test_pass = true;
                        while (testi != last && test_pass) {
                            if (!comp(*testi, *(testi - 1))) testi++;
                            else {
                                test_pass = false;
                                testi = first + 1;
                                test_reverse = true;
                                while (testi != last && test_reverse) {
                                    if (!comp(*(testi - 1), *testi)) testi++;
                                    else test_reverse = false;
                                }
                            }
                        }

                        if (test_reverse) {
                            i = first + 1;
                            twist_wait = 0;
                        }
                    }

                    if (i + 1 > last) {
                        i = first + 1;
                        testi = first + 1;
                        test_pass = true;
                        while (testi != last && test_pass) {
                            if (!comp(*testi, *(testi - 1))) testi++;
                            else {
                                test_pass = false;
                                testi = first + 1;
                                test_reverse = true;
                                while (testi != last && test_reverse) {
                                    if (!comp(*(testi - 1), *testi)) testi++;
                                    else test_reverse = false;
                                }
                            }
                        }

                        if (test_reverse) {
                            i = last - 1;
                            twist_wait = 0;
                        }
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void merry_go_round_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto im = first + 1;
                auto very = first + 1;
                auto dizzy = first + 1;

                bool can_someone_please = false;
                bool stop_the_ride = false;

                while (!stop_the_ride) {
                    im = very;
                    can_someone_please = false;
                    while (im + 1 <= last) {
                        if (comp(*im, *(im - 1))) {
                            can_someone_please = true;
                            dizzy = im;
                            while (dizzy + 1 <= last) {
                                std::iter_swap(dizzy - 1, dizzy);
                                dizzy += 2;
                            }
                            if (im > first + 1) im--;
                        }
                        else im += 2;
                    }

                    if (!can_someone_please) {
                        very = first + 1;
                        stop_the_ride = true;
                        while (very != last && stop_the_ride) {
                            if (!comp(*very, *(very - 1))) very++;
                            else stop_the_ride = false;
                        }
                    }
                }
            }
            template<class RandomAccessIterator, class Compare>
            void x_pattern_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int gap = length;
                auto i = first + 1;
                auto xleft = first + 1, xright = first + 1;
                bool any_swaps = false, test_pass = false;

                while (!test_pass) {
                    any_swaps = false;
                    i = first + 1;

                    while (i - 1 + gap < last) {
                        if (comp(*(i - 1 + gap), *(i - 1))) {
                            std::iter_swap(i - 1 + gap, i - 1);
                            any_swaps = true;
                            xleft = i + 1;
                            xright = i + gap - 1;
                            if (gap != 1) {
                                for (int r = 0; r < gap - 1; r++) {
                                    if (comp(*(xright - 1), *(xleft - 1))) std::iter_swap(xleft - 1, xright - 1);
                                    xleft++;
                                    xright--;
                                }
                            }
                        }
                        i++;
                    }

                    if (gap == 1 && !any_swaps) test_pass = true;
                    else if (gap != 1 && !any_swaps) gap--;
                }
            }
        }

        namespace playground_sorting {
            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator select_lowest(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto lowest_index = first;
                for (auto j = first + 1; j < last; j++)
                    if (comp(*j, *lowest_index)) lowest_index = j;
                return lowest_index;
            }

            template<class RandomAccessIterator, class Compare>
            RandomAccessIterator select_next(RandomAccessIterator first, RandomAccessIterator last, RandomAccessIterator target, Compare comp) {
                auto lowest_high = last;
                auto right = first;

                while (right < last) {
                    if (comp(*target, *right)) {
                        if (lowest_high == last) lowest_high = right;
                        else if (comp(*right, *lowest_high)) lowest_high = right;
                    }
                    right++;
                }
                return lowest_high;
            }

            template<class RandomAccessIterator>
            void chase(RandomAccessIterator first, RandomAccessIterator item, RandomAccessIterator target) {
                int dir = 0;
                auto chase = first;
                if (std::abs(target - item) != 1) {
                    if (target - item > 0) dir = 1;
                    else dir = -1;
                    chase = item;
                    while (std::abs(target - chase) != 1) {
                        std::iter_swap(chase, chase + dir);
                        chase += dir;
                    }
                }
            }

            template<class RandomAccessIterator>
            void quit(RandomAccessIterator bound, RandomAccessIterator item) {
                auto pull = item;
                while (pull + 1 < bound) {
                    std::iter_swap(pull, pull + 1);
                    pull++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void playground_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto bound = last, last_target = first, target = first;
                while (bound > first + 1) {
                    last_target = select_lowest(first, bound, comp);
                    target = last_target;

                    while (target != bound) {
                        target = select_next(first, bound, last_target, comp);
                        if (target != bound) {
                            chase(first, last_target, target);
                            last_target = target;
                        }
                        else quit(bound, last_target);
                    }

                    bound--;
                }
            }
        }

        namespace gnome_weave_sorting {
            template<class RandomAccessIterator, class Compare>
            void gnome_weave_sort(RandomAccessIterator first, RandomAccessIterator last, bool high, Compare comp) {
                int length = last - first;
                int gap = length;

                auto icheck = first + 1, i = first + 1, bound = first + 1;

                int prime_test_i = 2;
                int prime_test_running = 1;
                bool prime_test = false, final_gap = false;

                while (!final_gap) {
                    i = icheck;
                    bound = icheck;

                    while (i - 1 + gap < last) {
                        if (comp(*(i - 1 + gap), *(i - 1))) {
                            std::iter_swap(i - 1, i - 1 + gap);
                            if (i - gap > first) i -= gap;
                        }
                        else {
                            bound += gap;
                            i = bound;
                        }
                    }

                    if (gap == 1) final_gap = true;
                    if (icheck + 1 > gap + first && !final_gap) {
                        prime_test_running = gap;
                        while (high ? prime_test_running != 1 : prime_test_running == gap) {
                            prime_test = false;
                            prime_test_i = 2;

                            while (!prime_test) {
                                if (prime_test_running % prime_test_i == 0) {
                                    prime_test_running /= prime_test_i;
                                    prime_test = true;
                                }
                                else prime_test_i++;
                            }
                        }
                        gap /= prime_test_i;
                        icheck = first + 1;
                    }
                    else icheck++;
                }
            }
        }

        namespace odd_even_weave_sorting {
            template<class RandomAccessIterator, class Compare>
            void odd_even_weave_sort(RandomAccessIterator first, RandomAccessIterator last, bool high, Compare comp) {
                int length = last - first;
                auto check = first + 1, last_bound = first + 2, i = first + 1, boundi = first + 1;
                int last_move = length, move = length, no_swaps_for = 0;
                int gap = length, prime_testi = 2, prime_test_running = 1;
                bool prime_test = false, any_swaps = false, test_pass = false;
                bool visual_aesthetic = true;

                while (!test_pass) {
                    i = check;
                    any_swaps = false;

                    while (i - 1 + gap < last) {
                        if (comp(*(i - 1 + gap), *(i - 1))) {
                            std::iter_swap(i - 1 + gap, i - 1);
                            any_swaps = true;
                        }
                        i += move;
                    }

                    if (!any_swaps && gap != 1) {
                        no_swaps_for++;
                        if (no_swaps_for == move) {
                            no_swaps_for = 0;
                            last_move = move;
                            prime_test_running = move;

                            if (!visual_aesthetic) {
                                while (prime_test_running != 1) {
                                    prime_test = false;
                                    prime_testi = 2;
                                    while (!prime_test) {
                                        if (prime_test_running % prime_testi == 0) {
                                            prime_test_running /= prime_testi;
                                            prime_test = true;
                                        }
                                        else prime_testi++;
                                    }
                                }
                                move /= prime_testi;
                            }
                            visual_aesthetic = false;

                            if (move != 1) {
                                prime_test_running = move;
                                while (high ? prime_test_running != 1 : prime_test_running == move) {
                                    prime_test = false;
                                    prime_testi = 2;
                                    while (!prime_test) {
                                        if (prime_test_running % prime_testi == 0) {
                                            prime_test_running /= prime_testi;
                                            prime_test = true;
                                        }
                                        else prime_testi++;
                                    }
                                }
                                gap = move / prime_testi;
                            }
                            else {
                                move = last_move;
                                gap = 1;
                            }
                            check = first;
                        }
                    }
                    else no_swaps_for = 0;

                    if (gap == 1) {
                        if (last_bound > first + 1) boundi = last_bound - 1;
                        else boundi = first + 1;

                        test_pass = true;
                        while (boundi < last && test_pass) {
                            if (comp(*boundi, *(boundi - 1))) {
                                test_pass = false;
                                last_bound = boundi;
                                check = boundi;
                            }
                            else boundi++;
                        }
                    }
                    else check = first + (check - first) % move + 1;
                }
            }
        }

        namespace zipper_sorting {
            int _log2(int x) {
                int n = 1;
                while (1 << n < x) n++;
                if (1 << n > x) n--;
                return n;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator bin_search(RandomAccessIterator a, RandomAccessIterator b, const T& value, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(value, *m)) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class RandomAccessIterator, class Compare>
            void ending(RandomAccessIterator first, RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                auto right = start;
                if (right < first + 1) right++;
                while (right < end) {
                    if (comp(*right, *(right - 1))) {
                        auto left = bin_search(first, right - 1, *right, comp);
                        while (left < right) {
                            std::iter_swap(left, right);
                            left++;
                        }
                    }
                    right++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_zipper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                auto i = first, head = first;
                int gap = length;
                int threshold = std::max(2 * _log2(length), (int)std::sqrt(length));

                while (gap > threshold) {
                    gap = 1;
                    i = (head > first + 1) ? head - 1 : first;
                    while (i + gap < last) {
                        if (comp(*(i + gap), *i)) {
                            std::iter_swap(i, i + gap);
                            if (gap == 1) head = i;
                            gap++;
                        }
                        else i++;
                    }
                }

                if (gap != 1) ending(first, head, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void zipper_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                auto i = first;
                int gap = 2;
                auto head = first;
                while (gap > 1) {
                    gap = 1;
                    i = head > first + 1 ? head - 1 : first;
                    while (i + gap < last) {
                        if (comp(*(i + gap), *i)) {
                            std::iter_swap(i, i + gap);
                            if (gap == 1) head = i;
                            gap++;
                        }
                        else i++;
                    }
                }
            }
        }

        namespace queue_sorting {
            template<class Container, class Compare>
            struct GComperator {
                Compare comp_;
                GComperator() {}
                GComperator(const Compare& comp) : comp_(comp) {}
                bool operator()(const Container& a, const Container& b) const {
                    return comp_(b.front(), a.front());
                }
            };

            template<class Container, class Compare>
            struct LComperator {
                Compare comp_;
                LComperator() {}
                LComperator(const Compare& comp) : comp_(comp) {}
                bool operator()(const Container& a, const Container& b) const {
                    return comp_(a.front(), b.front());
                }
            };

            template<class RandomAccessIterator, class Compare>
            void deque_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                using Deque = std::deque<T>;

                auto vleft = *first, vright = *first;

                std::vector<Deque> queues;
                queues.emplace_back();
                queues.back().push_back(vright);

                for (auto i = first + 1; i < last; i++) {
                    if (!comp(*i, vright)) {
                        vright = *i;
                        queues.back().push_back(*i);
                    }
                    else if (!comp(vleft, *i)) {
                        vleft = *i;
                        queues.back().push_front(*i);
                    }
                    else {
                        vleft = vright = *i;
                        queues.emplace_back();
                        queues.back().push_back(*i);
                    }
                }

                GComperator<Deque, Compare> comperator(comp);
                std::priority_queue<Deque, std::vector<Deque>, GComperator<Deque, Compare>> pq(comperator);
                for (auto& q : queues) pq.push(std::move(q));

                auto j = first;
                while (!pq.empty()) {
                    Deque cur = pq.top();
                    pq.pop();

                    *j = cur.front();
                    cur.pop_front();

                    if (!cur.empty()) pq.push(std::move(cur));
                    j++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void queue_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                using Queue = std::queue<T>;

                auto val = *first;

                std::vector<Queue> queues;
                queues.emplace_back();
                queues.back().push(val);

                for (auto i = first + 1; i < last; i++) {
                    if (!comp(*i, *(i - 1))) queues.back().push(*i);
                    else {
                        val = *i;
                        queues.emplace_back();
                        queues.back().push(*i);
                    }
                }

                GComperator<Queue, Compare> comperator(comp);
                std::priority_queue<Queue, std::vector<Queue>, GComperator<Queue, Compare>> pq(comperator);
                for (auto& q : queues) pq.push(std::move(q));

                auto j = first;
                while (!pq.empty()) {
                    Queue cur = pq.top();
                    pq.pop();

                    *j = cur.front();
                    cur.pop();

                    if (!cur.empty()) pq.push(std::move(cur));
                    j++;
                }
            }

        }

        namespace patience_sorting {
            template<class Container, class Compare>
            struct GComperator {
                Compare comp_;
                GComperator() {}
                GComperator(const Compare& comp) : comp_(comp) {}
                bool operator()(const Container& a, const Container& b) const {
                    return comp_(b.top(), a.top());
                }
            };

            template<class Container, class Compare>
            struct LComperator {
                Compare comp_;
                LComperator() {}
                LComperator(const Compare& comp) : comp_(comp) {}
                bool operator()(const Container& a, const Container& b) const {
                    return comp_(a.top(), b.top());
                }
            };

            template<class RandomAccessIterator, class Compare>
            void patience_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                using Pile = std::stack<T>;
                using sortcpp::_internal::equals;
                using sortcpp::_internal::binary_find;

                LComperator<Pile, Compare> less_comp(comp);
                GComperator<Pile, Compare> greater_comp(comp);

                std::vector<Pile> piles;
                for (auto x = first; x < last; x++) {
                    Pile new_pile;
                    new_pile.push(*x);

                    auto i = binary_find(piles.begin(), piles.end(), new_pile, less_comp);
                    if (i == piles.end()) piles.push_back(std::move(new_pile));
                    else (*i).push(*x);
                }


                std::priority_queue<Pile, std::vector<Pile>, GComperator<Pile, Compare>> heap(greater_comp);
                for (auto& pile : piles) heap.push(std::move(pile));
                for (auto c = first; c < last; c++) {
                    auto small_pile = heap.top();
                    heap.pop();

                    *c = small_pile.top();
                    small_pile.pop();

                    if (!small_pile.empty()) heap.push(std::move(small_pile));
                }
            }
        }


        namespace library_sorting {
            using sortcpp::_internal::equals;
            using sortcpp::_internal::bogo_sorting::rd;
            using sortcpp::_internal::bogo_sorting::randint;

            const int G = 15, R = 4;

            template<class RandomAccessIterator, class T>
            void shift_ext(RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, const T& max) {
                auto m1 = a + std::min(m - a, b - m);
                while (m > a) {
                    b--;
                    m--;
                    *b = *m;
                }
                while (a < m1) {
                    *a = max;
                    a++;
                }
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator left_block_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                int s = G + 1;

                while (a < b) {
                    auto m = a + (((b - a) / s) / 2) * s;
                    if (comp(*m, val)) a = m + s;
                    else b = m;
                }
                return a;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator right_block_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                int s = G + 1;

                while (a < b) {
                    auto m = a + (((b - a) / s) / 2) * s;
                    if (comp(val, *m)) b = m;
                    else a = m + s;
                }
                return a;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator loc_search(RandomAccessIterator a, RandomAccessIterator b, const T& max, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(*m, max)) a = m + 1;
                    else b = m;
                }
                return a;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator right_bin_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(val, *m)) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class RandomAccessIterator>
            void insert_to(RandomAccessIterator a, RandomAccessIterator b) {
                auto temp = *a;
                while (a > b) {
                    *a = *(a - 1);
                    a--;
                }
                *b = temp;
            }

            template<class RandomAccessIterator, class Compare>
            void binary_insertion(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                for (auto i = a + 1; i < b; i++) insert_to(i, right_bin_search(a, i, *i, comp));
            }

            template<class RandomAccessIterator, class Buffer, class T, class Compare>
            void retrieve(RandomAccessIterator i, Buffer tmp, int p_end, const T& max, Compare comp) {
                auto loc = i - 1;

                for (auto k = p_end - (G + 1) + tmp; k > G + tmp; ) {
                    auto m = loc_search(k - G, k, max, comp) - 1;
                    k -= G + 1;

                    while (m >= k) {
                        *loc = *m;
                        *m = max;
                        loc--;
                        m--;
                    }
                }

                auto m = loc_search(tmp, tmp + G, max, comp) - 1;
                while (m >= tmp) {
                    *loc = *m;
                    *m = max;
                    loc--;
                    m--;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void library_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                std::mt19937 gen(rd());

                const int G = 15, R = 4;
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                auto max = *first;
                for (auto i = first + 1; i < last; i++)
                    if (comp(*i, max)) max = *i;

                int length = last - first;
                std::vector<T> tmp(length * (G + 1) - 1, max);

                int s = length;
                while (s >= 32) s = (s - 1) / R + 1;

                auto i = first + s;
                int j = R * s, p_end = (s + 1) * (G + 1) + G;
                binary_insertion(first, i, comp);

                for (int k = 0; k < s; k++) tmp[k * (G + 1) + G] = *(first + k);
                for (; i < last; i++) {
                    if (i - first == j) {
                        retrieve(i, tmp.begin(), p_end, max, comp);

                        s = i - first;
                        p_end = (s + 1) * (G + 1) + G;
                        j *= R;

                        for (int k = 0; k < s; k++) tmp[k * (G + 1) + G] = *(first + k);
                    }

                    auto b_loc = left_block_search(tmp.begin() + G, tmp.begin() + p_end - (G + 1), *i, comp);

                    if (equals(*i, *b_loc, comp)) {
                        auto eq_end = right_block_search(b_loc + (G + 1), tmp.begin() + p_end - (G + 1), *i, comp);
                        b_loc += randint(0, (int)((eq_end - b_loc) / (G + 1)), gen) * (G + 1);
                    }

                    auto loc = loc_search(b_loc - G, b_loc, max, comp);

                    if (loc == b_loc) {
                        do b_loc += G + 1;
                        while (b_loc < tmp.begin() + p_end && loc_search(b_loc - G, b_loc, max, comp) == b_loc);

                        if (b_loc == tmp.begin() + p_end) {
                            retrieve(i, tmp.begin(), p_end, max, comp);

                            s = i - first;
                            p_end = (s + 1) * (G + 1) + G;
                            j = R * s;

                            for (int k = 0; k < s; k++) tmp[k * (G + 1) + G] = *(first + k);
                        }
                        else {
                            auto rot_p = loc_search(b_loc - G, b_loc, max, comp);
                            auto rot_s = b_loc - (std::max(rot_p, b_loc - G / 2) - tmp.begin());
                            shift_ext(loc - rot_s + tmp.begin(), b_loc - rot_s + tmp.begin(), b_loc, max);
                        }
                        i--;
                    }
                    else {
                        *loc = *i;
                        insert_to(loc, right_bin_search(b_loc - G, loc, *loc, comp));
                    }
                }
                retrieve(last, tmp.begin(), p_end, max, comp);
            }
        }

        namespace simplified_library_sorting {
            using sortcpp::_internal::binary_insertion_sorting::binary_insertion_sort;
            const int R = 4;
            
            int get_min_level(int n) {
                while (n >= 32) n = (n - 1) / R + 1;
                return n;
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator _binary_search(RandomAccessIterator a, RandomAccessIterator b, const T& val, Compare comp) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (comp(val, *m)) b = m;
                    else a = m + 1;
                }
                return a;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void rebalance(RandomAccessIterator array, int m, int b, 
                           std::vector<T> &temp, std::vector<int> &cnts, std::vector<int> &locs, Compare comp) {
                for (int i = 0; i < m; i++) cnts[i + 1] += cnts[i] + 1;

                for (int i = m, j = 0; i < b; i++, j++) {
                    temp[cnts[locs[j]]] = *(array + i);
                    cnts[locs[j]]++;
                }

                for (int i = 0; i < m; i++) {
                    temp[cnts[i]] = *(array + i);
                    cnts[i]++;
                }

                for (int i = 0; i < b; i++) *(array + i) = temp[i];
                binary_insertion_sort(array, array + cnts[0] - 1, comp);
                for (int i = 0; i < m - 1; i++) binary_insertion_sort(array + cnts[i], array + cnts[i + 1] - 1, comp);
                binary_insertion_sort(array + cnts[m - 1], array + cnts[m], comp);

                for (int i = 0; i < m + 2; i++) cnts[i] = 0;
            }

            template<class RandomAccessIterator, class Compare>
            void library_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

                int length = last - first;
                int j = get_min_level(length);
                binary_insertion_sort(first, first + j, comp);

                int max_level = j;
                for (; max_level * R < length; max_level *= R);

                std::vector<T> temp(length);
                std::vector<int> cnts(max_level + 2), locs(length - max_level);

                for (int i = j, k = 0; i < length; i++) {
                    if (R * j == i) {
                        rebalance<T>(first, j, i, temp, cnts, locs, comp);
                        j = i;
                        k = 0;
                    }

                    int loc = _binary_search(first, first + j, *(first + i), comp) - first;
                    cnts[loc + 1]++;
                    locs[k++] = loc;
                }
                rebalance<T>(first, j, length, temp, cnts, locs, comp);
            }
        }

        namespace concurrent_sorting {
            using sortcpp::_internal::insertion_sorting::insertion_sort;
            using sortcpp::_internal::circle_sorting::circle_sort_routine;

            template<class RandomAccessIterator, class Compare>
            void bitonic_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int i, j, k;

                for (k = 2; k < length * 2; k *= 2) {
                    bool m = (((length + (k - 1)) / k) % 2) != 0;

                    for (j = k >> 1; j > 0; j >>= 1)
                        for (i = 0; i < length; i++) {
                            int ij = i ^ j;
                            if (ij > i && ij < length) {
                                if ((((i & k) == 0) == m) && comp(*(first + ij), *(first + i)))
                                    std::iter_swap(first + i, first + ij);
                                if ((((i & k) != 0) == m) && comp(*(first + i), *(first + ij)))
                                    std::iter_swap(first + i, first + ij);
                            }
                        }

                }
            }

            int hyperfloor(int n) {
                int k = 1;
                while (k < n) k <<= 1;
                return k >> 1;
            }

            template<class RandomAccessIterator, class Compare>
            void compare(RandomAccessIterator i, RandomAccessIterator j, bool dir, Compare comp) {
                if (comp(*j, *i) == dir) std::iter_swap(i, j);
            }

            template<class RandomAccessIterator, class Compare>
            void bitonic_merge(RandomAccessIterator lo, int n, bool dir, Compare comp) {
                if (n > 1) {
                    int m = hyperfloor(n);

                    for (auto i = lo; i < lo + n - m; i++) compare(i, i + m, dir, comp);

                    bitonic_merge(lo, m, dir, comp);
                    bitonic_merge(lo + m, n - m, dir, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void bitonic_sort_rec(RandomAccessIterator lo, int n, bool dir, Compare comp) {
                if (n > 1) {
                    int m = n / 2;
                    bitonic_sort_rec(lo, m, !dir, comp);
                    bitonic_sort_rec(lo + m, n - m, dir, comp);
                    bitonic_merge(lo, n, dir, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void bitonic_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bitonic_sort_rec(first, last - first, true, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void comp_swap(RandomAccessIterator a, RandomAccessIterator b, RandomAccessIterator last, Compare comp) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            template<class RandomAccessIterator, class Compare>
            void range_comp(RandomAccessIterator a, RandomAccessIterator b, RandomAccessIterator last, int offset, Compare comp) {
                int half = (b - a) / 2;
                auto m = a + half;
                a += offset;

                for (int i = 0; i < half - offset; i++)
                    if ((i & ~offset) == i) comp_swap(a + i, m + i, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void bose_nelson_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int size = 1 << (int)(std::ceil(std::log2(last - first)));
                for (int k = 2; k <= size; k *= 2)
                    for (int j = 0; j < k / 2; j++)
                        for (auto i = first; i + j < last; i += k)
                            range_comp(i, i + k, last, j, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void bose_nelson_merge(RandomAccessIterator start1, int len1, 
                                   RandomAccessIterator start2, int len2, RandomAccessIterator last, Compare comp) {
                if (len1 == 1 && len2 == 1) comp_swap(start1, start2, last, comp);
                else if (len1 == 1 && len2 == 2) {
                    comp_swap(start1, start2 + 1, last, comp);
                    comp_swap(start1, start2, last, comp);
                }
                else if (len1 == 2 && len2 == 1) {
                    comp_swap(start1, start2, last, comp);
                    comp_swap(start1 + 1, start2, last, comp);
                }
                else {
                    int mid1 = len1 / 2;
                    int mid2 = (len1 & 1) ? (len2 / 2) : ((len2 + 1) / 2);
                    bose_nelson_merge(start1, mid1, start2, mid2, last, comp);
                    bose_nelson_merge(start1 + mid1, len1 - mid1, start2 + mid2, len2 - mid2, last, comp);
                    bose_nelson_merge(start1 + mid1, len1 - mid1, start2, mid2, last, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void bose_nelson_rec(RandomAccessIterator start, int length, RandomAccessIterator last, Compare comp) {
                if (length > 1) {
                    int mid = length / 2;
                    bose_nelson_rec(start, mid, last, comp);
                    bose_nelson_rec(start + mid, length - mid, last, comp);
                    bose_nelson_merge(start, mid, start + mid, length - mid, last, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void bose_nelson_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                bose_nelson_rec(first, last - first, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void crease_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int max = 1;
                for (; max * 2 < length; max *= 2);

                int next = max;
                while (next > 0) {
                    for (auto i = first; i + 1 < last; i += 2) comp_swap(i, i + 1, last, comp);
                    for (int j = max; j >= next && j > 1; j /= 2)
                        for (auto i = first + 1; i + j - 1 < last; i += 2)
                            comp_swap(i, i + j - 1, last, comp);
                    next /= 2;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void diamond_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int n = 1;
                for (; n < length; n *= 2);

                int m = 4;
                for (; m <= n; m *= 2) {
                    for (int k = 0; k < m / 2; k++) {
                        int cnt = (k <= m / 4) ? k : m / 2 - k;
                        for (auto j = first; j < last; j += m)
                            if (j + cnt + 1 < last)
                                for (auto i = j + cnt; i + 1 < std::min(last, j + m - cnt); i += 2)
                                    comp_swap(i, i + 1, last, comp);
                    }
                }

                m /= 2;
                for (int k = 0; k <= m / 2; k++)
                    for (auto i = first + k; i + 1 < std::min(last, first + m - k); i += 2)
                        comp_swap(i, i + 1, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void diamond_sort_rec(RandomAccessIterator start, RandomAccessIterator stop, bool merge, Compare comp) {
                if (stop - start == 2) {
                    if (comp(*(stop - 1), *start)) std::iter_swap(start, stop - 1);
                    return;
                }
                if (stop - start >= 3) {
                    int len = stop - start;
                    int p = len / 4, q = (int)(len * 0.75);

                    auto mid = (stop - start) / 2 + start;
                    if (merge) {
                        diamond_sort_rec(start, mid, true, comp);
                        diamond_sort_rec(mid, stop, true, comp);
                    }
                    diamond_sort_rec(p + start, q + start, false, comp);
                    diamond_sort_rec(start, mid, false, comp);
                    diamond_sort_rec(mid, stop, false, comp);
                    diamond_sort_rec(p + start, q + start, false, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void diamond_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                diamond_sort_rec(first, last, true, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void halver(RandomAccessIterator low, RandomAccessIterator high, RandomAccessIterator last, Compare comp) {
                while (low < high) {
                    comp_swap(low, high, last, comp);
                    low++;
                    high--;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void fold_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int ceil_log = 1;
                for (; (1 << ceil_log) < length; ceil_log++);

                int size = 1 << ceil_log;
                for (int k = size >> 1; k > 0; k >>= 1)
                    for (int i = size; i >= k; i >>= 1)
                        for (auto j = first; j < last; j += i)
                            halver(j, j + i - 1, last, comp);
            }

            struct matrix_shape {
                int width;
                bool unbalanced, insert_last;
                matrix_shape() {}
                matrix_shape(int width, int height, bool insert_last) {
                    this->width = width;
                    this->unbalanced = (width == 1) ^ (height == 1);
                    this->insert_last = insert_last || unbalanced;
                }
            };

            template<class RandomAccessIterator>
            void gap_reverse(RandomAccessIterator start, RandomAccessIterator end, int gap) {
                for (auto i = start, j = end; i < j; i += gap, j -= gap) std::iter_swap(i, j - gap);
            }

            template<class T, class Compare>
            bool vcompare(const T& a, const T& b, bool dir, Compare comp) {
                return dir ? comp(a, b) : comp(b, a);
            }

            template<class RandomAccessIterator, class Compare>
            bool insert_last(RandomAccessIterator a, RandomAccessIterator b, int gap, bool dir, Compare comp) {
                bool did = false;
                auto key = *b;
                auto j = b - gap;

                while (j >= a && vcompare(key, *j, dir, comp)) {
                    *(j + gap) = *j;
                    did = true;
                    j -= gap;
                }
                *(j + gap) = key;
                return did;
            }


            matrix_shape get_matrix_dims(int len) {
                int dim = (int)std::sqrt(len);
                bool insert_last = false;
                if (dim * dim == len - 1) insert_last = true;
                for (; len % dim != 0; dim--);
                return matrix_shape(dim, len / dim, insert_last);
            }

            template<class RandomAccessIterator, class Compare>
            bool matrix_sort_rec(RandomAccessIterator start, RandomAccessIterator end, int gap, bool dir, Compare comp) {
                bool did = false;
                int length = (end - start) / gap;
                if (length < 2) return false;
                else if (length <= 16) {
                    did = false;
                    for (auto i = start; i < end; i += gap)
                        did = insert_last(start, i, gap, dir, comp) || did;
                }
                else {
                    bool newdid;
                    auto shape = get_matrix_dims(length);

                    if (shape.insert_last) {
                        bool did1 = matrix_sort_rec(start, end - gap, gap, dir, comp);
                        bool did2 = insert_last(start, end - gap, gap, dir, comp);
                        return did1 || did2;
                    }

                    int width = shape.width;
                    for (auto i = start + width * gap; i < end; i += 2 * width * gap)
                        gap_reverse(i, i + width * gap, gap);

                    did = false;
                    do {
                        newdid = false;

                        bool curdir = dir;
                        for (auto i = start; i < end; i += width * gap) {
                            newdid = matrix_sort_rec(i, i + width * gap, gap, curdir, comp) || newdid;
                            did = did || newdid;
                            curdir = !curdir;
                        }

                        newdid = false;
                        for (int i = 0; i < width; i++) {
                            newdid = matrix_sort_rec(start + i * gap, end + i * gap, gap * width, dir, comp) || newdid;
                            did = did || newdid;
                        }

                    } while (newdid);

                    for (auto i = start + width * gap; i < end; i += 2 * width * gap)
                        gap_reverse(i, i + width * gap, gap);
                }

                return did;
            }

            template<class RandomAccessIterator, class Compare>
            void matrix_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                matrix_sort_rec(first, last, 1, true, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void odd_even_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;

                for (int p = 1; p < length; p += p)
                    for (int k = p; k > 0; k /= 2)
                        for (int j = k % p; j + k < length; j += k + k)
                            for (int i = 0; i < k; i++)
                                if (i + j + k < length)
                                    if ((i + j) / (p + p) == (i + j + k) / (p + p))
                                        if (comp(*(first + i + j + k), *(first + i + j)))
                                            std::iter_swap(first + i + j, first + i + j + k);
            }

            template<class RandomAccessIterator, class Compare>
            void oddeven_merge(RandomAccessIterator lo, int m2, int n, int r, Compare comp) {
                int m = r * 2;
                if (m < n) {
                    if ((n / r) % 2 != 0) {
                        oddeven_merge(lo, (m2 + 1) / 2, n + r, m, comp);
                        oddeven_merge(lo + r, m2 / 2, n - r, m, comp);
                    }
                    else {
                        oddeven_merge(lo, (m2 + 1) / 2, n, m, comp);
                        oddeven_merge(lo + r, m2 / 2, n, m, comp);
                    }

                    if (m2 % 2 != 0) {
                        for (auto i = lo; i + r < lo + n; i += m) compare(i, i + r, true, comp);
                    }
                    else {
                        for (auto i = lo + r; i + r < lo + n; i += m) compare(i, i + r, true, comp);
                    }
                }
                else if (n > r) compare(lo, lo + r, true, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void oddeven_merge_sort_rec(RandomAccessIterator lo, int n, Compare comp) {
                if (n > 1) {
                    int m = n / 2;
                    oddeven_merge_sort_rec(lo, m, comp);
                    oddeven_merge_sort_rec(lo + m, n - m, comp);
                    oddeven_merge(lo, m, n, 1, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void odd_even_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                oddeven_merge_sort_rec(first, last - first, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int n = 1;
                for (; n < length; n <<= 1);

                for (int k = n >> 1; k > 0; k >>= 1)
                    for (auto j = first; j < last; j += k << 1)
                        for (int i = 0; i < k; i++)
                            comp_swap(j + i, j + i + k, last, comp);

                for (int k = 2; k < n; k <<= 1)
                    for (int m = k >> 1; m > 0; m >>= 1)
                        for (auto j = first; j < last; j += k << 1)
                            for (int p = m; p < ((k - m) << 1); p += m << 1)
                                for (int i = 0; i < m; i++)
                                    comp_swap(j + p + i, j + p + m + i, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_merge(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator b, RandomAccessIterator last, Compare comp) {
                auto m = a + (b - a) / 2;
                auto m1 = a + (m - a) / 2;
                int g = m - m1;

                for (int i = 0; m1 + i < m; i++)
                    for (int j = m1 - first, k = g; k > 0; k >>= 1, j -= k - (i & k))
                        comp_swap(first + j + i, first + j + i + k, last, comp);

                if (b - a > 4) pairwise_merge(first, m, b, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_merge_sort_rec(RandomAccessIterator first, RandomAccessIterator a, 
                                         RandomAccessIterator b, RandomAccessIterator last, Compare comp) {
                auto m = a + (b - a) / 2;
                for (auto i = a, j = m; i < m; i++, j++) comp_swap(i, j, last, comp);

                if (b - a > 2) {
                    pairwise_merge_sort_rec(first, a, m, last, comp);
                    pairwise_merge_sort_rec(first, m, b, last, comp);
                    pairwise_merge(first, a, b, last, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int n = 1;
                for (; n < length; n <<= 1);
                pairwise_merge_sort_rec(first, first, first + n, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void iterative_pairwise_loop(RandomAccessIterator start, RandomAccessIterator end, Compare comp) {
                int length = end - start;
                int a = 1, b = 0, c = 0, d = 0, e = 0;
                while (a < length) {
                    b = a;
                    c = 0;
                    while (b < length) {
                        comp_swap(start + b - a, start + b, end, comp);
                        c = (c + 1) % a;
                        b++;
                        if (c == 0) b += a;
                    }
                    a *= 2;
                }

                a /= 4;
                e = 1;

                while (a > 0) {
                    d = e;
                    while (d > 0) {
                        b = (d + 1) * a;
                        c = 0;
                        while (b < length) {
                            comp_swap(start + b - d * a, start + b, end, comp);
                            c = (c + 1) % a;
                            b++;
                            if (c == 0) b += a;
                        }
                        d /= 2;
                    }
                    a /= 2;
                    e = e * 2 + 1;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                iterative_pairwise_loop(first, last, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_sort_rec(RandomAccessIterator start, RandomAccessIterator end, RandomAccessIterator last, int gap, Compare comp) {
                if (start == end - gap) return;

                auto b = start + gap;
                while (b < end) {
                    comp_swap(b - gap, b, last, comp);
                    b += 2 * gap;
                }

                if (((end - start) / gap) % 2 == 0) {
                    pairwise_sort_rec(start, end, last, gap * 2, comp);
                    pairwise_sort_rec(start + gap, end + gap, last, gap * 2, comp);
                }
                else {
                    pairwise_sort_rec(start, end + gap, last, gap * 2, comp);
                    pairwise_sort_rec(start + gap, end, last, gap * 2, comp);
                }

                int a = 1;
                while (a < ((end - start) / gap)) a = a * 2 + 1;

                b = start + gap;
                while (b + gap < end) {
                    int c = a;
                    while (c > 1) {
                        c /= 2;
                        comp_swap(b, b + c * gap, last, comp);
                    }
                    b += 2 * gap;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void pairwise_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                pairwise_sort_rec(first, last, last, 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void weave_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int n = 1;
                for (; n < length; n *= 2);

                for (int i = 1; i < n; i *= 2)
                    for (int j = 1; j <= i; j *= 2)
                        for (int k = 0; k < n; k += n / j)
                            for (int d = n / i / 2, m = 0, l = n / j - d; l >= n / j / 2; l -= d)
                                for (int p = 0; p < d; p++, m++)
                                    comp_swap(first + k + m, first + k + l + p, last, comp);

            }

            template<class RandomAccessIterator, class Compare>
            void circle(RandomAccessIterator pos, int len, RandomAccessIterator last, int gap, Compare comp) {
                if (len < 2) return;

                for (int i = 0; 2 * i < (len - 1) * gap; i += gap)
                    comp_swap(pos + i, pos + (len - 1) * gap - i, last, comp);

                circle(pos, len / 2, last, gap, comp);
                if (pos + len * gap / 2 < last) circle(pos + len * gap / 2, len / 2, last, gap, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void weave_circle(RandomAccessIterator pos, int len, RandomAccessIterator last, int gap, Compare comp) {
                if (len < 2) return;

                weave_circle(pos, len / 2, last, 2 * gap, comp);
                weave_circle(pos + gap, len / 2, last, 2 * gap, comp);
                circle(pos, len, last, gap, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void weave_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int n = 1;
                for (; n < length; n *= 2);
                weave_circle(first, n, last, 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void apollyon_merge(RandomAccessIterator lo, int n, bool dir, Compare comp) {
                if (n > 1) {
                    int m = hyperfloor(n);

                    for (auto i = lo; i < lo + n - m; i++) compare(i, i + m, dir, comp);

                    apollyon_merge(lo, m, dir, comp);
                    apollyon_merge(lo + m, n - m, dir, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void apollyon_sort_rec(RandomAccessIterator lo, int n, bool dir, Compare comp) {
                if (n > 1) {
                    int m = n / 2;
                    apollyon_sort_rec(lo, m, !dir, comp);
                    apollyon_merge(lo, n, dir, comp);
                }
            }

            template<class RandomAccessIterator, class Compare>
            void apollyon_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                int threshold = 0, n = 1;
                for (; n < length; n *= 2, threshold++);

                threshold /= 2;
                int iterations = 0;

                apollyon_sort_rec(first, length, true, comp);

                while (circle_sort_routine(first, last - 1, last, false, comp)) {
                    iterations++;
                    if (iterations >= threshold) {
                        insertion_sort(first, last, comp);
                        break;
                    }
                }
            }

            template<class RandomAccessIterator, class Compare>
            void comp_swap(RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            }

            template<class RandomAccessIterator, class Compare>
            void comp_range(RandomAccessIterator a, RandomAccessIterator m, int s, Compare comp) {
                for (int i = s; a + i < m; i++) comp_swap(a + i, m + i, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void comp_range_extd(RandomAccessIterator a, RandomAccessIterator m, int p, Compare comp) {
                int l = m - a;

                if (l > p) {
                    auto i = a;
                    int j, d = l - p;

                    for (j = 0; j < d; j++, i++) comp_swap(i, i + p, comp);
                    for (j = 0; j < p - d; j++, i++) comp_swap(i, i + l, comp);
                    for (j = 0; j < d; j++, i++) comp_swap(i + d, i + l, comp);
                }
                else comp_range(a, m, 0, comp);
            }
            

            template<class RandomAccessIterator, class Compare>
            void optimized_oddeven_merge(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                int s = (b - a) & 1;
                a -= s;
                auto m = a + (b - a) / 2;
                comp_range(a, m, s, comp);

                int l = b - a;
                if (l < 4) return;

                int p;
                for (p = 1; 2 * p < l; p *= 2);

                while (p > 0) {
                    auto i = a + p;

                    while (i + 2 * p <= m) {
                        comp_range(i, i + p, 0, comp);
                        i += 2 * p;
                    }

                    comp_range_extd(i, m, p, comp);
                    i = first + 2 * (m - first) - (i - first);

                    while (i < b - p) {
                        comp_range(i, i + p, 0, comp);
                        i += 2 * p;
                    }
                    p /= 2;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_oddeven_merge_sort_loop(RandomAccessIterator first, RandomAccessIterator a, RandomAccessIterator b, Compare comp) {
                auto m = a + (b - a) / 2;

                if (m - a > 1) optimized_oddeven_merge_sort_loop(first, a, m, comp);
                if (b - m > 1) optimized_oddeven_merge_sort_loop(first, m, b, comp);

                optimized_oddeven_merge(first, a, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void optimized_oddeven_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                optimized_oddeven_merge_sort_loop(first, first, last, comp);
            }
        }

        namespace twin_sorting {
            template<class RandomAccessIterator, class Compare>
            int twin_swap(RandomAccessIterator left, int nmemb, Compare comp) {
                RandomAccessIterator index, start, end;

                index = left;
                end = left + nmemb - 2;

                while (index <= end) {
                    if (!comp(*(index + 1), *index)) {
                        index += 2;
                        continue;
                    }

                    start = index;
                    index += 2;

                    while (true) {
                        if (index > end) {
                            if (start == left) {
                                if (nmemb % 2 == 0 || comp(*index, *(index - 1))) {
                                    end = left + nmemb - 1;

                                    while (start < end)
                                        std::iter_swap(start++, end--);
                                    return 1;
                                }
                            }
                            break;
                        }

                        if (comp(*(index + 1), *index)) {
                            if (comp(*index, *(index - 1))) {
                                index += 2;
                                continue;
                            }

                            std::iter_swap(index, index + 1);
                        }
                        break;
                    }

                    end = index - 1;

                    while (start < end)
                        std::iter_swap(start++, end--);

                    end = left + nmemb - 2;

                    index += 2;
                }
                return 0;
            }

            template<class T, class RandomAccessIterator, class Compare>
            void tail_merge(RandomAccessIterator left, std::vector<T> &swap, int nmemb, int block, Compare comp) {
                int offset, s, c, c_max;
                RandomAccessIterator a, d, d_max, e;

                s = 0;

                while (block < nmemb) {
                    for (offset = 0; offset + block < nmemb; offset += block * 2) {
                        a = left + offset;
                        e = a + block - 1;

                        if (!comp(*(e + 1), *e)) continue;

                        if (offset + block * 2 <= nmemb) {
                            c_max = s + block;
                            d_max = a + block * 2;
                        }
                        else {
                            c_max = s + nmemb - (offset + block);
                            d_max = left + nmemb;
                        }

                        d = d_max - 1;

                        while (!comp(*d, *e)) {
                            d_max--;
                            d--;
                            c_max--;
                        }

                        c = s;
                        d = a + block;

                        while (c < c_max) {
                            swap[c++] = *d;
                            d++;
                        }
                        c--;

                        d = a + block - 1;
                        e = d_max - 1;

                        if (!comp(*(a + block), *a)) {
                            *(e--) = *(d--);

                            while (c >= s) {
                                while (comp(swap[c], *d)) *(e--) = *(d--);
                                *(e--) = swap[c--];
                            }
                        }

                        else {
                            *(e--) = *(d--);

                            while (d >= a) {
                                while (!comp(swap[c], *d)) *(e--) = swap[c--];
                                *(e--) = *(d--);
                            }

                            while (c >= s) *(e--) = swap[c--];
                        }
                    }
                    block *= 2;
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void twinsort_swap(RandomAccessIterator start, std::vector<T>& swap, int nmemb, Compare comp) {
                if (twin_swap(start, nmemb, comp) == 0) tail_merge(start, swap, nmemb, 2, comp);
            }

            template<class T, class RandomAccessIterator, class Compare>
            void _twinsort(RandomAccessIterator first, int nmemb, Compare comp) {
                if (twin_swap(first, nmemb, comp) == 0) {
                    std::vector<T> swap(nmemb / 2);
                    tail_merge(first, swap, nmemb, 2, comp);
                }
            }

            template<class T, class RandomAccessIterator, class Compare>
            void _tailsort(RandomAccessIterator first, int nmemb, Compare comp) {
                if (nmemb < 2) return;

                std::vector<T> swap(nmemb / 2);
                tail_merge(first, swap, nmemb, 1, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void twin_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                _twinsort<T>(first, last - first, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void tail_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
                _tailsort<T>(first, last - first, comp);
            }
        }

        namespace reverse_lazy_stable_sorting {
            template<class RandomAccessIterator>
            void rotate_left(RandomAccessIterator start, RandomAccessIterator dest, int size);

            template<class RandomAccessIterator>
            int rotate_right(RandomAccessIterator start, RandomAccessIterator dest, int size);

            template<class RandomAccessIterator>
            void rotate_smart(RandomAccessIterator start, RandomAccessIterator dest, int size);


            template<class RandomAccessIterator>
            void rotate_left(RandomAccessIterator start, RandomAccessIterator dest, int size) {
                int amount = start - dest;
                if (size > 1) {
                    while (amount >= size) {
                        for (auto i = start; i > start - size; i--) std::iter_swap(i - 1, i + size - 1);
                        start -= size;
                        amount -= size;
                    }
                    if (amount > 0) {
                        rotate_smart(start, dest, size / 2);
                        rotate_smart(start + size / 2, dest + size / 2, size - (size / 2));
                    }
                }
                else {
                    auto tmp = *start;
                    for (auto i = start; i > dest; i--) *i = *(i - 1);
                    *dest = tmp;
                }
            }

            template<class RandomAccessIterator>
            int rotate_right(RandomAccessIterator start, RandomAccessIterator dest, int size) {
                int amount = dest - start;
                int moved = 0;
                if (size > 1) {
                    while (amount >= size) {
                        for (auto i = start; i < start + size; i++) std::iter_swap(i, i + size);
                        start += size;
                        amount -= size;
                        moved += size;
                    }
                }
                else {
                    auto tmp = *start;
                    for (auto i = start; i < dest; i++) *i = *(i + 1);
                    *dest = tmp;
                    moved += dest - start;
                }
                return moved;
            }

            template<class RandomAccessIterator>
            void rotate_smart(RandomAccessIterator start, RandomAccessIterator dest, int size) {
                if (size > start - dest) {
                    int start_dest = start - dest;
                    int moved = rotate_right(dest, start + size - start_dest, start_dest);
                    size -= moved;
                    dest = dest + moved;
                    start = dest + start_dest;
                }
                if (size > 0) {
                    rotate_left(start, dest, size);
                }
            }

            template<class RandomAccessIterator>
            void rotate_common(RandomAccessIterator start, RandomAccessIterator dest, int size) {
                if (start > dest) rotate_smart(start, dest, size);
                else rotate_smart(start + size, start, dest - (start + size));
            }

            template<class RandomAccessIterator, class T, class Compare>
            RandomAccessIterator bin_search(RandomAccessIterator start, RandomAccessIterator i, const T& num, Compare comp) {
                auto lo = start, hi = i;

                while (lo < hi) {
                    auto mid = lo + (hi - lo) / 2;

                    if (comp(num, *mid)) hi = mid;
                    else lo = mid + 1;
                }
                return lo;
            }

            template<class RandomAccessIterator, class Compare>
            void lazy_merge(RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end, Compare comp) {
                int threshold = 2;
                for (int len = mid - start; len >= 1; len /= 2) threshold++;

                while (start < mid && mid < end) {
                    if (comp(*start, *mid)) {
                        start++;
                        int i;
                        for (i = 0; i < threshold; i++) {
                            if (comp(*start, *mid)) start++;
                            else break;
                        }
                        if (i == threshold) start = bin_search(start, mid, *mid, comp);
                    }
                    if (start >= mid) break;

                    int size = bin_search(mid, end, *start, comp) - mid;
                    rotate_smart(mid, start, size);

                    start += size + 1;
                    mid += size;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void lazy_stable_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;

                for (auto i = first; i < last - 1; i += 2)
                    if (comp(*(i + 1), *i)) std::iter_swap(i, i + 1);

                int gap;
                for (gap = 4; gap <= length; gap *= 2)
                    for (auto i = first; i + gap <= last; i += gap)
                        lazy_merge(i, i + gap / 2, i + gap, comp);

                if (length - gap / 2 > 0)
                    lazy_merge(first, first + gap / 2, last, comp);
            }
        }

        namespace index_merge_sorting {
            template<class RandomAccessIterator>
            void index_sort(RandomAccessIterator first, std::vector<int>& idx, int a, int b) {
                while (a < b) {
                    if (a != idx[a]) {
                        auto t = *(first + a);
                        int i = a, nxt = idx[a];

                        do {
                            *(first + i) = *(first + nxt);
                            idx[i] = i;

                            i = nxt;
                            nxt = idx[nxt];
                        } while (nxt != a);

                        *(first + i) = t;
                        idx[i] = i;
                    }
                    a++;
                }
            }

            template<class RandomAccessIterator, class Compare>
            void index_merge(RandomAccessIterator first, std::vector<int>& idx, int a, int m, int b, Compare comp) {
                int i = a, j = m, c = a;

                while (i < m && j < b) {
                    if (!comp(*(first + j), *(first + i))) idx[c++] = i++;
                    else idx[c++] = j++;
                }

                while (i < m) idx[c++] = i++;
                while (j < b) idx[c++] = j++;
                index_sort(first, idx, a, b);
            }

            template<class RandomAccessIterator, class Compare>
            void _index_merge_sort(RandomAccessIterator first, std::vector<int>& idx, int a, int b, Compare comp) {
                if (b - a < 2) return;

                int m = (a + b) / 2;
                _index_merge_sort(first, idx, a, m, comp);
                _index_merge_sort(first, idx, m, b, comp);
                index_merge(first, idx, a, m, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void index_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> idx(length);
                std::iota(idx.begin(), idx.end(), 0);

                _index_merge_sort(first, idx, 0, length, comp);
            }
        }

        namespace index_quick_sorting {
            using sortcpp::_internal::index_merge_sorting::index_sort;
            using sortcpp::_internal::bogo_sorting::rd;
            using sortcpp::_internal::bogo_sorting::randint;

            template<class RandomAccessIterator, class Compare>
            int index_partition(RandomAccessIterator first, std::vector<int>& idx, int a, int b, Compare comp) {
                int c0 = a, c1 = c0, i;

                std::mt19937 gen(rd());
                int m = randint(a, b, gen);

                for (i = a; i < m; i++)
                    if (!comp(*(first + m), *(first + i))) c1++;
                
                i++;
                c1++;
                for (; i < b; i++)
                    if (comp(*(first + i), *(first + m))) c1++;
                
                int p = c1 - 1;
                for (i = a; i < m; i++) {
                    if (!comp(*(first + m), *(first + i))) idx[c0++] = i;
                    else idx[c1++] = i;
                }
                idx[p] = i++;

                for (; i < b; i++) {
                    if (comp(*(first + i), *(first + m))) idx[c0++] = i;
                    else idx[c1++] = i;
                }

                index_sort(first, idx, a, b);
                return p;
            }
            
            template<class RandomAccessIterator, class Compare>
            void _index_quick_sort(RandomAccessIterator first, std::vector<int>& idx, int a, int b, Compare comp) {
                if (b - a < 2) return;

                int p = index_partition(first, idx, a, b, comp);
                _index_quick_sort(first, idx, a, p, comp);
                _index_quick_sort(first, idx, p + 1, b, comp);
            }

            template<class RandomAccessIterator, class Compare>
            void index_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
                int length = last - first;
                std::vector<int> idx(length);
                std::iota(idx.begin(), idx.end(), 0);
                _index_quick_sort(first, idx, 0, length, comp);
            }
        }

        
    }
}