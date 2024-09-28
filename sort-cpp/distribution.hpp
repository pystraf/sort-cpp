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
            auto max_number_of_digits = [&](RandomAccessIterator first, RandomAccessIterator last) -> int {
                int max = 0;
                for (auto i = first; i < last; i++) {
                    auto val = key(*i);
                    int tmp = (int)(std::log(val) / std::log(buckets)) + 1;
                    if (tmp > max) max = tmp;
                }
                return max;
            };

            auto get_digit = [&](int integer, int divisor) -> int {
                return (integer / divisor) % buckets;
            };

            auto _sort = [&](auto &&self, RandomAccessIterator start, RandomAccessIterator end, int divisor) -> void {
                std::vector<int> count(buckets, 0);
                std::vector<RandomAccessIterator> offset(buckets);
                int digit = 0;

                for (auto i = start; i < end; i++) {
                    auto d = key(*i);
                    digit = get_digit(d, divisor);
                    count[digit]++;
                }

                offset[0] = start;
                for (int i = 1; i < buckets; i++) offset[i] = offset[i - 1] + count[i - 1];

                for (int b = 0; b < buckets; b++) {
                    while (count[b] > 0) {
                        auto origin = offset[b];
                        auto from = origin;
                        auto num = *from;

                        do {
                            digit = get_digit(key(num), divisor);
                            auto to = offset[digit];

                            offset[digit]++;
                            count[digit]--;

                            auto temp = *to;
                            *to = num;

                            num = temp;
                            from = to;
                        } while (from != origin);
                    }
                }

                if (divisor > 1) {
                    for (int i = 0; i < buckets; i++) {
                        auto bg = (i > 0) ? offset[i - 1] : start;
                        auto ed = offset[i];
                        if (ed - bg > 1) self(self, bg, ed, divisor / buckets);
                    }
                }
            };

            int number_of_digits = max_number_of_digits(first, last);
            int max = 1;
            for (int i = 0; i < number_of_digits - 1; i++) max *= buckets;
            _sort(_sort, first, last, max);
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
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            
            auto max = key(*first);
            for (auto i = first + 1; i < last; i++) {
                if (key(*i) > max) max = key(*i);
            }
            
            int length = last - first;
            std::vector<int> counts(max + 1, 0);
            std::vector<T> output(length);
            
            for (auto i = first; i < last; i++) counts[key(*i)]++;
            for (int i = 1; i < counts.size(); i++) counts[i] += counts[i - 1];

            for (auto i = last - 1; i >= first; i--) {
                output[counts[key(*i)] - 1] = *i;
                counts[key(*i)]--;
            }

            for (int i = 0; i < length; i++) *(first + i) = output[i];
        }


        // Flash Sort
        template<class RandomAccessIterator, class Key>
        void flash_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::insertion_sorting::insertion_sort;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            
            if (last - first <= 1) return;

            int length = last - first;
            int m = (int)((0.2 * length) + 2);

            auto min = *first, max = *first;
            auto max_index = first;

            for (auto i = first + 1; i < last - 1; i += 2) {
                T small, big;
                RandomAccessIterator big_index;

                if (key(*i) < key(*(i + 1))) {
                    small = *i;
                    big = *(i + 1);
                    big_index = i + 1;
                }
                else {
                    big = *i;
                    big_index = i;
                    small = *(i + 1);
                }

                if (key(big) > key(max)) {
                    max = big;
                    max_index = big_index;
                }
                if (key(small) < key(min)) min = small;
            }
            
            if (key(*(last - 1)) < key(min)) min = *(last - 1);
            if (key(*(last - 1)) > key(max)) {
                max = *(last - 1);
                max_index = last - 1;
            }
            if (key(min) == key(max)) return;

            std::vector<int> L(m + 1, 0);
            double c = (m - 1.0) / (key(max) - key(min));
            int K;
            for (auto h = first; h < last; h++) {
                K = (int)((key(*h) - key(min)) * c) + 1;
                L[K]++;
            }
            for (K = 2; K <= m; K++) L[K] += L[K - 1];

            std::iter_swap(first, max_index);

            auto j = first;
            K = m;
            int num_moves = 0;

            while (num_moves < length) {
                while (j - first >= L[K]) {
                    j++;
                    K = (int)((key(*j) - key(min)) * c) + 1;
                }

                auto evicted = *j;

                while (j - first < L[K]) {
                    K = (int)((key(evicted) - key(min)) * c) + 1;

                    auto loc = first + (L[K] - 1);
                    auto temp = *loc;
                    *loc = evicted;
                    evicted = temp;

                    L[K]--;
                    num_moves++;
                }
            }

            int threshold = (int)(1.25 * ((length / m) + 1));
            int min_elements = 30;

            for (K = m - 1; K >= 1; K--) {
                int class_size = L[K + 1] - L[K];
                if (class_size > threshold && class_size > min_elements)
                    flash_sort(first + L[K], first + L[K + 1], key);
            }
            insertion_sort(first, last, [&](auto a, auto b) { return key(a) < key(b); });
        }


        // In-Place LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void inplace_lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::get_digit;
            
            auto pos = first;
            int length = last - first;

            std::vector<int> vregs;

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++) 
                if (key(*i) > max) max = key(*i);
            int maxpower = (int)(std::log(max) / std::log(base));

            for (int p = 0; p <= maxpower; p++) {
                vregs.assign(base - 1, length - 1);
                pos = first;

                for (int i = 0; i < length; i++) {
                    auto digit = get_digit(key(*pos), p, base);
                    if (digit == 0) pos++;
                    else {
                        auto end = first + vregs[digit - 1];
                        for (auto i = pos; i < end; i++) std::iter_swap(i, i + 1);
                        for (int j = digit - 1; j > 0; j--) vregs[j - 1]--;
                    }
                }
            }
        }


        // LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::get_digit;
            using sortcpp::_internal::fancy_transcribe;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++)
                if (key(*i) > max) max = key(*i);
            int maxpower = (int)(std::log(max) / std::log(base));

            std::vector<std::vector<T>> registers(base);

            for (int p = 0; p <= maxpower; p++) {
                for (auto i = first; i < last; i++) {
                    auto digit = get_digit(key(*i), p, base);
                    registers[digit].push_back(*i);
                }

                fancy_transcribe<T>(first, last, registers);
            }
        }

        // MSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void msd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::get_digit;
            using sortcpp::_internal::transcribe_msd;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++)
                if (key(*i) > max) max = key(*i);
            int maxpower = (int)(std::log(max) / std::log(base));

            auto radix_msd = [&](auto&& self, RandomAccessIterator min, RandomAccessIterator max, int radix, int pow) -> void {
                if (min >= max || pow < 0) return;

                std::vector<std::vector<T>> registers(radix);

                for (auto i = min; i < max; i++) {
                    auto digit = get_digit(key(*i), pow, radix);
                    registers[digit].push_back(*i);
                }

                transcribe_msd(min, registers);

                int sum = 0;
                for (int i = 0; i < registers.size(); i++) {
                    self(self, min + sum, min + sum + registers[i].size(), radix, pow - 1);
                    sum += registers[i].size();
                }
            };
            radix_msd(radix_msd, first, last, base, maxpower);
        }


        // Rotate LSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void rotate_lsd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::get_digit;

            auto multi_swap = [&](RandomAccessIterator a, RandomAccessIterator b, int len) {
                for (int i = 0; i < len; i++) std::iter_swap(a + i, b + i);
            };

            auto rotate = [&](RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) {
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
            };

            auto bin_search = [&](RandomAccessIterator a, RandomAccessIterator b, const auto& d, int p) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (get_digit(key(*m), p, base) < d) a = m + 1;
                    else b = m;
                }
                return a;
            };

            auto merge = [&](auto&& self, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, int da, int db, int p) -> void {
                if (b - a < 2 || db - da < 2) return;

                int dm = da + (db - da) / 2;
                auto m1 = bin_search(a, m, dm, p);
                auto m2 = bin_search(m, b, dm, p);

                rotate(m1, m, m2);
                m = m1 + (m2 - m);

                self(self, m, m2, b, dm, db, p);
                self(self, a, m1, m, da, dm, p);
            };

            auto merge_sort = [&](auto &&self, RandomAccessIterator a, RandomAccessIterator b, int p) -> void {
                if (b - a < 2) return;

                auto m = a + (b - a) / 2;

                self(self, a, m, p);
                self(self, m, b, p);
                merge(merge, a, m, b, 0, base, p);
            };

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++)
                if (key(*i) > max) max = key(*i);
            int maxpower = (int)(std::log(max) / std::log(base));

            for (int p = 0; p <= maxpower; p++) merge_sort(merge_sort, first, last, p);
        }


        // Rotate MSD Radix Sort
        template<class RandomAccessIterator, class Key>
        void rotate_msd_radix_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int base = 10) {
            using sortcpp::_internal::get_digit;

            auto shift = [&](auto n, int q) {
                while (q > 0) {
                    n /= base;
                    q--;
                }
                return n;
            };

            auto multi_swap = [&](RandomAccessIterator a, RandomAccessIterator b, int len) {
                for (int i = 0; i < len; i++) std::iter_swap(a + i, b + i);
            };

            auto rotate = [&](RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) {
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
            };

            auto bin_search = [&](RandomAccessIterator a, RandomAccessIterator b, const auto& d, int p) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    if (get_digit(key(*m), p, base) < d) a = m + 1;
                    else b = m;
                }
                return a;
            };

            auto merge = [&](auto&& self, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b, int da, int db, int p) -> void {
                if (b - a < 2 || db - da < 2) return;

                int dm = da + (db - da) / 2;
                auto m1 = bin_search(a, m, dm, p);
                auto m2 = bin_search(m, b, dm, p);

                rotate(m1, m, m2);
                m = m1 + (m2 - m);

                self(self, m, m2, b, dm, db, p);
                self(self, a, m1, m, da, dm, p);
            };

            auto merge_sort = [&](auto&& self, RandomAccessIterator a, RandomAccessIterator b, int p) -> void {
                if (b - a < 2) return;

                auto m = a + (b - a) / 2;

                self(self, a, m, p);
                self(self, m, b, p);
                merge(merge, a, m, b, 0, base, p);
            };

            auto dist = [&](RandomAccessIterator a, RandomAccessIterator b, int p) {
                merge_sort(merge_sort, a, b, p);
                return bin_search(a, b, 1, p);
            };

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++)
                if (key(*i) > max) max = key(*i);
            int q = (int)(std::log(max) / std::log(base));

            auto i = first, b = last;
            int m = 0;

            while (i < last) {
                auto p = (b - i < 1) ? i : dist(i, b, q);

                if (q == 0) {
                    m += base;
                    int t = m / base;
                    while (t % base == 0) {
                        t /= base;
                        q++;
                    }

                    i = b;
                    while (b < last && shift(key(*b), q + 1) == shift(m, q + 1)) b++;
                }
                else {
                    b = p;
                    q--;
                }
            }
        }

        
        // Stackless American Flag Sort
        template<class RandomAccessIterator, class Key>
        void stackless_american_flag_sort(RandomAccessIterator first, RandomAccessIterator last, Key key, int buckets = 128) {
            using sortcpp::_internal::get_digit;
            
            auto shift = [&](int n, int q, int r) {
                while (q > 0) {
                    n /= r;
                    q--;
                }
                return n;
            };

            std::vector<int> cnts, offs;

            auto dist = [&](RandomAccessIterator a, RandomAccessIterator b, int q, int r) {
                for (int i = 1; i < r; i++) {
                    cnts[i] += cnts[i - 1];
                    offs[i] = cnts[i - 1];
                }

                for (int i = 0; i < r - 1; i++) {
                    auto pos = a + offs[i];

                    if (cnts[i] > offs[i]) {
                        auto t = *pos;

                        do {
                            int digit = get_digit(key(t), q, r);
                            cnts[digit]--;

                            auto t1 = *(a + cnts[digit]);
                            *(a + cnts[digit]) = t;
                            t = t1;
                        } while (cnts[i] > offs[i]);
                    }
                }

                auto p = a + offs[1];
                for (int i = 0; i < r; i++) offs[i] = cnts[i] = 0;
                return p;
            };

            int r = buckets;
            auto i = first, b = last;

            auto max = key(*first);
            for (auto i = first + 1; i < last; i++)
                if (key(*i) > max) max = key(*i);
            int q = (int)(std::log(max) / std::log(r)), m = 0;

            cnts.resize(r);
            offs.resize(r);

            for (auto j = i; j < b; j++) {
                int digit = get_digit(key(*j), q, r);
                cnts[digit]++;
            }

            while (i < last) {
                auto p = (b - i < 1) ? i : dist(i, b, q, r);

                if (q == 0) {
                    m += r;
                    int t = m / r;

                    while (t % r == 0) {
                        t /= r;
                        q++;
                    }

                    i = b;
                    while (b < last && shift(key(*b), q + 1, r) == shift(m, q + 1, r)) {
                        int digit = get_digit(key(*b), q, r);
                        cnts[digit]++;
                        b++;
                    }
                }
                else {
                    b = p;
                    q--;

                    for (auto j = i; j < b; j++) {
                        int digit = get_digit(key(*j), q, r);
                        cnts[digit]++;
                    }
                }
            }
        }


        // Stackless Binary Quick Sort
        template<class RandomAccessIterator, class Key>
        void stackless_binary_quick_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::get_bit;
            using sortcpp::_internal::binary_quick_sorting::get_max_bit;
            
            auto partition = [&](RandomAccessIterator a, RandomAccessIterator b, int bit) {
                auto i = a - 1, j = b;

                while (true) {
                    do i++;
                    while (i < j && !get_bit(key(*i), bit));

                    do j--;
                    while (j > i && get_bit(key(*j), bit));

                    if (i < j) std::iter_swap(i, j);
                    else return i;
                }
            };

            auto i = first, b = last;
            int q = get_max_bit(first, last, key) - 1, m = 0;

            while (i < last) {
                auto p = (b - i < 1) ? i : partition(i, b, q);

                if (q == 0) {
                    m += 2;
                    while (!get_bit(m, q + 1)) q++;

                    i = b;
                    while (b < last && (key(*b) >> (q + 1)) == (m >> (q + 1))) b++;
                }
                else {
                    b = p;
                    q--;
                }
            }
        }


        // Static Sort
        template<class RandomAccessIterator, class Key>
        void static_sort(RandomAccessIterator first, RandomAccessIterator last, Key key) {
            using sortcpp::_internal::insertion_sorting::insertion_sort;
            using sortcpp::_internal::heap_sorting::heap_sort;

            auto key_comp = [&](const auto& a, const auto& b) {
                return key(a) < key(b);
            };

            auto _static_sort = [&](RandomAccessIterator a, RandomAccessIterator b) {
                auto min = *a, max = *b;
                for (auto i = a + 1; i < b; i++) {
                    if (key(*i) < min) min = *i;
                    if (key(*i) > max) max = *i;
                }

                int aux_len = b - a;
                std::vector<int> count(aux_len + 1);
                std::vector<RandomAccessIterator> offset(aux_len + 1);

                double CONST = (double)aux_len / (key(max) - key(min) + 1);

                int idx;
                for (auto i = a; i < b; i++) {
                    idx = (int)((key(*i) - key(min)) * CONST);
                    count[idx]++;
                }

                offset[0] = a;
                for (int i = 1; i < aux_len; i++) offset[i] = offset[i - 1] + count[i - 1];

                for (int v = 0; v < aux_len; v++) {
                    while (count[v] > 0) {
                        auto origin = offset[v];
                        auto from = origin;
                        auto num = *from;

                        do {
                            idx = (int)((key(num) - key(min)) * CONST);
                            auto to = offset[idx];

                            offset[idx]++;
                            count[idx]--;

                            auto temp = *to;
                            *to = num;

                            num = temp;
                            from = to;
                        } while (from != origin);
                    }
                }

                for (auto i = 0; i < aux_len; i++) {
                    auto s = (i > 1)? offset[i - 1] : a;
                    auto e = offset[i];

                    if (e - s <= 1) continue;
                    if (e - s > 16) heap_sort(s, e, true, key_comp);
                    else insertion_sort(s, e, key_comp);
                }
            };
            _static_sort(first, last);
        }
    }
}