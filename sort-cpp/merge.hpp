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
        // Block Swap Merge Sort
        template<class RandomAccessIterator, class Compare>
        void block_swap_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            
            auto multi_swap = [&](RandomAccessIterator a, RandomAccessIterator b, int len) -> void {
                for (int i = 0; i < len; i++) std::iter_swap(a + i, b + i);
            };

            auto binary_search_mid = [&](RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> int {
                int a = 0, b = std::min(mid - start, end - mid), m = a + (b - a) / 2;

                while (b > a) {
                    if (comp(*(mid + m), *(mid - m - 1))) a = m + 1;
                    else b = m;
                    m = a + (b - a) / 2;
                }
                return m;
            };

            auto multi_swap_merge = [&](auto &&self, RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> void {
                int m = binary_search_mid(start, mid, end);
                while (m > 0) {
                    multi_swap(mid - m, mid, m);
                    self(self, mid, mid + m, end);

                    end = mid;
                    mid -= m;

                    m = binary_search_mid(start, mid, end);
                }
            };

            auto multi_swap_merge_sort = [&](auto &&self, RandomAccessIterator a, RandomAccessIterator b) -> void {
                int len = b - a;
                RandomAccessIterator i;

                for (int j = 1; j < len; j *= 2) {
                    for (i = a; i + 2 * j <= b; i += 2 * j) multi_swap_merge(multi_swap_merge, i, i + j, i + 2 * j);
                    if (i + j < b) multi_swap_merge(multi_swap_merge, i, i + j, b);
                }
            };
            multi_swap_merge_sort(multi_swap_merge_sort, first, last);
        }

        template<class RandomAccessIterator>
        void block_swap_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            block_swap_merge_sort(first, last, std::less<void>());
        }


        // Bottom-Up Merge Sort
        template<class RandomAccessIterator, class Compare>
        void bottom_up_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            
            int length = last - first;
            std::vector<T> scratch;
            int copy_length;
            
            auto merge = [&](RandomAccessIterator first, RandomAccessIterator last, RandomAccessIterator index, int merge_size) -> void {
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
                } else copy_length = left - first;
            };

            scratch.resize(length);
            int merge_size = 2;

            while (merge_size <= length) {
                copy_length = length;
                for (auto i = first; i < last; i += merge_size) merge(first, last, i, merge_size);
                for (int i = 0; i < copy_length; i++) *(first + i) = scratch[i];
                merge_size *= 2;
            }

            if ((merge_size / 2) != length) {
                merge(first, last, first, merge_size);
                for (int i = 0; i < length; i++) *(first + i) = scratch[i];
            }
        }

        template<class RandomAccessIterator>
        void bottom_up_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            bottom_up_merge_sort(first, last, std::less<void>());
        }


        // Buffered Stooge Sort
        template<class RandomAccessIterator, class Compare>
        void buffered_stooge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto wrapper = [&](auto&& self, RandomAccessIterator start, RandomAccessIterator stop) -> void {
                if (stop - start > 1) {
                    if (stop - start == 2 && comp(*(stop - 1), *start)) std::iter_swap(start, stop - 1);
                    if (stop - start > 2) {
                        auto third = start + std::ceil((stop - start) / 3.0);
                        auto two_third = start + std::ceil((stop - start) / 3.0 * 2.0);

                        if (two_third - third < third - start) two_third--;
                        if ((stop - start - 2) % 3 == 0) two_third--;

                        self(self, third, two_third);
                        self(self, two_third, stop);

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

                        self(self, two_third, stop);
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
            wrapper(wrapper, first, last);
        }

        template<class RandomAccessIterator>
        void buffered_stooge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            buffered_stooge_sort(first, last, std::less<void>());
        }


        // Improved In-Place Merge Sort
        template<class RandomAccessIterator, class Compare>
        void improved_inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto push = [&](RandomAccessIterator p, RandomAccessIterator a, RandomAccessIterator b) -> void {
                if (a == b) return;

                auto temp = *p;
                *p = *a;

                for (auto i = a + 1; i < b; i++) *(i - 1) = *i;
                *(b - 1) = temp;
            };

            auto merge = [&](RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) -> void {
                auto i = a, j = m;

                while (i < m && j < b) {
                    if (comp(*j, *i)) j++;
                    else push(i++, m, j);
                }

                while (i < m) push(i++, m, b);
            };

            auto merge_sort = [&](auto &&self, RandomAccessIterator a, RandomAccessIterator b) -> void {
                auto m = a + (b - a) / 2;
                
                if (b - a > 2) {
                    if (b - a > 3) self(self, a, m);
                    self(self, m, b);
                }
                merge(a, m, b);
            };
            merge_sort(merge_sort, first, last);
        }

        template<class RandomAccessIterator>
        void improved_inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            improved_inplace_merge_sort(first, last, std::less<void>());
        }


        // In-Place Merge Sort
        template<class RandomAccessIterator, class Compare>
        void inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto push = [&](RandomAccessIterator low, RandomAccessIterator high) -> void {
                for (auto i = low; i < high; i++) 
                    if (comp(*(i + 1), *i)) std::iter_swap(i, i + 1);
            };

            auto merge = [&](RandomAccessIterator min, RandomAccessIterator max, RandomAccessIterator mid) -> void {
                auto i = min;
                while (i <= mid) {
                    if (comp(*(mid + 1), *i)) {
                        std::iter_swap(i, mid + 1);
                        push(mid + 1, max);
                    }
                    i++;
                }
            };

            auto merge_sort = [&](auto&& self, RandomAccessIterator min, RandomAccessIterator max) -> void {
                if (min == max) return;
                if (max - min == 1) {
                    if (comp(*max, *min)) std::iter_swap(min, max);
                    return;
                }

                auto mid = min + (max - min) / 2;
                self(self, min, mid);
                self(self, mid + 1, max);
                merge(min, max, mid);
            };
            merge_sort(merge_sort, first, last - 1);
        }

        template<class RandomAccessIterator>
        void inplace_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            inplace_merge_sort(first, last, std::less<void>());
        }


        // Iterative Top-Down Merge Sort
        template<class RandomAccessIterator, class Compare>
        void iterative_top_down_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::ceil_pow2;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;

            int length = last - first;
            std::vector<T> tmp;

            auto merge = [&](RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> void {
                auto low = start, high = mid;
                int nxt = start - first;

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

                for (auto i = start; i < end; i++) *i = tmp[i - first];
            };

            auto run_sort_small = [&](RandomAccessIterator first, RandomAccessIterator last) -> void {
                int length = last - first;
                for (int subarray_count = ceil_pow2(length); subarray_count > 1; subarray_count >>= 1)
                    for (int i = 0; i < subarray_count; i += 2)
                        merge(first + length * i / subarray_count, first + length * (i + 1) / subarray_count, first + length * (i + 2) / subarray_count);
            };

            auto run_sort_large = [&](RandomAccessIterator first, RandomAccessIterator last) -> void {
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

                        merge(start, mid, whole);
                    }

                    subarray_count >>= 1;
                    whole_i <<= 1;
                    if (frac_i >= subarray_count) {
                        whole_i++;
                        frac_i -= subarray_count;
                    }
                }
            };

            tmp.resize(length);
            
            if (length < 1 << 15) run_sort_small(first, last);
            else run_sort_large(first, last);
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
            using sortcpp::_internal::ceil_pow2;
            using sortcpp::_internal::equals;
            
            auto rotate_equal = [&](RandomAccessIterator a, RandomAccessIterator b, int size) -> void {
                for (int i = 0; i < size; i++) std::iter_swap(a + i, b + i);
            };

            auto rotate = [&](RandomAccessIterator mid, RandomAccessIterator a, RandomAccessIterator b) -> void {
                while (a > first && b > first) {
                    if (a > b) {
                        rotate_equal(mid - (b - first), mid, b - first);
                        mid -= b - first;
                        a -= b - first;
                    }
                    else {
                        rotate_equal(mid - (a - first), mid, a - first);
                        mid += a - first;
                        b -= a - first;
                    }
                }
            };

            auto shuffle_easy = [&](RandomAccessIterator start, int size) -> void {
                for (int i = 1; i < size; i *= 3) {
                    auto val = *(start + i - 1);
                    for (int j = i * 2 % size; j != i; j = j * 2 % size) {
                        auto nval = *(start + j - 1);
                        *(start + j - 1) = val;
                        val = nval;
                    }
                    *(start + i - 1) = val;
                }
            };

            auto shuffle = [&](RandomAccessIterator start, RandomAccessIterator end) -> void {
                while (end - start > 1) {
                    int n = (end - start) / 2;
                    int l = 1;
                    while (l * 3 - 1 <= 2 * n) l *= 3;
                    int m = (l - 1) / 2;

                    rotate(start + n, first + n - m, first + m);
                    shuffle_easy(start, l);
                    start += l - 1;
                }
            };

            auto rotate_shuffled_equal = [&](RandomAccessIterator a, RandomAccessIterator b, int size) -> void {
                for (int i = 0; i < size; i += 2) std::iter_swap(a + i, b + i);
            };

            auto rotate_shuffled = [&](RandomAccessIterator mid, RandomAccessIterator a, RandomAccessIterator b) -> void {
                while (a > first && b > first) {
                    if (a > b) {
                        rotate_shuffled_equal(mid - (b - first), mid, b - first);
                        mid -= b - first;
                        a -= b - first;
                    }
                    else {
                        rotate_shuffled_equal(mid - (a - first), mid, a - first);
                        mid += a - first;
                        b -= a - first;
                    }
                }
            };

            auto rotate_shuffled_outer = [&](RandomAccessIterator mid, RandomAccessIterator a, RandomAccessIterator b) -> void {
                if (a > b) {
                    rotate_shuffled_equal(mid - (b - first), mid + 1, b - first);
                    mid -= b - first;
                    a -= b - first;
                    rotate_shuffled(mid, a, b);
                }
                else {
                    rotate_shuffled_equal(mid - (a - first), mid + 1, a - first);
                    mid += a - first + 1;
                    b -= a - first;
                    rotate_shuffled(mid, a, b);
                }
            };

            auto unshuffle_easy = [&](RandomAccessIterator start, int size) -> void {
                for (int i = 1; i < size; i *= 3) {
                    int prev = i;
                    auto val = *(start + i - 1);
                    for (int j = i * 2 % size; j != i; j = j * 2 % size) {
                        *(start + prev - 1) = *(start + j - 1);
                        prev = j;
                    }
                    *(start + prev - 1) = val;
                }
            };

            auto unshuffle = [&](RandomAccessIterator start, RandomAccessIterator end) -> void {
                while (end - start > 1) {
                    int n = (end - start) / 2;
                    int l = 1;
                    while (l * 3 - 1 <= 2 * n) l *= 3;
                    int m = (l - 1) / 2;
                    
                    rotate_shuffled_outer(start + 2 * m, first + 2 * m, first + 2 * n - 2 * m);
                    unshuffle_easy(start, l);
                    start += l - 1;
                }
            };

            auto merge_up = [&](RandomAccessIterator start, RandomAccessIterator end, bool type) -> void {
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
                        rotate(j, j - i + first, first + 1);
                        break;
                    }
                    else {
                        int r = 0;
                        if (type) while (j + 2 * r < end && !comp(*i, *(j + 2 * r))) r++;
                        else while (j + 2 * r < end && comp(*(j + 2 * r), *i)) r++;
                        j--;
                        unshuffle(j, j + 2 * r);
                        rotate(j, j - i + first, first + r);
                        i += r + 1;
                        j += 2 * r + 1;
                    }
                }
            };

            auto merge = [&](RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> void {
                if (mid - start <= end - mid) {
                    shuffle(start, end);
                    merge_up(start, end, true);
                }
                else {
                    shuffle(start + 1, end);
                    merge_up(start, end, false);
                }
            };

            auto run_sort_small = [&](RandomAccessIterator first, RandomAccessIterator last) -> void {
                int length = last - first;
                for (int subarray_count = ceil_pow2(length); subarray_count > 1; subarray_count >>= 1)
                    for (int i = 0; i < subarray_count; i += 2)
                        merge(first + length * i / subarray_count, first + length * (i + 1) / subarray_count, first + length * (i + 2) / subarray_count);
            };

            auto run_sort_large = [&](RandomAccessIterator first, RandomAccessIterator last) -> void {
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

                        merge(start, mid, whole);
                    }

                    subarray_count >>= 1;
                    whole_i <<= 1;
                    if (frac_i >= subarray_count) {
                        whole_i++;
                        frac_i -= subarray_count;
                    }
                }
            };

            int length = last - first;
            if (length < 1 << 15) run_sort_small(first, last);
            else run_sort_large(first, last);
        }

        template<class RandomAccessIterator>
        void new_shuffle_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            new_shuffle_merge_sort(first, last, std::less<void>());
        }


        // Pattern-Defeating Merge Sort (Simplified Tim Sort)
        template<class RandomAccessIterator, class Compare>
        void pattern_defeating_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::reversal;
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            
            std::vector<T> copied;
            int run_count;

            auto merge_up = [&](RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> void {
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
            };

            auto merge_down = [&](RandomAccessIterator start, RandomAccessIterator mid, RandomAccessIterator end) -> void {
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
            };

            auto merge = [&](RandomAccessIterator left_start, RandomAccessIterator right_start, RandomAccessIterator end) -> void {
                if (end - right_start < right_start - left_start)
                    merge_down(left_start, right_start, end);
                else 
                    merge_up(left_start, right_start, end);
            };

            auto compare = [&](RandomAccessIterator a, RandomAccessIterator b) -> bool {
                return !comp(*b, *a);
            };

            auto identify_run = [&](RandomAccessIterator index, RandomAccessIterator max_index) -> RandomAccessIterator {
                auto start_index = index;
                if (index >= max_index) return max_index + 1;

                bool cmp = compare(index, index + 1);
                index++;

                while (index < max_index) {
                    bool check_cmp = compare(index, index + 1);
                    if (check_cmp != cmp) break;
                    index++;
                }

                if (!cmp) reversal(start_index, index + 1);
                if (index >= max_index) return max_index + 1;
                return index + 1;
            };

            auto find_runs = [&](RandomAccessIterator start, RandomAccessIterator max_index) {
                std::vector<RandomAccessIterator> runs((max_index - start) / 2 + 2);
                run_count = 0;

                auto last_run = start;
                while (last_run != max_index + 1) {
                    runs[run_count++] = last_run;
                    auto new_run = identify_run(last_run, max_index);
                    last_run = new_run;
                }
                return runs;
            };

            int length = last - first;
            auto runs = find_runs(first, last - 1);
            copied.resize(length / 2);

            while (run_count > 1) {
                for (int i = 0; i < run_count - 1; i += 2) {
                    auto end = i + 2 >= run_count ? last : runs[i + 2];
                    merge(runs[i], runs[i + 1], end);
                }
                for (int i = 1, j = 2; i < run_count; i++, j += 2, run_count--) 
                    runs[i] = runs[j];
            }
        }

        template<class RandomAccessIterator>
        void pattern_defeating_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            pattern_defeating_merge_sort(first, last, std::less<void>());
        }

        
        // TODO: Quad Sort


        // Rotate Merge Sort
        template<class RandomAccessIterator, class Compare>
        void rotate_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            
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

            auto binary_search = [&](RandomAccessIterator a, RandomAccessIterator b, const T& value, bool left) {
                while (a < b) {
                    auto m = a + (b - a) / 2;
                    bool cmp = left ? (!comp(*m, value)) : comp(value, *m);
                    if (cmp) b = m;
                    else a = m + 1;
                }
                return a;
            };

            auto rotate_merge = [&](auto &&self, RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) -> void {
                RandomAccessIterator m1, m2, m3;

                if (m - a >= b - m) {
                    m1 = a + (m - a) / 2;
                    m2 = binary_search(m, b, *m1, true);
                    m3 = m1 + (m2 - m);
                }
                else {
                    m2 = m + (b - m) / 2;
                    m1 = binary_search(a, m, *m2, false);
                    m3 = (m2++) - (m - m1);
                }
                rotate(m1, m, m2);

                if (m2 - (m3 + 1) > 0 && b - m2 > 0) self(self, m3 + 1, m2, b);
                if (m1 - a > 0 && m3 - m1 > 0) self(self, a, m1, m3);
            };

            auto rotate_merge_sort = [&](RandomAccessIterator a, RandomAccessIterator b) {
                int len = b - a;
                RandomAccessIterator i;

                for (int j = 1; j < len; j *= 2) {
                    for (i = a; i + 2 * j <= b; i += 2 * j) rotate_merge(rotate_merge, i, i + j, i + 2 * j);
                    if (i + j < b) rotate_merge(rotate_merge, i, i + j, b);
                }
            };

            rotate_merge_sort(first, last);
        }

        template<class RandomAccessIterator>
        void rotate_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            rotate_merge_sort(first, last, std::less<void>());
        }

        
        // Strand Sort
        template<class RandomAccessIterator, class Compare>
        void strand_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            std::vector<T> sublist;

            auto merge_to = [&](RandomAccessIterator a, RandomAccessIterator m, RandomAccessIterator b) -> void {
                int i = 0, s = m - a;

                while (i < s && m < b) {
                    if (comp(sublist[i], *m)) *(a++) = sublist[i++];
                    else *(a++) = *(m++);
                }

                while (i < s) *(a++) = sublist[i++];
            };

            int length = last - first;
            sublist.resize(length);

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

                merge_to(k, j, last);
                j = k;
            }
        }

        template<class RandomAccessIterator>
        void strand_sort(RandomAccessIterator first, RandomAccessIterator last) {
            strand_sort(first, last, std::less<void>());
        }


        // Weaved Merge Sort
        template<class RandomAccessIterator, class Compare>
        void weaved_merge_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using T = typename std::iterator_traits<RandomAccessIterator>::value_type;
            using sortcpp::_internal::equals;

            std::vector<T> tmp;

            auto _merge = [&](auto&& self, RandomAccessIterator residue, RandomAccessIterator end, int modulus) -> void {
                if (residue + modulus >= end) return;

                auto low = residue, high = residue + modulus;
                int dmodulus = modulus << 1;

                self(self, low, end, dmodulus);
                self(self, high, end, dmodulus);

                int nxt = residue - first;
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

                for (auto i = residue; i < end; i += modulus) *i = tmp[i - first];
            };

            int length = last - first;
            tmp.resize(length);
            _merge(_merge, first, last, 1);
        }

        template<class RandomAccessIterator>
        void weaved_merge_sort(RandomAccessIterator first, RandomAccessIterator last) {
            weaved_merge_sort(first, last, std::less<void>());
        }
    }
}