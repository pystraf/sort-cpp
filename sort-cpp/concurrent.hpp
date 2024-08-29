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
    namespace concurrent {
        // Bitonic Sort (Iterative)
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

        template<class RandomAccessIterator>
        void bitonic_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            bitonic_sort_iterative(first, last, std::less<void>());
        }


        // Bitonic Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void bitonic_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto hyperfloor = [&](int n) {
                int k = 1;
                while (k < n) k <<= 1;
                return k >> 1;
            };

            auto compare = [&](RandomAccessIterator i, RandomAccessIterator j, bool dir) {
                if (comp(*j, *i) == dir) std::iter_swap(i, j);
            };

            auto bitonic_merge = [&](auto&& self, RandomAccessIterator lo, int n, bool dir) -> void {
                if (n > 1) {
                    int m = hyperfloor(n);

                    for (auto i = lo; i < lo + n - m; i++) compare(i, i + m, dir);

                    self(self, lo, m, dir);
                    self(self, lo + m, n - m, dir);
                }
            };

            auto bitonic_sort = [&](auto&& self, RandomAccessIterator lo, int n, bool dir) -> void {
                if (n > 1) {
                    int m = n / 2;
                    self(self, lo, m, !dir);
                    self(self, lo + m, n - m, dir);
                    bitonic_merge(bitonic_merge, lo, n, dir);
                }
            };

            bitonic_sort(bitonic_sort, first, last - first, true);
        }

        template<class RandomAccessIterator>
        void bitonic_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            bitonic_sort_recursive(first, last, std::less<void>());
        }


        // Bose Nelson Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void bose_nelson_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto range_comp = [&](RandomAccessIterator a, RandomAccessIterator b, int offset) {
                int half = (b - a) / 2;
                auto m = a + half;
                a += offset;

                for (int i = 0; i < half - offset; i++)
                    if ((i & ~offset) == i) comp_swap(a + i, m + i);
            };

            int size = 1 << (int)(std::ceil(std::log2(last - first)));
            for (int k = 2; k <= size; k *= 2)
                for (int j = 0; j < k / 2; j++)
                    for (auto i = first; i + j < last; i += k)
                        range_comp(i, i + k, j);
        }

        template<class RandomAccessIterator>
        void bose_nelson_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            bose_nelson_sort_iterative(first, last, std::less<void>());
        }


        // Bose Nelson Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void bose_nelson_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto compare_swap = [&](RandomAccessIterator start, RandomAccessIterator end) {
                if (comp(*end, *start)) std::iter_swap(start, end);
            };

            auto bose_nelson_merge = [&](auto&& self, RandomAccessIterator start1, int len1, RandomAccessIterator start2, int len2) -> void {
                if (len1 == 1 && len2 == 1) compare_swap(start1, start2);
                else if (len1 == 1 && len2 == 2) {
                    compare_swap(start1, start2 + 1);
                    compare_swap(start1, start2);
                }
                else if (len1 == 2 && len2 == 1) {
                    compare_swap(start1, start2);
                    compare_swap(start1 + 1, start2);
                }
                else {
                    int mid1 = len1 / 2;
                    int mid2 = (len1 & 1) ? (len2 / 2) : ((len2 + 1) / 2);
                    self(self, start1, mid1, start2, mid2);
                    self(self, start1 + mid1, len1 - mid1, start2 + mid2, len2 - mid2);
                    self(self, start1 + mid1, len1 - mid1, start2, mid2);
                }
            };

            auto bose_nelson = [&](auto&& self, RandomAccessIterator start, int length) -> void {
                if (length > 1) {
                    int mid = length / 2;
                    self(self, start, mid);
                    self(self, start + mid, length - mid);
                    bose_nelson_merge(bose_nelson_merge, start, mid, start + mid, length - mid);
                }
            };
            bose_nelson(bose_nelson, first, last - first);
        }

        template<class RandomAccessIterator>
        void bose_nelson_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            bose_nelson_sort_recursive(first, last, std::less<void>());
        }


        // Crease Sort
        template<class RandomAccessIterator, class Compare>
        void crease_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            int length = last - first;
            int max = 1;
            for (; max * 2 < length; max *= 2);

            int next = max;
            while (next > 0) {
                for (auto i = first; i + 1 < last; i++) comp_swap(i, i + 1);
                for (int j = max; j >= next && j > 1; j /= 2)
                    for (auto i = first + 1; i + j - 1 < last; i += 2)
                        comp_swap(i, i + j - 1);
                next /= 2;
            }
        }

        template<class RandomAccessIterator>
        void crease_sort(RandomAccessIterator first, RandomAccessIterator last) {
            crease_sort(first, last, std::less<void>());
        }


        // Diamond Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void diamond_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

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
                                comp_swap(i, i + 1);
                }
            }

            m /= 2;
            for (int k = 0; k <= m / 2; k++)
                for (auto i = first + k; i + 1 < std::min(last, first + m - k); i += 2)
                    comp_swap(i, i + 1);
        }

        template<class RandomAccessIterator>
        void diamond_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            diamond_sort_iterative(first, last, std::less<void>());
        }


        // Diamond Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void diamond_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto _sort = [&](auto&& self, RandomAccessIterator start, RandomAccessIterator stop, bool merge) -> void {
                if (stop - start == 2) {
                    if (comp(*(stop - 1), *start)) std::iter_swap(start, stop - 1);
                    return;
                }
                if (stop - start >= 3) {
                    double div = (stop - start) / 4.0;
                    auto mid = (stop - start) / 2 + start;
                    if (merge) {
                        self(self, start, mid, true);
                        self(self, mid, stop, true);
                    }
                    self(self, (int)div + start, (int)(div * 3) + start, false);
                    self(self, start, mid, false);
                    self(self, mid, stop, false);
                    self(self, (int)div + start, (int)(div * 3) + start, false);
                }
            };
            _sort(_sort, first, last, true);
        }

        template<class RandomAccessIterator>
        void diamond_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            diamond_sort_recursive(first, last, std::less<void>());
        }


        // Fold Sort
        template<class RandomAccessIterator, class Compare>
        void fold_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto halver = [&](RandomAccessIterator low, RandomAccessIterator high) {
                while (low < high) {
                    comp_swap(low, high);
                    low++;
                    high--;
                }
            };

            int length = last - first;
            int ceil_log = 1;
            for (; (1 << ceil_log) < length; ceil_log++);

            int size = 1 << ceil_log;
            for (int k = size >> 1; k > 0; k >>= 1)
                for (int i = size; i >= k; i >>= 1)
                    for (auto j = first; j < last; j += i)
                        halver(j, j + i - 1);
        }

        template<class RandomAccessIterator>
        void fold_sort(RandomAccessIterator first, RandomAccessIterator last) {
            fold_sort(first, last, std::less<void>());
        }


        // Matrix Sort
        template<class RandomAccessIterator, class Compare>
        void matrix_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            using sortcpp::_internal::matrix_shape;

            auto gap_reverse = [&](RandomAccessIterator start, RandomAccessIterator end, int gap) {
                for (auto i = start, j = end; i < j; i += gap, j -= gap)
                    std::iter_swap(i, j - gap);
            };

            auto compare = [&](const auto& a, const auto& b, bool dir) {
                return dir ? comp(a, b) : comp(b, a);
            };

            auto insert_last = [&](RandomAccessIterator a, RandomAccessIterator b, int gap, bool dir) {
                bool did = false;
                auto key = *b;
                auto j = b - gap;

                while (j >= a && compare(key, *j, dir)) {
                    *(j + gap) = *j;
                    did = true;
                    j -= gap;
                }
                *(j + gap) = key;
                return did;
            };

            auto get_matrix_dims = [&](int len) {
                int dim = (int)std::sqrt(len);
                bool insert_last = false;
                if (dim * dim == len - 1) insert_last = true;
                for (; len % dim != 0; dim--);
                return matrix_shape(dim, len / dim, insert_last);
            };

            auto _matrix_sort = [&](auto&& self, RandomAccessIterator start, RandomAccessIterator end, int gap, bool dir) -> bool {
                bool did = false;
                int length = (end - start) / gap;
                if (length < 2) return false;
                else if (length <= 16) {
                    did = false;
                    for (auto i = start; i < end; i += gap)
                        did = insert_last(start, i, gap, dir) || did;
                }
                else {
                    bool newdid;
                    auto shape = get_matrix_dims(length);

                    if (shape.insert_last) {
                        bool did1 = self(self, start, end - gap, gap, dir);
                        bool did2 = insert_last(start, end - gap, gap, dir);
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
                            newdid = self(self, i, i + width * gap, gap, curdir) || newdid;
                            did = did || newdid;
                            curdir = !curdir;
                        }

                        newdid = false;
                        for (int i = 0; i < width; i++) {
                            newdid = self(self, start + i * gap, end + i * gap, gap * width, dir) || newdid;
                            did = did || newdid;
                        }

                    } while (newdid);

                    for (auto i = start + width * gap; i < end; i += 2 * width * gap)
                        gap_reverse(i, i + width * gap, gap);
                }

                return did;
            };
            _matrix_sort(_matrix_sort, first, last, 1, true);
        }

        template<class RandomAccessIterator>
        void matrix_sort(RandomAccessIterator first, RandomAccessIterator last) {
            matrix_sort(first, last, std::less<void>());
        }


        // Odd-Even Merge Sort (Iterative)
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

        template<class RandomAccessIterator>
        void odd_even_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_merge_sort_iterative(first, last, std::less<void>());
        }


        // Odd-Even Merge Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void odd_even_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto compare = [&](RandomAccessIterator i, RandomAccessIterator j) {
                if (comp(*j, *i)) std::iter_swap(i, j);
            };

            auto oddeven_merge = [&](auto&& self, RandomAccessIterator lo, int m2, int n, int r) -> void {
                int m = r * 2;
                if (m < n) {
                    if ((n / r) % 2 != 0) {
                        self(self, lo, (m2 + 1) / 2, n + r, m);
                        self(self, lo + r, m2 / 2, n - r, m);
                    }
                    else {
                        self(self, lo, (m2 + 1) / 2, n, m);
                        self(self, lo + r, m2 / 2, n, m);
                    }

                    if (m2 % 2 != 0) {
                        for (auto i = lo; i + r < lo + n; i += m) compare(i, i + r);
                    }
                    else {
                        for (auto i = lo + r; i + r < lo + n; i += m) compare(i, i + r);
                    }
                }
                else if (n > r) compare(lo, lo + r);
            };

            auto oddeven_merge_sort = [&](auto&& self, RandomAccessIterator lo, int n) -> void {
                if (n > 1) {
                    int m = n / 2;
                    self(self, lo, m);
                    self(self, lo + m, n - m);
                    oddeven_merge(oddeven_merge, lo, m, n, 1);
                }
            };
            oddeven_merge_sort(oddeven_merge_sort, first, last - first);
        }

        template<class RandomAccessIterator>
        void odd_even_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            odd_even_merge_sort_recursive(first, last, std::less<void>());
        }


        // Pairwise Merge Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void pairwise_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };
            
            int length = last - first;
            int n = 1;
            for (; n < length; n <<= 1);

            for (int k = n >> 1; k > 0; k >>= 1)
                for (auto j = first; j < last; j += k << 1)
                    for (int i = 0; i < k; i++)
                        comp_swap(j + i, j + i + k);

            for (int k = 2; k < n; k <<= 1)
                for (int m = k >> 1; m > 0; m >>= 1)
                    for (auto j = first; j < last; j += k << 1)
                        for (int p = m; p < ((k - m) << 1); p += m << 1)
                            for (int i = 0; i < m; i++)
                                comp_swap(j + p + i, j + p + m + i);
        }

        template<class RandomAccessIterator>
        void pairwise_merge_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_merge_sort_iterative(first, last, std::less<void>());
        }


        // Pairwise Merge Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void pairwise_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto pairwise_merge = [&](auto&& self, RandomAccessIterator a, RandomAccessIterator b) -> void {
                auto m = a + (b - a) / 2;
                auto m1 = a + (m - a) / 2;
                int g = m - m1;

                for (int i = 0; m1 + i < m; i++)
                    for (int j = m1 - first, k = g; k > 0; k >>= 1, j -= k - (i & k))
                        comp_swap(first + j + i, first + j + i + k);
                
                if (b - a > 4) self(self, m, b);
            };

            auto pairwise_merge_sort = [&](auto&& self, RandomAccessIterator a, RandomAccessIterator b) -> void {
                auto m = a + (b - a) / 2;
                for (auto i = a, j = m; i < m; i++, j++) comp_swap(i, j);
                
                if (b - a > 2) {
                    self(self, a, m);
                    self(self, m, b);
                    pairwise_merge(pairwise_merge, a, b);
                }
            };

            int length = last - first;
            int n = 1;
            for (; n < length; n <<= 1);
            pairwise_merge_sort(pairwise_merge_sort, first, first + n);
        }

        template<class RandomAccessIterator>
        void pairwise_merge_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_merge_sort_recursive(first, last, std::less<void>());
        }


        // Pairwise Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void pairwise_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto iterative_pairwise = [&](RandomAccessIterator start, int length) {
                int a = 1, b = 0, c = 0, d = 0, e = 0;
                while (a < length) {
                    b = a;
                    c = 0;
                    while (b < length) {
                        comp_swap(start + b - a, start + b);
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
                            comp_swap(start + b - d * a, start + b);
                            c = (c + 1) % a;
                            b++;
                            if (c == 0) b += a;
                        }
                        d /= 2;
                    }
                    a /= 2;
                    e = e * 2 + 1;
                }
            };

            iterative_pairwise(first, last - first);
        }

        template<class RandomAccessIterator>
        void pairwise_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_sort_iterative(first, last, std::less<void>());
        }


        // Pairwise Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void pairwise_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto pairwise_rec = [&](auto&& self, RandomAccessIterator start, RandomAccessIterator end, int gap) -> void {
                if (start == end - gap) return;

                auto b = start + gap;
                while (b < end) {
                    comp_swap(b - gap, b);
                    b += 2 * gap;
                }

                if (((end - start) / gap) % 2 == 0) {
                    self(self, start, end, gap * 2);
                    self(self, start + gap, end + gap, gap * 2);
                }
                else {
                    self(self, start, end + gap, gap * 2);
                    self(self, start + gap, end, gap * 2);
                }

                int a = 1;
                while (a < ((end - start) / gap)) a = a * 2 + 1;

                b = start + gap;
                while (b + gap < end) {
                    int c = a;
                    while (c > 1) {
                        c /= 2;
                        comp_swap(b, b + c * gap);
                    }
                    b += 2 * gap;
                }
            };

            pairwise_rec(pairwise_rec, first, last, 1);
        }

        template<class RandomAccessIterator>
        void pairwise_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            pairwise_sort_recursive(first, last, std::less<void>());
        }

        
        // Weave Sort (Iterative)
        template<class RandomAccessIterator, class Compare>
        void weave_sort_iterative(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };
            
            int length = last - first;
            int n = 1;
            for (; n < length; n *= 2);

            for (int i = 1; i < n; i *= 2)
                for (int j = 1; j <= i; j *= 2)
                    for (int k = 0; k < n; k += n / j)
                        for (int d = n / i / 2, m = 0, l = n / j - d; l >= n / j / 2; l -= d)
                            for (int p = 0; p < d; p++, m++)
                                comp_swap(first + k + m, first + k + l + p);

        }

        template<class RandomAccessIterator>
        void weave_sort_iterative(RandomAccessIterator first, RandomAccessIterator last) {
            weave_sort_iterative(first, last, std::less<void>());
        }


        // Weave Sort (Recursive)
        template<class RandomAccessIterator, class Compare>
        void weave_sort_recursive(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
            auto comp_swap = [&](RandomAccessIterator a, RandomAccessIterator b) {
                if (b >= last) return;
                if (comp(*b, *a)) std::iter_swap(a, b);
            };

            auto circle = [&](auto&& self, RandomAccessIterator pos, int len, int gap) -> void {
                if (len < 2) return;

                for (int i = 0; 2 * i < (len - 1) * gap; i += gap)
                    comp_swap(pos + i, pos + (len - 1) * gap - i);

                self(self, pos, len / 2, gap);
                if (pos + len * gap / 2 < last)
                    self(self, pos + len * gap / 2, len / 2, gap);
            };

            auto weave_circle = [&](auto&& self, RandomAccessIterator pos, int len, int gap) -> void {
                if (len < 2) return;

                self(self, pos, len / 2, 2 * gap);
                self(self, pos + gap, len / 2, 2 * gap);
                circle(circle, pos, len, gap);
            };

            int length = last - first;
            int n = 1;
            for (; n < length; n *= 2);
            weave_circle(weave_circle, first, n, 1);
        }

        template<class RandomAccessIterator>
        void weave_sort_recursive(RandomAccessIterator first, RandomAccessIterator last) {
            weave_sort_recursive(first, last, std::less<void>());
        }
    }
}