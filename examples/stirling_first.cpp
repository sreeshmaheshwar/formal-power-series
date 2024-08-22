/**
 * Solution to:
 * https://judge.yosupo.jp/problem/stirling_number_of_the_first_kind
 *
 * Computing S(n, k) for 0 <= k <= n, where S(n, k) are the signed Stirling
 * numbers of the first kind. S(n, k) is (or rather, can be) defined as the
 * coefficient of x^k in the polynomial (x)_n, where (x)_n is the falling
 * factorial, (x)_n = x(x - 1)(x - 2)...(x - n + 1).
 *
 * There is an O(N log N) solution to this task. However, we can in fact
 * compute the above falling factorial product of linear polynomials directly
 * through a simple yet elegant divide and conquer (D&C) approach in O(N log^2
 * N): recurse, then multiply the left half with the right half.
 *
 * Concretely, split the product expansion into two equal halves, compute the
 * product of each and then convolve them together (using NTT). Since the
 * product of the first half has degree equal to the size of the first half, and
 * similarly for the second half, the convolution operates on polynomials of
 * degree N/2 in O(N log N) time (where N is the size of the current D&C range).
 * Thus, the recurrence is T(N) = 2T(N/2) + O(N log N), which solves to O(N
 * log^2 N).
 *
 * In fact, the above D&C idea can be used to compute the product of several
 * large integers much more efficiently than naive left to right multiplication.
 */

#include "../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

PowerSeries falling_factorial_range(int l, int r) {
  if (l > r) {
    return {1};
  }
  if (l == r) {
    return {mint(-l), 1};
  }
  int m = (l + r) / 2;
  return falling_factorial_range(l, m) * falling_factorial_range(m + 1, r);
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  int n;
  std::cin >> n;
  for (const auto &x : falling_factorial_range(0, n - 1)) {
    std::cout << x.val() << ' ';
  }
}
