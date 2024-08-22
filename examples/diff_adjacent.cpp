/**
 * Solution to:
 * https://atcoder.jp/contests/abc297/tasks/abc297_h
 */

#include "../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

int main() {
  int n;
  std::cin >> n;
  PowerSeries p(n + 1), q(n + 1);
  for (int i = 1; i <= n; ++i) {
    for (int j = i, k = 1, l = 1; j <= n; j += i, ++k, l = -l) {
      p[j] += l, q[j] += l * k;
    }
  }
  std::cout << (q * (PowerSeries{1} - p).pow(2, n + 1).inverse(n + 1))[n].val();
}
