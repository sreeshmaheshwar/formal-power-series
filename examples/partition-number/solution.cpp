#include "../../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

// Computes, up to x^n, the product (1 - x^i) from i = 1 to inf.
// https://en.wikipedia.org/wiki/Pentagonal_number_theorem
PowerSeries pentagonal_series(std::size_t n) {
  PowerSeries pentagonal(n + 1);
  pentagonal[0] = 1;

  for (int i = 1;; ++i) {
    int64_t j1 = 1LL * i * (3 * i - 1) / 2;
    int64_t j2 = 1LL * (-i) * (3 * (-i) - 1) / 2;

    if (j1 > n && j2 > n) {
      break;
    }
    if (j1 <= n) {
      pentagonal[j1] = (i & 1) ? -1 : 1;
    }
    if (j2 <= n) {
      pentagonal[j2] = (i & 1) ? -1 : 1;
    }
  }

  return pentagonal;
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n;
  std::cin >> n;
  for (const auto &x : pentagonal_series(n).inverse(n + 1)) {
    std::cout << x.val() << ' ';
  }
}