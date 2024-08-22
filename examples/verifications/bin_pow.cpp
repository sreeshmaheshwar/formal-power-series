// https://judge.yosupo.jp/problem/pow_of_formal_power_series

#include "../../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <cstdint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n;
  std::int64_t m;
  std::cin >> n >> m;

  PowerSeries a(n);
  for (int i = 0; i < n; ++i) {
    int x;
    std::cin >> x;
    a[i] = x;
  }

  for (const auto &x : a.bin_pow(m, n)) {
    std::cout << x.val() << ' ';
  }
}