// https://judge.yosupo.jp/problem/exp_of_formal_power_series

#include "../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n;
  std::cin >> n;

  PowerSeries a(n);
  for (int i = 0; i < n; ++i) {
    int x;
    std::cin >> x;
    a[i] = x;
  }

  for (const auto &x : a.exp(n)) {
    std::cout << x.val() << ' ';
  }
}