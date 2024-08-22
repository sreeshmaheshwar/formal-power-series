#include "../../FormalPowerSeries.h"
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
