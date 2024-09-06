#include "../../FormalPowerSeries.h"
#include "../../ModCombinatorics.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>
#include <vector>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

int main() {
  int n, k;
  std::cin >> n >> k;

  ModCombinatorics<mint> combinatorics(n - 2);

  PowerSeries p(n - 1);
  for (int i = 0, s; i < k; ++i) {
    std::cin >> s;
    s -= 1;
    p[s] = combinatorics.inverse_facts[s];
  }

  std::cout << (combinatorics.facts[n - 2] * p.pow(n, n - 1)[n - 2]).val();
}
