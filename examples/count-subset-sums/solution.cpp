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
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n, t;
  std::cin >> n >> t;

  ModCombinatorics<mint> combinatorics(t);

  std::vector<int> freq(t + 1);
  for (int i = 0, s; i < n; ++i) {
    std::cin >> s;
    freq[s] += 1;
  }

  PowerSeries p(t + 1);
  for (int i = 1; i <= t; ++i) {
    for (int j = i, k = 1, l = 1; j <= t; j += i, ++k, l = -l) {
      p[j] += freq[i] * l * combinatorics.inverses[k];
    }
  }

  p = p.exp(t + 1);
  for (int i = 1; i <= t; ++i) {
    std::cout << p[i].val() << " ";
  }
}
