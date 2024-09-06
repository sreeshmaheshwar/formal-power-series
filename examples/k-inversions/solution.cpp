// TODO: Problem link?

#include "../../FormalPowerSeries.h"
#include "../../ModCombinatorics.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

int main() {
  int n, k;
  std::cin >> n >> k;

  PowerSeries p(k + 1);
  ModCombinatorics<mint> combinatorics(k + 1);
  for (int i = 1; i <= k; ++i) {
    p[i] += n * combinatorics.inverses[i];
    for (int j = i, c = 1; i <= n && j <= k; j += i, ++c) {
      p[j] -= combinatorics.inverses[c];
    }
  }

  std::cout << p.exp(k + 1)[k].val() << std::endl;
}
