/**
 * Solution to:
 * https://atcoder.jp/contests/abc303/tasks/abc303_h
 */

#include "../FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <cstddef>
#include <iostream>
#include <vector>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

struct Combinatorics {
  std::size_t n;
  std::vector<mint> facts;         // Factorials.
  std::vector<mint> inverse_facts; // Multiplicative inverses of factorials.
  std::vector<mint> inverses;      // Multiplicative inverses.

  Combinatorics(std::size_t maximum)
      : n(maximum + 1), facts(n), inverse_facts(n), inverses(n) {
    facts[0] = 1;
    for (std::size_t i = 1; i < n; ++i) {
      facts[i] = facts[i - 1] * i;
    }
    inverse_facts[n - 1] = facts[n - 1].inv();
    for (std::size_t i = n - 1; i > 0; --i) {
      inverse_facts[i - 1] = inverse_facts[i] * i;
      inverses[i] = facts[i - 1] * inverse_facts[i];
    }
  }
};

int main() {
  int n, k;
  std::cin >> n >> k;

  Combinatorics combinatorics(n - 2);

  PowerSeries p(n - 1);
  for (int i = 0, s; i < k; ++i) {
    std::cin >> s;
    s -= 1;
    p[s] = combinatorics.inverse_facts[s];
  }

  std::cout << (combinatorics.facts[n - 2] * p.pow(n, n - 1)[n - 2]).val();
}
