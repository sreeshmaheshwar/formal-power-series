/**
 * Solution to:
 * https://judge.yosupo.jp/problem/sharp_p_subset_sum
 *
 * Given a sequence S, count its subsequences with sum t, for all 1 <= t <= T
 * in O(T log T) (and mod 998244353).
 *
 * The ordinary generating function of the answers is given by the product of (1
 * + x^{S_j}), since for each element of S, we can either include it in the
 * subsequence, increasing the sum by S_j, or not include it, leaving the sum
 * unchanged.
 *
 * To compute this efficiently, let c_i be the count of i in S, then we can
 * group terms together to get the product of (1 + x^i)^{c_i} instead.
 *
 * To compute this new product efficiently, we first take its (natural)
 * logarithm to go from a product to a sum, then raise e to the power of the
 * resultant formal power series.
 *
 * We therefore sum log((1 + x^i)^{c_i}) = c_i log(1 + x^i) over i. Note that
 *
 *    c_i log(1 + x^i) = c_i x^i - c_i x^{2i} / 2 + c_i x^{3i} / 3 - ...
 *
 * so we can simply add these coefficients to their corresponding indices of
 * a "running sum" polynomial, after precomputing multiplicative modular
 * inverses.
 *
 * This iteration will be O(T log T) overall due to us enumerating over
 * multiples of i for all i till T (the "harmonic series trick").
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
    facts[0] = inverse_facts[0] = 1;
    for (std::size_t i = 1; i < n; ++i) {
      facts[i] = facts[i - 1] * i;
    }
    inverse_facts[n - 1] = facts[n - 1].inv();
    for (std::size_t i = n - 2; i > 0; --i) {
      inverse_facts[i] = inverse_facts[i + 1] * (i + 1);
      inverses[i] = facts[i - 1] * inverse_facts[i];
    }
  }
};

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n, t;
  std::cin >> n >> t;

  Combinatorics combinatorics(t + 1);

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
