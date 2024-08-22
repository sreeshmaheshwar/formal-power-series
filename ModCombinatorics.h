#include <cstddef>
#include <vector>

/// Precomputed modular combinatorial quantities.
template <typename ModInt> struct ModCombinatorics {
  std::size_t n;
  std::vector<ModInt> facts;         // Factorials.
  std::vector<ModInt> inverse_facts; // Multiplicative inverses of factorials.
  std::vector<ModInt> inverses;      // Multiplicative inverses.

  /// Computes multiplicative modualr factorials, inverse factorials, and
  /// inverses up to and including `maximum` in 'linear' time (excluding the
  /// cost of computing a stand-alone multiplicative modular inverse directly
  /// via the implementation of `ModInt`).
  explicit constexpr ModCombinatorics(std::size_t maximum)
      : n(maximum + 1), facts(n), inverse_facts(n), inverses(n) {
    facts[0] = 1;
    for (std::size_t i = 1; i < n; ++i) {
      facts[i] = facts[i - 1] * i;
    }
    inverse_facts[n - 1] = ModInt(1) / facts[n - 1]; // ModInt::inv().
    for (std::size_t i = n - 1; i > 0; --i) {
      inverse_facts[i - 1] = inverse_facts[i] * i;
      inverses[i] = facts[i - 1] * inverse_facts[i];
    }
  }
};