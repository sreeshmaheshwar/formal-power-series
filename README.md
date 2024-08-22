# Formal Power Series

Lightweight C++ implementation of various (non-exhaustive) operations on truncated formal power series for use in programming contests, and an idiomatic interface to them. Centred around interoperability with and delegation to any suitable existing competitive programming library, such as [AtCoder](https://atcoder.jp/)'s popular [AC Library](https://github.com/atcoder/ac-library).

## Usage

The `FormalPowerSeries` class provides a `std::vector` interface and requires only a modular integer implementation and a (corresponding) convolution operation (typically, [NTT](https://mathworld.wolfram.com/NumberTheoreticTransform.html)), both common among libraries. The polynomial multiplication implementation underlying the formal power series operations may therefore be switched out at will.

The solution to an [AtCoder problem](https://atcoder.jp/contests/abc297/tasks/abc297_h), using this library backed by the [AC Library](https://github.com/atcoder/ac-library), is shown below.

```cpp
#include "FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <iostream>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution<mint>(a, b);
}>;

int main() {
  int n;
  std::cin >> n;

  PowerSeries p(n + 1), q(n + 1);
  for (int i = 1; i <= n; ++i) {
    for (int j = i, k = 1, l = 1; j <= n; j += i, ++k, l = -l) {
      p[j] += l, q[j] += l * k;
    }
  }

  std::cout << (q * (PowerSeries{1} - p).pow(2, n + 1).inverse(n + 1))[n].val();
}

```
