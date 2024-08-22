# Formal Power Series

Lightweight C++ implementation of various (non-exhaustive) operations on truncated [formal power series](https://en.wikipedia.org/wiki/Formal_power_series) for use in programming contests, an idiomatic interface to them, and some example problems with write-ups/explanations. Centred around interoperability with and delegation to any suitable existing competitive programming library, such as [AtCoder](https://atcoder.jp/)'s popular [AC Library](https://github.com/atcoder/ac-library).

## Usage

The `FormalPowerSeries` class provides a `std::vector`-like interface and requires only a modular integer implementation and a (corresponding) convolution operation (typically, [NTT](https://mathworld.wolfram.com/NumberTheoreticTransform.html)), both common among libraries. The polynomial multiplication implementation underlying the formal power series operations may therefore be switched out at will.

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

## Examples

The `examples` directory contains subdirectories corresponding to example competitive programming problems that can be solved with this library - write-ups/explanations are in each `README.md` and source codes in each `solution.cpp`. These problems were chosen for simple implementations that highlight the library's usage.

Each `.cpp` solution can be compiled from the `examples` directory via `make`, specifically `make single file=<file-path>` which generates the `<file-path>.out` executable. For instance,

```sh
❯ make single file=partition-number/solution.cpp && echo "10" | partition-number/solution.out # First 11 partition numbers
g++ -std=c++20 -Wall -Wextra -Wpedantic -I ../ac-library partition-number/solution.cpp -o partition-number/solution.out
1 1 2 3 5 7 11 15 22 30 42
```
