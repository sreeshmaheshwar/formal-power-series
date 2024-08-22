# Formal Power Series

Lightweight C++ implementation of various (non-exhaustive) operations on truncated [formal power series](https://en.wikipedia.org/wiki/Formal_power_series) for use in programming contests, an idiomatic interface to them, and some example problems with write-ups/explanations. Centred around interoperability with and delegation to any suitable existing competitive programming library, such as [AtCoder](https://atcoder.jp/)'s popular [AC Library](https://github.com/atcoder/ac-library) (ACL).

## Usage

The `FormalPowerSeries` class provides a `std::vector`-like interface and requires only a modular integer implementation and a (corresponding) convolution operation (typically, [NTT](https://mathworld.wolfram.com/NumberTheoreticTransform.html)), both common among libraries. The polynomial multiplication implementation underlying the formal power series operations may therefore be switched out at will.

The solution to an [AtCoder problem](https://atcoder.jp/contests/abc297/tasks/abc297_h), using this library backed by ACL, is shown below.

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

The `examples` directory contains subdirectories corresponding to example competitive programming problems that can be solved with this library. These tasks were chosen for simple implementations that highlight the library's usage.

Apart from `verifications`, that contains [Library Checker](https://judge.yosupo.jp/) verification submissions for individual formal power series operations, each such subdirectory has the following structure:

- `README.md`: a write-up / solution explanation of the problem.
-  `solution.cpp`: source code of the solution to the problem, using this library together with ACL.

Before compiling any examples, remember to instantiate submodules via the one-off 

```
git submodule init 
git submodule update
```

Each solution can be then compiled from the `examples` directory via `make`. Specifically running `make single file=<file-path>` generates the `<file-path>.out` executable. For instance,

```sh
❯ make single file=partition-number/solution.cpp && echo "10" | partition-number/solution.out # First 11 partition numbers
g++ -std=c++20 -Wall -Wextra -Wpedantic -I ../ac-library partition-number/solution.cpp -o partition-number/solution.out
1 1 2 3 5 7 11 15 22 30 42
```

## Submission

In competitive programming, a single, self-contained source file is typically submitted to the judge. Bundling tools such as [OJ-Bundle](https://github.com/online-judge-tools/verification-helper) are therefore commonly used to _"expand"_ out `#include`s of a source file where relevant, producing a single, submission-ready output. This tool is compatible with this library's headers. ACL's [expander.py](https://github.com/atcoder/ac-library/blob/master/expander.py) provides similar functionality but for ACL headers.

Install the `oj-bundle` CLI tool can be installed with `pip3 install online-judge-verify-helper` (for details, see the [repository](https://github.com/online-judge-tools/verification-helper)) and be sure that it is in your `PATH` (append the installation location to it if not). This should then allow code that uses this library's headers to be expanded.

For instance, an `examples` solution that uses this library with ACL can be expanded into the submission-ready `combined.cpp` file by running (from the `examples` directory):

```sh
❯ oj-bundle partition-number/solution.cpp > combined.cpp && python3 ../ac-library/expander.py --lib=../ac-library combined.cpp
```

## Notes

- There are formal power series operations required by some competitive programming problems that are not yet supported - for example, finding the square root of a formal power series or composing two together. Moreover, *sparse* variants (meaning, on large polynomials with comparatively few non-zero coefficients) of the operations that _are_ supported have not yet been implemented.
- Inspecting [Library Checker](https://judge.yosupo.jp/) submissions shows other implementations of operations being faster in practice. We rely on Newton's method for efficient - $O(N \log N)$, assuming $O(N \log N)$ convolution - yet simple implementations, but it would seem that other methods occasionally have better constant factors. Note that in some cases, different NTT performance is the culprit.
