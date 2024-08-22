# Formal Power Series for Competitive Programming

Lightweight implementations of various (non-exhaustive) operations on truncated formal power series for use in programming contests, and an idiomatic interface to them. Centred around interoperability with and delegation to any suitable existing competitive programming library, such as [AtCoder's](https://atcoder.jp/) popular [AC Library](https://github.com/atcoder/ac-library).

Concretely, the templated `FormalPowerSeries` class, that provides a `std::vector`-like programmatic interface, requires only implementations of modular integers and the corresponding convolution operation - typically, [NTT](https://mathworld.wolfram.com/NumberTheoreticTransform.html).
