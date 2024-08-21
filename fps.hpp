#include <algorithm>
#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <type_traits>
#include <vector>

template <typename F, typename T>
concept ConvolutionFunction = requires(F f, const std::vector<T> &a,
                                       const std::vector<T> &b) {
  { f(a, b) } -> std::same_as<std::vector<T>>;
};

/// Formal Power Series operations that rely on a provided convolution function
/// to multiply polynomials. A polynomial of degree n is represented as a
/// std::vector of coefficients of size (n + 1), where the i-th element is the
/// coefficient of x^i. For example, {1, 2, 0, 4} represents the polynomial 1 +
/// 2x + 4x^3.
template <typename ModInt, ConvolutionFunction<ModInt> auto Convolution>
class FormalPowerSeries : public std::vector<ModInt> {
public:
  using Base = std::vector<ModInt>;

  constexpr FormalPowerSeries() noexcept = default;

  explicit constexpr FormalPowerSeries(size_t n) : Base(n) {}

  explicit constexpr FormalPowerSeries(size_t n, const ModInt &value)
      : Base(n, value) {}

  explicit constexpr FormalPowerSeries(const std::vector<ModInt> &vec)
      : Base(vec) {}

  constexpr FormalPowerSeries(const std::initializer_list<ModInt> &list)
      : Base(list) {}

  constexpr FormalPowerSeries(std::vector<ModInt> &&vec) noexcept
      : Base(std::move(vec)) {}

  template <typename Iter>
  constexpr FormalPowerSeries(Iter first, Iter last) : Base(first, last) {}

  constexpr FormalPowerSeries(const FormalPowerSeries &) = default;

  constexpr FormalPowerSeries(FormalPowerSeries &&) noexcept = default;

  constexpr FormalPowerSeries &operator=(const FormalPowerSeries &) = default;

  constexpr FormalPowerSeries &
  operator=(FormalPowerSeries &&) noexcept = default;

  constexpr FormalPowerSeries operator*=(const ModInt &scalar) {
    for (auto &x : *this) {
      x *= scalar;
    }
    return *this;
  }

  constexpr friend FormalPowerSeries operator*(const FormalPowerSeries &fps,
                                               const ModInt &scalar) {
    return FormalPowerSeries(fps) *= scalar;
  }

  constexpr friend FormalPowerSeries operator*(const ModInt &scalar,
                                               const FormalPowerSeries &fps) {
    return fps * scalar;
  }

  constexpr FormalPowerSeries operator*(const FormalPowerSeries &other) const {
    return FormalPowerSeries(Convolution(*this, other));
  }

  constexpr FormalPowerSeries operator+(const FormalPowerSeries &other) const {
    const auto n = std::max(this->size(), other.size());
    FormalPowerSeries result(n);
    for (size_t i = 0; i < n; ++i) {
      if (i < this->size()) {
        result[i] += (*this)[i];
      }
      if (i < other.size()) {
        result[i] += other[i];
      }
    }
    return result;
  }

  constexpr FormalPowerSeries operator-(const FormalPowerSeries &other) const {
    const auto n = std::max(this->size(), other.size());
    FormalPowerSeries result(n);
    for (size_t i = 0; i < n; ++i) {
      if (i < this->size()) {
        result[i] += (*this)[i];
      }
      if (i < other.size()) {
        result[i] -= other[i];
      }
    }
    return result;
  }

  /// Returns the formal power series consisting of the first `size` terms of
  /// this formal power series.
  constexpr FormalPowerSeries take(size_t size) const {
    FormalPowerSeries result(*this);
    result.resize(size);
    return result;
  }

  /// Returns the derivative of this formal power series.
  constexpr FormalPowerSeries derivative() const {
    if (this->empty()) {
      return *this;
    }
    FormalPowerSeries result(this->size() - 1);
    for (size_t i = 1; i < this->size(); ++i) {
      result[i - 1] = (*this)[i] * ModInt(i);
    }
    return result;
  }

  /// Returns the anti-derivative of this formal power series.
  constexpr FormalPowerSeries antiderivative() const {
    FormalPowerSeries result(this->size() + 1);
    for (size_t i = 1; i < result.size(); ++i) {
      result[i] = (*this)[i - 1] / ModInt(i);
    }
    return result;
  }

  /// Returns the first `size` terms of the formal power series that is the
  /// natural logarithm of this formal power series.
  constexpr FormalPowerSeries log(size_t size) const {
    assert(!this->empty() && this->front() == ModInt(1));
    // d/dx (ln P(x)) = P'(x) / P(x).
    return (derivative() * inverse(size)).antiderivative().take(size);
  }

  /// Returns the first `size` terms of the formal power series that is the
  /// multiplicative inverse of this formal power series.
  constexpr FormalPowerSeries inverse(size_t size) const {
    assert(!this->empty() && this->front() != ModInt(0));
    // Newton's Method: Q_{k+1} = Q_k - F(Q_k) / F'(Q_k) (mod x^{2^{k+1}}).
    //
    // Since Q(x), the true inverse of our formal power series P(x), satisfies
    // Q(x) = 1 / P(x) and so P(x) = 1 / Q(x), we take F(Q) = 1 / Q - P = 0.
    //
    // Thus, F'(Q) = -1 / Q^2, and the Newton iteration becomes Q_{k+1} = Q_k *
    // (2 - P * Q_k) (mod x^{2^{k+1}}).
    //
    // As we assert a non-zero constant term of P(x), we can take the
    // multiplicative inverse of the constant term of P(x) as the initial Q_0
    // since it is the constant term of P(x)^{-1}.
    FormalPowerSeries res = {ModInt(1) / this->front()};
    while (res.size() != size) {
      const auto next_size = std::min(res.size() * 2, size);
      res = (res * (FormalPowerSeries{ModInt(2)} - take(next_size) * res))
                .take(next_size);
    }
    return size ? res : FormalPowerSeries();
  }

  /// Returns the first `size` terms of the formal power series that is e raised
  /// to the power of this formal power series.
  constexpr FormalPowerSeries exp(size_t size) const {
    assert(!this->empty() && this->front() == ModInt(0));
    // Newton's Method: Q_{k+1} = Q_k - F(Q_k) / F'(Q_k) (mod x^{2^{k+1}}).
    //
    // Since Q(x), the true exponential, satisfies Q(x) = e^{P(x)} and so P(x)
    // = ln Q(x), we take F(Q) = ln(Q) - P = 0.
    //
    // Thus, F'(Q) = 1 / Q, and the Newton iteration becomes Q_{k+1} = Q_k * (1
    // + P - ln(Q_k)) (mod x^{2^{k+1}}).
    //
    // As we assert a zero constant term of P(x), we can take 1 as the initial
    // Q_0 since it is the constant term of e^{P(x)}.
    FormalPowerSeries res = {ModInt(1)};
    while (res.size() != size) {
      const auto next_size = std::min(res.size() * 2, size);
      res = (res * (FormalPowerSeries{ModInt(1)} + take(next_size) -
                    res.log(next_size)))
                .take(next_size);
    }
    return size ? res : FormalPowerSeries();
  }

  /// Returns the first `size` terms of the formal power series that is this
  /// formal power series raised to the power of `k`, where `k` is a
  /// non-negative integer, in O(C(size)) time (excluding the cost of
  /// computing a multiplicative modular inverse), where C(N) is the time
  /// complexity of convolution.
  constexpr FormalPowerSeries pow(std::uint64_t k, size_t size) const {
    // We make no assumptions about the polynomial, unlike in other methods.
    //
    // If its constant term were 1, (indeed, as `FormalPowerSeries::log`
    // requires), we may delegate to `log` and `exp`, seeing as P^k(x) = exp(k *
    // ln(P(x))).
    //
    // So, we instead extract out a factor to exponentiate a polynomial with a
    // constant term of 1. Concretely, we find the first non-zero coefficient
    // index, in other words the first i such that a = [x^i]P(x) is non-zero,
    // then P(x) = a * x^i * Q(x) for some Q(x) with Q(0) = 1. Then, P^k(x) =
    // a^k * x^{ik} * Q^k(x).
    if (k == 0) {
      return FormalPowerSeries::mult_identity(size);
    }

    size_t i = 0;
    while (i < this->size() && (*this)[i] == ModInt(0)) {
      ++i;
    }
    // As our answer has a factor of x^{ik}, its first i * k coefficients are
    // zero. We can therefore exit early if i * k >= size.
    if (i == this->size() /* <-> P(x) = 0 */ ||
        i >= (size + k - 1) / k /* <-> i * k >= size */) {
      return FormalPowerSeries(size);
    } // i * k < size.

    const auto a = (*this)[i];
    // Q(x) = (P(x) / a) / x^i. Note that dividing by x^i is a left shift.
    auto q = *this * (ModInt(1) / a);
    q.erase(q.begin(), q.begin() + i); // Left shift.

    // P^k(x) = a^k * x^{ik} * Q^k(x). Multiplication by x^{ik} is a shift right
    // by i * k, meaning we only need to compute the first (size - i * k) terms
    // of Q^k(x).
    size_t n = size - i * k;
    q = (q.log(n) * k).exp(n) * a.pow(k);
    q.insert(q.begin(), i * k, ModInt(0)); // Right shift, pad with zeros.
    assert(q.size() == size);

    return q;
  }

  /// Returns the first `size` terms of the formal power series that is this
  /// formal power series raised to the power of `k`, where `k` is a
  /// non-negative integer, using naive binary exponentiation in
  /// O(C(size) * log K) time, where C(N) is the time complexity of convolution.
  /// Asymptotically slower than `FormalPowerSeries::pow`.
  constexpr FormalPowerSeries bin_pow(std::uint64_t k, size_t size) const {
    FormalPowerSeries result = FormalPowerSeries::mult_identity(size);
    FormalPowerSeries power = this->take(size);
    while (k > 0) {
      if (k & 1) {
        result = (result * power).take(size);
      }
      power = (power * power).take(size);
      k >>= 1;
    }
    return result;
  }

  /// Returns the first `size` terms of the formal power series P(x) = 1.
  static constexpr FormalPowerSeries mult_identity(size_t size) {
    if (!size) {
      return {};
    }
    FormalPowerSeries result(size);
    result[0] = ModInt(1);
    return result;
  }
};