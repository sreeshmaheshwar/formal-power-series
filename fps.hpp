#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <type_traits>
#include <vector>

template <typename F, typename T>
concept ConvolutionFunction = requires(F f, const std::vector<T> &a,
                                       const std::vector<T> &b) {
  { f(a, b) } -> std::same_as<std::vector<T>>;
};

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
    // We take the multiplicative inverse of the constant term of P(x) as the
    // initial Q_0.
    FormalPowerSeries res = {ModInt(1) / this->front()};
    while (res.size() != size) {
      const auto next_size = std::min(res.size() * 2, size);
      res = (res * (FormalPowerSeries{2} - take(next_size) * res))
                .take(next_size);
    }
    return res;
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
    // Since the exponential series has a constant term of 1, we take the
    // initial Q_0 to be 1.
    FormalPowerSeries res = {ModInt(1)};
    while (res.size() != size) {
      const auto next_size = std::min(res.size() * 2, size);
      res =
          (res * (FormalPowerSeries{1} + take(next_size) - res.log(next_size)))
              .take(next_size);
    }
    return res;
  }
};