#pragma once

#include <cstddef>
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
/// std::vector of coefficients of size (n + 1) whose i-th element is the
/// coefficient of x^i. For example, {1, 2, 0, 4} represents the polynomial 1 +
/// 2x + 4x^3.
template <typename ModInt, ConvolutionFunction<ModInt> auto Convolution>
class FormalPowerSeries : public std::vector<ModInt> {
public:
  using Base = std::vector<ModInt>;

  constexpr FormalPowerSeries() noexcept = default;

  explicit constexpr FormalPowerSeries(std::size_t n);

  explicit constexpr FormalPowerSeries(std::size_t n, const ModInt &value);

  explicit constexpr FormalPowerSeries(const std::vector<ModInt> &vec);

  constexpr FormalPowerSeries(const std::initializer_list<ModInt> &list);

  constexpr FormalPowerSeries(std::vector<ModInt> &&vec) noexcept;

  template <typename Iter> constexpr FormalPowerSeries(Iter first, Iter last);

  constexpr FormalPowerSeries(const FormalPowerSeries &) = default;

  constexpr FormalPowerSeries(FormalPowerSeries &&) noexcept = default;

  constexpr FormalPowerSeries &operator=(const FormalPowerSeries &) = default;

  constexpr FormalPowerSeries &
  operator=(FormalPowerSeries &&) noexcept = default;

  constexpr FormalPowerSeries operator+(const FormalPowerSeries &) const;

  constexpr FormalPowerSeries operator-(const FormalPowerSeries &) const;

  constexpr FormalPowerSeries operator*(const FormalPowerSeries &) const;

  constexpr FormalPowerSeries &operator*=(const ModInt &scalar);

  /// Returns the formal power series consisting of the first `size` terms of
  /// this formal power series.
  [[nodiscard]] constexpr FormalPowerSeries take(std::size_t size) const;

  /// Returns the derivative of this formal power series.
  [[nodiscard]] constexpr FormalPowerSeries derivative() const;

  /// Returns the anti-derivative of this formal power series.
  [[nodiscard]] constexpr FormalPowerSeries antiderivative() const;

  /// Returns the first `size` terms of the formal power series that is the
  /// natural logarithm of this formal power series.
  /// Precondition: this polynomial is non-empty with constant term of one.
  [[nodiscard]] constexpr FormalPowerSeries log(std::size_t size) const;

  /// Returns the first `size` terms of the formal power series that is the
  /// multiplicative inverse of this formal power series.
  /// Precondition: this polynomial is non-empty with a non-zero constant term.
  [[nodiscard]] constexpr FormalPowerSeries inverse(std::size_t size) const;

  /// Returns the first `size` terms of the formal power series that is e raised
  /// to the power of this formal power series.
  /// Precondition: this polynomial is non-empty with a zero constant term.
  [[nodiscard]] constexpr FormalPowerSeries exp(std::size_t size) const;

  /// Returns the first `size` terms of the formal power series that is this
  /// formal power series raised to the power of `k`, where `k` is a
  /// non-negative integer.
  [[nodiscard]] constexpr FormalPowerSeries pow(std::uint64_t k,
                                                std::size_t size) const;

  /// Returns the first `size` terms of the formal power series that is this
  /// formal power series raised to the power of `k`, where `k` is a
  /// non-negative integer, using naive binary exponentiation in
  /// O(C(size) * log K) time, where C(N) is the time complexity of convolution.
  /// Generally slower than `FormalPowerSeries::pow` when C(N) is O(N log N).
  [[nodiscard]] constexpr FormalPowerSeries bin_pow(std::uint64_t k,
                                                    std::size_t size) const;

  /// Returns the first `size` terms of the formal power series P(x) = 1.
  [[nodiscard]] static constexpr FormalPowerSeries
  mult_identity(std::size_t size);

  template <typename ModIntT, ConvolutionFunction<ModIntT> auto ConvolutionT>
  constexpr friend FormalPowerSeries<ModIntT, ConvolutionT>
  operator*(const FormalPowerSeries<ModIntT, ConvolutionT> &fps,
            const ModIntT &scalar) {
    return FormalPowerSeries(fps) *= scalar;
  }

  template <typename ModIntT, ConvolutionFunction<ModIntT> auto ConvolutionT>
  constexpr friend FormalPowerSeries<ModIntT, ConvolutionT>
  operator*(const ModIntT &scalar,
            const FormalPowerSeries<ModIntT, ConvolutionT> &fps) {
    return fps * scalar;
  }
};

#include "FormalPowerSeries.cpp" // Templated class, so include implementation.
