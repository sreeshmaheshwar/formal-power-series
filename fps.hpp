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

  constexpr FormalPowerSeries take(size_t size) const {
    FormalPowerSeries result(*this);
    result.resize(size);
    return result;
  }

  constexpr FormalPowerSeries inverse(size_t size) const {
    assert(!this->empty() && this->front() != ModInt(0));
    FormalPowerSeries res = {ModInt(1) / this->front()};
    while (res.size() != size) {
      const auto next_size = std::min(res.size() * 2, size);
      res = (res * (FormalPowerSeries{2} - take(next_size) * res))
                .take(next_size);
    }
    return res;
  }
};