#include "FormalPowerSeries.h"
#include <atcoder/convolution>
#include <atcoder/modint>
#include <cstddef>
#include <gtest/gtest.h>
#include <ostream>
#include <vector>

using mint = atcoder::modint998244353;
using PowerSeries = FormalPowerSeries<mint, [](const auto &a, const auto &b) {
  return atcoder::convolution(a, b);
}>;

class FormalPowerSeriesTest : public ::testing::Test {
protected:
  void check_content(const PowerSeries &p, const std::vector<mint> &expected) {
    ASSERT_EQ(p.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      EXPECT_EQ(p[i], expected[i]);
    }
  }
};

TEST_F(FormalPowerSeriesTest, Constructors) {
  const std::vector<mint> example = {mint(1), mint(2), mint(3)};
  const std::size_t example_size = 3;

  check_content(PowerSeries{}, {});
  check_content(PowerSeries(example_size),
                std::vector<mint>(example_size, mint(0)));
  check_content(PowerSeries(example_size, mint(5)),
                std::vector<mint>(example_size, mint(5)));
  check_content(PowerSeries(example), example);
  check_content(PowerSeries{mint(1), mint(2), mint(3)}, example);
  check_content(PowerSeries(example.begin(), example.end()), example);

  PowerSeries q{mint(1), mint(2), mint(3)};
  check_content(PowerSeries(q), example);
  check_content(PowerSeries(std::move(PowerSeries(example))), example);

  PowerSeries p;
  p = q;
  check_content(p, example);

  p = PowerSeries(example);
  check_content(p, example);
}

TEST_F(FormalPowerSeriesTest, Take) {
  PowerSeries p{1, 2, 3, 4, 5};

  check_content(p.take(3), {1, 2, 3});
  check_content(p.take(5), {1, 2, 3, 4, 5});
  check_content(p.take(0), {});
  check_content(p.take(7), {1, 2, 3, 4, 5, 0, 0});
}

TEST_F(FormalPowerSeriesTest, Derivative) {
  PowerSeries p{1, 2, 3, 4, 5};
  check_content(p.derivative(), {2, 6, 12, 20});

  PowerSeries empty;
  check_content(empty.derivative(), {});
}

TEST_F(FormalPowerSeriesTest, Antiderivative) {
  PowerSeries p{1, 2, 3, 4};
  check_content(p.antiderivative(), {0, 1, 1, 1, 1});

  PowerSeries empty;
  check_content(empty.antiderivative(), {0});
}

TEST_F(FormalPowerSeriesTest, BasicArithmetic) {
  PowerSeries p{1, 2, 3};
  PowerSeries q{4, 5, 6, 7};

  check_content(p + q, {5, 7, 9, 7});
  check_content(p - q, {-3, -3, -3, -7});

  mint scalar = 2;
  check_content(p * scalar, {2, 4, 6});
  check_content(scalar * p, {2, 4, 6});
}

TEST_F(FormalPowerSeriesTest, Multiplication) {
  PowerSeries p{1, 2};
  PowerSeries q{3, 4, 5};

  check_content(p * q, {3, 10, 13, 10});
}

TEST_F(FormalPowerSeriesTest, LogPrecondition) {
  PowerSeries valid{1, 2, 3};
  EXPECT_NO_THROW(valid.log(3));

  PowerSeries invalid_empty;
  EXPECT_DEATH(invalid_empty.log(3), "");

  PowerSeries invalid_constant{2, 3, 4};
  EXPECT_DEATH(invalid_constant.log(3), "");
}

TEST_F(FormalPowerSeriesTest, InversePrecondition) {
  PowerSeries valid{1, 2, 3};
  EXPECT_NO_THROW(valid.inverse(3));

  PowerSeries invalid_empty;
  EXPECT_DEATH(invalid_empty.inverse(3), "");

  PowerSeries invalid_zero_constant{0, 1, 2};
  EXPECT_DEATH(invalid_zero_constant.inverse(3), "");
}

TEST_F(FormalPowerSeriesTest, ExpPrecondition) {
  PowerSeries valid{0, 1, 2};
  EXPECT_NO_THROW(valid.exp(3));

  PowerSeries invalid_empty;
  EXPECT_DEATH(invalid_empty.exp(3), "");

  PowerSeries invalid_nonzero_constant{1, 2, 3};
  EXPECT_DEATH(invalid_nonzero_constant.exp(3), "");
}

TEST_F(FormalPowerSeriesTest, MultIdentity) {
  check_content(PowerSeries::mult_identity(0), {});
  check_content(PowerSeries::mult_identity(1), {1});
  check_content(PowerSeries::mult_identity(3), {1, 0, 0});
}

TEST_F(FormalPowerSeriesTest, InverseSamples) {
  PowerSeries p{5, 4, 3, 2, 1};
  check_content(p.inverse(5),
                {598946612, 718735934, 862483121, 635682004, 163871793});

  // Empty result edge case:
  check_content(p.inverse(0), {});
}

TEST_F(FormalPowerSeriesTest, ExpSamples) {
  PowerSeries p{0, 1, 2, 3, 4};
  check_content(p.exp(5), {1, 1, 499122179, 166374064, 291154613});

  // Empty result edge case:
  check_content(p.exp(0), {});
}

TEST_F(FormalPowerSeriesTest, LogSamples) {
  PowerSeries p{1, 1, 499122179, 166374064, 291154613};
  check_content(p.log(5), {0, 1, 2, 3, 4});

  // Empty result edge case:
  check_content(p.log(0), {});
}

// To reduce duplication between testing `FormalPowerSeries::pow` and
// `FormalPowerSeries::bin_pow`, we use a value-parameterized test suite.
struct PowerMethodParam {
  using power_func_t = std::function<PowerSeries(const PowerSeries &,
                                                 std::uint64_t, std::size_t)>;

  power_func_t power;
  std::string name;

  // Required for readable instantiated test names.
  friend std::ostream &operator<<(std::ostream &os,
                                  const PowerMethodParam &param) {
    return os << param.name;
  }
};

class FormalPowerSeriesPowerTest
    : public FormalPowerSeriesTest,
      public ::testing::WithParamInterface<PowerMethodParam> {
protected:
  PowerMethodParam::power_func_t power;

public:
  void SetUp() override { power = GetParam().power; }
};

TEST_P(FormalPowerSeriesPowerTest, BinomialSamples) {
  PowerSeries p{1, 1};
  check_content(power(p, 2, 3), {1, 2, 1});
  check_content(power(p, 5, 3), {1, 5, 10});
  check_content(power(p, 2, 0), {});
}

TEST_P(FormalPowerSeriesPowerTest, ZeroPower) {
  PowerSeries p{1, 2, 3};
  check_content(power(p, 0, 3), {1, 0, 0});
}

TEST_P(FormalPowerSeriesPowerTest, LeadingZeroes) {
  PowerSeries p{0, 0, 9, 12};
  check_content(power(p, 1, 5), {0, 0, 9, 12, 0});
  check_content(power(p, 3, 4), {0, 0, 0, 0});
}

TEST_P(FormalPowerSeriesPowerTest, ZeroBaseEdgeCases) {
  PowerSeries p{2, 0}, empty;
  check_content(power(p, 0, 3), {1, 0, 0});
  check_content(power(empty, 0, 3), {1, 0, 0});
  check_content(power(empty, 2, 3), {0, 0, 0});
}

INSTANTIATE_TEST_SUITE_P(
    PowerMethods, FormalPowerSeriesPowerTest,
    ::testing::Values(
        PowerMethodParam{[](const PowerSeries &p, std::uint64_t n,
                            std::size_t deg) { return p.pow(n, deg); },
                         "Pow"},
        PowerMethodParam{[](const PowerSeries &p, std::uint64_t n,
                            std::size_t deg) { return p.bin_pow(n, deg); },
                         "BinPow"}));