# Counting Subset Sums

[Problem Link](https://judge.yosupo.jp/problem/sharp_p_subset_sum).

Given a sequence $S$, we aim to count its subsequences with sum $t$, for all $1 \le t \le T$ in $O(T \log T)$ time (and modulo $998244353$).

The ordinary generating function of the answers is given by the product of $(1 + x^{S_i})$, since for each element of $S$, we can either include it in the subsequence, increasing its sum by $S_i$, or not include it, leaving the sum unchanged.

To compute this efficiently, let $c_i$ be the count of $i$ in $S$, then we can group terms together to express our answer as the product of $(1 + x^i)^{c_i}$ instead.

To compute this new product efficiently, we first take its natural
logarithm to go from a product to a sum, then raise $e$ to the power of the
resultant formal power series.

We therefore sum $\ln((1 + x^i)^{c_i}) = c_i \ln(1 + x^i)$ over $i$. Note that

$$ c_i \log(1 + x^i) = c_i x^i - \frac{c_i}{2} x^{2i} + \frac{c_i}{3} x^{3i} - \dots, $$

so we can simply add these coefficients to their corresponding indices of a *running sum* polynomial, after precomputing multiplicative modular inverses.

This iteration will be $O(T \log T)$ overall, since for all $i$ until $T$, we iterate through all multiples of $i$ that are at most $T$ (recall that we only need to compute the terms of this formal power series up to $x^T$) - the well-known ["harmonic series trick"](https://discuss.codechef.com/t/more-intuitive-explanation-for-the-harmonic-seriess-sum/67287). The consequent `FormalPowerSeries::exp` invocation is also $O(T \log T)$, and so this is our complexity.

See the source code of `./solution.cpp` for more implementation details.
