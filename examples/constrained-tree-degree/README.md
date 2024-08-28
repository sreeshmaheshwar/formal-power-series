# Constrained Tree Degree

[Problem Link](https://atcoder.jp/contests/abc303/tasks/abc303_h).

We consider the *exponential* generating function of the quantity in question - in other words, the formal power series $g$ such that

$$ (N - 2)! [x^{N - 2}] g(x) $$

is our answer. As explained in the [official editorial](https://atcoder.jp/contests/abc303/editorial/6456), $g$ is given by $f^N$ where 

$$ f(x) = \sum_{i = 0}^{n - 2} b_i x^i  $$

and $b_i$ is $1 / i!$ if $(i + 1)$ is present in $S$ and $0$ otherwise.

The idea here is that the product of $N$ copies of $f$, from $i = 1$ to $N$, accumulates the number of elements chosen as the power of $x$ and the product of factorials of the counts of each $i$ as the reciprocal of the corresponding coefficient.

Our desired resultant power of $x$ is therefore $(N - 2)$ and we can multiply $(N - 2)!$ by its coefficient to produce the desired sum of all multinomial coefficients, which is our answer. 

This stems from the known result that the number of trees on $N$ vertices with degrees $d_i$ is given by

$$ \binom{N-2}{(d_1 - 1),\dots,(d_N - 1)}. $$

As for the implementation, we find the coefficients $b_i$ directly (having precomputed inverse factorials) and raise the resultant truncated formal power series to the power $N$. This yields an $O(N \log N)$ solution using `FormalPowerSeries::pow`.

See the source code of `./solution.cpp` for more implementation details.
