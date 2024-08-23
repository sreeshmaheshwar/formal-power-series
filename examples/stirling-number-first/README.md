# Signed Stirling Numbers of the First Kind

[Problem Link](https://judge.yosupo.jp/problem/stirling_number_of_the_first_kind).

We are tasked with computing $S(n, k)$ for $0 \leq k \leq n$, where $S(n, k)$ are the signed Stirling numbers of the first kind.

$S(n, k)$ can be defined as the coefficient of $x^k$ in the polynomial $(x)_n$, where $(x)_n$ is the falling factorial:

$$ (x)_n = x(x - 1)(x - 2)\cdots(x - n + 1). $$

**Note** that there is an $O(N \log N)$ solution to this task that we do not describe here.

Instead, we can compute the above falling factorial product of linear polynomials directly through a simple yet elegant divide-and-conquer approach in $O(N \log^2 N)$: recurse, then multiply the left half with the right half.

Concretely, split the product expansion into two equal halves, compute the product of each, and then convolve them together (using NTT). Since the product of the first half has a degree equal to the size of the first half, and similarly for the second half, the convolution operates on polynomials of degree $N/2$ in $O(N \log N)$ time (where $N$ is the size of the current divide-and-conquer range).

Thus, the recurrence is:

$$ T(N) = 2T(N/2) + O(N \log N), $$

which solves to $O(N \log^2 N)$.

In fact, the above divide-and-conquer idea can be used to compute the product of several large integers much more efficiently than naive left-to-right multiplication.

See the source code of `./solution.cpp` for implementation details. 