# Diff Adjacent

[Problem Link](https://atcoder.jp/contests/abc297/tasks/abc297_h).

As explained thoroughly in the [official editorial](https://atcoder.jp/contests/abc297/editorial/6199), the ordinary generating function of the quantity in question can be expressed as 

$$ \frac{Q(x)}{1 - P(x)^2}, $$

where $Q$ and $P$ are formal power series given by

$$ Q(x) = \sum_{i = 1}^\infty \frac {x^i} {(1 + x^i)^2} $$

and

$$ P(x) = \sum_{i = 1}^\infty \frac {x^i} {1 + x^i}. $$

Computing the first $(N + 1)$ terms of $Q$ and $P$, taking the inverse of the square of $1$ minus latter, and multiplying by the former then yields the result.

The problem is therefore reduced to finding the truncated $P$ and $Q$. The trick with both tasks is to expand each summand into a formal power series itself.

The expansion of the $i$'th term of $P$'s summation relies on a simple geometric series expansion of $1/(1+x^i)$, multiplied by $x^i$:

$$ \frac {x^i} {1 + x^i} = x^i \cdot \frac {1} {1 - (-x^i)} = \sum_{k = 1}^\infty (-1)^{k+1} x^{ik}.$$

Now for $Q$, differentiating the standard geometric series result yields a well-known expansion for $- 1 / (1 + x^i)^2$ given by $-1 + 2x^i - 3x^{2i} + \cdots$, and multiplying by $-x^i$ yields the following expansion of the $i$'th term of $Q$'s summation:

$$ \frac {x^i} {(1 + x^i)^2} = \sum_{k = 1}^\infty (-1)^{k+1} k x^{ik}. $$

Therefore,

$$ P(x) =  \sum_{i = 1}^\infty \sum_{k = 1}^\infty (-1)^{k+1} x^{ik}, $$

and 

$$ Q(x) = \sum_{i = 1}^\infty \sum_{k = 1}^\infty (-1)^{k+1} k x^{ik}. $$

At this point, the expressions above can be ported to code via nested loops. Within the inner most $k$ loop, we iterate through all multiples of $i$ that are at most $N$ (recall that we need to compute $P$ and $Q$ only up to the $N$'th power of $x$), and doing so for all $i$ till $N$ yields an $O(N \log N)$ time complexity (the well-known ["harmonic series trick"](https://discuss.codechef.com/t/more-intuitive-explanation-for-the-harmonic-seriess-sum/67287)).

See the source code of `./solution.cpp` for implementation details. 
