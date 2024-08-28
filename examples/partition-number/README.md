# Partition Numbers

[Problem Link](https://judge.yosupo.jp/problem/partition_function).

Our task is to find the $i$'th [partition number](https://en.wikipedia.org/wiki/Partition_function_(number_theory)) for all $0 \le i \le N$, modulo $998244353$.

The key is to consider, for each $i = 1, 2, \dots$, how many times we pick $i$ for our partition. Picking, for example, $3$ zero times keeps our sum the same, picking $3$ once increases our sum by $3$, picking it twice increases our sum by $6$, picking it three times does so by $9$, and so on. Thus, the formal power series factor corresponding to $i = 3$ is

$$(1 + x^3 + x^6 + \dots) = \frac {1} {1 - x^3},$$

with the desired ordinary generating function thus being 

$$\frac{1}{\prod_{i = 1}^\infty (1 - x^i)}.$$

The formal power series corresponding to the denominator, up to $x^N$, can be computed in linear time via the [Pentagonal Number Theorem](https://en.wikipedia.org/wiki/Pentagonal\_number\_theorem). Inverting it in $O(N \log N)$ time then solves our problem with this time complexity.

See the source code of `./solution.cpp` for implementation details.
