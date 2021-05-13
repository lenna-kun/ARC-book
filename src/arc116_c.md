# ARC116 C - Multiple Sequences

[https://atcoder.jp/contests/arc116/tasks/arc116_c](https://atcoder.jp/contests/arc116/tasks/arc116_c)

## 問題文

整数\\(N \\)， \\(M \\)が与えられる．長さ \\(N \\)の数列 \\(A \\)であって，以下の条件を満たすものの数を答えよ．

- \\(1 \leq A_i \leq M \\)（\\(i=1,2,...,N \\)）
- \\(A_{i+1} \\)は\\(A_i \\)の倍数（\\(i=1,2,...,N \\)）

ただし，答えは非常に大きくなる場合があるので，\\(998244353\\)で割ったあまりを答えよ．

## 制約

- 入力は全て整数
- \\(1 \leq N \leq 2 \times {10}^{5} \\)
- \\(1 \leq M \leq 2 \times {10}^{5} \\)

## 解法

数列\\(A \\)の各項で，前の項に何を掛けるかに注目する．わかりやすくするため，数列\\(A \\)の先頭に\\(A_{0} = 1 \\)を加える．

この問題が，約数と密接に関係があることは明白．約数を考える時は，素因数ごとに独立して考えるとわかりやすい．そこで，次のような数列\\(A \\)を考えてみる．

| \\( A_{0}\\) | \\( A_{1}\\) | \\( A_{2}\\) | \\( A_{3}\\) | \\( A_{4}\\) | \\( A_{5}\\) |
|:-:|:-:|:-:|:-:|:-:|:-:|
| \\(1 \\) | \\(1 \\) | \\(2 \\) | \\(4 \\) | \\(4 \\) | \\(16 \\)|

前述の通り，各項が前の項の何倍になっているかをみてみる．

| \\( A_{1} \div A_{0}\\) | \\( A_{2} \div A_{1}\\) | \\( A_{3} \div A_{2}\\) | \\( A_{4} \div A_{3}\\) | \\( A_{5} \div A_{4}\\) |
|:-:|:-:|:-:|:-:|:-:|
| \\( 2^0\\) | \\( 2^1\\) | \\( 2^1\\) | \\( 2^0\\) | \\( 2^2\\) |

\\(A_{N} = 16 \\)である時，\\( 16 = 2^4\\)であるため，最終的には\\(2 \\)は\\( 4\\)回掛けられる．どこで何回\\(2 \\)を掛けるかということについて，何通りあるかを考えればいいことがわかる．これは，

\\[
    \begin{align}
        {}_{n} \mathrm{H} _{4}
    \end{align}
\\]

で求めることができる．以上で，素因数が1種類のときの求め方が分かった．素因数の種類が複数のときも上と同様に，以下の式が成り立つ．

\\[
    \begin{align}
        \prod_{i=1}^{N} (A_{i} \div A_{i-1}) = A_{N}
    \end{align}
\\]

よって，素因数ごとに独立して考え，それらを最後に掛け合わせれることで答えを得ることができる．

最後に計算量について考える．まず，\\(M \\)以下の全ての数においてそれぞれ素因数分解をする必要がある．前準備としてエラトステネスの篩\\(O(M \log \log M) \\)をする際，最小の素因数を記憶しておくことで，素因数分解を\\(O(\log M) \\)を行うことができる．それぞれにおける組合せの計算は，前準備（\\(O(N) \\)）をしておくことで\\(O(1) \\)で行うことができる．よって，上記のアルゴリズムにより，\\(O(N + M \log M) \\)でこの問題を解くことができる．

## 実装

```py
import math
MOD = 998244353

# 組合せ計算用前準備
def init(n): # O(n)
  fact = [1]*(n+1)
  rfact = [1]*(n+1)
  r = 1
  for i in range(1, n+1):
    fact[i] = r = r * i % MOD
  rfact[n] = r = pow(fact[n], MOD-2, MOD)
  for i in range(n, 0, -1):
    rfact[i-1] = r = r * i % MOD
  return fact, rfact

# nCk (mod MOD)
def comb(fact, rfact, n, k):
  return fact[n] * rfact[k] * rfact[n-k] % MOD

# nHk (mod MOD)
def homb(fact, rfact, n, k):
  return comb(fact, rfact, n+k-1, k)

# 素因数分解の前準備
def eratosthenes(n): # O(nlog(log(n)))
  limit = math.sqrt(n)
  data = list(range(2, n+1))
  ret = list(range(n+1))
  while True:
    p = data[0]
    if limit < p:
      return ret
    tmp = []
    for e in data:
      if e % p == 0:
        ret[e] = p
      else:
        tmp.append(e)
    data = tmp

# 素因数分解
def prime_factorization(data, n): # O(nlog(n))
  ret = {}
  while True:
    if data[n] not in ret:
      ret[data[n]] = 0
    ret[data[n]] += 1
    if data[n] == n:
      return ret
    n //= data[n]

n, m = map(int, input().split())

fact, rfact = init(10**6)
if m > 1:
  data = eratosthenes(m)
ans = 1 # 全て1のとき
for i in range(2, m+1):
  p = prime_factorization(data, i)
  tmp = 1
  for v in p.values():
    tmp = (tmp * homb(fact, rfact, n, v))%MOD
  ans = (ans + tmp)%MOD

print(ans)
```