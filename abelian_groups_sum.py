import numpy as np
from math import isqrt

LIMIT = 3 * 10**18


def partitions(m):
    p = [0] * (m + 1)
    p[0] = 1
    for k in range(1, m + 1):
        for i in range(k, m + 1):
            p[i] += p[i - k]
    return p


def icbrt(y):
    r = int(round(y ** (1.0 / 3.0)))
    while r * r * r > y:
        r -= 1
    while (r + 1) ** 3 <= y:
        r += 1
    return r


def primes_upto(m):
    if m < 2:
        return []
    s = np.ones(m + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(m) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].tolist()


def floor_sqrt(v):
    r = np.sqrt(v.astype(np.float64)).astype(np.int64)
    r -= (r * r > v).astype(np.int64)
    r += ((r + 1) * (r + 1) <= v).astype(np.int64)
    return r


def square_floor_sum(y):
    a_max = icbrt(y)
    a = np.arange(1, a_max + 1, dtype=np.int64)
    total = int((np.int64(y) // (a * a)).sum())
    j_max = y // ((a_max + 1) ** 2)
    if j_max:
        j = np.arange(1, j_max + 1, dtype=np.int64)
        total += int((floor_sqrt(np.int64(y) // j) - a_max).sum())
    return total


def square_floor_table(size):
    f = np.zeros(size + 1, dtype=np.int64)
    for a in range(1, isqrt(size) + 1):
        f[a * a::a * a] += 1
    return np.cumsum(f)


def cubefull_terms(n, weights):
    primes = primes_upto(icbrt(n))
    values = [1]
    mult = [1]

    def walk(start, t, w):
        for k in range(start, len(primes)):
            p = primes[k]
            q = t * p * p * p
            if q > n:
                break
            e = 3
            while q <= n:
                nw = w * weights[e]
                values.append(q)
                mult.append(nw)
                walk(k + 1, q, nw)
                q *= p
                e += 1

    walk(0, 1, 1)
    return np.array(values, dtype=np.int64), np.array(mult, dtype=np.int64)


def cubefull_weights(emax):
    p = partitions(emax)
    w = [0] * (emax + 1)
    for e in range(emax + 1):
        w[e] = p[e] - (p[e - 1] if e >= 1 else 0) - (p[e - 2] if e >= 2 else 0) + (p[e - 3] if e >= 3 else 0)
    return w


def count(n, table_size=None):
    if n < 1:
        return 0
    assert n <= LIMIT
    weights = cubefull_weights(64)
    t, w = cubefull_terms(n, weights)
    if table_size is None:
        table_size = min(max(icbrt(n) * 8, 1000), 2 * 10**7, n)
    table = square_floor_table(table_size)
    y = np.int64(n) // t
    near = y <= table_size
    total = int((w[near] * table[y[near]]).sum())
    for ti, wi in zip(t[~near].tolist(), w[~near].tolist()):
        total += wi * square_floor_sum(n // ti)
    return total


def powerful_walk(n):
    p = partitions(64)
    q = [p[e] - (p[e - 1] if e else 0) for e in range(65)]
    primes = primes_upto(isqrt(n))
    total = 0

    def walk(start, m, h):
        nonlocal total
        total += h * (n // m)
        for k in range(start, len(primes)):
            pr = primes[k]
            if m * pr * pr > n:
                break
            e = 2
            x = m * pr * pr
            while x <= n:
                walk(k + 1, x, h * q[e])
                x *= pr
                e += 1

    walk(0, 1, 1)
    return total


def sieve_prefix(n):
    p = partitions(64)
    f = np.ones(n + 1, dtype=np.int64)
    f[0] = 0
    for pr in primes_upto(isqrt(n)):
        step = pr * pr
        sub = f[step::step]
        ex = np.full(len(sub), 2, dtype=np.int64)
        e = 3
        pe = step * pr
        while pe <= n:
            s = pe // step
            ex[s - 1::s] += 1
            pe *= pr
            e += 1
        sub *= np.array(p, dtype=np.int64)[ex]
    return np.cumsum(f)
