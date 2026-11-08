import sys
import time
from abelian_groups_sum import count, powerful_walk, sieve_prefix


def timed(f, *a):
    t = time.perf_counter()
    r = f(*a)
    return r, time.perf_counter() - t


print("n  sieve  powerful  cubefull  value")
for e in range(6, 19):
    n = 10**e
    s = pw = None
    if e <= 8:
        s = timed(lambda m: int(sieve_prefix(m)[m]), n)[1]
    if e <= 12:
        pw = timed(powerful_walk, n)[1]
    v, t = timed(count, n)
    print(e, None if s is None else round(s, 3), None if pw is None else round(pw, 3), round(t, 3), v, flush=True)
