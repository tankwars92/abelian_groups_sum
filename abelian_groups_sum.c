#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

static u64 n;
static u64 table_size;
static u64 *table;
static uint32_t *primes;
static size_t prime_count;
static u64 weights[65];
static u128 total;

static u64 isqrt64(u64 v)
{
    u64 r = (u64)sqrtl((long double)v);
    while (r * r > v)
        r--;
    while ((r + 1) * (r + 1) <= v)
        r++;
    return r;
}

static u64 icbrt64(u64 v)
{
    u64 r = (u64)cbrtl((long double)v);
    while ((u128)r * r * r > v)
        r--;
    while ((u128)(r + 1) * (r + 1) * (r + 1) <= v)
        r++;
    return r;
}

static u64 square_floor_sum(u64 y)
{
    u64 a_max = icbrt64(y);
    u64 s = 0;
    for (u64 a = 1; a <= a_max; a++)
        s += y / (a * a);
    u64 j_max = y / ((a_max + 1) * (a_max + 1));
    for (u64 j = 1; j <= j_max; j++)
        s += isqrt64(y / j) - a_max;
    return s;
}

static void walk(size_t start, u64 t, u64 w)
{
    u64 y = n / t;
    total += (u128)w * (y <= table_size ? table[y] : square_floor_sum(y));
    for (size_t k = start; k < prime_count; k++) {
        u64 p = primes[k];
        u128 q = (u128)t * p * p * p;
        if (q > n)
            break;
        int e = 3;
        while (q <= n) {
            walk(k + 1, (u64)q, w * weights[e]);
            q *= p;
            e++;
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
        return 1;
    n = strtoull(argv[1], 0, 10);
    if (n > 11000000000000000000ULL)
        return 3;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    u64 part[65] = {1};
    for (int k = 1; k <= 64; k++)
        for (int i = k; i <= 64; i++)
            part[i] += part[i - k];
    for (int e = 0; e <= 64; e++)
        weights[e] = part[e] - (e >= 1 ? part[e - 1] : 0) - (e >= 2 ? part[e - 2] : 0) + (e >= 3 ? part[e - 3] : 0);
    u64 root = icbrt64(n);
    char *sieve = calloc(root + 2, 1);
    primes = malloc((root / 2 + 16) * sizeof(uint32_t));
    for (u64 i = 2; i <= root; i++) {
        if (sieve[i])
            continue;
        primes[prime_count++] = (uint32_t)i;
        for (u64 j = i * i; j <= root; j += i)
            sieve[j] = 1;
    }
    table_size = 8 * root;
    if (table_size < 1000)
        table_size = 1000;
    if (table_size > (1u << 23))
        table_size = 1u << 23;
    if (table_size > n)
        table_size = n;
    table = calloc(table_size + 1, sizeof(u64));
    for (u64 a = 1; a * a <= table_size; a++)
        for (u64 m = a * a; m <= table_size; m += a * a)
            table[m]++;
    for (u64 i = 1; i <= table_size; i++)
        table[i] += table[i - 1];
    walk(0, 1, 1);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double dt = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    char buf[64];
    int len = 0;
    u128 v = total;
    if (!v)
        buf[len++] = '0';
    while (v) {
        buf[len++] = '0' + (int)(v % 10);
        v /= 10;
    }
    printf("n=%llu time=%.3fs value=", (unsigned long long)n, dt);
    while (len)
        putchar(buf[--len]);
    putchar('\n');
    return 0;
}
