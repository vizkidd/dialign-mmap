#!/usr/bin/env python3
"""Regenerates the synthetic DNA/protein test sets in tests/data (seeded)."""
import random, os
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DNA, PROT = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'

def mut(s, alpha, r, rate):
    o = []
    for c in s:
        x = r.random()
        if x < rate: o.append(r.choice(alpha))
        elif x < rate * 1.3: continue
        elif x < rate * 1.6: o.append(c); o.append(r.choice(alpha))
        else: o.append(c)
    return ''.join(o)

def gen(name, n, L, alpha, rate, seed, wrap=60, names=None):
    r = random.Random(seed); root = ''.join(r.choice(alpha) for _ in range(L))
    with open(os.path.join(here, name + '.fa'), 'w') as f:
        for i in range(n):
            s = mut(root, alpha, r, rate)
            s = ''.join(r.choice(alpha) for _ in range(r.randint(0, 15))) + s + \
                ''.join(r.choice(alpha) for _ in range(r.randint(0, 15)))
            f.write('>' + (names(i) if names else 'seq%d' % (i + 1)) + '\n')
            for k in range(0, len(s), wrap): f.write(s[k:k + wrap] + '\n')

gen('dna8', 8, 150, DNA, 0.15, 1)
gen('dna20', 20, 120, DNA, 0.2, 2, names=lambda i: 'sample_%d_long_name_x desc' % (i + 1))
gen('dna36', 36, 80, DNA, 0.15, 3)
gen('prot6', 6, 120, PROT, 0.25, 4)
gen('prot15', 15, 90, PROT, 0.3, 5)
gen('prot40', 40, 60, PROT, 0.3, 6)
gen('dna2', 2, 300, DNA, 0.2, 7)
gen('dna_orf5', 5, 210, DNA, 0.1, 8)
gen('dna_nowrap', 5, 400, DNA, 0.2, 9, wrap=100000)
