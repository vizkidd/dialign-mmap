#!/usr/bin/env python3
"""Differential fuzzer.  fuzz.py <reference-binary> <new-binary> [iterations] [seed]
Random multi-FASTA inputs x random option sets are run with both programs;
every output file (and stdout / exit code) must be byte-identical.
The reference is the original DIALIGN 2.2.1 plus the minimal bug fixes that
make its undefined behaviour defined (see CHANGES.md, tests/make_reference.sh)."""
import os, sys, random, subprocess, tempfile, shutil, filecmp

ref, new = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
iters = int(sys.argv[3]) if len(sys.argv) > 3 else 100
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
here = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, DIALIGN2_DIR=os.path.join(here, '..', 'dialign2_dir'))
R = random.Random(seed)
DNA, PROT = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'

def mutate(s, alpha, rate):
    o = []
    for c in s:
        x = R.random()
        if x < rate: o.append(R.choice(alpha))
        elif x < rate * 1.25: continue
        elif x < rate * 1.5: o += [c, R.choice(alpha)]
        else: o.append(c)
    return ''.join(o)

def gen():
    dna = R.random() < 0.6
    alpha = DNA if dna else PROT
    n = R.choice([2, 2, 3, 4, 5, 6, 8, 10, 14, 20, 30, 36, 40])
    L = R.choice([6, 10, 20, 40, 80, 120, 200])
    root = ''.join(R.choice(alpha) for _ in range(L))
    rate = R.choice([0.02, 0.1, 0.2, 0.35, 0.6, 1.0])
    seqs = []
    for i in range(n):
        s = mutate(root, alpha, rate) if rate < 1.0 else ''.join(R.choice(alpha) for _ in range(R.randint(max(3, L // 2), L)))
        if R.random() < 0.3: s = ''.join(R.choice(alpha) for _ in range(R.randint(0, 20))) + s
        if R.random() < 0.1: s = s.lower()
        if not s: s = R.choice(alpha) * 3
        seqs.append(s)
    if R.random() < 0.15: seqs[-1] = seqs[0]                      # duplicate
    wrap = R.choice([60, 70, 10, 100000])
    txt = ''
    for i, s in enumerate(seqs):
        txt += '>%s\n' % (R.choice(['s%d', 'seq_%d', 'a_long_sequence_name_%d x']) % (i + 1))
        txt += ''.join(s[k:k + wrap] + '\n' for k in range(0, len(s), wrap))
    if R.random() < 0.1: txt = txt.rstrip('\n')
    return dna, n, txt

def options(dna, n):
    o = []
    mode = R.choice(['-n', '-nt', '-ma', ''] if dna else ['', '', '-mat'])
    if not dna and mode == '-mat': o += ['-mat', '-mat_thr', str(R.randint(0, 3))]
    elif mode: o.append(mode)
    pool = [['-o'], ['-ds'], ['-cs'], ['-thr', str(R.randint(0, 4))], ['-lmax', str(R.choice([10, 20, 30, 40]))],
            ['-smin', str(R.randint(1, 10))], ['-it'], ['-istep', str(R.randint(1, 4))], ['-stars', str(R.randint(1, 6))],
            ['-lo'], ['-ff'], ['-fop'], ['-fsm'], ['-fsmv'], ['-afc'], ['-afc_v'], ['-msf'], ['-cw'], ['-fn', 'res'],
            ['-max_link'], ['-min_link'], ['-mask'], ['-nta'], ['-ow'], ['-iw'], ['-pst'], ['-csc'], ['-ref_seq'], ['-nas'],
            ['-lgs'], ['-lgs_t'], ['-lgsx'], ['-wtp'], ['-pand'], ['-pamnd'], ['-pao'], ['-d1w'], ['-cd_gobics']]
    for _ in range(R.choice([0, 0, 1, 1, 2, 3, 5])):
        o += R.choice(pool)
    return o

fails = 0; done = 0
for it in range(iters):
    dna, n, txt = gen()
    opts = options(dna, n)
    w = tempfile.mkdtemp()
    outs = {}
    for tag, b in (('r', ref), ('r2', ref), ('n', new), ('n2', new)):
        d = os.path.join(w, tag); os.mkdir(d)
        # the original miscounts the last sequence when the file has no final newline
        # (it reads one byte too many); the reference is fed the terminated file
        open(os.path.join(d, 'in.fa'), 'w').write(txt + ('\n' if tag == 'r' and not txt.endswith('\n') else ''))
        shutil.copy(b, os.path.join(d, 'dialign2-2'))
        try:
            p = subprocess.run(['./dialign2-2'] + opts + ['in.fa'], cwd=d, env=env, capture_output=True, timeout=60)
            outs[tag] = p.returncode
            open(os.path.join(d, 'stdout'), 'wb').write(p.stdout)
        except subprocess.TimeoutExpired:
            outs[tag] = 'timeout'
    ok = outs['r'] == outs['n']
    for f in sorted(os.listdir(os.path.join(w, 'r'))):
        if f in ('dialign2-2', 'in.fa'): continue
        if not os.path.exists(os.path.join(w, 'n', f)) or not filecmp.cmp(os.path.join(w, 'r', f), os.path.join(w, 'n', f), False):
            ok = False
    for f in os.listdir(os.path.join(w, 'n')):
        if f not in os.listdir(os.path.join(w, 'r')): ok = False
    done += 1
    if not ok:
        fails += 1
        keep = os.path.join(os.environ.get('FUZZ_KEEP', '/tmp'), 'fuzz_fail_%d_%d' % (seed, it)); shutil.copytree(w, keep)
        print('MISMATCH it=%d opts=%s rc=%s  -> %s' % (it, ' '.join(opts), outs, keep))
    shutil.rmtree(w, ignore_errors=True)
print('fuzz: %d runs, %d mismatches (seed %d)' % (done, fails, seed))
sys.exit(1 if fails else 0)
