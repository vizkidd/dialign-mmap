import random, sys, os
N = int(sys.argv[1])
outdir = sys.argv[2]
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 4242
random.seed(seed)
bases = "ACGT"

ACCEPTOR_STEM = "GCGGATTT"
D_ARM = "AGCTCAGTTGGGAGAGCGTTAGACTGAAGATCTAAAGGTC"
ANTICODON_ARM_5 = "CCTTAG"
ANTICODON_LOOP = "TTGAAGAA"
ANTICODON_ARM_3 = "CTAAGG"
VARIABLE_LOOP = "A"
TPSI_ARM = "TGGCGGATGTAGCCAAGTGGATCAAGGCAGT"
ACCEPTOR_STEM_3 = "GGATTCG"
CCA = "CCA"
base_seq = (ACCEPTOR_STEM + D_ARM + ANTICODON_ARM_5 + ANTICODON_LOOP + ANTICODON_ARM_3
            + VARIABLE_LOOP + TPSI_ARM + ACCEPTOR_STEM_3 + CCA)
anticodon_start = len(ACCEPTOR_STEM + D_ARM + ANTICODON_ARM_5)
anticodon_core = ANTICODON_LOOP[1:7]

def mutate(seq, rate, protect):
    s = list(seq)
    for i in range(len(s)):
        if protect[0] <= i < protect[1]:
            continue
        if random.random() < rate:
            s[i] = random.choice(bases)
    return "".join(s)

def maybe_intron(seq, protect, p=0.3, lo=20, hi=90):
    if random.random() > p:
        return seq, protect
    pos = random.choice([random.randint(0, max(0, protect[0]-5)),
                          random.randint(protect[1]+5, len(seq))])
    ins = "".join(random.choice(bases) for _ in range(random.randint(lo, hi)))
    new_seq = seq[:pos] + ins + seq[pos:]
    if pos <= protect[0]:
        protect = (protect[0] + len(ins), protect[1] + len(ins))
    return new_seq, protect

def indel(seq, protect, max_indels=4):
    for _ in range(random.randint(0, max_indels)):
        pos = random.randrange(len(seq))
        if protect[0] - 3 <= pos <= protect[1] + 3:
            continue
        if random.random() < 0.5 and len(seq) > 60:
            seq = seq[:pos] + seq[pos+1:]
            if pos < protect[0]:
                protect = (protect[0]-1, protect[1]-1)
        else:
            seq = seq[:pos] + random.choice(bases) + seq[pos:]
            if pos < protect[0]:
                protect = (protect[0]+1, protect[1]+1)
    return seq, protect

records = []
positions = {}
for i in range(N):
    protect = (anticodon_start, anticodon_start + len(ANTICODON_LOOP))
    seq = mutate(base_seq, rate=0.08, protect=protect)
    seq, protect = maybe_intron(seq, protect)
    seq, protect = indel(seq, protect)
    leader = "".join(random.choice(bases) for _ in range(random.randint(30, 220)))
    trailer = "".join(random.choice(bases) for _ in range(random.randint(30, 220)))
    full = leader + seq + trailer
    core_pos = len(leader) + protect[0] + 1
    found = full.find(anticodon_core, max(0, core_pos-5))
    assert found != -1
    name = f"tRNA_{i+1:04d}"
    records.append((name, full))
    positions[name] = (found, len(anticodon_core))

os.makedirs(outdir, exist_ok=True)
base = f"trna{N}"
with open(os.path.join(outdir, base+".fa"), "w") as f:
    for name, seq in records:
        f.write(f">{name}\n")
        for k in range(0, len(seq), 70):
            f.write(seq[k:k+70] + "\n")

names = [r[0] for r in records]
n_anc = 0
with open(os.path.join(outdir, base+".anc"), "w") as f:
    for i in range(N):
        pi, L = positions[names[i]]
        for j in range(i+1, N):
            pj, _ = positions[names[j]]
            f.write(f"{i+1} {j+1} {pi+1} {pj+1} {L} 10.0\n")
            n_anc += 1

lens = [len(s) for _, s in records]
print(f"N={N} anchors={n_anc} len_range={min(lens)}-{max(lens)}")
