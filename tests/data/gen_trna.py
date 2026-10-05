import random

random.seed(1234)
bases = "ACGU".replace("U","T")  # use DNA alphabet (T instead of U) since DIALIGN mostly targets DNA/protein; keep consistent

# Canonical tRNA cloverleaf-ish template (~76nt), loosely based on conserved tRNA structural regions.
# We build a "base" tRNA sequence and mutate it per-sequence, keeping a strictly-conserved
# 7-nt anticodon loop (positions 33-39, 0-indexed 32-38) with anticodon "GAA" (Phe) in the middle.
ACCEPTOR_STEM = "GCGGATTT"
D_ARM = "AGCTCAGTTGGGAGAGCGTTAGACTGAAGATCTAAAGGTC"
ANTICODON_ARM_5 = "CCTTAG"
ANTICODON_LOOP = "TTGAAGAA"   # 8nt, includes conserved central "GAA" anticodon at pos 3-5
ANTICODON_ARM_3 = "CTAAGG"
VARIABLE_LOOP = "A"
TPSI_ARM = "TGGCGGATGTAGCCAAGTGGATCAAGGCAGT"
ACCEPTOR_STEM_3 = "GGATTCG"
CCA = "CCA"

base_seq = (ACCEPTOR_STEM + D_ARM + ANTICODON_ARM_5 + ANTICODON_LOOP + ANTICODON_ARM_3
            + VARIABLE_LOOP + TPSI_ARM + ACCEPTOR_STEM_3 + CCA)

anticodon_start = len(ACCEPTOR_STEM + D_ARM + ANTICODON_ARM_5)  # 0-indexed start of the 8nt anticodon loop
anticodon_len = len(ANTICODON_LOOP)

print("base_seq len:", len(base_seq))
print("anticodon loop 0-indexed:", anticodon_start, "-", anticodon_start+anticodon_len, "=>", base_seq[anticodon_start:anticodon_start+anticodon_len])

N = 600  # number of sequences

def mutate(seq, rate, protect_range=None):
    s = list(seq)
    for i in range(len(s)):
        if protect_range and protect_range[0] <= i < protect_range[1]:
            continue
        if random.random() < rate:
            s[i] = random.choice(bases)
    return "".join(s)

def indel(seq, protect_range, max_indels=2):
    # occasionally insert/delete a base OUTSIDE the protected anticodon loop region,
    # to create realistic unequal-length sequences (stress-tests dialign's gap handling)
    s = seq
    for _ in range(random.randint(0, max_indels)):
        pos = random.randrange(len(s))
        if protect_range[0] - 2 <= pos <= protect_range[1] + 2:
            continue
        if random.random() < 0.5 and len(s) > 60:
            s = s[:pos] + s[pos+1:]  # deletion
        else:
            s = s[:pos] + random.choice(bases) + s[pos:]  # insertion
    return s

records = []
for i in range(N):
    protect = (anticodon_start, anticodon_start + anticodon_len)
    seq = mutate(base_seq, rate=0.06, protect_range=protect)
    seq = indel(seq, protect_range=protect, max_indels=3)
    records.append((f"tRNA_{i+1:04d}", seq))

with open("/tmp/trna_600.fa", "w") as f:
    for name, seq in records:
        f.write(f">{name}\n")
        for j in range(0, len(seq), 70):
            f.write(seq[j:j+70] + "\n")

print("wrote", len(records), "sequences")
print("length range:", min(len(s) for _,s in records), "-", max(len(s) for _,s in records))

# Save the anticodon-loop position per sequence (needed to build the .anc anchor file later,
# since indels can shift its absolute position per-sequence).
with open("/tmp/trna_600_anticodon_positions.txt", "w") as f:
    for name, seq in records:
        pos = seq.find(ANTICODON_LOOP[1:7])  # search for the stable 6nt core (TGAAGA) to relocate after indels
        f.write(f"{name}\t{pos}\t{len(seq)}\n")
