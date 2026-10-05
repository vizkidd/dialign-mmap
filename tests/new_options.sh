#!/bin/bash
# new_options.sh <new-binary>  -- checks of the options the original does not have
here=$(dirname "$(readlink -f "$0")"); N=$(readlink -f "$1"); D=$here/data; fail=0
export DIALIGN2_DIR=$here/../dialign2_dir
w=$(mktemp -d); cp $D/dna20.fa $D/prot15.fa $w/; cd $w
strip(){ grep -v "program call\|program parameters" "$1"; }
ok(){ if [ "$1" = 0 ]; then echo "ok   $2"; else echo "FAIL $2"; fail=$((fail+1)); fi; }
# reference output
$N -n dna20.fa >/dev/null 2>&1; strip dna20.ali > base_dna.txt
$N prot15.fa >/dev/null 2>&1; strip prot15.ali > base_prot.txt
# -nommap / -threads / -tmpdir give the same output
mkdir td
for o in "-nommap" "-threads 1" "-threads 5" "-tmpdir td" "-nommap -threads 3"; do
  $N -n $o dna20.fa >/dev/null 2>&1; strip dna20.ali | cmp -s - base_dna.txt; ok $? "dna20 $o"
  $N $o prot15.fa >/dev/null 2>&1;    strip prot15.ali | cmp -s - base_prot.txt; ok $? "prot15 $o"
done
[ -z "$(ls -A td)" ]; ok $? "-tmpdir directory left clean"
# unusable tmpdir: warning + fallback, same output
$N -n -tmpdir /nonexistent/dir dna20.fa >/dev/null 2>err; strip dna20.ali | cmp -s - base_dna.txt; ok $? "unusable -tmpdir falls back"
grep -q "falling back" err; ok $? "fallback warning printed"
# -notree: identical alignment, only the tree text is missing
$N -n -notree dna20.fa >/dev/null 2>&1; strip dna20.ali > nt.txt
[ "$(diff base_dna.txt nt.txt | grep -c '^[<>]')" = 2 ] && diff base_dna.txt nt.txt | grep -q '^< *(('; ok $? "-notree removes only the tree"
grep -q "sample_1_lon" nt.txt; ok $? "-notree keeps the alignment"
# bad arguments are rejected with rc != 0
$N -n -threads 0 dna20.fa >/dev/null 2>&1; [ $? -ne 0 ]; ok $? "-threads 0 rejected"
$N -n -tmpdir dna20.fa >/dev/null 2>&1; [ $? -ne 0 ]; ok $? "-tmpdir without argument rejected"
# more than 10000 sequences are accepted by the reader (no MAX_SEQNUM)
python3 - <<'PY'
import random
r=random.Random(3)
with open('many.fa','w') as f:
    for i in range(10500): f.write('>q%d\n%s\n'%(i,''.join(r.choice('ACGT') for _ in range(9))))
PY
$N -n -thr 9 -lmax 5 -stdo -notree -ref_seq many.fa > many.out 2>&1; rc=$?
[ $rc = 0 ] && grep -q "q10499" many.out; ok $? "10500 sequences (-ref_seq): completes, last sequence present"
rm -rf $w; [ $fail = 0 ] && echo "NEW OPTIONS OK" || echo "$fail FAILED"; exit $fail
