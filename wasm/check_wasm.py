#!/usr/bin/env python3
"""wasm/check_wasm.py <file.wasm> [--min-functions N] [--quiet]

Inspects a compiled WebAssembly module (standard library only -- no wabt /
binaryen needed) and exits non-zero if it can't possibly be a working
dialign2-2 build:

  * no exported `main` / `__main_argc_argv`  -> the C sources were not
    linked in (or were garbage-collected away); the page would load, run
    nothing, print nothing, and report "done" with no output.
  * suspiciously few defined functions        -> same failure, softer form.

Why this exists: a build once produced a 9,713-byte module holding only
malloc/free/stack helpers (10 functions, no main). Nothing failed at build
time and nothing failed at load time -- it just silently did nothing, which
is a miserable thing to debug from a browser. This turns that into a loud
build failure instead.

Exit codes: 0 ok, 1 not a working build, 2 unreadable / not wasm.
"""
import sys

def leb(d, i):
    r = s = 0
    while True:
        b = d[i]; i += 1
        r |= (b & 0x7f) << s; s += 7
        if not b & 0x80:
            return r, i

def name(d, i):
    n, i = leb(d, i)
    return d[i:i + n].decode("utf-8", "replace"), i + n

def limits(d, i):
    flag = d[i]; i += 1
    _, i = leb(d, i)
    if flag & 1:
        _, i = leb(d, i)
    return i

KINDS = {0: "func", 1: "table", 2: "memory", 3: "global", 4: "tag"}
SECTIONS = {0: "custom", 1: "type", 2: "import", 3: "function", 4: "table",
            5: "memory", 6: "global", 7: "export", 8: "start", 9: "element",
            10: "code", 11: "data", 12: "datacount", 13: "tag"}

def inspect(data):
    if data[:4] != b"\0asm":
        raise ValueError("not a WebAssembly binary (bad magic)")
    i = 8
    secs = {}
    while i < len(data):
        sid = data[i]; i += 1
        sz, i = leb(data, i)
        secs.setdefault(sid, (i, sz))
        i += sz
    info = {"size": len(data), "sections": {SECTIONS.get(k, k): v[1] for k, v in secs.items()},
            "imports": [], "exports": [], "functions": 0}
    if 2 in secs:
        i, _ = secs[2]
        cnt, i = leb(data, i)
        for _ in range(cnt):
            mod, i = name(data, i); nm, i = name(data, i)
            kind = data[i]; i += 1
            if kind == 0:   _, i = leb(data, i)
            elif kind == 1: i += 1; i = limits(data, i)
            elif kind == 2: i = limits(data, i)
            elif kind == 3: i += 2
            elif kind == 4: i += 1; _, i = leb(data, i)
            info["imports"].append(f"{mod}.{nm}")
    if 3 in secs:
        i, _ = secs[3]
        info["functions"], _ = leb(data, i)
    if 7 in secs:
        i, _ = secs[7]
        cnt, i = leb(data, i)
        for _ in range(cnt):
            nm, i = name(data, i); kind = data[i]; i += 1
            _, i = leb(data, i)
            info["exports"].append((nm, KINDS.get(kind, kind)))
    return info

def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    quiet = "--quiet" in argv
    min_fn = 50
    if "--min-functions" in argv:
        min_fn = int(argv[argv.index("--min-functions") + 1])
        args = [a for a in args if a != str(min_fn)]
    if not args:
        print(__doc__); return 2
    try:
        info = inspect(open(args[0], "rb").read())
    except Exception as e:
        print(f"check_wasm: cannot read {args[0]}: {e}", file=sys.stderr)
        return 2
    exports = [n for n, _ in info["exports"]]
    has_main = "main" in exports or "__main_argc_argv" in exports
    problems = []
    if not has_main:
        problems.append("no exported main()/__main_argc_argv -- the C sources are not in this module")
    if info["functions"] < min_fn:
        problems.append(f"only {info['functions']} functions defined (expected at least {min_fn} "
                        f"for a real dialign2-2 build)")
    if not quiet or problems:
        print(f"check_wasm: {args[0]}")
        print(f"  size        : {info['size']} bytes")
        print(f"  functions   : {info['functions']} defined, {len(info['imports'])} imported")
        print(f"  exports     : {', '.join(exports) if exports else '(none)'}")
        print(f"  has main    : {'yes' if has_main else 'NO'}")
    if problems:
        for p in problems:
            print(f"check_wasm: FAIL: {p}")
        return 1
    if not quiet:
        print("check_wasm: OK")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
