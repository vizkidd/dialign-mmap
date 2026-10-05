/* wasm/wasm_inspect.js
 *
 * inspectWasm(Uint8Array) -> { size, functions, imports[], exports[{name,kind}],
 *                              hasMain, sections{} }
 *
 * Twin of wasm/check_wasm.py (kept deliberately in sync; both are tested
 * against the same binaries). The page runs this before instantiating the
 * module so that "the module has no main()" -- which otherwise looks like a
 * successful run that silently produced nothing -- becomes a loud, specific
 * error instead.
 */
(function (root) {
  function leb(d, i) {
    let r = 0, s = 0;
    for (;;) {
      const b = d[i++];
      r += (b & 0x7f) * Math.pow(2, s); s += 7;
      if (!(b & 0x80)) return [r, i];
    }
  }
  function str(d, i) {
    let n; [n, i] = leb(d, i);
    return [new TextDecoder().decode(d.subarray(i, i + n)), i + n];
  }
  function limits(d, i) {
    const flag = d[i++]; let _;
    [_, i] = leb(d, i);
    if (flag & 1) [_, i] = leb(d, i);
    return i;
  }
  const KINDS = { 0: 'func', 1: 'table', 2: 'memory', 3: 'global', 4: 'tag' };
  const SECTIONS = { 0: 'custom', 1: 'type', 2: 'import', 3: 'function', 4: 'table', 5: 'memory',
                     6: 'global', 7: 'export', 8: 'start', 9: 'element', 10: 'code', 11: 'data',
                     12: 'datacount', 13: 'tag' };

  function inspectWasm(d) {
    if (!(d[0] === 0 && d[1] === 0x61 && d[2] === 0x73 && d[3] === 0x6d))
      throw new Error('not a WebAssembly binary (bad magic) -- the server may be returning an HTML error page instead of the .wasm file');
    let i = 8; const secs = {};
    while (i < d.length) {
      const id = d[i++]; let sz; [sz, i] = leb(d, i);
      if (!(id in secs)) secs[id] = [i, sz];
      i += sz;
    }
    const info = { size: d.length, sections: {}, imports: [], exports: [], functions: 0, hasMain: false };
    for (const k in secs) info.sections[SECTIONS[k] || k] = secs[k][1];
    let cnt, _;
    if (2 in secs) {
      i = secs[2][0]; [cnt, i] = leb(d, i);
      for (let n = 0; n < cnt; n++) {
        let mod, nm; [mod, i] = str(d, i); [nm, i] = str(d, i);
        const kind = d[i++];
        if (kind === 0) [_, i] = leb(d, i);
        else if (kind === 1) { i += 1; i = limits(d, i); }
        else if (kind === 2) i = limits(d, i);
        else if (kind === 3) i += 2;
        else if (kind === 4) { i += 1; [_, i] = leb(d, i); }
        info.imports.push(mod + '.' + nm);
      }
    }
    if (3 in secs) { i = secs[3][0]; [info.functions] = leb(d, i); }
    if (7 in secs) {
      i = secs[7][0]; [cnt, i] = leb(d, i);
      for (let n = 0; n < cnt; n++) {
        let nm; [nm, i] = str(d, i);
        const kind = d[i++]; [_, i] = leb(d, i);
        info.exports.push({ name: nm, kind: KINDS[kind] || String(kind) });
      }
    }
    info.hasMain = info.exports.some(e => e.name === 'main' || e.name === '__main_argc_argv');
    return info;
  }

  root.inspectWasm = inspectWasm;
  if (typeof module !== 'undefined' && module.exports) module.exports = { inspectWasm };
})(typeof window !== 'undefined' ? window : globalThis);
