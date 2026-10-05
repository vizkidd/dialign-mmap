#!/usr/bin/env node
/* wasm/smoke_test.js [--dist DIR] [--module FILE.js]
 *
 * Loads the built WASM module under Node (the build uses
 * -sENVIRONMENT=web,node precisely so this is possible), runs dialign2-2 on
 * small inputs through the same mechanism the web page uses (config
 * `arguments` + `preRun` writing into the virtual filesystem), and requires
 * the alignment to be byte-identical to the native build's output stored in
 * wasm/golden/ (volatile "program call" lines excluded -- the same
 * convention tests/compare.sh uses).
 *
 * This is what stops a silently-broken build (e.g. one containing no main(),
 * which loads fine, runs nothing, prints nothing, and reports success) from
 * ever being deployed.
 *
 * Exit 0 = all cases identical to native; 1 = any failure (with a log tail).
 */
'use strict';
const fs = require('fs');
const path = require('path');

function arg(name, def) {
  const i = process.argv.indexOf(name);
  return i >= 0 ? process.argv[i + 1] : def;
}
const dist = path.resolve(arg('--dist', path.join(__dirname, 'dist')));
const modFile = path.resolve(arg('--module', path.join(dist, 'dialign2-2.js')));
const dataDir = path.join(__dirname, '..', 'tests', 'data');
const goldenDir = path.join(__dirname, 'golden');

// golden files were produced natively with: dialign2-2 [args] -nommap input.fa
const CASES = [
  { name: 'prot6', args: [] },
  { name: 'dna8',  args: ['-n'] },
];
const strip = (s) => s.split('\n').filter((l) => !/program call|program parameters/.test(l)).join('\n');

async function runCase(factory, c) {
  const log = [];
  const fasta = fs.readFileSync(path.join(dataDir, c.name + '.fa'));
  const cfg = {
    print: (t) => log.push(t),
    printErr: (t) => log.push('[stderr] ' + t),
    arguments: [...c.args, '-nommap', 'input.fa'],
    locateFile: (p) => path.join(dist, p),
    preRun: [function () {
      // Same access pattern as wasm/index.html: the config object becomes the
      // Module object, so FS is reachable through it by the time preRun runs.
      const fsx = cfg.FS || (typeof FS !== 'undefined' ? FS : null);
      if (!fsx) throw new Error('Emscripten FS is not exposed (need -sFORCE_FILESYSTEM=1 -sEXPORTED_RUNTIME_METHODS=["FS"])');
      fsx.writeFile('/input.fa', fasta);
    }],
  };
  let M;
  try {
    M = await factory(cfg);
  } catch (e) {
    return { ok: false, why: 'module threw: ' + ((e && e.stack) || e), log };
  }
  let ali;
  try {
    ali = Buffer.from((cfg.FS || (M && M.FS)).readFile('/input.ali')).toString('utf8');
  } catch (e) {
    return { ok: false, why: 'no /input.ali was produced (main() never ran, or failed before writing output)', log };
  }
  const want = fs.readFileSync(path.join(goldenDir, c.name + '.ali'), 'utf8');
  const got = strip(ali);
  if (got !== want) {
    const g = got.split('\n'), w = want.split('\n');
    let k = 0;
    while (k < g.length && k < w.length && g[k] === w[k]) k++;
    return { ok: false, log,
      why: `alignment differs from native golden output at line ${k + 1}:\n    wasm  : ${JSON.stringify(g[k])}\n    native: ${JSON.stringify(w[k])}` };
  }
  return { ok: true, log };
}

(async () => {
  console.log(`smoke_test: module ${modFile}`);
  let factory;
  try {
    factory = require(modFile);
  } catch (e) {
    console.error('smoke_test: FAIL: cannot load module: ' + e.message);
    process.exit(1);
  }
  let failed = 0;
  for (const c of CASES) {
    const r = await runCase(factory, c);
    if (r.ok) {
      console.log(`  ok    ${c.name} ${c.args.join(' ')}  (identical to native)`);
    } else {
      failed++;
      console.log(`  FAIL  ${c.name} ${c.args.join(' ')}: ${r.why}`);
      const tail = r.log.slice(-25);
      console.log('  --- last ' + tail.length + ' lines of module output ---');
      tail.forEach((l) => console.log('  | ' + l));
    }
  }
  if (failed) {
    console.error(`smoke_test: FAIL (${failed} of ${CASES.length})`);
    process.exit(1);
  }
  console.log('smoke_test: OK -- WASM output identical to native for all cases');
  // Emscripten's node runtime can leave handles open / set exitCode; be explicit.
  process.exit(0);
})();
