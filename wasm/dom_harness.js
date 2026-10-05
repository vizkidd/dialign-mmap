'use strict';
const fs = require('fs');
const path = require('path');
const vm = require('vm');

function makeEl(overrides) {
  const el = {
    value: '', checked: false, files: [], textContent: '', innerHTML: '',
    disabled: false, style: {}, listeners: {},
    addEventListener(ev, fn) { el.listeners[ev] = fn; },
    appendChild(child) { el.children = (el.children || []); el.children.push(child); },
  };
  return Object.assign(el, overrides);
}

const INDEX_HTML = process.env.DIALIGN_INDEX_HTML || path.join(__dirname, 'index.html');
const SCRIPT = (() => {
  const html = fs.readFileSync(INDEX_HTML, 'utf8');
  const m = html.match(/<script>([\s\S]*?)<\/script>/);
  if (!m) throw new Error('no inline <script> found in ' + INDEX_HTML);
  return m[1];
})();

async function run(scenario, { useAnc = true, extra = '', mode = '' } = {}) {
  const elements = {
    fasta: makeEl({ files: [{ name: 'input.fa' }] }),
    useAnc: makeEl({ checked: useAnc }),
    ancFile: makeEl({ files: useAnc ? [{ name: 'input.anc' }] : [] }),
    ancFileWrap: makeEl(),
    debugMode: makeEl({ checked: mode === 'debug' }),
    opt_cs: makeEl(), opt_lo: makeEl(), opt_ma: makeEl(), opt_nas: makeEl(),
    extraOpts: makeEl({ value: extra }),
    runBtn: makeEl(),
    status: makeEl(),
    output: makeEl(),
    downloadWrap: makeEl(),
  };
  const seqtypeRadio = { value: '' };
  const logs = { console: [] };
  const fakeConsole = {
    log: (...a) => logs.console.push(['log', a.join(' ')]),
    warn: (...a) => logs.console.push(['warn', a.join(' ')]),
    error: (...a) => logs.console.push(['error', a.join(' ')]),
  };

  const wasmScenarios = {
    faithful: { hasMain: true, factory: async (cfg) => {
      cfg.FS = { writeFile() {}, readFile: () => new TextEncoder().encode('FAKE ALIGNMENT OUTPUT'),
                 readdir: () => ['.', '..', 'input.fa', 'input.ali'], stat: () => ({ size: 42 }) };
      (cfg.preRun || []).forEach((f) => f());
      cfg.print('DIALIGN 2.2.1');
      return cfg;
    } },
    noop: { hasMain: true, factory: async (cfg) => {
      cfg.FS = { writeFile() {}, readFile: () => { throw new Error('ENOENT'); },
                 readdir: () => ['.', '..', 'input.fa'], stat: () => ({ size: 10 }) };
      (cfg.preRun || []).forEach((f) => f());
      return cfg;
    } },
    'no-main': { hasMain: false, factory: async () => { throw new Error('should never be called'); } },
    // Emscripten's documented behaviour for a C exit(1): models the case
    // where it DOES throw. A real build with EXIT_RUNTIME=0 was observed
    // to instead resolve normally with no output -- that's the 'noop'
    // scenario above, which is the one actually exercised in practice.
    // Keeping both: the catch branch this models is cheap to keep and
    // might still fire on some Emscripten version/config.
    'throws-exit': { hasMain: true, factory: async () => { throw Object.assign(new Error('exit'), { name: 'ExitStatus', status: 1 }); } },
    'throws-other': { hasMain: true, factory: async () => { throw new TypeError('something unrelated broke'); } },
  };
  const scn = wasmScenarios[scenario];
  if (!scn) throw new Error('unknown scenario ' + scenario);

  const sandbox = {
    document: {
      getElementById: (id) => elements[id],
      querySelector: () => seqtypeRadio,
      createElement: (tag) => makeEl({ tag }),
    },
    console: fakeConsole,
    alert: (msg) => { throw new Error('alert() called: ' + msg); },
    fetch: async () => ({ arrayBuffer: async () => new Uint8Array([0, 0x61, 0x73, 0x6d]).buffer }),
    inspectWasm: () => ({ size: 123, functions: scn.hasMain ? 100 : 10,
                           exports: scn.hasMain ? [{ name: 'main' }] : [], hasMain: scn.hasMain }),
    DialignModule: scn.factory,
    FileReader: function () {
      this.readAsArrayBuffer = (file) => { this.result = new TextEncoder().encode('>' + file.name).buffer; this.onload && this.onload(); };
    },
    Blob: function (parts) { this.parts = parts; },
    URL: { createObjectURL: () => 'blob:fake' },
    TextDecoder, TextEncoder, setTimeout, clearTimeout,
  };
  sandbox.window = sandbox;
  sandbox.globalThis = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(SCRIPT, sandbox, { filename: 'index.html-inline-script.js' });

  if (!elements.runBtn.listeners.click) throw new Error('runBtn click handler was never registered');
  await elements.runBtn.listeners.click();
  return { elements, logs };
}

module.exports = { run };
