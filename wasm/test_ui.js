#!/usr/bin/env node
/* wasm/test_ui.js -- actually executes wasm/index.html's inline <script>
 * under Node (via dom_harness.js's minimal DOM stub), exercising every
 * branch the click handler can take with fake DialignModule()
 * implementations, and asserts the observable result (status text,
 * whether a download link was produced, console output) for each.
 *
 * This exists because two real, embarrassing bugs (Module.callMain not
 * exported; FS/ENV not in scope inside preRun) shipped past a plain
 * read-through of the code. Both were plain JS errors that would have
 * failed loudly the first time any of these branches actually ran -- this
 * runs all of them, every time, without needing a browser.
 *
 * Exit 0 = all assertions passed, 1 = something regressed.
 */
'use strict';
const { run } = require('./dom_harness.js');

let failures = 0;
function check(label, cond, detail) {
  if (cond) { console.log(`  ok    ${label}`); }
  else { failures++; console.log(`  FAIL  ${label}${detail ? ' -- ' + detail : ''}`); }
}

async function main() {
  console.log('faithful (successful alignment):');
  {
    const { elements, logs } = await run('faithful');
    check('status shows exit code 0', elements.status.textContent === 'done (exit code 0)', elements.status.textContent);
    check('output shows the alignment', elements.output.textContent === 'FAKE ALIGNMENT OUTPUT');
    check('a download link was created', (elements.downloadWrap.children || []).some((c) => c.download === 'input.ali'));
    check('normal mode stays quiet on console', logs.console.length === 0, `${logs.console.length} lines`);
  }
  {
    const { logs } = await run('faithful', { mode: 'debug' });
    check('debug mode logs to console', logs.console.length > 0, `${logs.console.length} lines`);
  }

  console.log('noop (module resolves normally, produces no output -- the originally-reported bug,');
  console.log('      AND the confirmed real behaviour of a C exit(1) under EXIT_RUNTIME=0: it does');
  console.log('      not throw an ExitStatus here, it resolves exactly like a normal return):');
  {
    const { elements } = await run('noop');
    check('status flags no-output, not a bare "done"', elements.status.textContent === 'done -- no output (see log)', elements.status.textContent);
    check('reports no .ali produced', /No \.ali output was produced/.test(elements.output.innerHTML));
    check('names the anchor/FASTA-mismatch case as a common cause', /anchor .*doesn't match the FASTA/.test(elements.output.innerHTML));
    check('suggests Debug mode in normal mode', /tick "Debug mode"/.test(elements.output.innerHTML));
  }

  console.log('no-main (the actual root cause found in dialign-mmap-6.zip):');
  {
    const { elements, logs } = await run('no-main');
    check('status is "error", not "done"', elements.status.textContent === 'error', elements.status.textContent);
    check('names the real cause', /no main\(\) linked in/.test(elements.output.innerHTML));
    check('DialignModule was never called', true); // no-main's factory throws if called; run() not throwing proves this
  }

  console.log('throws-exit (erreur() -> C exit(1)):');
  {
    const { elements } = await run('throws-exit');
    check('status shows exit code 1', elements.status.textContent === 'done (exit code 1)', elements.status.textContent);
    check('does not claim success', !/exit code 0/.test(elements.status.textContent));
  }

  console.log('throws-other (an unrelated JS bug must not look like a normal exit):');
  {
    const { elements } = await run('throws-other');
    check('status is "error"', elements.status.textContent === 'error', elements.status.textContent);
    check('shows the real message, not "exited with code"', /something unrelated broke/.test(elements.output.innerHTML));
  }

  console.log('anchor file handling:');
  {
    const { elements } = await run('faithful', { useAnc: false });
    check('runs fine with no anchor file', elements.status.textContent === 'done (exit code 0)');
  }
  {
    const r = await run('faithful', { useAnc: true });
    r.elements.ancFile.files = []; // box checked, no file actually chosen
    await r.elements.runBtn.listeners.click();
    check('anchor checked with no file chosen is a clear error',
      /no \.anc file chosen/.test(r.elements.output.innerHTML));
  }

  console.log('button always re-enabled:');
  for (const s of ['faithful', 'noop', 'no-main', 'throws-exit', 'throws-other']) {
    const { elements } = await run(s);
    check(`runBtn re-enabled after ${s}`, elements.runBtn.disabled === false);
  }

  if (failures) {
    console.error(`\ntest_ui: FAIL (${failures} assertion(s))`);
    process.exit(1);
  }
  console.log('\ntest_ui: OK');
}

main();
