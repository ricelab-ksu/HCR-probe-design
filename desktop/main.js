/**
 * HCR Probe Designer - Electron main process.
 * Spawns an embedded R runtime (bundled R.framework on macOS, R-x.y on
 * Windows, or system Rscript in --dev mode), runs runner.R which sources
 * the HCR app script, picks a free port, and prints HCR_PORT=<n> on stdout.
 * The port handshake lets the Electron main process open a BrowserWindow
 * pointed at the running Shiny server.
 */

const { app, BrowserWindow, dialog, Menu, Notification, session, shell, ipcMain } = require('electron');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

const DEV_MODE = process.argv.includes('--dev') || !app.isPackaged;
const APP_VERSION = app.isPackaged ? app.getVersion() : require('./package.json').version;

// ── Resource-path resolution ──────────────────────────────────────────────

function resolveAppDir() {
  if (DEV_MODE) return path.join(__dirname, 'resources', 'app');
  if (process.resourcesPath) {
    const p = path.join(process.resourcesPath, 'app');
    if (fs.existsSync(p)) return p;
  }
  return __dirname;
}

function resolveReferenceDir() {
  const candidates = [];
  if (DEV_MODE) {
    // Dev: resources/reference_files is only a stub now (genomes are fetched on
    // demand), so prefer the lab's live reference store next to desktop/.
    candidates.push(path.join(__dirname, '..', 'reference_files'));
    candidates.push(path.join(app.getPath('home'), 'Downloads', 'reference_files'));
  } else if (process.resourcesPath) {
    candidates.push(path.join(process.resourcesPath, 'reference_files'));
  }
  candidates.push(path.join(__dirname, 'resources', 'reference_files'));
  for (const c of candidates) if (fs.existsSync(c)) return c;
  return null;
}

function resolveRRuntime() {
  if (DEV_MODE) {
    const override = process.env.HCR_RSCRIPT;
    if (override && fs.existsSync(override)) return { kind: 'system', executable: override, args: [] };
    return { kind: 'system', executable: process.env.RSCRIPT || 'Rscript', args: [] };
  }
  // macOS: bundled R.framework
  const resourcesR = path.join(process.resourcesPath, 'r');
  const versionsDir = path.join(resourcesR, 'R.framework', 'Versions');
  if (fs.existsSync(versionsDir)) {
    const vers = fs.readdirSync(versionsDir).filter((v) => /^\d/.test(v)).sort().reverse();
    if (vers.length) {
      const res = path.join(versionsDir, vers[0], 'Resources');
      const execR = path.join(res, 'bin', 'exec', 'R');
      if (fs.existsSync(execR)) {
        return { kind: 'framework', executable: execR, args: [], env: {
          R_HOME: res,
          R_SHARE_DIR: path.join(res, 'share'),
          R_INCLUDE_DIR: path.join(res, 'include'),
          R_DOC_DIR: path.join(res, 'doc'),
          R_LIBS: path.join(res, 'library'),
          R_LIBS_USER: path.join(res, 'library'),
          R_ENVIRON_USER: '', R_PROFILE_USER: '',
        }};
      }
    }
  }
  // Windows: bundled R-x.y
  if (fs.existsSync(resourcesR)) {
    const sub = fs.readdirSync(resourcesR).find((d) => /^R-\d/.test(d));
    if (sub) {
      const rHome = path.join(resourcesR, sub);
      const script = path.join(rHome, 'bin', 'x64', 'Rscript.exe');
      return { kind: 'windows', executable: script, args: [], env: {
        R_HOME: rHome, R_LIBS: path.join(resourcesR, 'library'),
        R_LIBS_USER: path.join(resourcesR, 'library'),
        R_ENVIRON_USER: '', R_PROFILE_USER: '',
      }};
    }
  }
  return null;
}

// ── R child process + port handshake ─────────────────────────────────────

let rProc = null;
let rBufferedError = '';
let rPortResolve = null;

function spawnR(appDir, refDir) {
  const rt = resolveRRuntime();
  if (!rt) {
    throw new Error('No bundled R runtime found and not running in --dev mode. ' +
      'Run `npm run bundle:mac` (macOS) or `npm run bundle:win` (Windows) first.');
  }
  const runner = path.join(appDir, 'runner.R');
  // exec/R framework takes -f <file> --args ...; Rscript takes the script positionally.
  const scriptArgs = rt.kind === 'framework'
    ? ['-f', runner, '--args', appDir, refDir]
    : [runner, appDir, refDir];
  const args = ['--no-save', '--no-restore', '--no-environ', '--no-site-file',
    '--no-init-file', '--slave', ...scriptArgs];

  const env = {
    ...process.env,
    HCR_REFERENCE_DIR: refDir,
    R_DEFAULT_PACKAGES: 'methods,stats,graphics,grDevices,utils,datasets',
    R_PROFILE_USER: '', R_ENVIRON_USER: '', R_LIBS_USER: '',
  };
  if (rt.env) Object.assign(env, rt.env);

  console.log(`[HCR] spawning R: ${rt.executable} ${args.join(' ')}`);
  rProc = spawn(rt.executable, args, { env, stdio: ['ignore', 'pipe', 'pipe'] });

  let stdoutBuf = '';
  rProc.stdout.on('data', (chunk) => {
    stdoutBuf += chunk.toString();
    const m = stdoutBuf.match(/HCR_PORT=(\d+)/);
    if (m && rPortResolve) { rPortResolve(parseInt(m[1], 10)); rPortResolve = null; }
  });
  rProc.stderr.on('data', (chunk) => {
    const s = chunk.toString();
    rBufferedError += s;
    if (rBufferedError.length > 200000) rBufferedError = rBufferedError.slice(-200000);
    console.error('[HCR] R >>', s.trimEnd().split('\n').join('\n[HCR] R >> '));
  });
  rProc.on('error', (err) => {
    if (rPortResolve) { rPortResolve(null); rPortResolve = null; }
    console.error('[HCR] failed to launch R:', err.message);
  });
  rProc.on('exit', (code, signal) => {
    console.log(`[HCR] R exited code=${code} signal=${signal}`);
    if (rPortResolve) { rPortResolve(null); rPortResolve = null; }
  });
  return rProc;
}

// ── Window management ────────────────────────────────────────────────────

let mainWindow = null;

function createLoadingWindow() {
  mainWindow = new BrowserWindow({
    width: 1280, height: 830, minWidth: 1000, minHeight: 700,
    title: 'HCR Probe Designer',
    backgroundColor: '#1e1e2e',
    show: false,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true, nodeIntegration: false, sandbox: true,
    },
  });
  mainWindow.on('closed', () => { mainWindow = null; });
  mainWindow.webContents.setWindowOpenHandler(({ url }) => {
    if (url.startsWith('http')) shell.openExternal(url);
    return { action: 'deny' };
  });
  const bootstrapHtml = `<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<style>
  body{margin:0;background:#1e1e2e;color:#cdd6f4;font-family:-apple-system,Segoe UI,sans-serif;
    display:flex;align-items:center;justify-content:center;height:100vh;flex-direction:column;gap:18px;}
  .logo{font-size:30px;letter-spacing:1px;font-weight:600;color:#a6e3a1;}
  .sub{color:#9399b2;font-size:14px;}
  .bar{width:260px;height:4px;background:#313244;border-radius:2px;overflow:hidden;}
  .bar i{display:block;width:40%;height:100%;background:#89b4fa;border-radius:2px;
    animation:move 1.1s ease-in-out infinite;}
  @keyframes move{0%{transform:translateX(-100%)}100%{transform:translateX(650%)}}
  #err{display:none;color:#f38ba8;font-size:13px;max-width:620px;text-align:center;
    white-space:pre-wrap;font-family:ui-monospace,Menlo,monospace;}
</style></head><body>
  <div class="logo">HCR Probe Designer</div>
  <div class="sub">Starting the design engine (R)…</div>
  <div class="bar"><i></i></div>
  <div id="err"></div>
  <script>
    window.hcr.onError((msg) => {
      const el = document.getElementById('err');
      el.style.display = 'block';
      el.textContent = msg;
    });
  </script></body></html>`;
  mainWindow.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(bootstrapHtml));
  mainWindow.show();
  return mainWindow;
}

function onShinyReady(port) {
  if (!mainWindow || mainWindow.isDestroyed()) return;
  const url = `http://127.0.0.1:${port}`;
  let attempts = 0;

  // runner.R prints HCR_PORT just before runApp() binds the socket, so the
  // first load may hit ERR_CONNECTION_REFUSED. Retry until the server is up.
  const load = () => {
    attempts += 1;
    mainWindow.loadURL(url).catch((err) => {
      if (attempts < 120 && err && err.errno === 'ERR_CONNECTION_REFUSED') {
        setTimeout(load, 400);
      }
    });
  };
  mainWindow.webContents.on('did-fail-load', (_e, code, desc, _url, isMainFrame) => {
    if (isMainFrame && code === -102 && attempts < 120) {
      setTimeout(load, 400);
    }
  });
  load();

  mainWindow.once('page-title-updated', (event, title) => {
    event.preventDefault();
    mainWindow.setTitle(`HCR Probe Designer - ${title}`);
  });
}

function showRFatalError(tail) {
  if (!mainWindow || mainWindow.isDestroyed()) return;
  const summary = tail.length > 8000 ? tail.slice(-8000) : tail;
  mainWindow.webContents.send('hcr:error',
    `The R design engine could not start.\n\n${summary}\n\n` +
    'Make sure the app resources are intact, then relaunch.');
}

// ── App lifecycle ─────────────────────────────────────────────────────────

async function boot() {
  Menu.setApplicationMenu(null);
  createLoadingWindow();

  // Route downloads to the user's Downloads folder (needs app-ready session).
  session.defaultSession.on('will-download', (event, item) => {
    const downloads = app.getPath('downloads');
    item.setSavePath(path.join(downloads, item.getFilename()));
    item.once('done', (_e, state) => {
      if (state === 'completed' && Notification.isSupported()) {
        new Notification({ title: 'Download complete', body: item.getFilename() }).show();
      }
    });
  });

  const appDir = resolveAppDir();
  let refDir = resolveReferenceDir();

  // Packaged: stage a writable copy of reference data outside the .app bundle.
  // The installer no longer ships genome payloads (melanogaster and the other
  // species are downloaded from NCBI on first use), so the staged folder is
  // created even when the bundle carries nothing to copy. Existing users keep
  // their previously-staged reference set untouched.
  if (!DEV_MODE) {
    try {
      const dest = path.join(app.getPath('userData'), 'reference_files');
      const done = path.join(app.getPath('userData'), '.references_staged');
      if (refDir && fs.existsSync(refDir) && !fs.existsSync(done)) {
        fs.rmSync(dest, { recursive: true, force: true });
        fs.mkdirSync(path.dirname(dest), { recursive: true });
        fs.cpSync(refDir, dest, { recursive: true });
        fs.writeFileSync(done, String(Date.now()));
      }
      if (!fs.existsSync(dest)) fs.mkdirSync(dest, { recursive: true });
      refDir = dest;
    } catch (err) {
      console.error('[HCR] staging failed', err);
    }
  }

  console.log(`[HCR] appDir=${appDir}  refDir=${refDir}  dev=${DEV_MODE}`);

  if (!fs.existsSync(path.join(appDir, 'HCR_probe_design_v41.R'))) {
    dialog.showErrorBox('HCR Probe Designer', 'App payload missing: HCR_probe_design_v41.R');
    app.quit();
    return;
  }

  const port = await new Promise((resolve) => {
    rPortResolve = resolve;
    try { spawnR(appDir, refDir); } catch (err) { resolve(null); }
  });

  if (port) {
    onShinyReady(port);
  } else {
    const tail = rBufferedError || '(no output)';
    console.error('[HCR] R failed to report a port.\n', tail);
    showRFatalError(tail);
  }
}

app.whenReady().then(boot);

app.on('window-all-closed', () => app.quit());

app.on('before-quit', () => {
  if (rProc && !rProc.killed) {
    try { rProc.kill('SIGTERM'); } catch (e) { /* ignore */ }
    setTimeout(() => {
      if (rProc && !rProc.killed) try { rProc.kill('SIGKILL'); } catch (e2) {}
    }, 1500);
  }
});

// ── IPC ──────────────────────────────────────────────────────────────────

ipcMain.on('hcr:quit', () => app.quit());
ipcMain.handle('hcr:version', () => ({
  app: APP_VERSION, electron: process.versions.electron, node: process.versions.node,
}));
ipcMain.on('hcr:open-reference-folder', () => {
  let ref = null;
  if (!DEV_MODE) {
    const staged = path.join(app.getPath('userData'), 'reference_files');
    if (fs.existsSync(staged)) ref = staged;
  }
  if (!ref) ref = resolveReferenceDir();
  if (ref && fs.existsSync(ref)) shell.openPath(ref);
});
ipcMain.on('hcr:open-outputs-folder', () => {
  const out = path.join(app.getPath('documents'), 'HCR_Probe_Designs');
  fs.mkdirSync(out, { recursive: true });
  shell.openPath(out);
});