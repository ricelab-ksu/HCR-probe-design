/**
 * HCR Probe Designer — preload script.
 * Exposes a minimal, safe IPC surface to the Shiny renderer.
 */

const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('hcr', {
  /** Open the writable reference-data folder (species FASTA store). */
  openReferenceFolder: () => ipcRenderer.send('hcr:open-reference-folder'),
  /** Open the folder where probe design CSV outputs are written. */
  openOutputsFolder: () => ipcRenderer.send('hcr:open-outputs-folder'),
  /** Close the app. */
  quit: () => ipcRenderer.send('hcr:quit'),
  /** Version info for the splash screen. */
  version: () => ipcRenderer.invoke('hcr:version'),
  /** Subscribe to fatal bootstrap errors. */
  onError: (cb) => {
    ipcRenderer.on('hcr:error', (_e, msg) => cb(msg));
  },
});