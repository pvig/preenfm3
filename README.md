# preenfm3

Here are the sources for the preenfm3 bootloader and firmware.

Binaries can be found in the [Release section](https://github.com/gresade/preenfm3/releases).

For the documentation, go to [the Wiki](https://github.com/Ixox/preenfm3/wiki).

## Branch Notes (`fixes` vs fork base)

Baseline used for comparison: merge-base with `origin/master` (`65cb963`).

Functional additions and fixes in this branch include:

- Per-operator start phase offset support in the synth engine, with UI visualization of the phase marker in the oscilloscope view and encoder/editor wiring.
- Per-operator **Warp** parameter: asymmetrically stretches the first and second halves of the wavetable cycle so the waveform leans towards attack or decay. Range ±4 (values beyond ±1 invert one half of the waveform). Warp is applied in all oscillator render paths (`getNextSample`, `getNextBlock`, `getNextBlockHQ`) including the decimated sub-paths. A deadband near zero preserves the fast no-warp render path under small modulation.
- Both `Phase` and `Warp` for operators 1–6 are exposed as modulation-matrix destinations (`o1Ph`–`o6Ph`, `o1Wr`–`o6Wr`), wired through the matrix update loop in `Voice::prepareMatrixForNewBlock`. Warp updates are applied every block; Phase is applied at note-on.
- `Warp` is mapped to MIDI CC per operator and accessible via the encoder row in the editor UI.
- Patch format bumped to **v1.5**: Warp defaults to 0 on load of older patches; out-of-range Warp values are clamped on load.
- DX7 SysEx import hardening and mapping fixes (bank validation, vibrato/AMS handling, fixed-frequency preservation).
- FM decimation control added to engine parameters with preset/file persistence support.
- FM decimation modes are capped at 19-bit precision (plus `Full`), with legacy preset values above 19-bit clamped on load for compatibility.
- FM decimation path is simplified for CPU efficiency: output sample quantization is retained, phase-accumulator quantization removed, and temporal decimation (`D=2` sample-and-hold) is applied in safe oscillator render paths.
- Feedback oscillator rendering remains full-rate for stability while still using decimation output quantization.
- High-quality wavetable interpolation is bypassed when decimation is enabled.
- LFO Sync selector expanded with one-shot to eight-shot modes for both internal and external sync clocks: `1Si..8Si` and `1Se..8Se`, in addition to `Int` and `Ext`.
- LFO one-shot mode and KSyn are now independent so KSyn can be combined with shot modes.
- LFO shape list expanded in the editor/oscilloscope: `SawD`, `DExp`, `DLog`, `RExp`, `RLog`, `AD`, `AHD`, `SDec`, `Plng`, `Plg2`, `SnSq`, `Sn0`, `Sn+`, and `Usr1..Usr6`.
- LFO Phase encoder now doubles as a startup delay: negative values set a delay of 0–4000 ms (displayed as integer ms); positive values set a phase offset of 0–1. Encoder stepping uses adaptive resolution — 1 ms steps below 50 ms, up to 50 ms steps near 4 s — with a float32 forward-progress guard that prevents the encoder from getting stuck at quantisation boundaries (notably the 100 ms transition between step sizes).
- Build and release workflow improvements for VS Code/CLI headless builds (`scripts/build_cli.sh`), including automatic release artifact refresh and checksum regeneration.
- New subrelease packaging updated to firmware `v1.06f` with `bl1.09` bundle naming.
- Oscillator CPU hot-path optimised: `quantizeOscOutputBeforeEnvelope` is now skipped when wave-decimation is disabled (the common case), removing a branch and two multiplies per sample across all 32-sample block-render loops in `getNextBlock`, `getNextBlockHQ`, and `getNextSample`. Decimation-enabled state is now cached once per block call rather than read per sample. Hot DSP translation units (`Osc.cpp`, `Voice.cpp`, `FxBus.cpp`, `TimbreFx.cpp`, `SimpleComp.cpp`, `SimpleEnvelope.cpp`) use selective `#pragma GCC optimize("Ofast","fast-math")` in Release builds to enable aggressive floating-point and speed optimisations without changing Debug behaviour.

Notes:

- Branch content includes both committed changes and additional in-branch FM decimation updates currently present in the workspace.
- For reproducible release output, run a clean release build with `CLEAN_BUILD=1 BUILD_CONFIG=Release ./scripts/build_cli.sh`.

## Build in VS Code (STM32CubeIDE Headless CLI)

This repository can be compiled from the VS Code terminal by using STM32CubeIDE in headless mode.

### Prerequisites

- STM32CubeIDE installed (default macOS path: `/Applications/STM32CubeIDE.app/Contents/MacOS/STM32CubeIDE`)
- ARM toolchain in PATH (`arm-none-eabi-gcc`, `arm-none-eabi-objcopy`)

### Clean Release Build

From the repository root, run:

```bash
CLEAN_BUILD=1 BUILD_CONFIG=Release ./scripts/build_cli.sh
```

The script builds in this order:

1. `preenfm3lib`
2. `preenfm3` (firmware)
3. `preenfm3 bootloader`

This order is required so firmware links against the freshly built library.

### Build Artifacts

- Firmware ELF: `firmware/Release/preenfm3.elf`
- Firmware BIN: `firmware/Release/preenfm3.bin`
- Bootloader ELF: `bootloader/Release/preenfm3 bootloader.elf`
- Bootloader BIN: `bootloader/Release/preenfm3 bootloader.bin`
- Static library: `lib/Release/libpreenfm3lib.a`

### Logs

Build logs are written to `build-logs/` with a timestamp prefix, for example:

- `<timestamp>-import.log`
- `<timestamp>-lib-clean.log`, `<timestamp>-lib-build.log`
- `<timestamp>-firmware-clean.log`, `<timestamp>-firmware-build.log`
- `<timestamp>-bootloader-clean.log`, `<timestamp>-bootloader-build.log`
- `<timestamp>-firmware-bin.log`, `<timestamp>-bootloader-bin.log`

### Notes

- The generated `.bin` files are raw binaries and do not embed a flash address.
- The firmware image is linked for the application flash offset configured in the linker script (`0x08020000`).
- The bootloader build flow compiles cleanly in CLI/VS Code, but the bootloader binary has not been tested on hardware yet.

Xavier

---

Recent VS Code/CLI build workflow and warning-cleanup related changes were done by Hans-Henning Klos.
