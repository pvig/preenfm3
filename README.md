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
- Added a new `HQ` mode for oscillator wavetable lookup: linear interpolation enabled, decimation disabled.
- HQ warped paths now use continuous warped-phase interpolation (warp before lookup interpolation), reducing residual stair-step artifacts from integer warp remapping.
- `Full` mode remains non-interpolated and non-decimated.
- FM decimation modes are capped at 19-bit precision plus `Full` and `HQ`, with legacy preset values above the supported range clamped on load for compatibility.
- Compatibility defaults now force `Full` when loading older presets without decimation settings and when importing DX7 patches.
- FM decimation path is simplified for CPU efficiency: output sample quantization is retained, phase-accumulator quantization removed, and temporal decimation (`D=2` sample-and-hold) is applied in safe oscillator render paths.
- Feedback oscillator rendering remains full-rate for stability while still using decimation output quantization.
- High-quality wavetable interpolation is bypassed when decimation is enabled.
- LFO Sync selector expanded with one-shot to eight-shot modes for both internal and external sync clocks: `1Si..8Si` and `1Se..8Se`, in addition to `Int` and `Ext`.
- LFO one-shot mode and KSyn are now independent so KSyn can be combined with shot modes.
- LFO shape list expanded in the editor/oscilloscope: `SawD`, `DExp`, `DLog`, `RExp`, `RLog`, `AD`, `AHD`, `SDec`, `Plng`, `Plg2`, `SnSq`, `Sn0`, `Sn+`, and `Usr1..Usr6`.
- LFO Phase encoder now doubles as a startup delay: negative values set a delay of 0–4000 ms (displayed as integer ms); positive values set a phase offset of 0–1. Encoder stepping uses adaptive resolution — 1 ms steps below 50 ms, up to 50 ms steps near 4 s — with a float32 forward-progress guard that prevents the encoder from getting stuck at quantisation boundaries (notably the 100 ms transition between step sizes).
- Build and release workflow improvements for VS Code/CLI headless builds (`scripts/build_cli.sh`), including automatic release artifact refresh and checksum regeneration.
- New subrelease packaging updated to firmware `v1.06g` with `bl1.09` bundle naming.
- Oscillator CPU hot-path optimised: `quantizeOscOutputBeforeEnvelope` is now skipped when wave-decimation is disabled (the common case), removing a branch and two multiplies per sample across all 32-sample block-render loops in `getNextBlock`, `getNextBlockHQ`, and `getNextSample`. Decimation-enabled state is now cached once per block call rather than read per sample. Hot DSP translation units (`Osc.cpp`, `Voice.cpp`, `FxBus.cpp`, `TimbreFx.cpp`, `SimpleComp.cpp`, `SimpleEnvelope.cpp`) use selective `#pragma GCC optimize("Ofast","fast-math")` in Release builds to enable aggressive floating-point and speed optimisations without changing Debug behaviour.
- Preset randomizer (`SynthState::randomizePreset`) overhauled with coherent cross-parameter behaviour and a wider musical range:
  - **Cross-parameter coherence**: the four randomizer controls (`Oper`, `EnvT`, `IM`, `Modl`) now influence each other — `EnvT` biases the algorithm selection (percussive patches favour feedback/carrier-heavy algos; pad patches favour complex modulator trees), oscillator mix is scaled by the `IM` level so modulator activity feels consistent, and pad mode overrides wavetable shapes toward sine/triangle/saw for smoother FM spectra.
  - **Coherent effect randomization**: `effect1` (filter) and `effect2` (modulation) are now both randomized with mode-appropriate pools. Percussive patches draw from HP/BP/saturation/fold/LP filters with drum-appropriate open cutoffs; pad patches draw from LP/shelving/tilt/stereo filters with low resonance; random/`--` mode uses the full `FILTER_TYPE` and `FILTER2_TYPE` ranges.
  - **Percussive LP filters**: LP, LP2, and LP3 are included in the percussive filter pool with a higher minimum cutoff (`p1Min = 0.55`) so bass-drum bodies are preserved while high-frequency content is still trimmed.
  - **Modulation effect pool**: `effect2` is randomized per mode — pads get chorus/dimension/wide/diffuser/doubler at high probability when `Modl > 0`; percussive patches get a low-probability flange or grain hit; random mode draws from the full `FILTER2_TYPE` range.
  - **LFO coherence**: percussive patches assign one-shot LFO shapes (exponential/log decay, attack-decay, Buchla) at audio-rate-adjacent frequencies with phase stagger; pad patches use slow smooth shapes (sine/triangle/saw) with long keyboard ramps.
  - **Matrix destination intelligence**: safe and advanced destination pools are filtered at selection time — oscillator-specific destinations are skipped when the oscillator is silent (mix + weighted IM contribution < 0.03 threshold); `FILTER1_PARAM1/2/AMP` and `FILTER2_PARAM1/2/AMP` are skipped unless the corresponding effect slot is active.
  - **Unique destination selection**: a `destinationUsed` tracker prevents duplicate matrix destination assignments; a `rowUsed` tracker prevents reuse of the same modulation source-row across randomizer passes.
  - **Destination-aware modulation depth**: matrix multiply values are scaled by destination type — pitch/phase/warp destinations get narrower ranges; pan/mix wider; catch-all destinations wider still; negative modulation is allowed at high `Modl` except for percussive patches.
  - **Anti-clipping output normalization**: `effect1.param3` (the mixer gain applied to every voice sample) is now divided by the algorithm's carrier count (`algoInformation[algo].mix`), keeping the summed output at roughly the same peak level regardless of whether the algorithm has 1 or 6 carrier operators.
  - **Infinite-sustain prevention**: pad envelope times are capped (decay 0.5–3 s, sustain ramp 0.5–2 s, release 1–5 s, down from 5/5/8 s ceilings); `releaseLevel` is fixed at `0.0` for all envelope modes so voices always decay fully to silence. This also closes the modulator-loop edge case where `releaseLevel == 1.0` combined with `releaseTime == 0.0` could cause an envelope to loop indefinitely.
  - **LFO DC-offset elimination**: `osc->bias` is set to `0.0` for all three LFO modes. A non-zero bias was causing one-shot shapes to freeze at a non-zero amplitude (e.g. a percussive decay freezing at `terminal + bias` instead of at zero), and continuous LFOs routed to release or amplitude destinations to add a persistent DC offset that lengthened notes or prevented voices from going fully silent.
  - **Expanded LFO shape pools** (curated per mode): percussive mode now draws from 7 shapes (`DECAY_EXP`, `DECAY_LOG`, `DECAY_S`, `ATTACK_DECAY`, `ATTACK_HOLD_DECAY`, `BUCHLA_PLONG`, `BUCHLA_PLONG2`) — all with `oneShotTerminalShapeValue = −1` so they freeze cleanly at zero; pad mode draws from 7 smooth zero-mean shapes (`SIN`, `TRIANGLE`, `SAW`, `SAW_DOWN`, `BROWNIAN`, `WANDERING`, `FLOW`); random/`--` mode draws from 9 shapes adding `SAW_DOWN`, `BROWNIAN`, `WANDERING`, and `FLOW` to the previous set. Always-positive shapes (`SIN_POS`, `SIN_ZERO`, `SIN_SQUARE`) and positive-terminal one-shots (`RISE_EXP`, `RISE_LOG`) are excluded from all pools.
  - Hardware RNG (`HAL_RNG_GenerateRandomNumber`) with LCG fallback (`HAL_GetTick()`-seeded) for all random draws throughout the randomizer.

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
