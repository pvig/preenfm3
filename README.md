# preenfm3

Bootloader and firmware sources for the **preenfm3** FM synthesizer (STM32H7).

- Binaries: [Release section](https://github.com/pvig/preenfm3/releases)
- General documentation: [Official Wiki](https://github.com/Ixox/preenfm3/wiki)

---

## Features & Enhancements Developed in This Repository

This fork integrates the developments from [gresade/preenfm3](https://github.com/gresade/preenfm3) and adds a considerable amount of refactoring, bug fixes, and optimizations across almost every aspect of the preenfm3 (audio engine, filters, FX bus, sequencer, arpeggiator, modulation system, and web-editor integration).

### 1. FM Engine, Oscillators & Anti-Aliasing

- **Per-operator start phase (`o1Ph`–`o6Ph`)**: configurable phase offset at note-on, with visual marker display in the oscilloscope view and modulation-matrix routing (Gresade).
- **Per-operator `Warp` (`o1Wr`–`o6Wr`)**: asymmetrical stretching of the two half-cycles of the wavetable (range ±4), mapped to MIDI CC, exposed in the modulation matrix, and applied across all render paths (`getNextSample`, `getNextBlock`, `getNextBlockHQ`) (Gresade).
- **`FULL` and `HQx2` quality modes** (`HQ` is now legacy):
  - `FULL`: non-decimated rendering (classic behaviour; legacy `HQ` presets map accordingly).
  - `HQx2`: 2× oversampling to significantly reduce aliasing on rich FM modulations.
- **Windowed Sync**: lower-aliasing windowed oscillator sync.
- **FM Decimation**: output quantization control up to 19-bit precision (with backward-compatible defaults for legacy presets and DX7 imports) (Gresade).
- **DSP & memory optimizations**: selective `#pragma GCC optimize("Ofast", "fast-math")` on hot DSP translation units in Release builds, fast-path bypass when decimation is off, freed `DTCMRAM`, and compile-time `PFM3_FAST_FEEDBACK_INTERP` switch in `Osc.h`.

### 2. Per-Voice Filters (`VoiceFx`)

- **ZDF SVF topology & modulation anti-aliasing**: Zero-Delay Feedback state-variable filters with **per-sample coefficient interpolation** across the entire audio block (`BLOCK_SIZE = 32`), eliminating stepping and modulation aliasing.
- **Transition safety**: synchronous filter-type change detection inside the audio loop, strict parameter clamping, and non-recursive `a1` interpolation to prevent numerical blow-ups (`NaN`/`Inf`) when switching between heterogeneous effects (e.g. `Crsh` → `BP`).
- **ADAA saturators & wavefolders**: analytical Antiderivative Anti-Aliasing (ADAA) used across `Sat`, `Sigm`, `Fold`, `Wrap`, `LpWs`, `Oryx`, `Lad+`, and `SvfLH`.
- **Centered & calibrated `Wrap`**: symmetric phase wrapping centered at zero (with exact ADAA antiderivative), progressive convex drive curve around the `±1` threshold, and RMS level compensation.
- **`Ap4D` (True Dual Nested Allpass Phaser)**: genuine 4-pole cascaded dual nested allpass topology with soft-saturated resonant feedback, stereo cross-feedback, and per-sample interpolation without startup sweep artifacts.
- **TB-303 emulation (`TeeB`)**: high-accuracy Padé approximation for frequency-to-pole mapping and dynamic Accent circuit modelling (C13 "gummy capacitor" + MIDI velocity).
- **`Oryx` series (`Oryx`, `Oryx2`, `Oryx3`)**: tamed resonance peaks and distinct sonic personalities across variants.

### 3. Timbre Effects, FX Bus & Compressor

- **Timbre effects (`TimbreFx`)**: engine cleanup and Hermite delay-line interpolation.
- **FX Bus (`FxBus`)**: reworked warm soft-knee saturation, overdrive, and fuzz algorithms for improved clarity and musicality.
- **Compressor**: dedicated presets (`Punch`, `Glue`, `Drive`, etc.) with refined dynamic response.

### 4. Envelopes, LFOs, Modulation & Web Editor

- **Envelopes**: adjustable Attack and Release curves (`Attack & Release Curves`).
- **Expanded LFOs**:
  - Multi-shot sync modes `1Si..8Si` (internal) and `1Se..8Se` (external), independent from `KSyn`.
  - Additional shapes: `SawD`, `DExp`, `DLog`, `RExp`, `RLog`, `AD`, `AHD`, `SDec`, `Plng`, `Plg2`, `SnSq`, `Sn0`, `Sn+`, and `Usr1..Usr6`.
  - Startup delay (0–4000 ms) on negative `Phase` encoder values with adaptive step resolution.
- **Matrix & MIDI**: new destinations (`note+`, `o1Ph..o6Ph`, `o1Wr..o6Wr`, etc.) and refactored **MIDI CC / NRPN** mapping for web-editor adaptation.

### 5. Arpeggiator, Step Sequencer & Voice Management

- **Enhanced arpeggiator**: `Oct+` mode, **Ratchet**, **Strumming**, accent masks, Euclidean patterns, `Arp+`, `ArpVl`, `ArpGt`, `ArpDv` parameters, and preset caching.
- **Step sequencer**: configurable loop length, `Alt2` / `Alt3` step actions, LFO modes, and per-step ratchet.
- **Voice allocation**: **Bass-protect voice stealing** and configurable per-timbre voice priority.

### 6. Preset Randomizer & Patch Compatibility

- **Coherent randomizer (`SynthState::randomizePreset`)**: cross-parameter coupling across `Oper`, `EnvT`, `IM`, and `Modl`, mode-appropriate effect and LFO pools (percussive / pad / random), carrier-count output normalization, LFO DC-offset safety, infinite-sustain prevention, and STM32 hardware RNG (`HAL_RNG`).
- **Patch format `v1.6f`** (updated from `v1.5`): backward-compatible loading of older patches and hardened DX7 SysEx import.

---

## Building the Firmware

### Prerequisites

- **STM32CubeIDE 1.5.1** with its bundled **GNU Tools for STM32 (`7-2018-q2-update`)** toolchain (`arm-none-eabi-gcc` / `arm-none-eabi-g++` **7.3.1**, `arm-none-eabi-objcopy`) in `PATH`

### Build Order

Projects must be built in the following order so the firmware links against the freshly built static library:

1. `preenfm3lib` (`lib/`)
2. `preenfm3` (`firmware/`)
3. `preenfm3 bootloader` (`bootloader/`, optional)

You can build directly from **STM32CubeIDE** or by running `make` inside `lib/Release`, `firmware/Release`, and `bootloader/Release`.

*(An optional headless CLI script `CLEAN_BUILD=1 BUILD_CONFIG=Release ./scripts/build_cli.sh` is also available in `scripts/`.)*

### Build Artifacts

- Firmware ELF / BIN: `firmware/Release/preenfm3.elf`, `firmware/Release/preenfm3.bin`
- Bootloader ELF / BIN: `bootloader/Release/preenfm3 bootloader.elf`, `bootloader/Release/preenfm3 bootloader.bin`
- Static library: `lib/Release/libpreenfm3lib.a`

### Notes

- The generated `.bin` files are raw binaries and do not embed a flash address.
- The firmware image is linked for the application flash offset configured in the linker script (`0x08020000`).
- The bootloader build flow compiles cleanly in CLI/VS Code, but the bootloader binary has not been tested on hardware yet.

---
