# preenfm3

Here are the sources for the preenfm3 bootloader and firmware.

Binaries can be found in the [Release section](https://github.com/Ixox/preenfm3/releases).

For the documentation, go to [the Wiki](https://github.com/Ixox/preenfm3/wiki).

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
