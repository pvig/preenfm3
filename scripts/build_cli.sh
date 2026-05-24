#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CUBEIDE_BIN="${CUBEIDE_BIN:-/Applications/STM32CubeIDE.app/Contents/MacOS/STM32CubeIDE}"
BUILD_CONFIG="${BUILD_CONFIG:-DebugLQFP144}"
CLEAN_BUILD="${CLEAN_BUILD:-0}"
OBJCOPY_BIN="${OBJCOPY_BIN:-arm-none-eabi-objcopy}"

LOG_DIR="${1:-$ROOT_DIR/build-logs}"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
CUBEIDE_WORKSPACE="${CUBEIDE_WORKSPACE:-/tmp/pfm3-cubews-$TIMESTAMP}"

mkdir -p "$LOG_DIR"

if [[ ! -x "$CUBEIDE_BIN" ]]; then
  echo "STM32CubeIDE binary not found or not executable: $CUBEIDE_BIN" >&2
  echo "Set CUBEIDE_BIN to your STM32CubeIDE launcher path." >&2
  exit 1
fi

run_headless() {
  "$CUBEIDE_BIN" \
    -nosplash \
    -consolelog \
    -application org.eclipse.cdt.managedbuilder.core.headlessbuild \
    -data "$CUBEIDE_WORKSPACE" \
    "$@"
}

import_projects_once() {
  local log_file="$LOG_DIR/${TIMESTAMP}-import.log"

  echo "=== Importing projects into workspace ==="
  echo "Workspace: $CUBEIDE_WORKSPACE"
  echo "Log: $log_file"

  if run_headless \
      -import "$ROOT_DIR/lib" \
      -import "$ROOT_DIR/firmware" \
      -import "$ROOT_DIR/bootloader" \
      >"$log_file" 2>&1; then
    echo "OK: project import"
  else
    echo "FAILED: project import" >&2
    echo "Last lines from ${log_file}:" >&2
    tail -n 80 "$log_file" >&2 || true
    exit 1
  fi
}

run_build_step() {
  local project_config="$1"
  local log_name="$2"
  local clean_log="$LOG_DIR/${TIMESTAMP}-${log_name}-clean.log"
  local build_log="$LOG_DIR/${TIMESTAMP}-${log_name}-build.log"

  if [[ "$CLEAN_BUILD" == "1" ]]; then
    echo "=== ${project_config} (-cleanBuild) ==="
    echo "Log: $clean_log"

    if run_headless \
        -cleanBuild "$project_config" \
        >"$clean_log" 2>&1; then
      echo "OK: ${project_config} clean"
    else
      echo "FAILED: ${project_config} clean" >&2
      echo "Last lines from ${clean_log}:" >&2
      tail -n 80 "$clean_log" >&2 || true
      exit 1
    fi
  fi

  echo "=== ${project_config} (-build) ==="
  echo "Log: $build_log"

  if run_headless \
      -build "$project_config" \
      >"$build_log" 2>&1; then
    echo "OK: ${project_config} build"
  else
    echo "FAILED: ${project_config} build" >&2
    echo "Last lines from ${build_log}:" >&2
    tail -n 80 "$build_log" >&2 || true
    exit 1
  fi
}

convert_flash_bin() {
  local elf_file="$1"
  local bin_file="$2"
  local log_name="$3"
  local log_file="$LOG_DIR/${TIMESTAMP}-${log_name}.log"

  if [[ ! -f "$elf_file" ]]; then
    echo "ELF not found, skipping binary conversion: $elf_file"
    return
  fi

  echo "=== Converting $(basename "$elf_file") to $(basename "$bin_file") ==="
  echo "Log: $log_file"

  if "$OBJCOPY_BIN" \
      -R .ram_d2b \
      -R .ram_d2 \
      -R .ram_d1 \
      -R .ram_d3 \
      -R .instruction_ram \
      -O binary \
      "$elf_file" \
      "$bin_file" \
      >"$log_file" 2>&1; then
    echo "OK: $bin_file"
  else
    echo "FAILED: binary conversion for $elf_file" >&2
    echo "Last lines from ${log_file}:" >&2
    tail -n 80 "$log_file" >&2 || true
    exit 1
  fi
}

if ! command -v "$OBJCOPY_BIN" >/dev/null 2>&1; then
  echo "objcopy not found: $OBJCOPY_BIN" >&2
  echo "Set OBJCOPY_BIN to your arm-none-eabi-objcopy path." >&2
  exit 1
fi

import_projects_once

run_build_step "preenfm3lib/${BUILD_CONFIG}" "lib"
run_build_step "preenfm3/${BUILD_CONFIG}" "firmware"
run_build_step "preenfm3 bootloader/${BUILD_CONFIG}" "bootloader"

convert_flash_bin "$ROOT_DIR/firmware/${BUILD_CONFIG}/preenfm3.elf" "$ROOT_DIR/firmware/${BUILD_CONFIG}/preenfm3.bin" "firmware-bin"
convert_flash_bin "$ROOT_DIR/bootloader/${BUILD_CONFIG}/preenfm3 bootloader.elf" "$ROOT_DIR/bootloader/${BUILD_CONFIG}/preenfm3 bootloader.bin" "bootloader-bin"

echo ""
echo "All builds succeeded."
echo "Logs written to: $LOG_DIR"
