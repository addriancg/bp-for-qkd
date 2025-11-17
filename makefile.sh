#!/usr/bin/env bash
set -euo pipefail

: "${CC:=gcc}"
: "${CFLAGS:=-std=c11 -O3 -Wall -Wextra -Wpedantic}"
: "${OUT:=main}"

CPPFLAGS="-Isrc -Iinclude -I."

SRC=(
  "src/bp_decoder.c"
  "src/dvbs2ldpcShort.c"
  "src/matrixUtils.c"
  "src/getchecknodetable.c"
  "main.c"
)

echo "[1/2] Compilando → $OUT"
$CC $CFLAGS $CPPFLAGS -o "$OUT" "${SRC[@]}" -lm

echo "[2/2] Ejecutando → ./$OUT $*"
./"$OUT" "$@"

