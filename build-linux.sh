#!/bin/bash
# Build astro-catalog binaries for Linux using Docker

set -e

echo "Building astro-catalog for Linux (x86_64)..."

docker run --rm \
  -v "$(pwd):/workspace" \
  -w /workspace \
  rust:latest \
  bash -c "
    cargo build --bin build-catalog --release --target x86_64-unknown-linux-gnu
    cargo build --bin query-catalog --release --target x86_64-unknown-linux-gnu
  "

echo " Build complete!"
echo "Binaries located at:"
echo "  target/x86_64-unknown-linux-gnu/release/build-catalog"
echo "  target/x86_64-unknown-linux-gnu/release/query-catalog"
