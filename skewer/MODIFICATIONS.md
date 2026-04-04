# Skewer Local Modifications

Vendored from: [relipmoc/skewer](https://github.com/relipmoc/skewer) v0.2.2 (2016)
License: MIT

## Why vendored

ReSeq uses only the core bit-masked k-difference adapter matching algorithm
(~1,200 lines from `matrix.cpp/h` and `fastq.cpp/h`). The vendored copy
carries local modifications that are not suitable for upstreaming.

## Modifications

1. **Namespace wrapping** — All code wrapped in `namespace skewer {}` to
   avoid global namespace pollution when linking as a static library.

2. **Removed unused code** (~4,700 lines) ��� Deleted `main.cpp` (CLI tool),
   `parameter.cpp/h` (argument parsing), and all code not used by ReSeq's
   `AdapterStats::Detect()` pathway.

3. **CMakeLists.txt** — Created minimal CMake build config producing a
   static library (`skewer_matrix`) from `matrix.cpp` and `fastq.cpp`.

4. **GCC 15 compatibility** — Added `const` qualifier to
   `ElementComparator::operator()` for C++20 compliance.
