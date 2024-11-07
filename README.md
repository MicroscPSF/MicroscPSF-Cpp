# C++17 port of the MicroscPSF-Matlab project

## Status

Currently, this work is a *very* early work in progress. Refer to https://github.com/MicroscPSF/MicroscPSF-Matlab

## Quick start

### Ubuntu/Linux

Install the compiler toolchain:

```bash
sudo apt install build-essentials
```

Install Meson the build system:

```bash
cd MicroscPSF-Cpp/
python3 -m venv .venv/
.venv/bin/pip3 install meson ninja
```

Compile everything

```bash
cd MicroscPSF-Cpp/
meson setup build/
ninja -C build all
```

Test everything

```bash
cd MicroscPSF-Cpp/build/
ninja test
```

(Optional) Install the example app

```bash
ninja install
cd MicroscPSF-Cpp/build/
meson configure -Dbuild_examples=true
ninja all
sudo ninja install
```

## Appendix: Bessel function support

- The original ISO C++ proposal: https://wg21.link/p0226r1
- GCC compiler support status: https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#status.iso.2017
- Clang support status: https://libcxx.llvm.org/Status/Cxx17.html
- Clang support ticket: https://github.com/llvm/llvm-project/issues/99939