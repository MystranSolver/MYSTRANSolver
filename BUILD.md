# Building MYSTRAN from source

###### Last updated 2023-12-18.

## Setting up a build environment

In order to build (compile) MYSTRAN using CMake, you first have to set up a
proper build environment (i.e. toolchain and required programs/libraries).

You can skip this part if you've done it already (or if you really know what
you're doing).

### Steps for Windows (x86_64)

First, download and install MSYS2 from the
[official site](https://www.msys2.org/).

Open the MSYS2 terminal and run the following commands:

  1. **`pacman -Syu`**
This updates repository information and installed packages, and might require
you close and reopen MSYS2 terminals.
  1. **`pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja mingw-w64-x86_64-openblas git`**
This installs the required compilers (the GNU C and Fortran compilers), CMake
itself, the Ninja build tool, OpenBLAS (the recommended BLAS/LAPACK
provider), and `git`.
  1. **`export PATH="/mingw64/bin:$PATH"`**
This makes the MinGW toolchain programs (such as `make` and the compilers)
visible so CMake can find them more easily. Note that this command's effects
are lost when you reopen the terminal, so you might want to append it to your
`~/.bashrc` to save time.

### Steps for Linux (any)

Follow your distribution's steps to install the following programs/libraries:
  - **`gcc`**
  - **`g++`**
  - **`gfortran`**
  - **`ninja`** (recommended; package is usually `ninja-build`)
  - **`cmake`**
  - **`git`**
  - **`openblas`** (recommended BLAS/LAPACK provider; package is usually
    `libopenblas-dev` on Debian/Ubuntu or `openblas-devel` on Fedora/RHEL)

All of those are fairly common, so get in touch in the MYSTRAN Forums or
MYSTRAN Discord if you have trouble installing any of them. Also, note that
most distros have a "base" package group for developers (e.g. Arch's
`base-devel` or Ubuntu's `build-essential`) that includes necessary tooling
such as `gcc` and `make`. If that's the case, install it!

If your distribution doesn't ship CMake 3.18+ yet, check if your distro has a
some sort of testing/unstable channel before attempting to
[install it manually](https://cmake.org/install/).

For WSL (Linux for Windows)
===========================
Mystran won't work with Ubuntu 20.04, hasn't been tested on 22.04 and should work on 24.04 (what we're testing).

If you're upgrading your WSL, open PowerShell as Administrator and run:
```
wsl --update
wsl --install --distribution Ubuntu-24.04
```

Now that you've got into a modern version of Ubuntu
```
sudo apt update
sudo apt upgrade
sudo apt install gcc g++ gfortran ninja-build cmake git libopenblas-dev
```

---

## Building MYSTRAN

If your build environment is already set up, building MYSTRAN is quite
straightforward.

### Steps for Windows (any)

  1. Open the MSYS2 shell.
  2. Re-run step #3 of the previous section if needed.
  3. Fetch the source code if you haven't already. If you're using Git, you can
  clone the repo with
  **`git clone https://github.com/MYSTRANsolver/MYSTRAN.git`**.
  4. Move the terminal to the MYSTRAN folder. If you've just run `git clone`,
     just do a **`cd MYSTRAN`**.
  5. Generate the build scripts by running **`cmake -G Ninja .`**.
  6. Compile with **`cmake --build .`**. If you have an N-core processor,
  running **`cmake --build . -jN`** will probably be much faster (Ninja already
  parallelizes by default, but `-jN` lets you cap it). You can find the number
  of cores/threads with the `nproc` command.
  7. The executable will reside at **`Binaries/mystran.exe`**.

### Steps for Linux (any)

  1. Open a terminal.
  2. Fetch the source code if you haven't already. If you're using Git, you can
  clone the repo with
  **`git clone https://github.com/MYSTRANsolver/MYSTRAN.git`**.
  3. Move the terminal to the MYSTRAN folder. If you've just run `git clone`,
  just do a **`cd MYSTRAN`**.
  1. Generate the build scripts by running **`cmake -G Ninja .`**.
  2. Compile with **`cmake --build .`**. If you have an N-core processor,
  running **`cmake --build . -jN`** will probably be much faster (Ninja already
  parallelizes by default, but `-jN` lets you cap it). You can find the number
  of cores/threads with the `nproc` command (not all distros ship it
  out-of-the-box though).
  1. The executable will reside at **`Binaries/mystran`**.

---

### CPU architecture and tuning

All C, C++, and Fortran code compiled by this project defaults (on x86_64
targets) to `-march=x86-64 -mtune=generic`, including the bundled GKlib, METIS,
SuperLU, and embedded BLAS/LAPACK when applicable.

Use the CMake cache settings `MYSTRAN_MARCH` and `MYSTRAN_MTUNE` to change
the instruction set and tuning independently, e.g.:

```sh
cmake -G Ninja . -DCMAKE_BUILD_TYPE=Release \
  -DMYSTRAN_MARCH=x86-64-v3 -DMYSTRAN_MTUNE=generic
cmake --build . -j8
```

For a host-specific experiment, use `-DMYSTRAN_MARCH=native
-DMYSTRAN_MTUNE=native`. Such binaries can be somewhat faster, but may not run
on other (older) CPUs!

Restore the portable defaults with:

```sh
cmake . -DMYSTRAN_MARCH=x86-64 -DMYSTRAN_MTUNE=generic
cmake --build . -j8
```

CPU names must be accepted by all three compilers; unsupported settings
fail at configure time. Set either option to an empty string to omit that
flag and use the compiler/toolchain default. On non-x86 targets or
toolchains without these flags, select appropriate supported values or
disable both with `-DMYSTRAN_MARCH= -DMYSTRAN_MTUNE=`.
Avoid adding competing `-march`/`-mtune` flags through `CMAKE_*_FLAGS`.

---

## Troubleshooting

While this process is meant to be straightforward, here is a list of some of
the more common issues that can arise. Other issues users find might be added
here if they're not too specific.

If your issue isn't here, you can always ask for help at the
[MYSTRAN forums](https://www.mystran.com/forums/) or the
[Discord server](https://discord.gg/9k76SkHpHM)

---

### "I'm getting "file not found" errors when running the step #2 setup command!"

Run a **`pacman -Syyu`** (note the two 'y's) and try again.

---

### "CMake is complaining about not being able to find the toolchain or the Fortran compiler or the build tool!"

Try running the commands `ninja`, `gcc`, and `gfortran`. If any of these comes
up as a "command not found", make sure they've been installed. If you're
**sure** they are, they might not be in the PATH.

Windows users, have a look at step #3 of the setup. Linux users, check out your
distro documentation, because whatever's happening should not be happening at
all.

---

### "CMake complains about `ARCHIVE_EXTRACT`!"

Check out the output of `cmake --version`. You must have version 3.18 or newer.
If you don't, first ensure it's up to date -- perform a system-wide update.
Windows users should not find this issue relevant -- MSYS2 ships CMake 3.27.1
as of this writing. Linux users should use their own package manager.

If your system is up to date and you still run into this issue, that means your
distro ships CMake 3.17 or older. Bad luck there. Here's what you can do:

  1. Enable a testing/unstable package channel (not all distros have one)
  2. Install the latest CMake [manually](https://cmake.org/install/)
  (might piss off your package manager)
  1. Download and extract `libf2c.zip` yourself, and comment out the
  `ARCHIVE_EXTRACT` stuff in `CMakeLists.txt`.

---

### "I'm getting random SuperLU build errors!"

SuperLU is included as a submodule. A recent update to the submodule might
require a clean build. Run `cmake --build . --target clean` and delete the
`superlu` subdirectory and run the appropriate `cmake` command again.

---

### "I'm getting cryptic linker errors related to BLAS or LAPACK!"

MYSTRAN and SuperLU both need BLAS and LAPACK. The build system
picks one of three providers via the `MYSTRAN_BLAS_LAPACK` CMake
option:

  - **`AUTO`** (default): try to locate a system BLAS and LAPACK (we
    recommend OpenBLAS); if either isn't found, fall back to the
    bundled Reference-LAPACK submodule.
  - **`SYSTEM`**: require system BLAS and LAPACK; configuration
    fails with a clear error if either cannot be found.
  - **`EMBEDDED`**: ignore the system entirely and always build the
    bundled Reference-LAPACK submodule (`Source/lapack/`), which
    provides both BLAS and LAPACK.

If the auto-detection picks up libraries that do not ship a static
`.a` archive (a common Windows situation), re-run CMake with
`-DMYSTRAN_BLAS_LAPACK=EMBEDDED` to force the bundled fallback.

On Windows we ship fully static binaries, so when requesting
`SYSTEM` mode you must have a static OpenBLAS available
(MSYS2 / MinGW64: `pacman -S mingw-w64-x86_64-openblas`).

The legacy `-Denable_internal_blaslib=YES` and `-DMYSTRAN_BLAS=...`
flags still work; they are mapped to `MYSTRAN_BLAS_LAPACK` with a
deprecation warning.

The bundled Reference-LAPACK is considerably slower than a tuned
implementation like OpenBLAS or MKL. That can have a significant
impact on the time it takes to run larger models. A small set of
mathematically-deviated routines (kept in
`Source/Modules/MYSTRAN_LAPACK_EXT/`) is always compiled into
MYSTRAN regardless of which provider is selected.

---

### "I want to build offline, but the CMake script attempts to download stuff!"

Download the `superlu` submodule and `libf2c.zip` beforehand, and you should be
fine.

---

### "The terminal output is garbled during compilation!"

Multiple threads are printing to standard output simultaneously. Ninja
serializes per-job output by default, so this should not happen with the
recommended generator. If you're using the legacy `make`/`mingw32-make`
generators and see garbled output, switch to Ninja (`cmake -G Ninja .`).

If you *really* need to keep using `make` and want readable output, ensure it
only runs with one thread by passing the option `-j1`. This will make
compilation slower, but at least you'll be able to read the output.

And if it's errors you're looking for, you can build fast with `-j[number]`,
and then `-j1` just to see the error again.

---

If your issue isn't here, you can always ask for help at the
[MYSTRAN forums](https://www.mystran.com/forums/) or the
[Discord server](https://discord.gg/9k76SkHpHM)
