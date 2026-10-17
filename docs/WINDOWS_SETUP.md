# Windows with CLion's existing bundled MinGW

This setup uses an existing CLion installation and Git for Windows. It creates
the GMP C++ library and Boost headers under the ignored `build/deps` directory.
It also creates ignored local Release/Debug CMake presets when
`CMakeUserPresets.json` is absent. An existing presets file is retained.

The tested toolchain is CLion 2026.1's bundled GCC 15.2.0, targeting
`x86_64-w64-mingw32` with the MSVCRT runtime. Its bundled static GMP C library is
6.3.0, using 64-bit `unsigned long long` limbs without nails. Setup checks these
GMP version and ABI assumptions by compiling and running a native program before
downloading sources. An incompatible installation stops with an error.

## Prepare dependencies and local presets

Open PowerShell in the project root:

```powershell
pwsh -NoProfile -File .\tools\setup-clion-deps.ps1 -ClionRoot 'C:\Program Files\JetBrains\CLion 2026.1'
```

Use the installed PowerShell 7 (`pwsh`) for this command. Omit `-ClionRoot` to
search installed `CLion*` folders under
`C:\Program Files\JetBrains`. For a Toolbox or custom installation, pass its
installation folder explicitly. If Git for Windows is in a custom location,
also pass `-GitBash 'D:\Apps\Git\bin\bash.exe'`. The script uses the bundled
compiler and existing Git Bash; it installs no Windows packages and changes no
system PATH entries. `-NoPresets` prepares dependencies only.

The first run downloads about 150 MB of source archives and can take several
minutes while GMP's configuration checks run in Git Bash. Subsequent runs reuse
verified archives and a configuration cache. Successful setup prints both
`Bundled GMP ... ABI: PASS` and `GMP C++ streams and Boost ... backend: PASS`.
Generated code, logs, archives, libraries, and headers stay under `build/deps`.

If the current local presets point to another installation, move that file
aside and rerun setup to generate presets for the selected CLion installation.
The generated names are `windows-clion` and `windows-clion-debug`; the matching
build folders are `build/windows-clion` and `build/windows-clion-debug`.

## Open and build in CLion

1. Open the project folder containing `CMakeLists.txt`.
2. In **Settings > Build, Execution, Deployment > Toolchains**, select the
   bundled **MinGW** toolchain. The generated presets refer to the toolchain
   named **MinGW**.
3. In **Settings > Build, Execution, Deployment > CMake**, enable
   **windows-clion**, then reload CMake. CLion may import a preset but initially
   leave its profile disabled.
4. Build `secular` or the complete project. Use **windows-clion-debug** when
   stepping through the C++ code.

For command-line verification in PowerShell, use the same CLion installation:

```powershell
$clionRoot = 'C:\Program Files\JetBrains\CLion 2026.1'
$cmake = Join-Path $clionRoot 'bin\cmake\win\x64\bin\cmake.exe'
$ctest = Join-Path $clionRoot 'bin\cmake\win\x64\bin\ctest.exe'
$env:PATH = (Join-Path $clionRoot 'bin\mingw\bin') + ';' + $env:PATH
& $cmake --preset windows-clion
& $cmake --build --preset windows-clion --parallel
& $ctest --preset windows-clion --output-on-failure
```

Python 3 is optional for the solver build, and required for the CLI and matrix
integrity tests. CMake reports when these optional tests are unavailable.

## Run the solver

The compiler's `bin` folder must also be on the executable's run/debug PATH so
Windows can find its matching `libstdc++-6.dll` and `libgcc_s_seh-1.dll`. The
PowerShell example above sets PATH only for that shell and its child processes.
For CLion, check its toolchain/run environment or add
`Path=<CLion installation>\bin\mingw\bin;$Path$` to the run configuration's
environment. There is no GMP DLL: this setup links both GMP libraries statically.

Start with the shared **01 - Small SPD example** run configuration, or run:

```powershell
.\build\windows-clion\secular.exe solve --matrix matrices/examples/spd3.mtx --digits 50 --out results/spd3
```

The output directory must be new or empty. The README documents the full
oil-rig commands and the other shared run configurations.

## Sources and build details

The setup pins official source archives and checks SHA256 before extraction:

| Source | SHA256 |
| --- | --- |
| [GMP 6.3.0](https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz) | `a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898` |
| [Boost 1.90.0](https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.gz) | `5e93d582aff26868d581a52ae78c7d8edf3f3064742c6e77901a1f18a437eea9` |

Boost publishes its hash in the [archive metadata](https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.gz.json).
The GMP hash is pinned to the GNU HTTPS archive used for this setup; this script
does not claim detached-signature verification. Dependency license files are
copied under `build/deps/mingw/share/licenses`.

The existing bundled `gmp.h` and `libgmp.a` are copied into the local prefix.
Only the eleven upstream GMP C++ wrapper sources are compiled, using the same
CLion compiler. Generic configuration avoids an m4 requirement from GMP's
unused timing-tool assembly. The wrapper then uses the checked bundled GMP
header, its 8-byte limb configuration, and internal tables generated for
64-bit limbs. The script does not rebuild GMP's C library. A final native smoke
test checks large integers, rational reduction, GMP C++ streams, and Boost's
GMP floating-point backend. Project tests remain the full application check.
