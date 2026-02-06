# Pre-built binaries

The easiest way to get started with Delphy is to use pre-built binaries. No build tools or dependencies are required.

### GitHub Releases

Pre-built binaries are available as tarballs attached to [GitHub Releases](https://github.com/broadinstitute/delphy/releases) for the following platforms:

- `delphy-linux-x86_64.tar.gz` — Linux x86_64 (statically linked)
- `delphy-linux-arm64.tar.gz` — Linux ARM64 (statically linked)
- `delphy-macos-arm64.tar.gz` — macOS ARM64 (Apple Silicon)
- `delphy-wasm.tar.gz` — WebAssembly bundle

Each tarball contains `delphy`, `delphy_mcc`, `beast_trees_to_dphy`, and `delphy_ui`.

To install, download and extract the tarball for your platform:
```
tar -xzvf delphy-linux-x86_64.tar.gz
./delphy --version
```

### Docker

Multi-arch Docker images (linux/amd64 and linux/arm64) are available from GitHub Container Registry:

```
docker run --rm ghcr.io/broadinstitute/delphy:latest delphy --version
```

Tagged releases are available as `ghcr.io/broadinstitute/delphy:<version>`.

To run an inference on a local FASTA file:
```
docker run --rm -v "$PWD":/data ghcr.io/broadinstitute/delphy:latest \
    delphy --v0-in-fasta /data/input.fasta \
    --v0-init-heuristic \
    --v0-steps 10000000 \
    --v0-out-log-file /data/delphy.log \
    --v0-log-every 200000 \
    --v0-out-trees-file /data/delphy.trees \
    --v0-tree-every 2000000
```

### Google Colab tutorials

Ready-to-use Google Colab tutorials that download and run pre-compiled Delphy binaries are also available — see [Google Colab Tutorials](README.md#google-colab-tutorials) in the README for details.

---

# Building from source

The remainder of this document describes how to build Delphy from source.

# Compatibility

These build instructions were last tested on
- Ubuntu 24.04 x86_64 (on a fresh AWS c7a.2xlarge instance, after running `apt update && apt upgrade` and restarting)
- Ubuntu 22.04 x86_64 (on a fresh AWS c7a.2xlarge instance, after running `apt update && apt upgrade` and restarting)
- macos Sequoia 15.6.1 x86-64 (on a fresh AWS mac1.metal instance)
- macos Sequoia 15.6.1 ARM (on a fresh AWS mac2.metal instance)

# Dependencies

First, install build dependencies.

### On Ubuntu:
```
sudo apt install cmake python3-venv g++  # Ubuntu (will bring on many other build dependencies)
```
On Ubuntu 22.04, the distribution's cmake version is too old.  Use a more recent version as follows:
```
mkdir -p ~/tools
cd ~/tools/
wget https://github.com/Kitware/CMake/releases/download/v3.31.9/cmake-3.31.9-linux-x86_64.tar.gz
tar -xvf cmake-3.31.9-linux-x86_64.tar.gz
export "PATH=${HOME}/tools/cmake-3.31.9-linux-x86_64/bin:${PATH}"
```

### On Mac
On macos, we only support installations that leverage Homebrew.
```
brew install cmake emscripten              # macos (AWS macos AMIs have lots of build tooling preinstalled, but not cmake)
```

### Other notes
Some third-party libraries (e.g., abseil) are linked directly in the source tree as git submodules, whereas more complex dependencies use Conan (see [conanfile.txt] for details).


# Building locally

Create and activate a python3 virtualenv to install packages
```
python3 -m venv delphy-venv      # Only needed once
source delphy-venv/bin/activate  # Needed every time you open up a new shell
```

Install Conan 2.8.1
```
pip3 install 'conan==2.8.1'   # See https://docs.conan.io/2/installation.html for details
```

Make some config adjustment to your Conan profile (see https://docs.conan.io/2/tutorial/consuming_packages/build_simple_cmake_project.html)
```
conan profile detect
```
On macos, you may have to adjust the resulting Conan profile in ~/.conan2/profiles/default to say "compiler.version=16", not "17" or later, for Conan 2.8.1 to recognize your compiler.

Make sure git submodules are checked out in the `third-party` directory:
```
git submodule update --init
```

Set up build directory and install dependencies
```
# For Debug binaries
conan install . --output-folder=build/debug --build=missing --settings=build_type=Debug
cmake --preset conan-debug

# For Release binaries
conan install . --output-folder=build/release --build=missing --settings=build_type=Release
cmake --preset conan-release
```

Then compile as follows:
```
cmake --build --preset conan-debug  # Or conan-release
```

This results in a program, `delphy`, that runs entirely on the command line.
It also builds a companion program, `delphy_ui`, with an OpenGL/GLUT-based visualization
of the tree inference process.

# Building for WebAssembly

NOTE: On macos systems with Apple Silicon (e.g., M1), conan will mistakenly try to
download Emscripten precompiled for x86-64
([issue](https://github.com/conan-io/conan-center-index/issues/27610)).  Instead, you need
to ask conan to rebuild emscripten locally by adding the option `'--build=emsdk/*'` (with
the single quotes) in both of the `conan install` commands below.

First set up the build directory and install dependencies using the `wasm.profile`:
```
# For Debug WASM blob
conan install . --profile=wasm.profile --output-folder=build/wasm-debug --build=missing --settings=build_type=Debug
cmake --preset conan-emscripten-debug

# For Release WASM blob
conan install . --profile=wasm.profile --output-folder=build/wasm-release --build=missing --settings=build_type=Release
cmake --preset conan-emscripten-release
```

Then compile as follows:
```
cmake --build --preset conan-emscripten-debug     # Or conan-emscripten-release
```

# Demo data

After succesful compilation, you should be able to run the command
```
./build/release/delphy --version   # or ./build/debug/delphy
```
and obtain a version string like the following:
```
Delphy Version 0.996 (build 2022, commit 33c06a8)
```

If that works, you can try to run a short inference to verify that everything is working.  We have prepared [several benchmarks](https://github.com/broadinstitute/delphy-2024-paper-data), any of which can be used for a test run after compilation.  In particular, the [SARS-CoV-2 benchmark there](https://github.com/broadinstitute/delphy-2024-paper-data/tree/main/sars-cov-2-lemieux/delphy_inputs/ma_sars_cov_2.fasta) is what is used in the `delphy-web` demo (also available [here](https://delphy.fathom.info/ma_sars_cov_2.fasta)).  The following command will execute a (too short) MCMC run using this data:
```
./build/release/delphy --v0-in-fasta ma_sars_cov_2.fasta \
   --v0-init-heuristic \
   --v0-steps 10000000 \
   --v0-out-log-file delphy.log \
   --v0-log-every 200000 \
   --v0-out-trees-file delphy.trees \
   --v0-tree-every 2000000
```
This should produces BEAST2-like `.log` and `.trees` files which can be analyzed accordingly.

The [benchmarks](https://github.com/broadinstitute/delphy-2024-paper-data) provide several complete, documented example of running Delphy end-to-end from the command line, including expected output.  See the files `07_plot_ess.py` in each benchmark for typical runtimes of each benchmark.

Run `delphy --help` for detailed command-line options.  Note: the CLI interface is not finalized and will change in the future.  The `--v0-` prefixes are there to emphasize this unstable nature.

You may also find it useful to run our development UI, where many internal parameters can be controlled interactively at a lower level than in the web interface.  The controls are undocumented and subject to change, but can be gleaned from studying the code for `keyboard_func` in [delphy_ui.cpp](tools/delphy_ui.cpp).