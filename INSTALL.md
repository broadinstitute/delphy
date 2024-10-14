# Dependencies

Some third-party libraries (e.g., abseil) are linked directly in the source tree as git submodules, whereas more complex dependencies use Conan (see [conanfile.txt] for details).

# Building locally

If necessary (e.g., Ubuntu >= 24.04), create and activate a python3 virtualenv to install packages
```
python3 -m venv delphy-venv
source delphy-venv/bin/activate
```

Install Conan 1.59.0 (we do not yet support using Conan 2.0):
```
pip3 install 'conan==1.59.0'   # See https://docs.conan.io/en/latest/installation.html for details
```

Make some config adjustment to your Conan profile (see https://docs.conan.io/en/latest/getting_started.html)
```
conan profile new default --detect
conan profile update settings.compiler.libcxx=libstdc++11 default
```

Make sure git submodules are checked out in the `third-party` directory:
```
git submodule update --init
```

Set up build directory and install dependencies
```
mkdir -p build/debug && cd build/debug
conan install ../..
```

Then compile as usual:
```
cmake ../.. -DCMAKE_BUILD_TYPE=Debug  # Or Release
make -j 6
```

This results in a program, `delphy`, that runs entirely on the command line.
It also builds a companion program, `delphy_ui`, with an OpenGL/GLUT-based visualization
of the tree inference process.

# Building for WebAssembly

First set up the build directory and install dependencies using the `wasm.profile`:
```
mkdir -p build/wasm-debug && cd build/wasm-debug
conan install ../.. -pr ../../wasm.profile
```

Then use the `emcmake` wrapper to cross-compile to WASM:
```
emcmake cmake ../.. -DCMAKE_BUILD_TYPE=Debug  # Or Release
make -j 6
```

# Demo data

After succesful compilation, you should be able to run the command
```
delphy --version
```
and obtain a version string like the following:
```
Delphy Version 0.996 (build 2022, commit 33c06a8)
```

If that works, you can try to run a short inference to verify that everything is working.  We have prepared [several benchmarks](https://github.com/broadinstitute/delphy-2024-paper-data), any of which can be used for a test run after compilation.  In particular, the [SARS-CoV-2 benchmark there](https://github.com/broadinstitute/delphy-2024-paper-data/tree/main/sars-cov-2-lemieux/delphy_inputs/ma_sars_cov_2.fasta) is what is used in the `delphy-web` demo (also available [here](https://delphy.fathom.info/ma_sars_cov_2.fasta)).  The following command will execute a (too short) MCMC run using this data:
```
delphy --v0-in-fasta ma_sars_cov_2.fasta \
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