# Building locally

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
