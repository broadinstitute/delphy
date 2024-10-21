include(default)

[settings]
os=Emscripten
arch=wasm
compiler=clang
compiler.version=18
compiler.libcxx=libc++
compiler.cppstd=20

[tool_requires]
emsdk/3.1.50

[conf]
tools.cmake.cmake_layout:build_folder_vars=['settings.os','settings.build_type']