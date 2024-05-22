Coding Style
============

The code mostly adheres to the [Google C++ Style guide](https://google.github.io/styleguide/cppguide.html), with a few
key exceptions:

* Compile against C++20.  Since we support both gcc and clang (Emscripten), some esoteric features of C++20 still
  can't be used yet.  We don't use modules (yet) because of poor tooling support, but we'd like to.

* Use exceptions instead of error codes or rogue return values.  It's ok to just crash when there's a problem.

* Use `snake_case` for names.  Type names (classes, class templates, aliases) use `Stroustrup_case`, so that we can
  immediately recognize types by the uppercase first letter.  Note: some variables use an uppercase letter when it's a
  useful math symbol or initials, e.g., `X`, `Ttwiddle`.

* Use (Almost-Always-Auto (AAA)
  style)[https://herbsutter.com/2013/08/12/gotw-94-solution-aaa-style-almost-always-auto/] for variable and function
  declarations, e.g., `auto double(int x) -> int;` and `auto m = node.mutations();`

* [TO ENFORCE] Prefer `not` to `!` in Boolean expressions.

* [TO ENFORCE] We use `log` to mean "log base e", never `ln`

* [TO ENFORCE] We use west-coast const style, and pointer and reference indicators (`*` and `&`) go with the types, not
  the names, e.g., `const char* str` instead of `char const *str`.

* [TO ENFORCE] Use std::ssize(c) (signed) instead of c.size() (unsigned) to avoid pitfalls of mixing signed and unsigned
  quantities in arithmetic

* [TO ENFORCE] Prefer absl::flat_hash_map and absl::flat_hash_set to std::unordered_map and stl::unordered_set

* [TO ENFORCE] Prefer structs with public fields to classes with private fields and accessors (see [C++ Core Guideline C.131](http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rh-get)).