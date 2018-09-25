# Dark Integers

`DarkIntegers.jl` is a lightweight library for working with unsigned integers and polynomials. It includes:

* Modulo arithmetic for primitive types;
* Multi-precision unsigned integers;
* Modulo arithmetic for multi-precision integers, with several representations optimized for different operations.
* (to be added) Modulo polynomials;
* (to be added) NTT.

It is intended as a partial replacement for [`Nemo.jl`](https://github.com/wbhart/Nemo.jl), which requires compilation of C libraries, and is often slower (in particular, because `DarkIntegers` does not need C calls and heap allocations for operations with multi-precision integers).
