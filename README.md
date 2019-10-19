# Dark Integers

Master branch: [![CircleCI](https://circleci.com/gh/nucypher/DarkIntegers.jl/tree/master.svg?style=svg)](https://circleci.com/gh/nucypher/DarkIntegers.jl/tree/master) [![codecov](https://codecov.io/gh/nucypher/DarkIntegers.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nucypher/DarkIntegers.jl)

`DarkIntegers.jl` is a lightweight library for working with unsigned integers and polynomials. It includes:

* Modulo arithmetic for primitive types;
* Regular and modulo arithmetic for multi-limb (fixed size) unsigned integers;
* Arithmetic for cyclic and negacyclic polynomials, including Karatsuba multiplication and multiplication via NTT

It is intended as a partial replacement for [`Nemo.jl`](https://github.com/wbhart/Nemo.jl), which requires compilation of C libraries, and is often slower (in particular, because `DarkIntegers` does not need C calls and heap allocations for operations with multi-limb integers).

**Note:** the library currently does not provide constant-time guarantees for arithmetic operations.
