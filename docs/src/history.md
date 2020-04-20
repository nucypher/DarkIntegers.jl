# Version history


## Current development version

* ADDED: `>>` and `<<` methods for `ModUInt`.


## v0.1.2 (19 February 2020)

* FIXED: a typo in `convert(::Type{MgModUInt{T, M}}, ::ModUInt{T, M})`
* FIXED: `known_isprime()` is now called even if it is defined in the user program.
* ADDED: `known_polynomial_mul_function()` to override the default choice.


## v0.1.1 (20 December 2019)

* FIXED: a bug in conversion between polynomials of different lengths
* FIXED: added a trivial `promote()` method to fix some promotion errors


## v0.1.0 (19 December 2019)

The first semi-stable release.


## v0.0.1

Initial version.
