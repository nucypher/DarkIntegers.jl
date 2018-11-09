# Dark Integers, an unsigned integer arithmetic toolbox

`DarkIntegers.jl` is a pure Julia implementation of unsigned integer modulo arithmetic for simple and multi-precision (arbitrary long) integers. Functionally, it is a subset of [Nemo](http://nemocas.org), but unlike it, `DarkIntegers` does not use external C libraries. In particular, it implements Montgomery representation of integers, allowing for efficient modulo multiplication. Also it has a slightly simpler interface, and does not show any messages on import.
