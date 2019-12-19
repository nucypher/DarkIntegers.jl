# Manual

In this section we assume that every example is preceded with
```julia
using DarkIntegers
```

```@meta
DocTestSetup = quote
    using DarkIntegers
end
```


## Multi-limb integers

The basic facility in `DarkIntegers` are multi-limb unsigned integers [`MLUInt{N, T}`](@ref MLUInt).
The type takes two parameters, `N` (an integer) being the number of limbs, or digits, and `T` (an unsigned integer type) is the type of a limb.
Its functionality is closer to that of `BigInt`, but [`MLUInt`](@ref) still has fixed length (although it can be much greater than what any Julia built-in type supports).
On the other hand, operations with it do not require dynamic memory allocation and are faster than the `BigInt` ones.

The [`MLUInt`](@ref) constructor takes any integer:
```@jldoctest mpnumber-basics
a = MLUInt{2, UInt8}(65534)
println(a)

# output

{(0xfe, 0xff)}
```
Note that the limbs are arranged in the big-endian order.

A multi-limb integer can be converted back to any integer:
```@jldoctest mpnumber-basics
b = convert(Int, a)
println(b)

# output

65534
```

Multi-precision numbers support arithmetic operations and comparisons and act identically to built-in unsigned integer types:
```@jldoctest mpnumber-arithmetic
a = MLUInt{2, UInt8}(65534)
b = MLUInt{2, UInt8}(65533)
println(a + b)
println(convert(Int, a + b))

# output
{(0xfb, 0xff)}
65531
```


## Modulo integers

The next level above multi-limb integers (and built-in unsigned integers) are unsigned integers with arithmetic operations performed modulo some number.
The object of type [`ModUInt{T, M}`](@ref ModUInt) is parametrized by the number type `T` (an unsigned integer) and the modulus `M` (which is a *value* of type `T`).

Similarly to [`MLUInt`](@ref) objects, [`ModUInt`](@ref) objects can be constructed out of an integer.
If the integer is greater than the modulus, it will be truncated (by applying `mod()`):
```jldoctest moduint-basics
modulus = UInt8(200)
a = ModUInt{UInt8, modulus}(250)
println(a)

# output
50RR
```

All arithmetic operations on [`ModUInt`](@ref) objects are performed modulo `modulus`.
Any regular integers in mixed expressions are promoted to [`ModUInt`](@ref) objects:
```jldoctest moduint-basics
println(a + 101)

# output
151RR
```

[`ModUInt`](@ref) objects can be converted back to integers after extracting its value with [`value`](@ref):
```jldoctest moduint-basics
println(convert(Int, value(a)))

# output
50
```


## Montgomery representation

Modulo integers can be converted to an alternate representation, which makes use of [Montgomery reduction](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication) for multiplication.
As a result, the multiplication becomes much faster, addition and subtraction stay the same, and division (and conversion to and from integers) become slower.
The representation implemented as the type [`MgModUInt{T, M}`](@ref MgModUInt), which is parametrized in the same way as [`ModUInt`](@ref).
Depending on the relative amount of different arithmetic operations one needs to perform, either [`ModUInt`](@ref) or [`MgModUInt`](@ref) may perform better.

The interface is the same as the one for [`ModUInt`](@ref), except for the restriction on `M` to be an odd number:
```jldoctest mgmoduint-basics
modulus = UInt8(201)
a = MgModUInt{UInt8, modulus}(250)
println(a)
println(convert(Int, value(a + 101)))

# output
49RRM
150
```


## Cyclic polynomials

Anything type supporting arithmetic operations (including [`MLUInt`](@ref), [`ModUInt`](@ref) and [`MgModUInt`](@ref)) can serve as the coefficient type in the [`Polynomial`](@ref) type.
`DarkIntegers` supports cyclic polynomials (with operations performed modulo `x^n-1`, where `n` is some non-negative integer called the length of the polynomial) and negacyclic ones (with operations performed modulo `x^n+1`).

Polynomials are created out of a coefficient array (where the `i`-th element corresponds to the `(i-1)`-th power of `x`) and negacyclic modulus (`x^N+1`):
```jldoctest polynomial-basics
p = Polynomial([1, 2, 3, 4], negacyclic_modulus) # creates a negacyclic polynomial
println(p)
println(p + 1)
println(p * 2)

# output
Polynomial{Int64,4}([1, 2, 3, 4], DarkIntegers.NegacyclicModulus(), DarkIntegers.karatsuba_mul, false)
Polynomial{Int64,4}([2, 2, 3, 4], DarkIntegers.NegacyclicModulus(), DarkIntegers.karatsuba_mul, false)
Polynomial{Int64,4}([2, 4, 6, 8], DarkIntegers.NegacyclicModulus(), DarkIntegers.karatsuba_mul, false)
```

The polynomial can be multiplied by a power of `x` using [`mul_by_monomial`](@ref).
Since the polynomial we created is cyclic, the coefficients with the powers greater or equal to `n` will reappear from the other side with the opposite sign (one can work it out by applying the `x^n+1` modulus manually):
```jldoctest polynomial-basics
println(mul_by_monomial(p, 2))

# output
Polynomial{Int64,4}([-3, -4, 1, 2], DarkIntegers.NegacyclicModulus(), DarkIntegers.karatsuba_mul, false)
```

Note the multiplication function that is the part of the [`Polynomial`](@ref) structure.
The default for the multiplication is Karatsuba algorithm; if possible a more faster NTT-based algorithm will be chosen.
It requires:
* the coefficient type to be a modulo integer ([`ModUInt`](@ref) or [`MgModUInt`](@ref)) with a prime modulus;
* the length of the polynomial to be a power of 2;
* the length of the polynomial to be a factor of `(modulus - 1)` for cyclic polynomials, and `(modulus - 1)/2` for negacyclic ones.
For example:
```jldoctest polynomial-mul
modulus = UInt8(241)
tp = ModUInt{UInt8, modulus}
p1 = Polynomial(tp[1, 2, 3, 4], negacyclic_modulus)
p2 = Polynomial(tp[1, 0, 1, 0], negacyclic_modulus)
println(p1)
println(p2)
println(p1 * p2)

# output
Polynomial{ModUInt{UInt8,0xf1},4}(ModUInt{UInt8,0xf1}[1RR, 2RR, 3RR, 4RR], DarkIntegers.NegacyclicModulus(), DarkIntegers.ntt_mul, false)
Polynomial{ModUInt{UInt8,0xf1},4}(ModUInt{UInt8,0xf1}[1RR, 0RR, 1RR, 0RR], DarkIntegers.NegacyclicModulus(), DarkIntegers.ntt_mul, false)
Polynomial{ModUInt{UInt8,0xf1},4}(ModUInt{UInt8,0xf1}[239RR, 239RR, 4RR, 6RR], DarkIntegers.NegacyclicModulus(), DarkIntegers.ntt_mul, false)
```
