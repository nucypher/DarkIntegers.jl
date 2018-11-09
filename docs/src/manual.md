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


## Multi-precision numbers

The basic facility in `DarkIntegers` are multi-precision unsigned integers [`MPNumber{N, T}`](@ref MPNumber).
The type takes two parameters, `N` (an integer) being the number of limbs, or digits, and `T` (an unsigned integer type) is the type of a limb.
Its functionality is closer to that of `BigInt`, but [`MPNumber`](@ref) still has fixed length (although it can be much greater than what any Julia built-in type supports).
On the other hand, operations with it do not require dynamic memory allocation and are faster than the `BigInt` ones.

The [`MPNumber`](@ref) constructor takes any integer:
```@jldoctest mpnumber-basics
a = MPNumber{2, UInt8}(65534)
println(a)

# output

{(0xfe, 0xff)}
```
Note that the limbs are arranged in the big-endian order.

A multi-precision number can be converted back to any integer:
```@jldoctest mpnumber-basics
b = convert(Int, a)
println(b)

# output

65534
```

Multi-precision numbers support arithmetic operations and comparisons and act identically to built-in unsigned integer types:
```@jldoctest mpnumber-arithmetic
a = MPNumber{2, UInt8}(65534)
b = MPNumber{2, UInt8}(65533)
println(a + b)
println(convert(Int, a + b))

# output
{(0xfb, 0xff)}
65531
```


## Residue ring elements

The next level above multi-precision integers (and built-in unsigned integers) are residue ring elements, that is, unsigned integers with the arithmetic operations performed modulo some number.
The residue ring element type [`RRElem{T, M}`](@ref RRElem) is parametrized by the number type `T` (an unsigned integer) and the modulus `M` (which is a *value* of type `T`).

Similarly to [`MPNumber`](@ref) objects, [`RRElem`](@ref) objects can be constructed out of an integer.
If the integer is greater than the modulus, it will be truncated (by applying `mod()`):
```jldoctest rrelem-basics
modulus = UInt8(200)
a = RRElem{UInt8, modulus}(250)
println(a)

# output
50RR
```

All arithmetic operations on [`RRElem`](@ref) objects are performed modulo `modulus`.
Any regular integers in mixed expressions are promoted to [`RRElem`](@ref) objects:
```jldoctest rrelem-basics
println(a + 101)

# output
151RR
```

[`RRElem`](@ref) objects can be converted back to integers:
```jldoctest rrelem-basics
println(convert(Int, a))

# output
50
```


## Montgomery representation

Residue ring elements can be converted to an alternate representation, which makes use of [Montgomery reduction](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication) for multiplication.
As a result, the multiplication becomes much faster, addition and subtraction stay the same, and division (and conversion to and from integers) become slower.
The representation implemented as the type [`RRElemMontgomery{T, M}`](@ref RRElemMontgomery), which is parametrized in the same way as [`RRElem`](@ref).
Depending on the relative amount of different arithmetic operations one needs to perform, either [`RRElem`](@ref) or [`RRElemMontgomery`](@ref) may perform better.

The interface is the same as the one for [`RRElem`](@ref), except for the restriction on `M` to be an odd number:
```jldoctest rrelemmontgomery-basics
modulus = UInt8(201)
a = RRElemMontgomery{UInt8, modulus}(250)
println(a)
println(convert(Int, a + 101))

# output
49RRM
150
```


## Cyclic polynomials

Anything type supporting arithmetic operations (including [`MPNumber`](@ref), [`RRElem`](@ref) and [`RRElemMontgomery`](@ref)) can serve as the coefficient type in the [`Polynomial`](@ref) type.
`DarkIntegers` supports cyclic polynomials (with operations performed modulo `x^n-1`, where `n` is some non-negative integer called the length of the polynomial) and negacyclic ones (with operations performed modulo `x^n+1`).

Polynomials are created out of a coefficient array (where the `i`-th element corresponds to the `(i-1)`-th power of `x`) and the negacyclicity flag:
```jldoctest polynomial-basics
p = Polynomial([1, 2, 3, 4], true) # creates a negacyclic polynomial
println(p)
println(p + 1)
println(p * 2)

# output
Polynomial{Int64}([1, 2, 3, 4], true, DarkIntegers.karatsuba_mul)
Polynomial{Int64}([2, 2, 3, 4], true, DarkIntegers.karatsuba_mul)
Polynomial{Int64}([2, 4, 6, 8], true, DarkIntegers.karatsuba_mul)
```

The polynomial can be multiplied by a power of `x` using [`shift_polynomial`](@ref).
Since the polynomial we created is cyclic, the coefficients with the powers greater or equal to `n` will reappear from the other side with the opposite sign (one can work it out by applying the `x^n+1` modulus manually):
```jldoctest polynomial-basics
println(shift_polynomial(p, 2))

# output
Polynomial{Int64}([-3, -4, 1, 2], true, DarkIntegers.karatsuba_mul)
```

Note the multiplication function that is the part of the [`Polynomial`](@ref) structure.
The default for the multiplication is Karatsuba algorithm; if possible a more faster NTT-based algorithm will be chosen.
It requires:
* the coefficient type to be a residue ring element ([`RRElem`](@ref) or [`RRElemMontgomery`](@ref)) with a prime modulus;
* the length of the polynomial to be a power of 2;
* the length of the polynomial to be a factor of `(modulus - 1)` for cyclic polynomials, and `(modulus - 1)/2` for negacyclic ones.
For example:
```jldoctest polynomial-mul
modulus = UInt8(241)
tp = RRElem{UInt8, modulus}
p1 = Polynomial(tp[1, 2, 3, 4], true)
p2 = Polynomial(tp[1, 0, 1, 0], true)
println(p1)
println(p2)
println(p1 * p2)

# output
Polynomial{RRElem{UInt8,0xf1}}(RRElem{UInt8,0xf1}[1RR, 2RR, 3RR, 4RR], true, DarkIntegers.ntt_mul)
Polynomial{RRElem{UInt8,0xf1}}(RRElem{UInt8,0xf1}[1RR, 0RR, 1RR, 0RR], true, DarkIntegers.ntt_mul)
Polynomial{RRElem{UInt8,0xf1}}(RRElem{UInt8,0xf1}[239RR, 239RR, 4RR, 6RR], true, DarkIntegers.ntt_mul)
```
