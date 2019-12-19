# API reference

```@meta
CurrentModule = DarkIntegers
```

## Utility functions

```@docs
bitsizeof
log_bitsizeof
num_bits
encompassing_type
as_builtin
```


## Single-limb arithmetic

```@docs
addhilo
mulhilo
divremhilo
divhilo
remhilo
```


## Single-limb modulo arithmetic

```@docs
addmod
submod
mulmod
powmod
```


## Multi-precision numbers

```@docs
MLUInt
MLInt
```


## Modulo integers

```@docs
ModUInt
```


## Modulo integers (Montgomery representation)

```@docs
MgModUInt
```


## Modulo integers helper functions

```@docs
value
raw_value
modulus
```


## (Nega)cyclic polynomials

```@docs
Polynomial
mul_by_monomial
cyclic_modulus
negacyclic_modulus
known_isprime
broadcast_into_polynomial
broadcast_into_polynomial!
with_modulus
resize
```

## NTT

```@docs
ntt
known_generator
```
