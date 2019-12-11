# API reference

```@meta
CurrentModule = DarkIntegers
```

## Utility functions

```@docs
bitsizeof
num_bits
encompassing_type
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
```


## Multi-precision numbers

```@docs
MLUint
```


## Residue ring elements

```@docs
RRElem
```


## Residue ring elements (Montgomery representation)

```@docs
RRElemMontgomery
```


## Cyclic polynomials

```@docs
Polynomial
shift_polynomial
```

## Partial modification of residue ring elements and polynomials

```@docs
change_representation
change_base_type
change_modulus
change_length
rescale
```
