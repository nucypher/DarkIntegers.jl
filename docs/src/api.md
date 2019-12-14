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
ModUInt
```


## Residue ring elements (Montgomery representation)

```@docs
MgModUInt
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
