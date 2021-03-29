# LiRA

Linear Response of (galactic) Arms

## Aim of the code

This code aims to compute the growth rate of the arms of a thin Kuzmin-Toomre galactic disc, influence by the presence of two Toomre sphere potential modeling a central bulge and a surrounding dark matter halo.

## Modeling

Following Aoki et al (1979), we consider thin disk made of polytropic gas with index 4/3, which gives analytical expression for the response matrix elements.
The implementation of the bulge component and of the self-gravity parameter slightly modify those matrix elements, which we implement in this semi-analytical code.

## Julia packages

Open the terminal in the folder `packages` and type

```
$ julia Install-pkg.jl
```

to install the following packages:

- `HDF5`
- `ArgParse`
- `Plots`
- `LinearAlgebra`
- `StaticArrays`
- `SphericalHarmonics`

### !! WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**
