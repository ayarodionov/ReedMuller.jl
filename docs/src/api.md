```@meta
CurrentModule = ReedMuller
```

# API Reference

```@docs
ReedMuller
```

## Codes

```@docs
RMCode
blocklength
dimension
minimum_distance
rate
generator_matrix
monomials
```

## Common interface

```@docs
AbstractEncoder
AbstractDecoder
encode
decode
basis
hard_llr
```

## Encoders

```@docs
MatrixEncoder
PlotkinEncoder
```

## Decoders

```@docs
ReedDecoder
FHTDecoder
DumerDecoder
DumerShabunovDecoder
SidelnikovPershakovDecoder
RPADecoder
BPDecoder
GLPDecoder
glp_permutations
MLDecoder
```

## Generic decoder wrappers

```@docs
AutomorphismEnsembleDecoder
ChaseDecoder
GMDDecoder
```

## Channels

```@docs
BSC
BIAWGN
BIAWGN_from_ebn0
transmit
```

## Simulation

```@docs
simulate
SimResult
```
