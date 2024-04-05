module SPGLM

using   CSV,
        DataFrames,
        Distributions,
        ForwardDiff,
        LinearAlgebra,
        LineSearches,
        MAT,
        MLBase,
        Optim,
        Parameters,
        Random,
        SpecialFunctions,
        StatsFuns
include("types.jl")
include("load.jl")
end
