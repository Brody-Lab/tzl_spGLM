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

export  GLM,
        Options

include("types.jl")
include("load.jl")
end
