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

import  StatsBase

export  GLM,
        Options

include("types.jl")
include("load.jl")
include("temporal_basis_functions.jl")
end
