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

include("types.jl")
include("load.jl")
include("basisfunctions.jl")
include("modelinputs.jl")
end
