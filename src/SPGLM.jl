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
include("basisfunctions.jl")
include("characterization.jl")
include("evidenceoptimization.jl")
include("inverselink.jl")
include("learnparameters.jl")
include("load.jl")
include("modelinputs.jl")
include("Poisson.jl")
end
