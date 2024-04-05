module SPGLM

using   Bootstrap,
        Distributed,
        Distributions,
        DSP,
        ForwardDiff,
        LinearAlgebra,
        LineSearches,
        MAT,
        Optim,
        Parameters,
        Random,
        SpecialFunctions,
        StatsFuns
import  CSV,
        DataFrames,
        MLBase
include("types.jl")
include("crossvalidation.jl")
include("evidenceoptimization.jl")
include("gaussianprior.jl")
include("loadmodel.jl")
include("maximumlikelihood.jl")
include("sampling.jl")
include("save.jl")
include("tests.jl")
