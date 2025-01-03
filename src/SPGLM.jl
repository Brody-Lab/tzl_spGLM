module SPGLM

using   Bootstrap,
        CSV,
        DataFrames,
        Distributed,
        Distributions,
        DSP,
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
include("baseline.jl")
include("basisfunctions.jl")
include("characterization.jl")
include("crossvalidation.jl")
include("evidenceoptimization.jl")
include("inverselink.jl")
include("learnparameters.jl")
include("load.jl")
include("modelinputs.jl")
include("perievent_time_histogram.jl")
include("Poisson.jl")
include("save.jl")
include("shortcuts.jl")
include(joinpath("projects", "auditory.jl"))
end
