function baseline(options::Options, trials::Vector{<:Trial}; kfold::Integer=5)
    𝐗 = read(matopen(options.baseline_spikecounts_path))["X"]
    𝐲 = collect(mean(trial.y) for trial in trials)./options.dt
    trialindices = collect(trial.trialindex for trial in trials)
    𝐗 = hcat(ones(length(trials)), 𝐗[trialindices,:])
    testindices, trainindices = SPGLM.cvpartition(kfold, length(trialindices))
    𝛌 = 10.0.^collect(range(log10(options.baseline_L2min), log10(options.baseline_L2max), length=options.baseline_L2n))
    mse = collect(SPGLM.gaussianmse(λ, testindices, trainindices, 𝐗, 𝐲) for λ in 𝛌)
    minindex = findmin(mse)[2]
    if (minindex==1) || (minindex == options.baseline_L2n)
        error("The range of L2 penalties is too limited. ")
    end
    λₒₚₜ = 𝛌[minindex]
    # λₒₚₜ = 1e6 # for testing
    𝐗*((𝐗'*𝐗 + λₒₚₜ*I) \ (𝐗'*𝐲))
end

"""
    gaussianmse(λ. testindices, trainindices, 𝐗, 𝐲)

Out-of-sample mean squared error, which proportional to the log-likelihood, of a gaussian linear model

ARGUMENT
-`λ`: L2 penalty
-`testindices`: nested vector whose each k-th element is a vector of integers indicating the indices of the samples used for testing in the k-th cross-validation fold
-`trainindices`: nested vector whose each k-th element is a vector of integers indicating the indices of the samples used for trainining in the k-th cross-validation fold
-`𝐗`: design matrix; i.e., the predictor. Each row corresponds to a sample, and each column a regressor
-`𝐲`: independent variable; i.e., the response variable.

RETURN
-a scalar value representing the out-of-sample mean-squared error
"""
function gaussianmse(λ::AbstractFloat, testindices::Vector{<:Vector{<:Integer}}, trainindices::Vector{<:Vector{<:Integer}}, 𝐗::Matrix{<:Real}, 𝐲::Vector{<:Real})
    mse = 0
    for (testindices, trainindices) in zip(testindices, trainindices)
        𝐗train = 𝐗[trainindices,:]
        𝐲train = 𝐲[trainindices]
        𝐰 = ((𝐗train'*𝐗train) + λ*I) \ (𝐗train'*𝐲train)
        𝐗test = 𝐗[testindices,:]
        𝐲test = 𝐲[testindices]
        𝛆 = (𝐲test - 𝐗test*𝐰)
        mse += (𝛆'*𝛆)*length(testindices)
    end
    mse/sum(length.(testindices))^2
end
