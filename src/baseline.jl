"""
    baselineweights(options,trials;kfold)

RETURN the weights of the pre-trial firing rates

ARGUMENT
-`options`: fixed hyperparameters
-`trials`: data

OPTIONAL ARGUMENT
-`kfold`: number of cross-validation fold
"""
function baselineweights(options::Options, trials::Vector{<:Trial})
    if !isempty(options.baseline_pretrial_spikecounts_path) || options.baseline_time_in_session
        𝐗 = baselineinput(options,trials)
        𝐲 = collect(mean(trial.y) for trial in trials)./options.dt
        if !isempty(options.baseline_pretrial_spikecounts_path)
            𝛌 = 10.0.^collect(range(log10(options.baseline_L2min), log10(options.baseline_L2max), length=options.baseline_L2n))
        else
            𝛌 = 10.0.^collect(-9:-1)
        end
        fit_linear_gaussian(𝐗,𝐲,𝛌)
    else
        zeros(0)
    end
end

"""
    fit_linear_gaussian(𝐗,𝐲,𝛌;kfold)

RETURN the parameters of a linear gaussian model

ARGUMENT
-`𝐗`: design matrix
-`𝐲`: response vector
-`𝛌`: vector of L2 regularization coefficients to try

OPTIONAL ARGUMENT
-`kfold`: number of cross-validation folds
"""
function fit_linear_gaussian(𝐗::Matrix{<:AbstractFloat}, 𝐲::Vector{<:AbstractFloat}, 𝛌::Vector{<:AbstractFloat}; kfold::Integer=5)
    testindices, trainindices = cvpartition(kfold, length(𝐲))
    mse = collect(gaussianmse(λ, testindices, trainindices, 𝐗, 𝐲) for λ in 𝛌)
    minindex = findmin(mse)[2]
    if (minindex==1) || (minindex == length(𝛌))
        println("minindex=", minindex, ": The range of L2 penalties is too limited.")
    end
    λₒₚₜ = 𝛌[minindex]
    println(λₒₚₜ)
    (𝐗'*𝐗 + λₒₚₜ*I) \ (𝐗'*𝐲)
end

"""
    baselineinput(options, trials)

RETURN a matrix parametrizing the baseline inputs
"""
function baselineinput(options::Options, trials::Vector{<:Trial})
    if !isempty(options.baseline_pretrial_spikecounts_path)
        baseline_spikecounts(options,trials)
    else
        baseline_time_in_session(options,trials)
    end
end

"""
    baseline_spikecounts(options,trials)

RETURN a matrix containing pre-trial spike counts for estimating baseline input
"""
function baseline_spikecounts(options::Options, trials::Vector{<:Trial})
    trialindices = collect(trial.trialindex for trial in trials)
    file = matopen(options.baseline_pretrial_spikecounts_path)
    𝐗 = hcat(ones(length(trials)), read(file,"X")[trialindices,:])
    close(file)
    return 𝐗
end

"""
    baseline_time_in_session(options, trials)

RETURN a matrix containing basis functions of time in session for estimating baseline input
"""
function baseline_time_in_session(options::Options, trials::Vector{<:Trial})
    N = ceil(Int,trials[1].last_reference_time_s-trials[1].first_reference_time_s)+1
    Φ = basisfunctions(N, options.baseline_time_in_session_D; begins0=false, ends0=false, η=0.0)
    reference_times_s = collect(ceil(Int, trial.reference_time_s-trial.first_reference_time_s)+1 for trial in trials)
    Φ[reference_times_s,:]
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
