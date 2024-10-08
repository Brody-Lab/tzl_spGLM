"""
    crosssvalidate(kfold,options,trials)

Assess how well the factorial hidden Markov drift-diffusion model generalizes to independent datasets

ARGUMENT
-`kfold`: number of cross-validation folds
-`model`: a structure containing the settings, data, and parameters of a factorial hidden-Markov drift-diffusion model

OPTIONAL ARGUMENT
-`choicesonly`: whether to train on only the behavioral choices and ignore the spike trains

OUTPUT
-an instance of `CVResults`
"""
crossvalidate(kfold::Integer, options::Options, trials::Vector{<:Trial}) = crossvalidate(CVIndices(kfold, trials),options,trials)

"""
    crosssvalidate(cvindices,options,trials)
"""
function crossvalidate(cvindices::Vector{<:CVIndices}, options::Options, trials::Vector{<:Trial})
	trainingmodels = collect(fit(options,trials[cvindex.trainingtrials]) for cvindex in cvindices)
	testmodels = collect(test(trials[cvindex.testingtrials], trainingmodel) for (cvindex, trainingmodel) in zip(cvindices, trainingmodels))
	characterization = Characterization(cvindices, testmodels, trainingmodels)
	peths = perievent_time_histograms(testmodels[1].basissets, characterization.inferredrate, options, trials)
    CVResults(cvindices=cvindices, characterization=characterization, peths=peths, trainingmodels=trainingmodels)
end

"""
    CVIndices(kfold, trials)

Create indices for cross-validation

ARGUMENT
-`model`: a structure containing the settings, data, and parameters of a factorial hidden-Markov drift-diffusion model
-`kfold`: number of cross-validation folds

OUTPUT
-a vector of instances of `CVIndices`
"""
CVIndices(kfold::Integer, trials::Vector{<:Trial}) = CVIndices(kfold::Integer, length(trials))

function CVIndices(kfold::Integer, ntrials::Integer)
	testingtrials, trainingtrials = cvpartition(kfold,ntrials)
    map(1:kfold) do k
        CVIndices(trainingtrials = trainingtrials[k],
                  testingtrials = testingtrials[k])
    end
end

"""
    cvpartition(kfold, nsamples)

Partition samples for cross-validation

ARGUMENT
-`kfold`: number of cross-validation folds
-`nsamples`: number of samples

RETURN
-`testing`: a nested array whose element `testing[k]` is a vector of integers indexing the samples to be used for testing for the k-th cross-validation fold.
-`training`: indices of the samples to be used for training, organized similarly to `testing`
"""
function cvpartition(kfold::Integer, nsamples::Integer)
	training = collect(convert.(Int,x) for x in MLBase.Kfold(nsamples,kfold))
	testing = collect(setdiff(1:nsamples, training) for training in training)
    return testing, training
end

"""
	test(testtrials, trainingmodel)

Construct a model with only the test data and the parameters learned from the training data

ARGUMENT
-`cvindices`: indices of the trials and timesteps used for training and testing in each fold
-`model`: structure containing both the test and training data
-`trainingmodel`: structure containing only the training data and the parameters learned for those data

OUTPUT
-`testmodel`: a structure containing only the test data and the parameters learned from the training data
"""
function test(testtrials::Vector{<:Trial}, trainingmodel::Model)
	testmodel = Model(trainingmodel.options, testtrials, trainingmodel.𝐰_baseline)
	testmodel.a[1] = trainingmodel.a[1]
	testmodel.𝐰 .= trainingmodel.𝐰
	testmodel
end

"""
	Characterization(cvindices, testmodels, trainingmodels)

Out-of-sample computation of quantities for characterizing the model

ARGUMENT
-`cvindices`: a vector of composites, each of which containing the indices of trials used for testing and training for a particular cross-validation fold
-`testmodels`: a vector of composites, each of which containing the held-out data for a cross-validation fold
"""
function Characterization(cvindices::Vector{<:CVIndices}, testmodels::Vector{<:Model}, trainingmodels::Vector{<:Model})
	characterization_each_fold = collect(Characterization(testmodel,trainingmodel) for (testmodel,trainingmodel) in zip(testmodels,trainingmodels))
	trialindices = sortperm(vcat((cvindex.testingtrials for cvindex in cvindices)...))
	values =
		map(fieldnames(Characterization)) do fieldname
				vcat((getfield(characterization, fieldname) for characterization in characterization_each_fold)...)[trialindices]
		end
	Characterization(values...)
end
