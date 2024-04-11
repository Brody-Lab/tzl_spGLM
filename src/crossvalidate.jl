"""
    crosssvalidate(model)

Assess how well the factorial hidden Markov drift-diffusion model generalizes to independent datasets

ARGUMENT
-`kfold`: number of cross-validation folds
-`model`: a structure containing the settings, data, and parameters of a factorial hidden-Markov drift-diffusion model

OPTIONAL ARGUMENT
-`choicesonly`: whether to train on only the behavioral choices and ignore the spike trains

OUTPUT
-an instance of `CVResults`
"""
function crossvalidate(kfold::Integer, options::Options, trials::Vector{<:Trial})
    cvindices = CVIndices(kfold, trials)
	trainingmodels = pmap(cvindices->train(cvindices, options, trials), cvindices)
	trainingsummaries = collect(ModelSummary(trainingmodel) for trainingmodel in trainingmodels)
	testmodels = collect(test(cvindex, model, trainingmodel) for (cvindex, trainingmodel) in zip(cvindices, trainingmodels))
	characterization = Characterization(cvindices, testmodels, trainingmodels)
	psthsets = poststereoclick_time_histogram_sets(characterization.expectedemissions, model)
    CVResults(cvindices=cvindices, characterization=characterization, psthsets=psthsets, trainingsummaries=trainingsummaries)
end
