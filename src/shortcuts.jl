"""
    fit(cvspath, csvrow)

Fit the model and save the results

ARGUMENT
-`cvspath`: absolute path to the comma-separated values (CSV) file containing information about the path of the data and output and the values of hyperparameters
-`csvrow`: the row in the CSV to consider
"""
function fit(csvpath::String, csvrow::Integer; save_characterization::Bool=false)
    options = Options(csvpath,csvrow)
    println(options.datapath)
    trials = loadtrials(options)
    model = Model(options, trials)
    if options.opt_method == "evidenceoptimization"
		eo = maximizeevidence!(model)
	    save(eo, model.options.outputpath)
	elseif options.opt_method == "gridsearch"
		gridsearch!(model)
	elseif options.opt_method == "maximumaposteriori"
		maximizeposterior!(model)
	end
    save(model)
    characterization = Characterization(model)
    peths = perievent_time_histograms(characterization.inferredrate,model)
    save(peths, model.options.outputpath)
    save_characterization && save(characterization, options.outputpath)
end

"""
    crossvalidate(cvspath, csvrow)

Fit the model under a cross-validation scheme and save the results

ARGUMENT
-`cvspath`: absolute path to the comma-separated values (CSV) file containing information about the path of the data and output and the values of hyperparameters
-`csvrow`: the row in the CSV to consider
-`kfold`: number of cross-validation folds
"""
function crossvalidate(csvpath::String, csvrow::Integer, kfold::Integer; save_characterization::Bool=false)
    options = Options(csvpath,csvrow)
    println(options.datapath)
    trials = loadtrials(options)
    cvresults = crossvalidate(kfold,options,trials)
    matwrite(joinpath(options.outputpath, "cvindices.mat"), Dict("cvindices"=>map(dictionary,cvresults.cvindices)))
    matwrite(joinpath(options.outputpath, "trainingmodels.mat"), Dict("models"=>map(dictionary, cvresults.trainingmodels)))
    save_characterization && save(cvresults.characterization, options.outputpath)
    save(cvresults.peths, options.outputpath)
    save(cvresults.characterization, trials, options.outputpath)
end
