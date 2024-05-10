"""
    fit(cvspath, csvrow)

Fit the model and save the results

ARGUMENT
-`cvspath`: absolute path to the comma-separated values (CSV) file containing information about the path of the data and output and the values of hyperparameters
-`csvrow`: the row in the CSV to consider
"""
function fit(csvpath::String, csvrow::Integer)
    options = Options(csvpath,csvrow)
    println(options.datapath)
    trials = loadtrials(options)
    model = Model(options, trials)
    eo = maximizeevidence!(model)
    characterization = Characterization(model)
    peths = perievent_time_histograms(characterization.inferredrate,model)
    save(model)
    save(eo, model.options.outputpath)
    save(characterization, model.options.outputpath)
    save(peths, model.options.outputpath)
end

"""
    crossvalidate(cvspath, csvrow)

Fit the model under a cross-validation scheme and save the results

ARGUMENT
-`cvspath`: absolute path to the comma-separated values (CSV) file containing information about the path of the data and output and the values of hyperparameters
-`csvrow`: the row in the CSV to consider
-`kfold`: number of cross-validation folds
"""
function crossvalidate(csvpath::String, csvrow::Integer, kfold::Integer)
    options = Options(csvpath,csvrow)
    println(options.datapath)
    trials = loadtrials(options)
    cvresults = crossvalidate(kfold,options,trials)
    matwrite(joinpath(options.outputpath, "cvindices.mat"), Dict("cvindices"=>map(dictionary,cvresults.cvindices)))
    matwrite(joinpath(options.outputpath, "trainingmodels.mat"), Dict("models"=>map(dictionary, cvresults.trainingmodels)))
    matwrite(joinpath(options.outputpath, "evidenceoptimizations.mat"), Dict("evidenceoptimizations"=>map(dictionary, cvresults.evidenceoptimizations)))
    save(cvresults.characterization, options.outputpath)
    save(cvresults.peths, options.outputpath)
    nats_per_spike = sum(sum.(cvresults.characterization.LL))/sum(sum.(cvresults.characterization.observed_spiketrains))
    matwrite(joinpath(options.outputpath, "nats_per_spike.mat"), Dict("nats_per_spike"=>nats_per_spike))
end
