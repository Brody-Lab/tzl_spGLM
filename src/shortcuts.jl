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

"""
	crossvalidate(settingspath, datapath, testingtrials, trainingtrials)

ARGUMENT
-`settingspath`: a string indicating the absolute to the path of the comma-separated values (CSV) file specifying the fixed hyperparameters (i.e., options) of the generalized linear model. Each column is an option.
-`datapath`: a string specifying the absolute path of the `MAT` file
-`testingtrials`: nested array whose element `testingtrials[k][i]` is an integer that specifies the index of the i-th trial among all trials used for testing in the k-th cross-validation fold
-`trainingtrials`: nested array whose element `trainingtrials[k][i]` is an integer that specifies the index of the i-th trial among all trials used for training in the k-th cross-validation fold

RETURN
-`testinginput`: a nested array whose element `testinginput[k][i][t]` corresponds to the k-th cross-validation fold, i-th trial used for testing in that fold, and t-th time step in that trial
-`traininginput`: a nested array whose element `traininginput[k][i][t]` corresponds to the k-th cross-validation fold, i-th trial used for training in that fold, and t-th time step in that trial

EXAMPLES
```python
import numpy as np
from sklearn.model_selection import KFold
from juliacall import Main as jl
jl.seval("using SPGLM")
kf = KFold(n_splits=5,shuffle=True)
X = np.random.rand(500,2)
trainingtrials = list((output[0]+1 for output in kf.split(X)))
testingtrials = list((output[1]+1 for output in kf.split(X)))
settingspath = "/mnt/cup/people/zhihaol/Documents/tzluo/analyses/analysis_2024_06_09a_pythoncall/options.csv"
datapath = "/mnt/cup/labs/brody/tzluo/analysis_data/analysis_2024_04_06a_separate_Cells/dataset/T176_2018_05_03_619040938_003.mat"
testinginputs, traininginputs = jl.SPGLM.crossvalidate(settingspath, datapath, testingtrials, trainingtrials)
```
"""
function crossvalidate(settingspath::String, datapath::String, testingtrials, trainingtrials)
	dict = dictionary(settingspath,1)
	dict["datapath"] = datapath
	options = Options(dict)
    trials = loadtrials(options)
	trainingmodels = collect(fit(options,trials[trainingtrials]) for trainingtrials in trainingtrials)
	testmodels = collect(test(trials[testingtrials], trainingmodel) for (testingtrials, trainingmodel) in zip(testingtrials, trainingmodels))
	testinginputs = collect(externalinput(testmodel) for testmodel in testmodels)
	traininginputs = collect(externalinput(trainingmodel) for trainingmodel in trainingmodels)
	return testinginputs, traininginputs
end
