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
