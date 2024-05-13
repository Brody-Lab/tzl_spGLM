"""
RETURN an instance of a type converted into an instance of a `Dict`
"""
dictionary(x) = Dict((String(fieldname)=>getfield(x,fieldname) for fieldname in fieldnames(typeof(x)))...)
function dictionary(eo::EvidenceOptimization)
    Dict("a"=>eo.a,
        "E"=>eo.𝐸,
        "hessian_loglikelihood"=>eo.hessian_loglikelihood,
        "MAP_g_residual"=>eo.MAP_g_residual,
        "w"=>eo.𝐰)
end

function dictionary(basisset::BasisFunctionSet)
    Dict("name"=>String(basisset.name),
          "timesteps_s"=>collect(basisset.timesteps_s),
          "Phi"=>basisset.Φ)
end

function dictionary(model::Model)
    Dict("a"=>model.a,
          "baseline"=>model.baseline,
          "kernels"=>map(dictionary, convolutionkernels(model)),
          "options"=>dictionary(model.options),
          "basissets"=>map(dictionary, model.basissets),
          "weightindices"=>Dict((String(fieldname)=>collect(getfield(model.weightindices,fieldname)) for fieldname in fieldnames(typeof(model.weightindices)))...),
          "w"=>model.𝐰)
end
"""
WRITE a MATLAB ".mat" file
"""
save(model::Model) = matwrite(joinpath(model.options.outputpath, "model.mat"), dictionary(model))
save(eo::EvidenceOptimization, outputpath::String) = matwrite(joinpath(outputpath, "evidenceoptimization.mat"), dictionary(eo))
save(peths::Vector{<:PerieventTimeHistogram}, outputpath::String) = matwrite(joinpath(outputpath, "peths.mat"), Dict("peths"=>map(dictionary, peths)))
function save(characterization::Characterization, outputpath::String)
 	map(fieldnames(SPGLM.Characterization)) do fieldname
		dict = Dict(string(fieldname)=>getfield(characterization,fieldname))
		path = joinpath(outputpath, string(fieldname)*".mat")
		matwrite(path, dict)
	end
end
"""
	convolutionkernels(model)

RETURN a vector whose each element is an instance of type `Kernel`
"""
function convolutionkernels(model::Model)
	inputnames = collect(fieldnames(typeof(model.weightindices)))
	inputnames = filter(!isequal(:baseline), inputnames)
	indices = collect(!isempty(getfield(model.weightindices, fieldname)) for fieldname in inputnames)
	inputnames = inputnames[indices]
	map(inputnames) do inputname
		basisset = filter((basis)->match_input_to_basis(inputname)==basis.name, model.basissets)[1]
		h = basisset.Φ*model.𝐰[getfield(model.weightindices, inputname)]
		Kernel(basisname=String(basisset.name),
				h=h,
				inputname=String(inputname),
				timesteps_s=collect(basisset.timesteps_s))
	end
end

"""
	save(characterization, trials, outputpath)

Save the observed and predicted autocorrelation functions and nats per spike
"""
function save(characterization::Characterization, trials::Vector{<:Trial}, outputpath::String)
	N = minimum(length(𝐫) for 𝐫 in characterization.autocorrelation)
	𝐑 = filter(x->!isnan(x[1]), characterization.autocorrelation)
	𝐫_pred = zeros(N)
	for 𝐫 in 𝐑
		for t = 1:N
			𝐫_pred[t] += 𝐫[t]/length(𝐑)
		end
	end
	𝐫_obsv = zeros(N)
	count = 0
	lags = collect(1:N)
	for trial in trials
		if sum(trial.y) > 0
			count += 1
			𝐫_obsv.+=StatsBase.autocor(trial.y,lags)
		end
	end
	𝐫_obsv ./= count
	Δt = trials[1].timesteps_s[2]-trials[1].timesteps_s[1]
	dict = Dict("autocor_obsv"=>𝐫_obsv,
				"autocor_pred"=>𝐫_pred,
				"lags"=>lags.*Δt)
	matwrite(joinpath(outputpath, "autocor.mat"), dict)
	nats_per_spike = sum(sum.(characterization.LL))/sum(sum.(characterization.observed_spiketrains))
	matwrite(joinpath(outputpath, "nats_per_spike.mat"), Dict("nats_per_spike"=>nats_per_spike))
end
