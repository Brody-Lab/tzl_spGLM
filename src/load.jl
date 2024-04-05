"""
	Options(csvpath, row)

RETURN a struct containing the fixed hyperparameters of the model

ARGUMENT
-`csvpath`: the absolute path to a comma-separated values (CSV) file
-`row`: the row of the CSV to be considered
"""
Options(csvpath::String) = Options(csvpath::String, 1)
Options(csvpath::String, row::Integer) = Options(DataFrames.DataFrame(CSV.File(csvpath)), row)
Options(df::DataFrames.DataFrame, row::Integer) = Options(df[row,:])
Options(dfrow::DataFrames.DataFrameRow) = Options(Dict((name=>dfrow[name] for name in names(dfrow))...))

"""
	Options(options::Dict)

RETURN a struct containing the fixed hyperparameters of the model

ARGUMENT
-`options`: a dictionary
"""
function Options(options::Dict)
	keyset = keys(options)
	defaults = Options()
	pathfields = (:inputpath, :outputfolder)
	entries = 	map(fieldnames(Options)) do fieldname
					defaultvalue = getfield(defaults,fieldname)
					if String(fieldname) ∈ keyset
						convert(typeof(defaultvalue), options[String(fieldname)])
					elseif fieldname ∈ pathfields
						convert(typeof(defaultvalue), options[String(fieldname)])
					else
						defaultvalue
					end
				end
	options = Options(entries...)
	for fieldname in pathfields
		@assert !isempty(getfield(options, fieldname))
	end
	@assert isfile(options.inputpath)
	!isdir(options.outputfolder) && mkpath(options.outputfolder)
	@assert isdir(options.outputfolder)
	return options
end

"""
	loadtrials(options)

RETURN a vector of objects of the composite type `trials`
"""
function loadtrials(inputpath::String)
	trials = read(matopen(inputpath))["trials"]
	if load_choice_as_stimulus
		stimulus = Float64.(2vec(trials["choice"])) .- 1.0
	else
		stimulus = haskey(trials, "stimulus") ? vec(trials["stimulus"]) : sign.(vec(trials["gamma"]))
	end
	map((spiketrain, stimulus)->Trial(spiketrain=convert.(UInt8, vec(spiketrain)), stimulus=stimulus), vec(trials["spiketrain"]), stimulus)
end
