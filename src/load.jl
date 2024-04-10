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
	pathfields = (:datapath, :outputpath)
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
	@assert isfile(options.datapath)
	!isdir(options.outputpath) && mkpath(options.outputpath)
	@assert isdir(options.outputpath)
	return options
end

"""
	loadtrials(options)

RETURN a vector of objects of the composite type `trials`
"""
function loadtrials(options::Options)
	matfile = read(matopen(options.datapath))
	Cell = matfile["Cell"]
	Trials = matfile["Trials"]
	trialindices = (.!Trials["violated"]) .& (Trials["trial_type"] .== "a") .& (Trials["stateTimes"]["cpoke_out"] .> Trials["stateTimes"]["clicks_on"])
	trialindices = vec(trialindices)
	ntimesteps = convert(Integer, (options.reference_end_s - options.reference_begin_s)/options.dt) # intentionally to trigger `InexactError` if `ntimesteps` is not an integer
	binedges_s = options.reference_begin_s:options.dt:(ntimesteps*options.dt)
	reference_times_s = vec(Trials["stateTimes"][options.reference_event][trialindices])
	spiketimes_s = vec(Cell["spiketimes_s"])
	𝐲 = map(reference_times_s) do t
			y = StatsBase.fit(Histogram, spiketimes_s, (t .+ binedges_s)).weights
			convert.(UInt8, y)
		end
	return 𝐲

	# map((spiketrain, stimulus)->Trial(spiketrain=convert.(UInt8, vec(spiketrain)), stimulus=stimulus), vec(trials["spiketrain"]), stimulus)

	map(findall(trialindices)) do i
		reference_time_s = Trials["stateTimes"][options.reference_event][i]

		y = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ binedges_s)).weights

		Trial(binedges_s=binedges_s,
			  e
			  reference_time_s=reference_time_s,
			  trialindex=i,
			  y=y)
	end
end

"""
    Clicks(a_latency_s, L, R, Δt, ntimesteps)

Create an instance of `Clicks` to compartmentalize variables related to the times of auditory clicks in one trial

The stereoclick is excluded.

ARGUMENT
-`a_latency_s`: latency of the accumulator with respect to the clicks
-`Δt`: duration, in seconds, of each time step
-`L`: a vector of floating-point numbers specifying the times of left clicks, in seconds. Does not need to be sorted.
-`ntimesteps`: number of time steps in the trial. Time is aligned to the stereoclick. The first time window is `[-Δt, 0.0)`, and the last time window is `[ntimesteps*Δt, (ntimesteps+1)*Δt)`, defined such that `tₘₒᵥₑ - (ntimesteps+1)*Δt < Δt`, where `tₘₒᵥₑ` is the time when movement away from the center port was first detected.
-`R`: a vector of floating-point numbers specifying the times of right clicks, in seconds. Does not need to be sorted.

RETURN
-an instance of the type `Clicks`
"""
function Clicks(a_latency_s::AbstractFloat,
				Δt::AbstractFloat,
                L::Vector{<:AbstractFloat},
                ntimesteps::Integer,
                R::Vector{<:AbstractFloat})
    L = L[.!isapprox.(L, 0.0)] #excluding the stereoclick
    R = R[.!isapprox.(R, 0.0)]
	L .+= a_latency_s
	R .+= a_latency_s
	rightmost_edge_s = (ntimesteps-1)*Δt
	L = L[L.<rightmost_edge_s]
	R = R[R.<rightmost_edge_s]
    clicktimes = [L;R]
    indices = sortperm(clicktimes)
    clicktimes = clicktimes[indices]
    isright = [falses(length(L)); trues(length(R))]
    isright = isright[indices]
    is_in_timestep =
        map(1:ntimesteps) do t
            ((t-2)*Δt .<= clicktimes) .& (clicktimes .< (t-1)*Δt) # the right edge of the first time step is defined as 0.0, the time of the stereoclick
        end
    right = map(is_in_timestep) do I
                findall(I .& isright)
            end
    isleft = .!isright
    left =  map(is_in_timestep) do I
                findall(I .& isleft)
            end
	inputtimesteps=findall(sum.(is_in_timestep).>0)
	inputindex = map(t->findall(inputtimesteps .== t), 1:ntimesteps)
    Clicks(time=clicktimes,
		   inputtimesteps=inputtimesteps,
		   inputindex=inputindex,
           source=isright,
           left=left,
           right=right)
end


"""

"""
function GLM(inputpath::String)

end
