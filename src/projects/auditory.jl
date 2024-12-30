"""
    auditory_project_trials_path

RETURN the absolute path to the "trials.csv" file for the "auditory" project.
"""
function auditory_project_trials_path()
    if Sys.iswindows()
        "x:\\tzluo\\auditory\\processed_data\\2024_10_29a_trials\\trials.csv"
    elseif Sys.isapple()
        error("not yet implemented")
    elseif Sys.islinux()
        "/mnt/cup/labs/brody/tzluo/auditory/processed_data/2024_10_29a_trials/trials.csv"
    else
        error("unknown operating system")
    end
end

"""
	loadtrials_auditory(options)

Load trials for the "auditory" project

RETURN a vector of objects of the composite type `Trial`

ARGUMENT
-`options`: a structure containing fixed hyperparameters
"""
function loadtrials_auditory(options::Options)
    file = matopen(options.datapath)
    Cell = read(file)
    close(file)
    trials = SPGLM.open_auditory_trials_csv(; recording_id = Cell["recording_id"])
    if !isempty(options.trial_indices_path)
        if options.trial_indices_path == "do_not_exclude"
            trialindices = collect(1:length(trials.pokedR))
        else
            trialindices = DataFrames.DataFrame(CSV.File(options.trial_indices_path; header=0))[:,1]
        end
    else
        trialindices = (trials.violated .== 0) .&
                        (trials.trial_type .== "a") .&
                        .!isnan.(trials.pokedR) .&
                        .!isnan.(trials.cpoke_in_s) .&
                        .!isnan.(trials.spoke_s) .&
                        (trials.cpoke_out_s .> trials.clicks_on_s)
        trialindices = findall(vec(trialindices))
    end
    Na = floor(Int, options.time_in_trial_begin_s/options.dt)
    Npre = ceil(Int, (options.bfs_postspike_end_s-options.bfs_postspike_begin_s)/options.dt)
    pre_binedges = -Npre*options.dt:options.dt:0
    spiketimes_s = vec(Cell["spiketimes_s"])
    col_reference_event = Symbol(options.reference_event*"_s")
    col_trim_after_event= Symbol(options.trim_after_event*"_s")
    if options.reference_event == "stereoclick"
        first_reference_time_s = trials.leftclick_s[trialindices[1]][1]
        last_reference_time_s = trials.leftclick_s[trialindices[end]][1]
    else
        first_reference_time_s = trials[trialindices[1], col_reference_event]
        last_reference_time_s = trials[trialindices[end], col_reference_event]
    end
    usepose = (options.pose_filepath != "")
    if usepose
        file = matopen(options.pose_filepath)
        pose = vec(vec.(read(file, "pose")))
        pose_time_s = vec(read(file, "pose_time_s"))
        close(file)
    end
    map(trialindices) do i
        if options.reference_event == "stereoclick"
            reference_time_s = trials.leftclick_s[i][1]
        else
            reference_time_s = trials[i, col_reference_event]
        end
        t₀ = reference_time_s + Na*options.dt
        if isnan(options.time_in_trial_end_s) && !isempty(options.trim_after_event) # the trial is truncated with respect to only `trim_after_event` and not `reference_event`
            Nb = floor(Int, (trials[i, col_trim_after_event]-t₀)/options.dt)
        elseif !isnan(options.time_in_trial_end_s) && !isempty(options.trim_after_event) # the trial is truncated with respect to both `trim_after_event` and `reference_event`
            lasttimestep1 = floor(Int, (trials[i, col_trim_after_event]-t₀)/options.dt)
            lasttimestep2 = ceil(Int, options.time_in_trial_end_s/options.dt)
            Nb = min(lasttimestep1,lasttimestep2)
            @assert Nb > 0
        else # the trial is truncated with respect to only `reference_event` and not `trim_after_event`
            Nb = ceil(Int, options.time_in_trial_end_s/options.dt)
        end
        binedges_s = (Na*options.dt):options.dt:(Nb*options.dt)
        hist = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ binedges_s), closed=:left)
        y = hist.weights
        ypre = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ pre_binedges), closed=:left).weights
        if usepose # assumes sampling rate of the pose signals is lower than that of the neural
            ntimesteps = length(binedges_s)-1
            pose_sorted = collect(fill(NaN,ntimesteps) for pose in pose)
            for j = 1:ntimesteps
                index = findfirst(pose_time_s .>= (reference_time_s+binedges_s[j]))
                for k in eachindex(pose)
                    pose_sorted[k][j] = sum(pose[k][index])
                end
            end
        else
            pose_sorted = collect(zeros(0) for i=1:length(trialindices))
        end
        @assert !isnan(trials[i,:cpoke_out_s])
        fixation_timestep = ceil(Int, (trials[i,:cpoke_in_s] - t₀)/options.dt)
        movement_timestep = ceil(Int, (trials[i,:cpoke_out_s] - t₀)/options.dt)
        leftclicks_s = trials[i,:leftclick_s] .- t₀
        rightclicks_s = trials[i,:rightclick_s] .- t₀
        @assert leftclicks_s[1] == rightclicks_s[1]
        stereoclick_timestep = ceil(Int, leftclicks_s[1]/options.dt)
        leftclicks_timestep = ceil.(Int, leftclicks_s[2:end]./options.dt)
        rightclicks_timestep = ceil.(Int, rightclicks_s[2:end]./options.dt)
        response_timestep = ceil(Int, (trials[i,:spoke_s] - t₀)/options.dt)
        clicks_timestep = [stereoclick_timestep; leftclicks_timestep; rightclicks_timestep]
        clicks_source = [2; zeros(Int, length(leftclicks_timestep)); ones(Int, length(rightclicks_timestep))]
        SPGLM.Trial(choice = trials[i,:pokedR] == 1,
                clicks_source=clicks_source,
                clicks_timestep=clicks_timestep,
                first_reference_time_s=first_reference_time_s,
                fixation_timestep = fixation_timestep,
                γ = trials[2,:gamma],
                last_reference_time_s=last_reference_time_s,
                movement_timestep = movement_timestep,
                pose=pose_sorted,
                reference_time_s=reference_time_s,
                response_timestep = response_timestep,
                timesteps_s=binedges_s[2:end],
                trialindex=i,
                y=y,
                ypre=ypre)
    end
end

"""
    open_auditory_trials_csv()

Open the "trials.csv" file for the "auditory" project and extract the left and right click times.

RETURN a data frame

OPTIONAL ARGUMENT
-`recording_id`: a string specifying the recording ID
"""
function open_auditory_trials_csv(; recording_id::String="")
    trialspath = SPGLM.auditory_project_trials_path()
    trials = DataFrames.DataFrame(CSV.File(trialspath))
    if recording_id != ""
        trials = trials[trials.recording_id .== recording_id,:]
    end
    ntrials = size(trials,1)
    leftclick_s = collect(zeros(0) for i in 1:ntrials)
    rightclick_s = collect(zeros(0) for i in 1:ntrials)
    for i = 1:ntrials
        if contains(trials.leftclick_s[i], "[")
            leftclick_s[i] = collect(parse(Float64,String(x)) for x in split(chop(trials.leftclick_s[i],head=1)))
        elseif contains(trials.leftclick_s[i], ".")
            leftclick_s[i] = [parse(Float64,trials.leftclick_s[i])]
        end
        if contains(trials.rightclick_s[i], "[")
            rightclick_s[i] = collect(parse(Float64,String(x)) for x in split(chop(trials.rightclick_s[i],head=1)))
        elseif contains(trials.rightclick_s[i], ".")
            rightclick_s[i] = [parse(Float64,trials.rightclick_s[i])]
        end
    end
    trials.leftclick_s = leftclick_s
    trials.rightclick_s = rightclick_s
    trials
end