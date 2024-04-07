"""
    Options

Model settings
"""
@with_kw struct Options{TB<:Bool, TF<:AbstractFloat, TI<:Integer, TS<:String, VSym<:Vector{<:Symbol}}
	"order by which inputs are concatenated"
	concatenationorder::VSym = [:click, :drift, :movement, :postspike, :stereoclick]
	"absolute path of the file containing the data"
    datapath::TS=""
	"duration of each timestep in seconds"
    Δt::TF=0.01
	"""
	the hyperparameters below specify the parametrization of the linear kernels
		- `begin_s`: time relative to the event at which the kernel begins
		- `end_s`: time relative to the event at which the kernel ends
		- `begins0`: whether the first value of the kernel begins at zero
		- `ends0`: whether the last value of the kernel begins at zero
		- `N`: whether the first value of the kernel begins at zero
		- `eta`: nonlinearity in the parametrization of the kernel
	"""
	"click-aligned linear filter"
	k_click_begin_s::TF=0.01
	k_click_end_s::TF=0.5
	k_click_begins0::TB=true
	k_click_ends0::TB=true
	k_click_N::TI=5
	k_click_eta::TF=1.0
	k_click_side_selective::TB=false
	"movement-aligned linear filter"
	k_movement_begin_s::TF= -1.0
	k_movement_end_s::TF=0.5
	k_movement_begins0::TB=true
	k_movement_ends0::TB=false
	k_movement_N::TI=5
	k_movement_eta::TF=0.2
	k_movement_choice_selective::TB=false
	"postspike filter"
	k_postspike_begin_s::TF=0.01
	k_postspike_end_s::TF=0.25
	k_postspike_begins0::TB=false
	k_postspike_ends0::TB=true
	k_postspike_N::TI=5
	k_postspike_eta::TF=1.0
	"kernel aligned to the stereoclick"
	k_stereoclick_begin_s::TF=0.01
	k_stereoclick_end_s::TF=0.5
	k_stereoclick_begins0::TB=false
	k_stereoclick_ends0::TB=false
	k_stereoclick_N::TI=5
	k_stereoclick_eta::TF=1.0
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
	"event on each trial aligned to which spikes are counted"
	reference_event::TS="cpoke_in"
	"Time, in second, after which spikes are included on each trial, aligned to the reference event on that trial."
	reference_before_s::TF=0.5
	"Time, in second, before which spikes are included on each trial, aligned to the reference event on that trial."
	reference_after_s::TF=2.0
end

"""
    Clicks

Information on the clicks delivered during one trial

The stereoclick is excluded.
"""
@with_kw struct Clicks{VF<:Vector{<:AbstractFloat},
                       BA1<:BitArray{1},
					   VI<:Vector{<:Integer},
                       VVI<:Vector{<:Vector{<:Integer}}}
    "A vector of floats indicating the time of each click. It is sorted in ascending order."
    time::VF
	"A vector of integers indicating the timesteps with clicks. Not the timestep of each click"
	inputtimesteps::VI
	"A vector of integers whose element `i = inputindex[t]` indicating the i-th input for that time step. `inputtimesteps[inputindex[t]] == t` and `inputindex[inputtimesteps[i]] == i`"
	inputindex::VVI
    "An one dimensional BitArray specifying the side of each click. The times of the right and left clicks are given by calling `time[source]` and `time[.!source]`, respectively. "
    source::BA1
    "A vector of integers indexing the left clicks that occured in each timestep. The times of the left clicks in the t-th timestep can be found by calling `time[left[t]]`."
    left::VVI
    "A vector of integers indexing the right clicks that occured in each timestep. The times of the right click in the t-th timestep can be found by calling `time[right[t]]`."
    right::VVI
end

"""
    Trial

Information on the sensory stimulus and behavior each trial

Spike trains are not included. In sampled data, the generatives values of the latent variables are stored.
"""
@with_kw struct Trial{TB<:Bool, TC<:Clicks, TF<:AbstractFloat, TI<:Integer, TVI<:Vector{<:Integer}}
    "information on the auditory clicks"
    clicks::TC
    "behavioral choice"
    choice::TB
	"log of the ratio of the generative right and left click rate"
	γ::TF
	"index of the trial in the trialset"
	index_in_trialset::TI
	"time of leaving the center port, relative to the time of the stereoclick, in seconds"
	movementtime_s::TF; @assert movementtime_s > 0
	"time of leaving the center port, relative to the time of the stereoclick, in time steps"
	movementtimestep::TI; @assert movementtimestep > 0
    "number of time steps in this trial. The duration of each trial is from the onset of the stereoclick to the end of the fixation period"
    ntimesteps::TI
	"a nested array whose element `spiketrains[n][t]` is the spike train response of the n-th neuron on the t-th time step of the trial"
	spiketrain::TVI
	"time of the stereoclick, in seconds, in the sessions"
	stereoclick_time_s::TF
	"number of timesteps in the trialset preceding this trial"
	τ₀::TI
end

"""
	GLMKernel
"""
@with_kw struct GLMKernel{TB<:Bool, TF<:AbstractFloat, TI<:Integer, TS<:Symbol, MF<:Matrix{<:AbstractFloat}}
	"name"
	name::TS
	"columns correspond to temporal basis functions, rows to the value of a function at each time step"
	Φ::MF
	"start time relative to the reference event"
	begin_s::TF
	"end time relative to the reference event"
	end_s::TF
	"is the value of the filter equal to 0 at the beginning"
	begin_0::TB
	"is the value of the filter equal to 0 at the end"
	end_0::TB
	"number of temporal basis functions"
	N::TI
	"stretch in the nonlinearity"
	η::TF
end

"""
	GaussianPrior

Information on the zero-meaned Gaussian prior distribution on the values of the parameters in real space
"""
@with_kw struct GaussianPrior{	VI<:Vector{<:Integer},
								VF<:Vector{<:AbstractFloat},
								VR<:Vector{<:Real},
								VS<:Vector{<:String},
								MR<:Matrix{<:Real},
								VVI<:Vector{<:Vector{<:Integer}},
								VMF<:Vector{<:Matrix{<:AbstractFloat}}}
	"L2 penalty matrices"
	𝐀::VMF
	"L2 penalty coefficients"
	𝛂::VR
	"minimum values of the L2 penalty coefficients"
	𝛂min::VF
	"maximum values of the L2 penalty coefficients"
	𝛂max::VF
	"Indices of the parameters related to each L2 penalty coefficient: element `index𝐀[i][j]` corresponds to the i-th group of parameters and the j-th parameter in that group"
	index𝐀::VVI
	"the precision matrix, i.e., inverse of the covariance matrix, of the gaussian prior on the model parameters"
	𝚲::MR
	"indices of the dimensions with finite variance"
	index𝚽::VI = sort(union(index𝐀...))
	"square submatrix of the precision matrix after deleting the columns and rows corresponding to the dimensions with infinite variance"
	𝚽::MR= 𝚲[index𝚽,index𝚽]
	"indices of 𝐀 within `index𝚽`"
	index𝐀_in_index𝚽::VVI = map(indexA->map(indexAᵢⱼ->findfirst(index𝚽.==indexAᵢⱼ), indexA), index𝐀)
	"the name of each L2 penalty"
	penaltynames::VS
end

"""
    GLM

Poisson generalized linear model
"""
@with_kw struct GLM{TO<:Options,
					VF<:Vector{<:AbstractFloat},
					MF<:Matrix{<:AbstractFloat},
					VK<:Vector{<:GLMKernel},
					VT<:Vector{<:Trial},
					GP<:GaussianPrior,
					VI<:Vector{<:UInt8}}
	"gaussian prior"
	gaussianprior::GP
    "fixed hyperparameters"
    options::TO
	"linear filters"
	glmkernels::VK
	"trials"
	trials::VT
	"concatenated weights"
	𝐰::VF
	"design matrix"
	𝐗::MF
	"Poisson observations"
	𝐲::VI
end
