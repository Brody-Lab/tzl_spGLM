"""
    Options

Model settings
"""
@with_kw struct Options{TB<:Bool, TS<:String, TF<:AbstractFloat, TI<:Integer, TVF<:Vector{<:AbstractFloat}}
	"response latency of the accumulator to the clicks"
    a_latency_s::TF=1e-2
	"value optimized when initializing the choice-related parameters"
	choiceobjective::TS="posterior"
	"Exponent used to compute the scale factor of the log-likelihood of the choices. The scale factor is computed by raising the product of the number of neurons and the average number of time steps in each trial to the exponent. An exponent of 0 means no scaling"
	choiceLL_scaling_exponent::TF=0.6
	"full path of the data"
    datapath::TS=""
	"duration of each timestep in seconds"
    Δt::TF=0.01
	"whether the transition probability of remaining in the first state is fitted"
	fit_Aᶜ₁₁::TB=false
	"whether the transition probability of remaining in the second state is fitted"
	fit_Aᶜ₂₂::TB=false
	"whether to fit the height of the sticky bounds"
	fit_B::TB=true
	"whether to fit the parameter for transforming the accumulator"
	fit_b::TB=false
	"whether to fit separate encoding weights for when the accumulator is at the bound"
	fit_𝛃::TB=true
	"whether to fit the exponential change rate of inter-click adaptation"
	fit_k::TB=false
	"whether to fit the parameter specifying leak or instability"
	fit_λ::TB=false
	"whether to fit the constant added to the mean of the distribution of the accumulator variable at the first time step"
	fit_μ₀::TB=true
	"whether to fit an overdispersion parameter for the model of each neuron's spike count response. If true, the count response model is negative binomial rather than Poisson"
	fit_overdispersion::TB=false
	"whether to fit the strength of inter-click adaptation and sign of the adaptation (facilitation vs. depression)"
	fit_ϕ::TB=false
	"whether the prior probability of the first state is fitted"
	fit_πᶜ₁::TB=false
	"whether to fit the behavioral lapse rate"
	fit_ψ::TB=false
	"whether to fit the variance of the Gaussian noise added at each time step"
	fit_σ²ₐ::TB=false
	"whether to fit the variance of the prior probability of the accumulator variable"
	fit_σ²ᵢ::TB=false
	"whether to fit the variance of Gaussian noise added as a result of the clicks"
	fit_σ²ₛ::TB=true
	"whether to fit the weight of the rewarded option of the previous trial on the mean of the accumulator at the first time step"
	fit_wₕ::TB=false
	"L2 norm of the gradient at which convergence of model's cost function is considered to have converged"
	g_tol::TF=1e-8
	"number of states of the coupling variable"
	K::TI = 1; @assert (K==1 || K == 2)
	"maximum and minimum of the L2 shrinkage penalty for each class of parameters. The penalty is initialized as (and if not being learned, set as) the geometric mean of the maximum and minimum."
	"accumulator transformation"
	L2_b_max::TF=1e1
	L2_b_min::TF=1e-3
	"latent variable parameter when fitting to only choices"
	L2_choices_max::TF=1e0
	L2_choices_min::TF=1e-4
	"gain"
	L2_gain_max::TF=1e-5
	L2_gain_min::TF=1e-7
	"spike history"
	L2_postspike_max::TF=1e-2
	L2_postspike_min::TF=1e-4
	"pre-movement"
	L2_premovement_max::TF=1e-2
	L2_premovement_min::TF=1e-4
	"post-photostimulus filter"
	L2_postphotostimulus_max::TF=1e-2
	L2_postphotostimulus_min::TF=1e-4
	"post-stereoclick filter"
	L2_poststereoclick_max::TF=1e-2
	L2_poststereoclick_min::TF=1e-4
	"latent variable parameter"
	L2_latent_max::TF=10.0
	L2_latent_min::TF=0.1
	"encoding of accumulated evidence"
	L2_accumulator_max::TF=1e-4
	L2_accumulator_min::TF=1e-6
	"Value in native space corresponding to the lower bound ('l'), zero-value in real space ('q'), and upper bound ('u')"
	"transition probability of the coupling variable to remain in the coupled state"
	Aᶜ₁₁_l::TF=1e-4
	Aᶜ₁₁_q::TF=0.5
	Aᶜ₁₁_u::TF= 1.0-1e-4
	"transition probability of the coupling variable to remain in the decoupled state"
	Aᶜ₂₂_l::TF=1e-4
	Aᶜ₂₂_q::TF=0.5
	Aᶜ₂₂_u::TF= 1.0-1e-4
	"bound height"
	B_l::TF=10.0
	B_q::TF=15.0
	B_u::TF=20.0
	"adaptation change rate"
	k_l::TF=10.0
	k_q::TF=100.0
	k_u::TF=1000.0
	"feedback"
	λ_l::TF = -1.0
	λ_q::TF = 0.0
	λ_u::TF = 2.0
	"bias"
	μ₀_l::TF=-5.0
	μ₀_q::TF= 0.0
	μ₀_u::TF= 5.0
	"adaptation strength"
	ϕ_l::TF=1e-4
	ϕ_q::TF=1e-3
	ϕ_u::TF=1.0-1e-4
	"prior probability of the coupled state"
	πᶜ₁_l::TF=1e-4
	πᶜ₁_q::TF=0.5
	πᶜ₁_u::TF=1-1e-4
	"behavioral lapse rate. A value of 0 will result in underflow"
	ψ_l::TF=1e-4
	ψ_q::TF=1e-2
	ψ_u::TF=1-1e-4
	"variance of the diffusion noise (i.e., zero-meaned, independent and identically distributed gaussian noise added at each time step)"
	σ²ₐ_l::TF=0.1
	σ²ₐ_q::TF=1.0
	σ²ₐ_u::TF=10.0
	"variance of the initial probability of the accumulator variable"
	σ²ᵢ_l::TF=0.1
	σ²ᵢ_q::TF=1.0
	σ²ᵢ_u::TF=10.0
	"variance of the variance of per-click noise"
	σ²ₛ_l::TF = 0.1
	σ²ₛ_q::TF = 5.0
	σ²ₛ_u::TF = 20.0
	"weight of previous answer"
	wₕ_l::TF = -5.0
	wₕ_q::TF = 0.0
	wₕ_u::TF = 5.0
	"maximum duration of each trial"
	maxduration_s::TF=1.0
	"minimum value of the prior and transition probabilities of the accumulator"
	minpa::TF=1e-8
	"value to maximized to learn the parameters"
	objective::TS="posterior"; @assert any(objective .== ["evidence", "posterior", "likelihood", "initialization"])
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
	"coefficient multiplied to any scale factor of the temporal basis functions of a Poisson mixture GLM"
	sf_tbf::TVF=[NaN]
    "scale factor of the conditional likelihood of the spiking of a neuron at a time step"
	sf_y::TF=1.2
	"whether the temporal basis functions parametrizing the weight of the accumulator is at the trough or at the peak in the beginning of the trial"
	tbf_accumulator_begins0::TB=false
	"whether the temporal basis functions parametrizing the weight of the accumulator is at the trough or at the peak in the end of the trial"
	tbf_accumulator_ends0::TB=false
	"number of temporal basis functions parametrizing the weight of the accumulator per second"
	tbf_accumulator_hz::TF=0.0
	"scale factor of the temporal basis functions"
	tbf_accumulator_scalefactor::TF=5.0
	"degree to which temporal basis functions centered at later times in the trial are stretched. Larger values indicates greater stretch. This value must be positive"
	tbf_accumulator_stretch::TF=0.2
	"scale factor of the gain parameter"
	tbf_gain_scalefactor::TF=5.0
	"maximum number of basis functions"
	tbf_gain_maxfunctions::TI=8
	"Options for the temporal basis function whose linear combination constitute the post-spike filter. The setting `tbf_postspike_dur_s` is the duration, in seconds, of the filter. The setting `tbf_postspike_linear` determines whether a linear function is included in the basis."
	tbf_postspike_begins0::TB=false
	tbf_postspike_dur_s::TF=0.25
	tbf_postspike_ends0::TB=true
	tbf_postspike_hz::TF=12.0
	tbf_postspike_scalefactor::TF=1.0
	tbf_postspike_stretch::TF=1.0
	"Options for the temporal basis associated with the pre-movement filter"
	tbf_premovement_begins0::TB=true
	tbf_premovement_dur_s::TF=0.6
	tbf_premovement_ends0::TB=false
	tbf_premovement_hz::TF=3.0
	tbf_premovement_scalefactor::TF=5.0
	tbf_premovement_stretch::TF=0.1
	"Options for the temporal basis associated with the post-photostimulus filter"
	tbf_postphotostimulus_begins0::TB=false
	tbf_postphotostimulus_ends0::TB=false
	tbf_postphotostimulus_hz::TF=NaN
	tbf_postphotostimulus_scalefactor::TF=5.0
	tbf_postphotostimulus_stretch::TF=1.0
	"Options for the temporal basis associated with the post-stereoclick filter"
	tbf_poststereoclick_begins0::TB=false
	tbf_poststereoclick_dur_s::TF=1.0
	tbf_poststereoclick_ends0::TB=true
	tbf_poststereoclick_hz::TF=5.0
	tbf_poststereoclick_scalefactor::TF=5.0
	tbf_poststereoclick_stretch::TF=0.2
	"scale factor for the accumulator transformation parameter"
	tbf_b_scalefactor::TF=1.0
    "number of states of the discrete accumulator variable"
    Ξ::TI=53; @assert isodd(Ξ) && Ξ > 1
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
@with_kw struct Trial{TB<:Bool, TC<:Clicks, TF<:AbstractFloat, TI<:Integer, TVVI<:Vector{<:Vector{<:Integer}}}
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
	"time when the offset ramp of the photostimulus began"
	photostimulus_decline_on_s::TF
	"time when the onset ramp of the photostimulus began"
	photostimulus_incline_on_s::TF
    "location of the reward baited in the previous trial (left:-1, right:1, no previous trial:0)"
    previousanswer::TI
	"a nested array whose element `spiketrains[n][t]` is the spike train response of the n-th neuron on the t-th time step of the trial"
	spiketrains::TVVI
	"time of the stereoclick, in seconds, in the sessions"
	stereoclick_time_s::TF
	"number of timesteps in the trialset preceding this trial"
	τ₀::TI
    "index of the trialset to which this trial belongs"
	trialsetindex::TI
end

"""
	Indices𝐮

Indices of the encoding weights of the temporal basis vectors of the filters that are independent of the accumulator
"""
@with_kw struct Indices𝐮{UI<:UnitRange{<:Integer}}
	gain::UI=1:1
	postspike::UI
	poststereoclick::UI
	premovement::UI
	postphotostimulus::UI
end

"""
	GLMθ

Parameters of a mixture of Poisson generalized linear model
"""
@with_kw struct GLMθ{B<:Bool, IU<:Indices𝐮, R<:Real, VR<:Vector{<:Real}, VS<:Vector{<:Symbol}, VVR<:Vector{<:Vector{<:Real}}}
	"overdispersion parameter in real space. It is mapped into a nonnegative value using the softplus function."
	a::VR=[-Inf]
    "nonlinearity in accumulator transformation"
	b::VR=[NaN]
	"scale factor for the nonlinearity of accumulator transformation"
	b_scalefactor::R
	"order by which parameters are concatenated"
	concatenationorder::VS = [:𝐮, :𝐯, :𝛃, :a, :b]
	"whether the nonlinearity parameter is fit"
	fit_b::B
	"whether to fit separate encoding weights for when the accumulator at the bound"
	fit_𝛃::B
	"whether to fit an overdispersion parameter. If so, the count response model is negative binomial rather than Poisson"
	fit_overdispersion::B
	"state-independent linear filter of inputs from the spike history and time in the trial"
    𝐮::VR
	"Indices of the encoding weights of the temporal basis vectors of the filters that are independent of the accumulator"
	indices𝐮::IU
    "state-dependent linear filters of the inputs from the accumulator "
    𝐯::VVR
	"state-dependent linear filters of the time-varying input from the transformed accumulated evidence"
	𝛃::VVR=deepcopy(𝐯)
end

"""
    MixturePoissonGLM

Mixture of Poisson generalized linear model
"""
@with_kw struct MixturePoissonGLM{TI<:Integer,
								  F<:AbstractFloat,
								  UI<:UnitRange{<:Integer},
                                  VF<:Vector{<:AbstractFloat},
								  VI<:Vector{<:UInt8},
								  Tθ<:GLMθ,
                                  MF<:Matrix{<:AbstractFloat}}
    "size of the time bin"
    Δt::F
	"Normalized values of the accumulator"
    d𝛏_dB::VF
	"scale factor multiplied to the likelihood to avoid underflow when computing the likelihood of the population response"
	likelihoodscalefactor::F
	"Values of the smooth temporal basis functions used to parametrize the time-varying weight of accumulator. Columns correspond to temporal basis functions, and rows correspond to time steps, concatenated across trials."
	Φaccumulator::MF
	"values of the basis functions parametrizing the slow drift in gain on each trial"
	Φgain::MF
	"Values of the smooth temporal basis functions used to parametrize the post-spike filter"
	Φpostspike::MF
	"Values of the smooth temporal basis functions used to parametrize the time-varying relationship between the timing of the animal leaving the center and the neuron's probability of spiking. The timing is represented by a delta function, and the delta function is convolved with a linear combination of the temporal basis functions to specify the filter, or the kernel, of the event. The columns correspond to temporal basis functions and rows correspond to time steps, concatenated across trials."
	Φpremovement::MF
	"temporal basis vectors for the photostimulus"
	Φpostphotostimulus::MF
	"time steps of the temporal basis vectors relative to the onset of the photostimulus"
	Φpostphotostimulus_timesteps::UI
	"Values of the smooth temporal basis functions used to parametrize the time-varying relationship between the timing of the stereoclick and the neuron's probability of spiking."
	Φpoststereoclick::MF
	"parameters"
	θ::Tθ
    "Input of the accumulator. The first column consists of ones. The subsequent columns, if any, correspond to the time-varying input of the accumulator. Element 𝐕[t,i] corresponds to the value of the i-th temporal basis function at the t-th time bin"
    𝐕::MF
	"design matrix. The first column are ones. The subsequent columns correspond to spike history-dependent inputs. These are followed by columns corresponding to the time-dependent input. The last set of columns are given by 𝐕"
	𝐗::MF
	"columns corresponding to the state-independent inputs"
	𝐗columns_𝐮::UI = 1:(size(𝐗,2)-size(𝐕,2))
	"columns corresponding to the input from the accumulator"
	𝐗columns_𝐯::UI = (size(𝐗,2)-size(𝐕,2)+1):size(𝐗,2)
	"number of accumulator states"
	Ξ::TI=length(d𝛏_dB)
	"Poisson observations"
    𝐲::VI
end

"""

Quantities and memory used for computing the conditional partial derivatives of spike count response generalized linear model
"""
@with_kw struct GLMDerivatives{B<:Bool, F<:AbstractFloat, VF<:Vector{<:AbstractFloat}, MF<:Matrix{<:AbstractFloat}}
	"overdispersion parameter"
	α::VF=fill(NaN,1)
	"time step duration, in seconds"
	Δt::F
	"log-likelihood"
	ℓ::VF=fill(NaN,1)
	"derivative of the overdispersion parameter with respect to its real-valued parameter"
	dα_da::VF=fill(NaN,1)
	"second derivative of the overdispersion parameter with respect to its real-valued parameter"
	d²α_da²::VF=fill(NaN,1)
	"first-order partial derivative of the log-likelihood with respect to the real-valued over-dispersion parameter"
	dℓ_da::VF=fill(NaN,1)
	"first-order partial derivative of the log-likelihood with respect to the linear predictor"
	dℓ_dL::VF=fill(NaN,1)
	"second-order partial derivative of the log-likelihood with respect to the real-valued over-dispersion parameter"
	d²ℓ_da²::VF=fill(NaN,1)
	"second-order partial derivative of the log-likelihood with respect to the real-valued over-dispersion parameter and the linear predictor"
	d²ℓ_dadL::VF=fill(NaN,1)
	"second-order partial derivative of the log-likelihood with respect to the linear predictor"
	d²ℓ_dL²::VF=fill(NaN,1)
	"whether the model is a gamma-poisson mixture or a poisson"
	fit_overdispersion::B
	"vector for in-place computation of first-order partial derivatives "
	g::VF=fill(NaN,2)
	"matrix for in-place computation of second-order partial derivatives"
	H::MF=fill(NaN,2,2)
end

"""
    Trialset

A group of trials in which a population of neurons were recorded simultaneously
"""
@with_kw struct Trialset{VM<:Vector{<:MixturePoissonGLM}, TI<:Integer, VT<:Vector{<:Trial}}
	"Mixture of Poisson GLM of each neuron in this trialset"
    mpGLMs::VM=MixturePoissonGLM[]
	"number of time steps summed across trials"
	ntimesteps::TI=size(mpGLMs[1].𝐗,1)
	"Information on the stimulus and behavior for each trial in this trial-set"
    trials::VT
	"Number of trials"
	ntrials::TI=length(trials)
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
	Model

A factorial hidden Markov drift-diffusion model
"""
@with_kw struct Model{Toptions<:Options,
					GP<:GaussianPrior,
					Tθ1<:Latentθ,
					Tθ2<:Latentθ,
					Tθ3<:Latentθ,
					VT<:Vector{<:Trialset}}
	"settings of the model"
	options::Toptions
	"Gaussian prior on the parameters"
	gaussianprior::GP
	"model parameters in their native space (the term 'space' is not meant to be mathematically rigorous. Except for the sticky bound `B`, the native space of all parameters are positive real numbers, which is a vector space. The native space of `B` is upper bounded because I am concerned a large value of `B` would result in a loss of precision in the discretization of the accumulator variable.)"
	θnative::Tθ1
	"model parameters in real vector space ℝ"
	θreal::Tθ2
	"initial values of the parameters in native space"
	θ₀native::Tθ3
	"data used to constrain the model"
	trialsets::VT
end

"""
    Indexθ

Index of each model parameter if all values that were being fitted were concatenated into a vector
"""
@with_kw struct Indexθ{L<:Latentθ, VVG<:Vector{<:Vector{<:GLMθ}}}
	"parameters specifying the mixture of Poisson generalized linear model"
	glmθ::VVG
	"parameters specifying the latent variables"
	latentθ::L
end

"""
	CVIndices

Indices of trials and timesteps used for training and testing
"""
@with_kw struct CVIndices{VVI<:Vector{<:Vector{<:Integer}}}
	"`testingtrials[i]` indexes the trials from the i-th trialset used for testing"
	testingtrials::VVI
	"`trainingtrials[i]` indexes the trials from the i-th trialset used for training"
	trainingtrials::VVI
	"`testingtimesteps[i]` indexes the time steps from the i-th trialset used for testing"
	testingtimesteps::VVI
	"`trainingtimesteps[i]` indexes the time steps from the i-th trialset used for training"
	trainingtimesteps::VVI
end

"""
	Adaptedclicks

The post-adaptation magnitude of each click and the first- and second-order partial derivatives of the post-adaptation magnitude
"""
@with_kw struct Adaptedclicks{TVR1<:Vector{<:Real}, TVR2<:Vector{<:Real}}
	"adapted strengths of the clicks"
	C::TVR1
	"derivative of adapted click strengths with respect to the adaptation change rate"
	dC_dk::TVR2=zeros(0)
	"derivative of adapted click strengths with respect to the adaptation strength"
	dC_dϕ::TVR2=zeros(0)
	"second derivative of adapted click strengths with respect to the adaptation change rate"
	d²C_dkdk::TVR2=zeros(0)
	"second derivative of adapted click strengths with respect to the adaptation change rate and adaptation strength"
	d²C_dkdϕ::TVR2=zeros(0)
	"second derivative of adapted click strengths with respect to the adaptation strength"
	d²C_dϕdϕ::TVR2=zeros(0)
end

"""
	TrialSample

Behavioral choice and neuronal spike trains simulated by running the model forward in time using the auditory click
"""
@with_kw struct TrialSample{B<:Bool, VF<:Vector{<:AbstractFloat}, VVF<:Vector{<:Vector{<:AbstractFloat}}, VVI<:Vector{<:Vector{<:Integer}}}
	"a nested array whose element `accumulator[i][m][t]` is the value of the accumulated evidence on the t-th time step of the m-th trial on the i-th trialset"
	accumulator::VF
	"`true` mean a right choice and `false` a left choice"
	choice::B
	"a nested array whose element `λ[n][t]` is the value of generative spikes per s of the n-th neuron on the t-th time step of the m-th trial on the i-th trialset"
	λ::VVF
	"a nested array whose element `spiketrains[n][t]` is the spike train response of the n-th neuron on the t-th time step of the m-th trial on the i-th trialset"
	spiketrains::VVI
end

"""
	TrialsetSample

Simulated behavioral choice and neuronal spike trains on every trial of a trialset
"""
@with_kw struct TrialsetSample{VS<:Vector{<:TrialSample}}
	trials::VS
end

"""
	Sample

Behavioral choice and neuronal spike trains simulated by running the model forward in time using the auditory clicks on every trial of every trialset
"""
@with_kw struct Sample{VS<:Vector{<:TrialsetSample}}
	trialsets::VS
end

"""
	Expectation of the choice and spike train on a single trial

The unconditioned spike train of the n-th neuron on this trial can be computed as `spiketrain_leftchoice[n].*(1-rightchoice) .+ spiketrain_rightchoice[n].*rightchoice`

If no left choice occurred during simulations, then `spiketrain_leftchoice[n][t]` is zero for all `n` and `t`.

If we have a vector `𝐄` whose each element is a composite of the type `ExpectedEmissions` corresponding to a single trial, then the left-choice conditioned peri-stimulus time historam of the n-th neuron can be computed as follows
```
"""
@with_kw struct ExpectedEmissions{F<:AbstractFloat, VVF<:Vector{<:Vector{<:AbstractFloat}}}
	"expectation of the choice random variable. This can be interpreted as the expected probability of a right choice"
	rightchoice::F
	"expectation of the spike trains, conditioned on a left choice. The element `spiketrain_leftchoice[n][t]` corresponds to the n-th neuron at the t-th time step."
	spiketrain_leftchoice::VVF
	"expectation of the spike trains, conditioned on right choice"
	spiketrain_rightchoice::VVF
end


"""
	Quantities helpful for understanding the model
"""
@with_kw struct Characterization{VVE<:Vector{<:Vector{<:ExpectedEmissions}},
							VVF<:Vector{<:Vector{<:AbstractFloat}},
							VVVVF<:Vector{<:Vector{<:Vector{<:Vector{<:AbstractFloat}}}}}
	"probability of the accumulator variable given the actual auditory clicks. The element `paccumulator[i][m][t][j]` corresponds to the probability of the accumulator in the j-th state during the t-th time step of the m-th trial in the i-th trialset"
	paccumulator::VVVVF
	"posterior probability of the accumulator variable conditioned on the clicks and the behavioral choice. The format is identical to that of `paccumulator`."
	paccumulator_choice::VVVVF
	"posterior probability of the accumulator variable conditioned on the clicks, behavioral choice, and spike trains. The format is identical to that of `paccumulator`."
	paccumulator_choicespikes::VVVVF
	"posterior probability of the accumulator variable conditioned on the clicks and spike trains. The format is identical to that of `paccumulator`."
	paccumulator_spikes::VVVVF
	"log(base 2)-likelihood of the emissions on each trial. The element `LL[i][m]` is the log-likelihood of the choice in the m-th trial in the i-th trialset."
	LL::VVF
	"log(base 2)-likelihood of the observed behavioral choice. The element `LLchoice[i][m]` is the log-likelihood of the choice in the m-th trial in the i-th trialset."
	LLchoice::VVF
	"log(base 2)-likelihood of the observed choice under a Bernoulli model. The format is identical to that of the field `LLchoice`."
	LLchoice_bernoulli::VVF
	"log(base 2)-likelihood of the observed choice conditioned on the spike trains. The format is identical to that of the field `LLchoice`."
	LLchoice_spikes::VVF
	"log(base 2)-likelihood of the spike count at each time step. Element `LLspikes[i][n][m][t]` is the log-likelihood of observed spike count at the t-time step of the m-th trial for the n-th neuron in the i-th trialset."
	LLspikes::VVVVF
	"log(base 2)-likelihood of the spike count at each time step under a homogeneous Poisson model. The format is identical to that of the field `LLspikes`."
	LLspikes_poisson::VVVVF
	"expectation of the choice and spike trains computed averaging across simulations of the model. The element `expectedemissions[i][m]` contains the expected choice and the choice-conditioned spike trains of each neuron on the m-th trial of the i-th trialset."
	expectedemissions::VVE
end

"""
	ModelSummary

Features of the model useful for analysis
"""
@with_kw struct ModelSummary{F<:AbstractFloat,
							LT<:Latentθ,
							MF<:Matrix{<:AbstractFloat},
							VF<:Vector{<:AbstractFloat},
							VVF<:Vector{<:Vector{<:AbstractFloat}},
							VMF<:Vector{<:Matrix{<:AbstractFloat}},
							VVVF<:Vector{<:Vector{<:Vector{<:AbstractFloat}}},
							VVMF<:Vector{<:Vector{<:Matrix{<:AbstractFloat}}},
							VS<:Vector{<:String},
							VVGT<:Vector{<:Vector{<:GLMθ}},
							VVI<:Vector{<:Vector{<:Integer}}}
	"Weighted inputs, except for those from the latent variables and the spike history, to each neuron on each time step in a trialset. The element `externalinputs[i][n][t]` corresponds to the input to the n-th neuron in the i-th trialset on the t-th time step in the trialset (time steps are across trials are concatenated)"
	externalinputs::VVVF
	"the log of the likelihood of the data given the parameters"
	loglikelihood::F
	"the log of the likelihood of the data in each trial given the parameters. Element `loglikelihood_each_trial[i][m]` corresponds to the m-th trial of the i-th trialset"
	loglikelihood_each_trial::VVF
	"the log of the posterior probability of the parameters"
	logposterior::F
	"values of the parameters of the latent variable in their native space"
	thetanative::LT
	"values of the parameters of the latent variable mapped to real space"
	thetareal::LT
	"initial values of parameters of the latent variable in their space"
	theta0native::LT
	"parameters of each neuron's GLM. The element `θglm[i][n]` corresponds to the n-th neuron in the i-th trialset"
	thetaglm::VVGT
	"temporal basis vectors for accumulator encoding"
	temporal_basis_vectors_accumulator::VMF
	"temporal basis vectors for the gain on each trial"
	temporal_basis_vectors_gain::VVMF
	"temporal basis vectors for the post-spike kernel"
	temporal_basis_vectors_postspike::VMF
	"temporal basis vectors for the post-stereoclick kernel"
	temporal_basis_vectors_poststereoclick::VMF
	"temporal basis vectors for the pre-movement kernel"
	temporal_basis_vectors_premovement::VMF
	"parameters concatenated into a vector"
	parametervalues::VF
	"name of each parameter"
	parameternames::VS
	"a vector of L2 penalty matrices"
	penaltymatrices::VMF
	"index of the parameters regularized by the L2 penalty matrices"
	penaltymatrixindices::VVI
	"cofficients of the penalty matrices"
	penaltycoefficients::VF
	"names of each L2 regularization penalty"
	penaltynames::VS
	"precision matrix of the gaussian prior on the parameters"
	precisionmatrix::MF
	"hessian of the log-likelihood function evaluated at the current parameters"
	hessian_loglikelihood::MF = fill(NaN, length(parametervalues), length(parametervalues))
	"hessian of the log-posterior function evaluated at the current parameters and hyperparameters"
	hessian_logposterior::MF = fill(NaN, length(parametervalues), length(parametervalues))
end

"""
	SpikeTrainLinearFilter

Linear filter used to smooth the spike train
"""
@with_kw struct SpikeTrainLinearFilter{VF<:Vector{<:AbstractFloat}}
	"the function with which the spike train is convolved"
	impulseresponse::VF
	"because of lack of spike train responses before the beginning of the trial, the initial response needs to reweighed"
	weights::VF
end

"""
	PerieventTimeHistogram

The mean across trials of a single condition (e.g. trials that ended with a left choice) of the filtered spike train of one neuron, and the estimated 95% confidence interval of the trial mean
"""
@with_kw struct PerieventTimeHistogram{CIM<:Bootstrap.ConfIntMethod,
									S<:String,
									STLF<:SpikeTrainLinearFilter,
									VF<:Vector{<:AbstractFloat}}
	"method used to compute the confidence interval. The default is the bias-corrected and accelerated confidence interval (Efron & Tibshirani, 1993) for a confidence level of 0.95. The confidence level is the fraction of time when a random confidence interval constructed using the method below contains the true peri-stimulus time histogram."
	confidence_interval_method::CIM = BCaConfInt(0.95)
	"the inear filter used to smooth the spike train"
	linearfilter::STLF
	"estimate of the lower limit of the confidence interval of the peri-event time histogram based on observed spike trains"
	lowerconfidencelimit::VF
	"estimate of the peri-event time histogram based on observed spike trains"
	observed::VF
	"estimate of the peri-event time histogram based on simulated spike trains"
	predicted::VF
	"An event in the trial (e.g. steroclick) corresponding to which the peri-stimulus time histogram is aligned, i.e., when time=0 is defined."
	referenceevent::S="stereoclick"
	"time, in seconds, from the reference event"
	time_s::VF
	"estimate of the upper limit of the confidence interval of the peri-event time histogram based on observed spike trains"
	upperconfidencelimit::VF
end

"""
	PETHSet

A set of peri-event time histogram of one neuron.
"""
@with_kw struct PETHSet{PETH<:PerieventTimeHistogram}
	"average across trials that ended in a left choice, including both correct and incorrect trials"
	leftchoice::PETH
	"average across trials on which the aggregate evidence favored left and the reward is baited on the left"
	leftevidence::PETH
	"average across trials on which the animal made a left choice and the evidence was strongly leftward. Therefore, only correct trials are included. Evidence strength is defined by the generative log-ratio of click rates: `γ` ≡ log(right click rate) - log(left click rate). Strong left evidence evidence include trials for which γ < -2.25"
	leftchoice_strong_leftevidence::PETH
	"average across trials on which the animal made a left choice and the generatative `γ` < -2.25 && `γ` < 0"
	leftchoice_weak_leftevidence::PETH
	rightchoice::PETH
	rightevidence::PETH
	rightchoice_strong_rightevidence::PETH
	rightchoice_weak_rightevidence::PETH
	"average across all trials"
	unconditioned::PETH
end

"""
	CVResults

Results of cross-validation
"""
@with_kw struct CVResults{C<:Characterization, VC<:Vector{<:CVIndices}, VS<:Vector{<:ModelSummary}, VVP<:Vector{<:Vector{<:PETHSet}}}
	"a composite containing quantities that are computed out-of-sample and used to characterize the model`"
	characterization::C
	"cvindices[k] indexes the trials and timesteps used for training and testing in the k-th resampling"
	cvindices::VC
	"post-stereoclick time histogram sets"
	psthsets::VVP
	"summaries of the training models"
	trainingsummaries::VS
end


"""
	PoissonGLM

Object containing the data and parameters of a Poisson GLM and also the quantities for estimating its parmaeters
"""
@with_kw struct PoissonGLM{R<:Real, VR<:Vector{<:Real}, MR<:Matrix{<:Real}, TI<:Integer, VI<:Vector{<:Integer}}
	"design matrix"
	𝐗::MR
	"Poisson observations"
	𝐲::VI
	"number of weights"
	n𝐰::TI = size(𝐗,2)
	"weights"
	𝐰::VR=rand(n𝐰)
	"time step"
	Δt::R=0.01
	"log-likelihood"
	ℓ::VR=fill(NaN,1)
	"gradient"
	∇ℓ::VR=fill(NaN,n𝐰)
	"hessian"
	∇∇ℓ::MR=fill(NaN,n𝐰,n𝐰)
end

"""
	GainMixtureGLM

Object containing the data and parameters of a mixture-of-gain Poisson GLM and also the quantities for estimating its parmaeters
"""
@with_kw struct GainMixtureGLM{R<:Real, VR<:Vector{<:Real}, MR<:Matrix{<:Real}, VVR<:Vector{<:Vector{<:Real}}, TI<:Integer, VI<:Vector{<:Integer}}
	"time step"
	Δt::R
	"design matrix. the first column corresponds to the gain"
	𝐗::MR
	"spike train"
	𝐲::VI
	"slowly-varying baseline"
	𝐆::VR = 𝐗[:,1]
	"weight of the baseline. Each element corresponds to the weight of the baseline at different states"
	𝐠::VR = fill(NaN,2)
	"other inputs"
	𝐔::MR = 𝐗[:,2:end]
	"weight of the other inputs"
	𝐮::VR = fill(NaN, size(𝐗,2)-1)
	"number of baseline-related parameters"
	n𝐠::TI = length(𝐠)
	"number of weights for other inputs "
	n𝐮::TI = length(𝐮)
	"total number of weights"
	nweights::TI=n𝐠+n𝐮
	"number of time steps"
	ntimesteps::TI=length(𝐲)
	"posterior probability of the baseline state"
	𝛄::VVR=collect(fill(NaN, ntimesteps) for k=1:n𝐠)
	"prior probability of the baseline state"
	π::VR=fill(NaN,1)
	"expectation of the conditional log-likelihood under the posterior probability of the baseline state"
	Q::VR=fill(NaN, 1)
	"gradient of the expectation of the conditional log-likelihood"
	∇Q::VR=fill(NaN, nweights)
	"hessian of the expectation of the conditional log-likelihood"
	∇∇Q::MR=fill(NaN, nweights, nweights)
	"product between the posterior probability in each state and the derivative of the conditional log-likelihood with respect to the conditional linear predictor. Element `γₖdℓ_dLₖ[k][t]` corresponds to the k-th state and the t-th time step"
	γₖdℓ_dLₖ::VVR=collect(fill(NaN, ntimesteps) for k=1:n𝐠)
	"product between the posterior probability in each state and the second derivative of the conditional log-likelihood with respect to the conditional linear predictor. Element `γₖd²ℓ_dLₖ²[k][t]` corresponds to the k-th state and the t-th time step"
	γₖd²ℓ_dLₖ²::VVR=collect(fill(NaN, ntimesteps) for k=1:n𝐠)
end
