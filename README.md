# tzl_spGLM
Generalized linear model of the spike train of a neuron recorded in the Poisson Clicks task

## example

### Julia REPL
The optimization procedure learns both the parameters (i.e. regressor weights) and a hyperparamer (precision of the Gaussian prior on the parameters, analogous to a L2 penalty coefficient):
```
julia> using SPGLM
julia> csvpath = "/mnt/cup/people/zhihaol/Documents/tzluo/analyses/analysis_2024_04_22a_test_SPGLM/models.csv";
julia> options = SPGLM.Options(csvpath);
julia> trials = SPGLM.loadtrials(options);
julia> model = SPGLM.Model(options, trials);
julia> eo = SPGLM.maximizeevidence!(model)

Iter     Function value   Gradient norm
     0     3.553781e+04     1.285527e+03
 * time: 6.103515625e-5
     1     3.340892e+04     1.124030e+03
 * time: 0.05215096473693848
     2     3.038869e+04     7.491107e+02
 * time: 0.10418510437011719
     3     2.782119e+04     3.308262e+02
 * time: 0.15614914894104004
     4     2.628437e+04     1.192378e+02
 * time: 0.21809101104736328
     5     2.546613e+04     2.909452e+01
 * time: 0.25608205795288086
     6     2.535228e+04     3.249330e+00
 * time: 0.29502296447753906
     7     2.535107e+04     4.872618e-02
 * time: 0.33403897285461426
     8     2.535107e+04     1.218098e-05
 * time: 0.38610410690307617
     9     2.535107e+04     1.225686e-12
 * time: 0.44210100173950195
Evidence optimization iteration: 1: precision (a) = 0.5908446386657102
Evidence optimization iteration: 1: norm of the residual gradient of the log-posterior (g_residual) = 1.2256862191861728e-12
Evidence optimization iteration: 1: MAP optimization converged
Evidence optimization iteration: 1: approximate log-evidence (𝐸) = -25382.187862691735
Iter     Function value   Gradient norm
     0     2.493416e+04     1.770758e+01
 * time: 6.198883056640625e-5
     1     2.491330e+04     1.655549e+01
 * time: 0.05198216438293457
     2     2.487574e+04     1.496841e+01
 * time: 0.1039891242980957
     3     2.481107e+04     1.249717e+01
 * time: 0.15512514114379883
     4     2.471561e+04     8.400673e+00
 * time: 0.20711398124694824
     5     2.462951e+04     2.547959e+00
 * time: 0.2629849910736084
     6     2.461760e+04     4.920018e-01
 * time: 0.30104517936706543
     7     2.461753e+04     2.734311e-03
 * time: 0.3390049934387207
     8     2.461753e+04     1.063443e-07
 * time: 0.369035005569458
     9     2.461753e+04     2.953931e-14
 * time: 0.4078941345214844
Evidence optimization iteration: 2: precision (a) = 0.006298120303099656
Evidence optimization iteration: 2: norm of the residual gradient of the log-posterior (g_residual) = 2.9539305029802065e-14
Evidence optimization iteration: 2: MAP optimization converged
Evidence optimization iteration: 2: approximate log-evidence (𝐸) = -24707.60295639283
Iter     Function value   Gradient norm
     0     2.461416e+04     8.121148e-02
 * time: 6.29425048828125e-5
     1     2.461414e+04     1.074387e-03
 * time: 0.05279088020324707
     2     2.461414e+04     1.031625e-08
 * time: 0.10295701026916504
     3     2.461414e+04     8.031770e-15
 * time: 0.1548769474029541
Evidence optimization iteration: 3: precision (a) = 0.0051224104379705695
Evidence optimization iteration: 3: norm of the residual gradient of the log-posterior (g_residual) = 8.031769693772617e-15
Evidence optimization iteration: 3: MAP optimization converged
Evidence optimization iteration: 3: approximate log-evidence (𝐸) = -24707.271175121816
SPGLM.EvidenceOptimization{Vector{Float64}, Vector{Vector{Float64}}}
  a: Array{Float64}((3,)) [0.5908446386657102, 0.006298120303099656, 0.0051224104379705695]
  𝐸: Array{Float64}((3,)) [-25382.187862691735, -24707.60295639283, -24707.271175121816]
  MAP_g_residual: Array{Float64}((3,)) [1.2256862191861728e-12, 2.9539305029802065e-14, 8.031769693772617e-15]
  𝐰: Array{Vector{Float64}}((3,))

julia> characterization = SPGLM.Characterization(model)
julia> SPGLM.save(model)
julia> SPGLM.save(eo, model.options.outputpath)
julia> SPGLM.save(characterization, model.options.outputpath)
```

### MATLAB command window
Let's check whether the model actually fit the data, using utilities in [+SPGLM](/src/+SPGLM/) and [+TZL](https://github.com/Brody-Lab/tzluo/tree/master/%2BTZL)
```
>> analysispath = 'V:\Documents\tzluo\analyses\analysis_2024_04_23a_test_SPGLM'
>> models = readtable(fullfile(analysispath, "models.csv"), 'delimiter', ',');
>> models.outputpath = TZL.cup2windows(string(models.outputpath));
>> load(fullfile(models.outputpath{1}, 'model.mat'))
>> load(fullfile(models.outputpath{1}, 'characterization.mat'))
>> peths = SPGLM.tabulatepeths(peths);
>> kernels = SPGLM.tabulatekernels(kernels);
```
Plotting the peri-event time histograms (PETH), aligned to the rat's movement onset, and conditioned on whether it made a left or a right choice:

```
>> reference_event = "movement";
>> conditions = ["leftchoice", "rightchoice"];
>> SPGLM.stylizeaxes
>> colors = get(gca, 'colororder');
>> h = [];
>> k = 0;
>> for i = 1:numel(conditions)
    index = peths.condition == conditions(i) & peths.reference_event == reference_event;
    k =k + 1;
    h(k) = plot(peths.timesteps_s{index}, peths.observed{index}, '-', 'color', colors(i,:));
    k =k + 1;
    h(k) = plot(peths.timesteps_s{index}, peths.predicted{index}, '-', ...
        'linewidth', 1, 'color', colors(i,:));
end
>> xlabel(['time from ' char(reference_event) ' (s)'])
>> ylabel('spikes/s')
>> ylim(ylim.*[0,1])
>> legend(h, {'left choice, observed', 'left choice, predicted', ...
    'right choice, observed', 'right choice, predicted'}, 'location', 'northwest')
```

<img src="/example/analysis_2024_04_22a_test_SPGLM/peth_perimovement.svg" width = "400">

The model captures the choice-dependent peri-movement neuronal activity by treating the moment of movement onset as a transient impulse and fits convolution filters aligned to the impulse. Separate convolution kernels are fit for different choices. 

<img src="/example/analysis_2024_04_22a_test_SPGLM/kernel_perimovement.svg" width = "400">

Using similar code as above to plot the PETH aligned to the onset of auditory click trains, which is always preceded by a left and right click being played simultaneously (i.e., a so-called "stereoclick"). The meaning of the color and the line thickness are as above

```
>> reference_event = "stereoclick";
```

<img src="/example/analysis_2024_04_22a_test_SPGLM/peth_stereoclick.svg" width = "400">

Or aligned to the onset of nose fixation

```
>> reference_event = "cpoke_in";
```

<img src="/example/analysis_2024_04_22a_test_SPGLM/peth_perifixation.svg" width = "400">

We can also check whether we capture the response aligned to either a left or a right click (also called the "click-triggered average").

<img src="/example/analysis_2024_04_22a_test_SPGLM/peth_click.svg" width = "400">