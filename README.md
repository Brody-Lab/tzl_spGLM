# tzl_spGLM
Generalized linear model of the spike train of a neuron recorded in the Poisson Clicks task

## example
```
csvpath =
"/mnt/cup/people/zhihaol/Documents/tzluo/analyses/analysis_2024_04_17a_test_SPGLM/models.csv"
options = SPGLM.Options(csvpath
trials = SPGLM.loadtrials(options)
model = SPGLM.Model(options, trials)
𝐰₀ = rand(length(model.𝐰))
a₀ = rand()
model.𝐰 .= 𝐰₀
model.a .= a₀
SPGLM.maximizeposterior!(model)
```