##############################################################################
##############################################################################
##############################################################################
##
## DESCRIPTION: This example shows how to calculate the Bayesian sensitivity for a given signal process
## using the SNSensitivityEstimate package. 
##############################################################################
##############################################################################
##############################################################################


import Pkg;
Pkg.activate("/sps/nemo/scratch/mpetro/Sensitivity_training/")
using SNSensitivityEstimate, UnROOT, FHist, CairoMakie, CSV, DataFramesMeta
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, SpecialFunctions, ValueShapes
using Distributions 

# STEP 1: Load the data

# Start is the same as in ex1_background.jl, we load data and create DataProcess objects. 
data_info = CSV.read("data/data_info.csv", DataFrame)
file_paths = joinpath.("data", data_info.file_name) 
data_files = [UnROOT.ROOTFile(file_path) for file_path in file_paths] 
tree_name = "tree" 
variables = keys(data_files[1][tree_name]) 
data_tables = [UnROOT.LazyTree(data_file, tree_name, variables) for data_file in data_files] 


# STEP 1.1: apply data-cuts to the data before creating DataProcess objects.
# In ex1, we did not impose any data-cuts on the data.
# If we want somewhat meaningful results on sensitivity, we must first filter data to some meaningful data-cuts:

# First we define variables we want to cut
var_names = [ 
    "sumE", 
    "dy", 
    "dz",
    "Pint",
    "Pext",
    ]

# Next we define the boundaries of the data-cuts, we use simple Esum >300keV, dy < 100mm, dz < 100mm, Pint > 0.01, Pext < 0.01
var_bounds = (
    sumE = (300, 3500),
    dy = (0, 100),
    dz = (0, 100),
    Pint = (0.01, 1),
    Pext = (0, 0.01),
)

# next we use filter_data! function from SNSensitivityEstimate to filter the data before loading it to DataProcess.
filtered_data_tables = [SNSensitivityEstimate.filter_data(table, var_names, var_bounds) for table in data_tables]


processes = [
    SNSensitivityEstimate.DataProcess(
        collect(filtered_data_tables[i].sumE), # Here we specify a single variable to be used. In later stages we will use N-dimensional data vectors.
        String(data_info.process[i]), 
        data_info.is_signal[i], 
        data_info.activity[i], 
        data_info.time_s[i], 
        data_info.n_sim[i], 
        0:100:3500, 
        data_info.amount[i], 
    ) for i in 1:length(filtered_data_tables)
]


######
# (Optional) you can redo the background example with this data-cuts, or try some of your own and see what changes. 
#####


# STEP 2: create the normalize histograms for the signal and background processes. These are inputs to the bayesian inference process.

# declare which process is signal and obtain histogram
signal = SNSensitivityEstimate.get_process("bb0nu_foil_bulk", processes)[1]
signal_hist = SNSensitivityEstimate.get_bkg_counts_1D(signal)


# declare background processes and obtain histograms
bkg_names = [ p.isotopeName for p in processes if !p.signal ]
background = [SNSensitivityEstimate.get_process(name, processes)[1] for name in bkg_names]
bkg_hist = [SNSensitivityEstimate.get_bkg_counts_1D(p) for p in background]

# we must normalize the histograms for the fitting procedure
bkg_hist_normed = FHist.normalize.(bkg_hist, width = true)
signal_hist_normed = FHist.normalize(signal_hist, width = true)

# STEP 2.1: create pseudo-data histogram which will be fitted in the procedure.

# First we sample pseudo-data from the background histograms by sampling with Poisson distribution with mean equal to the expected number of counts in each bin.
# This is done using the SNSensitivityEstimate.get_pseudo_spectrum function, which takes a histogram and returns a pseudo-data histogram.

# sample each process individually
sample_hists = [SNSensitivityEstimate.get_pseudo_spectrum(b) for b in bkg_hist] 
# merge them together into a single pseudo-experiment
data_hist = merge(sample_hists...)

let
    f = CairoMakie.Figure(size = (800, 600))
    a = CairoMakie.Axis(f[1, 1], title = "Pseudo-data vs background model", xlabel = "sum of energies", ylabel = "counts")
    CairoMakie.stephist!(a, sum(bkg_hist), label = "background model", color = :red)
    CairoMakie.stephist!(a, data_hist, label = "pseudo-data", color = :blue)
    CairoMakie.axislegend(a, position = :lt)
    f
end

# STEP 3: define prior, likelihood, and posterior distributions for the Bayesian inference process.

# STEP 3.1: prior
# uninformed prior for each activity
prior = ValueShapes.NamedTupleDist(
    As = Distributions.Uniform(1e-20, 1), # one prior for signal - flat prior between 0 and 1
    Ab = [Distributions.Uniform(1e-20,1) for _ in 1:length(bkg_hist)] # same with background priors
)   

# STEP 3.2: likelihood
# likelihood built from normalized histograms
# The likelihood is defined in section 5.4.3 of my phd thesis (will upload to docDB when possible)
# The fit parameters are the relative ratios of each process's contribution to the total data_hist
# A_s = relative amount of signal
# A_b = vector of relative amounts of background processes
likelihood = SNSensitivityEstimate.make_hist_likelihood_uniform(
    data_hist,
    signal_hist_normed,
    bkg_hist_normed
)

# STEP 3.3: posterior
# BAT provides a nice interface to define posterior, simply call BAT.PosteriorMeasure with the likelihood and prior as arguments.
posterior = BAT.PosteriorMeasure(likelihood, prior)



# STEP 3.4: BAT setup
# There's a few things BAT related that need to be set-up in the MCMC sampling process
# The documentation can be found in https://bat.github.io/BAT.jl/stable/
# Here we just use some basic defaults.
burnin = BAT.MCMCMultiCycleBurnin(max_ncycles = 30, nsteps_final=1000)
mcmcalgo = BAT.RandomWalk()

# STEP 4: run the MCMC sampling, obtain the posterior samples and visualize the results
# sampling is done using the BAT.bat_sample function, which takes the posterior and the MCMC sampling algorithm as arguments.
# the nsteps is the number of steps to run the MCMC sampling, and nchains is the number of chains to run in parallel.
samples, _ = BAT.bat_sample( posterior, BAT.MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 10_000, nchains = 4))

# We can print some basic results from the samples:
# For example mode of samples gives the most likely values for each fitted parameter. 
println("Mode: $(mode(samples))")
println("Mean: $(mean(samples))")
println("Stddev: $(std(samples))")

binned_unshaped_samples, _ = BAT.bat_transform(Vector, samples)

# We can look at the posterior distribution of each parameter by plotting the histogram of the samples.
let
    f = CairoMakie.Figure(size = (800, 1000))
    a1 = CairoMakie.Axis(f[1, 1], title = "Posterior distribution of signal", xlabel = "signal ratio", ylabel = "counts")
    d = [par[1] for par in binned_unshaped_samples.v]
    h = FHist.Hist1D(d; binedges= range(0, maximum(d), length = 100))
    CairoMakie.stephist!(a1, h, label = "A_s", color = :blue)
    CairoMakie.axislegend(a1, position = :ct)

    labels = data_info.process[.!data_info.is_signal]
    colors = Makie.wong_colors()
    for i in 1:length(bkg_hist)
        a2 = CairoMakie.Axis(f[i+1, 1], title = "Posterior distribution of $(labels[i])", xlabel = "background ratio", ylabel = "counts")       
        d = [par[i+1] for par in binned_unshaped_samples.v]
        h = FHist.Hist1D(d; binedges= range(0, maximum(d), length = 100))
        CairoMakie.stephist!(a2, h, label = "A_b$(i)", color = colors[i])
        CairoMakie.axislegend(a2, position = :ct)
    end
    save("data/out/ex2/posterior_distributions.png", f, px_per_unit = 2)
    f
end

# STEP 5: Calculate the sensitivity of single pseudo-experiment using the posterior samples.
# The sensitivity is defined as the 90% quantile of the posterior distribution of the signal
# The sensitivity is calculated using the SNSensitivityEstimate.calculate_sensitivity function, which takes the posterior samples and the signal histogram as arguments.
# The sensitivity is returned in terms of the expected number of signal events, which can be converted to a half-life using the formula:
# t1/2 = log(2) * (Nₐ * m * t * eff / W) / exp_mu_signal_90
# where Nₐ is Avogadro's number, m is the mass of the isotope, t is the live-time of the experiment, eff is the efficiency of the signal process, and W is the molar mass of the isotope.



nDataPoints = integral(data_hist)
muS = [par[1] for par in binned_unshaped_samples.v]
exp_mu_signal_90 = BAT.smallest_credible_intervals(muS, nsigma_equivalent=1.65)[1].right * nDataPoints # this gives the smallest credible interval

Na = SNparams["Nₐ"]
m = SNparams["foilMass"] * SNparams["a"]
t = SNparams["tYear"]
W = SNparams["W"]

ROI_a, ROI_b = 300, 3500-50
eff = lookup(signal, ROI_a, ROI_b) 
t12 = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

println("Sensitivity (90% CI) in terms of half-life: $(t12) years)")
println("for a signal efficiency of $(eff)")
println("in the ROI of $(ROI_a) to $(ROI_b) keV.")
println("equivalent to an expected number of signal events of $(round(exp_mu_signal_90, digits= 3))")


# (optional) STEP 6: plot the results of the fit to visually inspect the results

amps = mean(samples) # get the amplitudes of each process from the posterior samples, which are the mean of the samples for each parameter.

let 
    colors = Makie.wong_colors()
    bkg_amps = amps[2] ./ sum(amps[2]) # normalize the background amplitudes to sum to 1
    labels = data_info.process[.!data_info.is_signal]
    f = Figure(size = (1200,800), fontsize = 26)
    a = Axis(f[1,1], xlabel = "sum E", ylabel = "Counts / 17.5kg.yr", title = "Fitted spectrum of a single pseudo-experiment", limits = (300,3500, 0, nothing) )
    h_total = sum(bkg_hist)
    n_total = integral(h_total)
    CairoMakie.stephist!(a, h_total, color = :black, label = "Total pseudo-background", linewidth = 4)
    
    h_fit = [ a*normalize(h, width = false)*n_total for (h, a) in zip(bkg_hist, bkg_amps)]

    for i in 1:length(h_fit[1:end])
        CairoMakie.stephist!(a, h_fit[i], color = colors[i], label = labels[i], linewidth = 4)
    end
    CairoMakie.hist!(a, sum(h_fit), label = "Fitted total", color = (colors[1], 0.4), )
    
    axislegend(a; position = :rt, patchsize = (45,35), labelsize = 20)
    save("data/out/ex2/fit_result.png", f, px_per_unit = 2)
    f
end

# (For later) STEP 7: The provided process gives the sensitvity for a single pseudo-experiment. 
# Since MCMC sampling is a stochastic process, the sensitivity obtained from a single pseudo-experiment is not representative of the true sensitivity of the experiment.
# To obtain the sensitivity of the experiment which we can quote in publications, we must repeat the process for many pseudo-experiments 
# and take the median of the sensitivities obtained. 
# This can be done using the SNSensitivityEstimate.get_sens_bayes_uniform function
# Which encapsulates the entire process of generating pseudo-data, running MCMC sampling, and calculating the sensitivity for a given number of pseudo-experiments.
# We repeat this process many times and extract median sensitivity. 

# this will take a while...
t_halfs = []
n_max = 10 # number of pseudo-experiments to run
for i in 1:n_max

    try 
        println("Running pseudo-experiment $(i) of $(n_max)")
        sens = SNSensitivityEstimate.get_sens_bayes_uniform(
            bkg_hist, 
            signal, 
            prior; 
            ROI_a = ROI_a, ROI_b = ROI_b, nsteps = 5*10^4, nchains = 4)
        println(sens)
        push!(t_halfs, sens)
    catch
        @warn "failed fit" 
        continue
    end
end

# Let's plot the distribution and show median sensitivity
let
    f = CairoMakie.Figure(size = (800, 600))
    a = CairoMakie.Axis(f[1, 1], title = "Distribution of sensitivities from $(n_max) pseudo-experiments", xlabel = "sensitivity (years)", ylabel = "counts")
    h = FHist.Hist1D(t_halfs; binedges= range(0, maximum(t_halfs), length = n_max))
    CairoMakie.stephist!(a, h, label = "sensitivity distribution", color = :blue)

    median_sens = StatsBase.median(t_halfs)
    CairoMakie.vlines!(a, [median_sens], color = :red, label = "median sensitivity: $(round(median_sens, sigdigits= 3)) years")
    CairoMakie.axislegend(a, position = :ct)
    save("data/out/ex2/sensitivity_distribution.png", f, px_per_unit = 2)
    f
end

