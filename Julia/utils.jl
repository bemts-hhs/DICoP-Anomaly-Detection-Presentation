"""
    simulate_and_plot(dist_type::String; n=500, seed=10232015, xlab::String="Value", ylab::String="Density")

Generate simulated data from a specified probability distribution and produce
a kernel density plot with a filled region and a vertical line at the theoretical mean. The function is intended for illustrating distributional characteristics for anomaly detection and related epidemiologic applications.

# Arguments
- `dist_type::String`: One of `"normal"`, `"poisson"`, `"quasipoisson"`, `"nb"`.
- `n::Int`: Number of observations to simulate. Default is 500.
- `seed::Int`: Random seed for reproducibility. Default is 10232015.

# Details
The function generates synthetic data according to the selected distribution,
computes the distribution-defined mean, and overlays a vertical reference line
at that value. This assists with demonstrating how distributional assumptions
relate to expected behavior in count-based epidemiologic data.

# Output
Returns a `Plots.Plot` object and saves a `.png` file to:
`./plots/distribution_plots/`.

"""
function simulate_and_plot(dist_type::String; n=500, seed=10232015)

    # Reproducibility
    Random.seed!(seed)

    # Ensure output directory exists
    outdir = "./plots/distribution_plots"
    isdir(outdir) || mkpath(outdir)

###_____________________________________________________________________________
# Select distribution, simulate data, assign theoretical mean, and filename
###_____________________________________________________________________________
    
    # Define a constant for plot fill color
    DIST_COLORS = Dict(
    "normal" => :gold,
    "poisson" => :dodgerblue,
    "quasipoisson" => :forestgreen,
    "nb" => :crimson
    )

    # Subset the active color
    color = DIST_COLORS[dist_type]

    # Define the distribution names
    DIST_NAMES = Dict(
        "normal" => "Normal (μ, σ²) Distribution",
        "poisson" => "Poisson (μ = λ = σ²) Distribution",
        "quasipoisson" => "Gamma-Poisson Mixture (σ² > λ) Simulation",
        "nb" => "Negative Binomial (σ² > μ) Distribution"
    )

    # Subset distribution name
    dist_name = DIST_NAMES[dist_type]

    # Initiate the loop to create distributions, filenames
    if dist_type == "normal"
        mu = 0
        sigma = 1
        data = rand(Normal(mu, sigma), n)
        fname = "normal.png"

    elseif dist_type == "poisson"
        lambda = 25
        mu = lambda
        # Must simulate Poisson, not assign the distribution
        data = rand(Poisson(lambda), n)
        fname = "poisson.png"

    elseif dist_type == "quasipoisson"
        lambda = 25
        phi = 10.0
        # Strong overdispersion for visible right‑skew
        rate = rand(Gamma(lambda / phi, phi), n)
        data = rand.(Poisson.(rate))
        mu = lambda
        fname = "quasipoisson.png"

    elseif dist_type == "nb"
        mu = 25
        var = 300               # ensures strong tail
        r = mu^2 / (var - mu)
        p = mu / var
        data = rand(NegativeBinomial(r, p), n)
        fname = "negative_binomial.png"

    else
        error("Unknown distribution type: $dist_type")
    end

    
###_____________________________________________________________________________
# Create the histogram with filled area and vertical line at theoretical mean
###_____________________________________________________________________________
    
    plt_hist = histogram(
        data,
        fill=(0, 0.5, color),
        xlabel = "Value",
        ylabel = "Frequency",
        legend = false
    )
    
###_____________________________________________________________________________
# Create the density plot with filled area and vertical line at theoretical mean
###_____________________________________________________________________________

    plt_density = density(
        data,
        fill=(0, 0.5, color),
        xlabel = "Value",
        ylabel = "Density",
        legend = false
    )

    # Add vertical line at theoretical mean
    vline!([mu], linecolor=:black, linewidth=2)

###_____________________________________________________________________________
# Combine the histogram and density plots
###_____________________________________________________________________________

    hist_density = plot(
        plt_hist, plt_density, layout=2, size=(800, 800 / (16/9)),
        title = ["Histogram" "Kernel Density Plot"]
        )
    
    return hist_density
end