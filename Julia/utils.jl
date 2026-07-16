###_____________________________________________________________________________
# Custom helper functions
###_____________________________________________________________________________

# -------------------------------------------------------------------------
# nb_pois_pred_interval
#
# Computes a two‑tailed prediction interval for annual registry counts by
# selecting either a Negative Binomial distribution (for overdispersed data)
# or a Poisson distribution (for equidispersed or mildly underdispersed
# data). The model selection is driven entirely by the empirical mean and
# variance of the supplied counts vector.
#
# Arguments
# ----------
# counts_raw
#     Inbound object containing registry counts. TidierData.jl may supply
#     this as many different structures, including:
#         • Vector{<:Real}           – standard rowwise vector
#         • NamedTuple               – from across(), one value per year
#         • Number                   – degenerate scalar case
#         • Vector{Vector{<:Real}}   – entire column of rowwise vectors
#
# upper_prob :: Real
#     Upper cumulative probability for the prediction interval. The lower
#     bound is computed as (1 - upper_prob). Example: 0.9332 for the
#     1.5‑sigma equivalent.
#
# Returns
# ----------
# (lower, upper) :: Tuple{Float64,Float64} OR Vector of such tuples
#     Returns lower and upper prediction interval bounds. If counts_raw is
#     a column‑level Vector{Vector}, the function returns a vector of row‑
#     level interval tuples.
#
# Interval Logic
# ----------
# Negative Binomial parameters follow the NB2 mean‑variance inversion:
#     r = mu^2 / (var - mu)
#     p = mu / var
#
# If var <= mu, NB is non‑identifiable; Poisson(mu) is used instead.
#
# Integration Notes
# ----------
# For TidierData.jl, generate per‑row vectors of year counts using:
#
#     counts_vec = c(`2020`, `2021`, ..., `2026`)
#
# Then call:
#
#     pred_interval = Main.nb_pois_pred_interval(counts_vec, 0.9332)
#
# -------------------------------------------------------------------------

function nb_pois_pred_interval(counts_raw, upper_prob::Real)

    # ------------------------------------------------------------------
    # Coerce ANY inbound structure into a usable Vector{Float64} or,
    # in the column‑level case, produce a vector of interval tuples.
    #
    # The ternary structure below is explicit. Each condition converts
    # counts_raw into a consistent representation for interval logic.
    # ------------------------------------------------------------------

    counts =
        # Case 1: Standard rowwise vector of numeric values
        counts_raw isa AbstractVector{<:Real} ?
            Float64.(counts_raw) :

        # Case 2: NamedTuple from across(), e.g. (2020 = 123, ...)
        counts_raw isa NamedTuple ?
            Float64.(collect(values(counts_raw))) :

        # Case 3: Degenerate scalar case, wrap into a 1‑element vector
        counts_raw isa Number ?
            [Float64(counts_raw)] :

        # Case 4: Full column supplied as Vector{Vector}. Apply the
        # interval logic rowwise and return a vector of interval tuples.
        counts_raw isa AbstractVector{<:AbstractVector} ?
            [begin
                # Extract row vector and coerce its elements
                row = Float64.(counts_raw[i])

                # Compute summary statistics for this row
                mu = mean(row)
                var_row = Statistics.var(row)
                lower_prob = 1 - upper_prob

                # Select distribution based on dispersion pattern
                dist =
                    var_row > mu ?
                        Distributions.NegativeBinomial(
                            (mu^2) / (var_row - mu),
                            mu / var_row
                        ) :
                        Distributions.Poisson(mu)

                # Return rowwise interval tuple
                (
                    Statistics.quantile(dist, lower_prob),
                    Statistics.quantile(dist, upper_prob)
                )
            end for i in eachindex(counts_raw)] :

        # Case 5: Unsupported inbound type
        error("Unsupported counts structure: $(typeof(counts_raw))")

    # ------------------------------------------------------------------
    # At this point, counts is guaranteed to be Vector{Float64}. Compute
    # interval normally for the row‑level case (not the column‑level case,
    # which already returned inside the mapping block).
    # ------------------------------------------------------------------

    mu = Statistics.mean(counts)
    var_counts = Statistics.var(counts)
    lower_prob = 1 - upper_prob

    # Check dispersion pattern to select distribution
    if var_counts > mu

        # Compute NB parameters for overdispersed counts
        r = (mu^2) / (var_counts - mu)
        p = mu / var_counts
        dist = Distributions.NegativeBinomial(r, p)

    else

        # Use Poisson for equidispersed or underdispersed counts
        dist = Distributions.Poisson(mu)
    end

    # Return prediction interval tuple
    lower = Statistics.quantile(dist, lower_prob)
    upper = Statistics.quantile(dist, upper_prob)

    return lower, upper
end

###############################################################################
# simulate_and_plot
#
# Description:
#     Generates simulated data from a specified probability distribution and
#     produces a kernel density plot with a filled density region and a vertical
#     line at the theoretical mean of the distribution. The function is intended
#     for educational and demonstration purposes in the DICoP anomaly detection
#     project, where foundational distributional behavior must be illustrated.
#
# Arguments:
#     dist_type::String
#         A string indicating the target distribution. Supported values include:
#             • "normal"       – standard Gaussian distribution
#             • "poisson"      – equidispersed Poisson count model
#             • "quasipoisson" – overdispersed Poisson (gamma–Poisson mixture)
#             • "nb"           – negative binomial for strong overdispersion
#
# Optional Arguments:
#     n::Integer (default = 500)
#         Number of observations to simulate.
#
#     seed::Integer (default = 10232015)
#         Seed used for reproducibility.
#
# Outputs:
#     Returns the plot object and also writes a .png file to:
#         ./plots/distribution_plots/
#
# Notes:
#     • For quasi-Poisson, a gamma–Poisson mixture is used because quasi-Poisson
#       is a quasi-likelihood model and does not define a proper probability
#       distribution with a pmf or quantile function.
#     • For negative binomial, parameters are derived using the NB2 mean–variance
#       relationship appropriate for count data with substantial overdispersion.
###############################################################################

function simulate_and_plot(dist_type::String; n=500, seed=10232015)

    # Ensure reproducibility
    Random.seed!(seed)

    # Create output directory if missing
    outdir = "./plots/distribution_plots"
    isdir(outdir) || mkpath(outdir)

    # Select distribution and simulate data
    if dist_type == "normal"
        # Mean and standard deviation for Gaussian data
        mu = 0
        sigma = 1
        
        # Simulated observations
        data = rand(Normal(mu, sigma), n)
        
        # Output filename
        fname = "normal.png"

    elseif dist_type == "poisson"
        # Rate parameter for Poisson process
        lambda = 25
        mu = lambda    # Poisson mean equals lambda
        
        # Simulated observations
        data = rand(Poisson(lambda), n)
        
        # Output filename
        fname = "poisson.png"

    elseif dist_type == "quasipoisson"
        # Quasi-Poisson mean
        lambda = 25
        phi = 10.0      # Overdispersion factor
        
        # Gamma-Poisson mixture to emulate quasi-Poisson behavior
        # Gamma controls rate variability; Poisson applies to generated rates
        rate = rand(Gamma(lambda / phi, phi), n)
        data = rand.(Poisson.(rate))
        
        # Theoretical mean remains lambda
        mu = lambda
        
        # Output filename
        fname = "quasipoisson.png"

    elseif dist_type == "nb"
        # Negative binomial mean and variance for overdispersed counts
        mu = 25
        var = 300
        
        # NB2 parameterization
        r = mu^2 / (var - mu)   # size parameter
        p = mu / var            # success probability
        
        # Simulated observations
        data = rand(NegativeBinomial(r, p), n)
        
        # Output filename
        fname = "negative_binomial.png"

    else
        # Defensive coding for unsupported arguments
        error("Unknown distribution type: $dist_type")
    end

    # Kernel density plot with filled region
    plt = density(
        data,
        bw = 2.0,
        fill = :steelblue,   # fill opacity and color
        linecolor = :black,         # density curve color
        xlabel = "Value",
        ylabel = "Density",
        title = "Kernel Density: $(dist_type)"
    )

    # Add vertical line at theoretical mean
    vline!([mu], linecolor=:black, linewidth=2)

    # Save figure to output directory
    savefig(plt, joinpath(outdir, fname))

    # Return plot object to support interactive exploration if needed
    return plt
end