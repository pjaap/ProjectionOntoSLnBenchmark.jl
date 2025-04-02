module ProjectionOntoSLnBenchmark

using ProjectionOntoSLn

using Plots: plot, plot!, xlabel!, ylabel!, savefig
using StatsPlots
using LaTeXStrings
using Random: shuffle, seed!
using StaticArrays: SVector
using JLD2: save, jldopen
using LinearAlgebra: svd, norm, Diagonal
using DrWatson: projectdir

import BenchmarkTools

function tune_benchmark()
    # set these benchmark parameters to avoid expensive tuning of each single benchmark run
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 50
    BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
    BenchmarkTools.DEFAULT_PARAMETERS.evals_set = true
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = Inf
    # this removes very expensive garbage collections before each benchmark run (we do not need it)
    BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
    return nothing
end

function do_run(a, p, walltimes, iterations; stop_on_error, kwargs...)

    caught_exception = false

    methods = Dict(
        :cs => composite_step,
        :cn => constrained_newton,
        :rf => root_finding,
        :un => unconstrained_newton,
    )

    # track this run, append if all succeeded
    walltime = Dict()
    iteration = Dict()

    for key in keys(methods)
        try
            tune_benchmark()

            # @elapsed needs explicitly interpolated kwargs....
            time = BenchmarkTools.@belapsed project_onto_sln($(methods[key]), $a; tolerance = $(kwargs[:tolerance]), debug = $(kwargs[:debug]), maxIter = $(kwargs[:maxIter]))
            result = project_onto_sln(methods[key], a; kwargs...)
            p[key] = result.projection
            walltime[key] = time
            iteration[key] = result.iterations
        catch e
            println("caught an exception: ", e)
            caught_exception = true
            @show prod(a), a
            if stop_on_error
                rethrow()
            end
        end
    end

    if !caught_exception
        for key in keys(methods)
            push!(walltimes[key], walltime[key])
            push!(iterations[key], iteration[key])
        end
    end

    ### Debug output ####

    if !caught_exception &&
            (
            !testStationarity(a, p[:rf]) ||
                !testStationarity(a, p[:cn]) ||
                !testStationarity(a, p[:cs]) ||
                !testStationarity(a, p[:un])
        )

        dist_a(p) = norm(a - p)
        @show dist_a(p[:rf])
        @show dist_a(p[:cn])
        @show dist_a(p[:cs])
        @show dist_a(p[:un])
        @show testStationarity(a, p[:rf])
        @show testStationarity(a, p[:cn])
        @show testStationarity(a, p[:cs])
        @show testStationarity(a, p[:un])
        @show a
    end

    return !caught_exception
end

"""
    run benchmarks.

    dim: size of the matrices
    samples: how many different matrices should be computed?
    mode: (:convex, :concave, :boundary, :singular)
                convex: only matrices with det ∈ (0,1) in the order cone
                concave only matrices with det > 1
                boundary: matrices with det > 0, on the boundary of the order cone
                singular: matrices with det = 0; 1/3 of the singular values are zero
"""
function runBenchmark(; dim, samples, mode, kwargs...)

    n = dim

    # put everything in a dictionary, then we have named access with the
    # symbols :rf, etc.
    walltimes = Dict(
        :rf => [],
        :cs => [],
        :cn => [],
        :un => [],
        :svd => [],
    )

    # store the solutions for comparing
    p = Dict(
        :rf => [],
        :cs => [],
        :cn => [],
        :un => [],
    )


    # store the number of iterations
    iterations = Dict(
        :cs => [],
        :cn => [],
        :un => [],
        :rf => [],
    )

    successful_runs = 0
    total_runs = 0

    # track the determinants
    determinants = []


    while successful_runs < samples


        svd_cutoff_dim = Inf

        # for small matrices, we perform a full SVD on a dense random matrix
        # for larger matrices, we consider diagonal matrices only
        if dim ≤ svd_cutoff_dim
            A = generateRandomMatrix(n; diag_only = false, kwargs...)
        else
            A = generateRandomMatrix(n; diag_only = true, kwargs...)
        end

        tune_benchmark()

        # svd
        timeSVD = BenchmarkTools.@belapsed svd($A)
        U, a, V = svd(A)

        timeProduct = BenchmarkTools.@belapsed $U * Diagonal($a) * $V'

        time = timeSVD + timeProduct

        if mode == :singular
            # set zero in the last third of the entries
            m = n ÷ 3
            a[(end - m):end] .= 0
        end

        if mode == :boundary
            # pick 1/3 of the entries and use geometric mean
            # for neighboring selected entries
            m = max(1, (n ÷ 3))
            boundary_indices = sort(shuffle(1:(n - 1))[1:m])
            geo_mean(v) = prod(v)^(1 / length(v))

            i = 1
            while i ≤ m
                bi = boundary_indices[i]
                j = 1
                # find consecutive selected boundary_indices
                while i + j ≤ m && boundary_indices[i + j] == bi + j
                    j += 1
                end
                a[bi:(bi + j)] .= geo_mean(a[bi:(bi + j)])
                i += j + 1
            end


        end

        # make a static constant array
        a = SVector{n}(a)

        Π = prod(a)

        if mode == :concave && Π < 1
            continue
        end

        if mode == :convex && Π > 1
            continue
        end

        # proceed only if do_run reports success
        if do_run(a, p, walltimes, iterations; kwargs...)
            if samples ≥ 10 && successful_runs % div(samples, 10) == 0
                println(div(successful_runs * 100, samples), "% runs done")
            end

            if dim ≤ svd_cutoff_dim
                push!(walltimes[:svd], time)
            end

            push!(determinants, prod(a))

            successful_runs += 1
        end
        total_runs += 1
    end

    @info "$successful_runs successful runs, $total_runs runs in total."
    open("total_runs.log", "a") do file
        write(file, "mode = $mode, dim = $dim: $successful_runs successful runs, $total_runs runs in total\n")
    end

    return walltimes,
        determinants,
        iterations
end

function runBenchmarks(; dims, kwargs...)

    N = length(dims)

    # wall time array
    wt = Vector{Dict{Symbol, Vector{Float64}}}(undef, N)

    # iterations array
    it = Vector{Dict{Symbol, Vector{Int64}}}(undef, N)

    # determinants of all runs
    determinants = []

    # append the wall times for all runs to the corresponding lists (there is probably a better solution)
    for i in 1:N

        println("\nrun test with dim = ", dims[i])

        wt[i], new_det, it[i] = runBenchmark(; dim = dims[i], kwargs...)
        determinants = vcat(determinants, new_det)
    end

    return wt, determinants, it

end


function compute_min_max_average(vector_of_dicts::AbstractVector)

    mins = []
    maxs = []
    avgs = []

    for v in vector_of_dicts
        # create a Dict of floats with the same keys as in vector_of_dicts
        av = Dict(k => 0.0 for k in keys(v))
        mi = Dict(k => 0.0 for k in keys(v))
        ma = Dict(k => 0.0 for k in keys(v))

        for k in keys(v)

            s = sum(v[k])
            n = length(v[k])

            av[k] = s / n
            mi[k] = minimum(v[k])
            ma[k] = maximum(v[k])
        end

        push!(avgs, av)
        push!(mins, mi)
        push!(maxs, ma)
    end

    return mins, maxs, avgs
end

# mode ∈ { :convex, :concave, :singular, :boundary}
function makePlots(walltimes, determinants::AbstractArray, iterations; dims, maxIter, mode, kwargs...)

    # generate labels
    labels = Dict(
        :svd => "SVD",
        :rf => "Root Finding",
        :cn => "Newton (constrained)",
        :cs => "Composite Step",
        :un => "Newton (eliminated)",
    )

    min_wt, max_wt, avg_wt = compute_min_max_average(walltimes)

    # a small portion of the dims for the O(n) graphs
    O_dims = 2:7
    O_factor = 5.0e-6

    # convert dims to a vector of strings
    dim_string = string.(dims)

    bar_width = 0.7

    # plot wall times
    symbols = [:rf :cs :un  :cn]

    # convert to matrices of size N × M,
    # where M is the number of methods and N the number of dims
    N = length(dims)
    M = length(symbols)

    min_matrix = zeros(N, M)
    max_matrix = zeros(N, M)
    avg_matrix = zeros(N, M)

    for i in 1:N, (j, symbol) in enumerate(symbols)
        min_matrix[i, j] = min_wt[i][symbol]
        max_matrix[i, j] = max_wt[i][symbol]
        avg_matrix[i, j] = avg_wt[i][symbol]
    end

    pp = groupedbar(dim_string, max_matrix, fillto = min_matrix, bar_width = bar_width, label = [labels[s] for s in symbols], yscale = :log10, legend = :topleft)

    # plot small avg bars
    groupedbar!(pp, dim_string, avg_matrix .* 1.02, fillto = avg_matrix .* 0.98, bar_width = bar_width, label = nothing, color = :black)

    # plot average SVD time
    avg_svd_wt = [ avg_wt[i][:svd] for i in 1:N ]
    # mapping to dim_string does not work: hack it by using [ 0.5, ... , N-0.5] for the x-range
    plot!(pp, 0.5:(N - 0.5), avg_svd_wt, style = :dash, mark = :circ, label = "average SVD wall-time")

    # plot O(n)
    c = avg_wt[end][:rf] / (3dims[end])
    # mapping to dim_string does not work: hack it by using [ 0.5, ... , N-0.5] for the x-range
    plot!(pp, 3.5:(N - 0.5), c .* dims[4:end], style = :dash, label = L"\mathcal{O}(n)")

    # plot O(n^3)
    c = avg_wt[end][:un] / (4dims[end]^3)
    # mapping to dim_string does not work: hack it by using [ 0.5, ... , N-0.5] for the x-range
    plot!(pp, 3.5:(N - 0.5), c .* dims[4:end] .^ 3, style = :dash, label = L"\mathcal{O}(n^3)")

    yticks!(pp, [10.0^i for i in -10:0])
    ylims!(pp, 10.0^-7, 10.0^-0)

    xlabel!(pp, L"n")
    ylabel!(pp, "wall-time [s]")

    savefig(pp, projectdir("gfx","walltimes_$mode.pdf"))

    # det plots only for boundary
    if mode == :boundary

        n = length(determinants)
        m = length(dims)
        split = div(n, m)

        # no log-style density plot: plot log_10(determinants) and remap the ticks
        xticks = collect(-3:3)
        mapped_ticks = (xticks, [ L"10^{%$i}" for i in xticks ])

        pp = density(xticks = mapped_ticks, xlimits = (-4, 4), size = (500, 250))

        for i in eachindex(dims)
            # split range in split long parts
            range = (1:split) .+ (i - 1) * split
            density!(pp, log10.(determinants[range]), label = L"n=%$(dims[i])", linewidth = 2, bandwidth = 0.1)
        end


        xlabel!(pp, "determinant")
        ylabel!(pp, "relative frequency")

        savefig(pp, projectdir("gfx","determinants_$mode.pdf"))
    end

    # plot average iterations
    min_it, max_it, avg_it = compute_min_max_average(iterations)

    # convert to matrices of size N × M,
    # where M is the number of methods and N the number of dims
    N = length(dims)
    M = length(symbols)

    min_matrix = zeros(N, M)
    max_matrix = zeros(N, M)
    avg_matrix = zeros(N, M)

    for i in 1:N, (j, symbol) in enumerate(symbols)
        min_matrix[i, j] = min_it[i][symbol]
        max_matrix[i, j] = max_it[i][symbol]
        avg_matrix[i, j] = avg_it[i][symbol]
    end

    pp = groupedbar(dim_string, max_matrix, fillto = min_matrix, bar_width = bar_width, label = [labels[s] for s in symbols])

    # plot small avg bars
    groupedbar!(pp, dim_string, avg_matrix .+ 0.2, fillto = avg_matrix .- 0.2, bar_width = bar_width, label = nothing, color = :black)

    ylims!(pp, 0, maxIter)

    xlabel!(pp, L"n")
    ylabel!(pp, "number of iterations")

    savefig(pp, projectdir("gfx","iterations_$mode.pdf"))

    return nothing

end

# do both in one call
function runBenchmarksAndPlot(; mode, kwargs...)
    @info "runBenchmarksAndPlot with mode = $mode"

    result = runBenchmarks(; mode, kwargs...)
    makePlots(result...; mode, kwargs...)

    return result
end

# do all benchmarks for all modes
function do_everything(; random_seed, kwargs...)

    # remove log files
    rm(projectdir("total_runs.log"), force = true)
    rm(projectdir("non_convergence.log"), force = true)

    # create gfx folder if not already present
    if !isdir(projectdir("gfx"))
        mkdir(projectdir("gfx"))
    end

    # this makes the rand() calls a bit more reproducible
    seed!(random_seed)

    results = Dict()

    for mode in [:boundary, :convex, :concave, :singular]
        results[mode] = runBenchmarksAndPlot(; mode, kwargs...)
    end

    save(projectdir("gfx","Benchmark_results.jld2"), "results", results, "random_seed", random_seed, "kwargs", kwargs)

    return results
end

function plot_from_file(filepath)
    file = jldopen(filepath)

    results = file["results"]
    kwargs = file["kwargs"]

    for mode in keys(results)
        result = results[mode]
        makePlots(result...; mode, kwargs...)
    end
    return nothing
end

function main(;
        dims = [2, 3, 4, 8, 16, 32, 64],
        samples = 1000,
        tolerance = 1.0e-8,
        det_range = [0.01, 100.0],
        maxIter = 200,
        debug = false,
        stop_on_error = false,
        random_seed = 12345
    )

    return ProjectionOntoSLnBenchmark.do_everything(;
        dims,
        samples,
        tolerance,
        det_range,
        maxIter,
        debug,
        stop_on_error,
        random_seed
    )

end


end # module ProjectionOntoSLnBenchmark
