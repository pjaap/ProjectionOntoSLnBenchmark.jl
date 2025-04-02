# ProjectionOntoSLnBenchmark.jl

This is a companion package for [ProjectionOntoSLn.jl](https://github.com/pjaap/ProjectionOntoSLn.jl) containing the benchmark code for the paper
[How to project onto SL(n)](https://arxiv.org/abs/2501.19310) with the DOI [https://doi.org/10.48550/arXiv.2501.19310](
https://doi.org/10.48550/arXiv.2501.19310).

This packages also generates the images used in the paper.

## Usage

To reproduce the results from the paper, you will need a Julia installation with version ≥ 1.11.0

1. Clone this repo
```
git clone https://github.com/pjaap/ProjectionOntoSLnBenchmark.jl
cd ProjectionOntoSLnBenchmark.jl
```

2. Start Julia with this project
```
julia --project
```

3. Instantiate the project (with the exact dependency versions of the paper's benchmark) from the `Manifest.toml` (enter `]` to get into the `pkg` mode)
```
pkg> instantiate

julia> using ProjectionOntoSLnBenchmark
julia> ProjectionOntoSLnBenchmark.main()
```

The whole benchmark run takes a while.
For some quick testing, you can reduce, e.g.,
- the dimension list (default is `dims = [2, 3, 4, 8, 16, 32, 64]`)
- number of samples (default is `samples = 1000`)

A quick test may look like this:
```
julia> ProjectionOntoSLnBenchmark.main(dims = [2,3,4,8], samples = 100)
```

All available arguments are stated at the end of `src/ProjectionOntoSLnBenchmark.jl`.

The resulting images will be generated in the `gfx` folder as PDFs.



