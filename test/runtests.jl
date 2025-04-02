using ProjectionOntoSLnBenchmark
using Test
using DrWatson: projectdir

@testset "small benchmark" begin

    # make sure gfx folder is removed
    rm(projectdir("gfx"), force = true)

    results = ProjectionOntoSLnBenchmark.main(dims = [2, 3, 4, 8, 16], samples = 50)
    @test results !== nothing

    # are there any graphics?
    @test isdir(projectdir("gfx"))
    @test length(readdir(projectdir("gfx"))) > 0

end
