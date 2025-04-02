using Test
using ProjectionOntoSLnBenchmark

@testset "small benchmark" begin

    current_path = @__DIR__

    # make sure gfx folder is removed
    rm(joinpath(current_path, "../gfx"), force = true)

    results = ProjectionOntoSLnBenchmark.main(dims = [2, 3, 4, 8, 16], samples = 50)
    @test results !== nothing

    # are there any graphics?
    @test isdir(joinpath(current_path, "../gfx"))
    @test length(readdir(joinpath(current_path, "../gfx"))) > 0

end
