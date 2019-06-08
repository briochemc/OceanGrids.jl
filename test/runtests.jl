using OceanGrids, Test

@testset "Small grid tests" begin
    grid = OceanGrid(2,3,4)
    @test grid isa OceanRectilinearGrid
    @test grid isa OceanGrid
    @testset "OceanGridBox" for box in grid
        @test box isa OceanGridBox
    end
end
