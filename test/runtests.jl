using OceanGrids, Test, Unitful, SparseArrays

@testset "Small grid tests" begin
    nlat, nlon, ndepth = 2, 3, 4
    grd = OceanGrid(nlat, nlon, ndepth)
    @test grd isa OceanRectilinearGrid
    @test grd isa OceanGrid
    @test size(grd) == (nlat, nlon, ndepth)
    @test length(grd) == *(nlat, nlon, ndepth)
    @testset "OceanGridBox" for box in grd
        @test box isa OceanGridBox
        show(box)
        @test OceanGrids.area(box) isa Quantity
        @test unit(OceanGrids.area(box)) == u"km^2"
        @test unit(OceanGrids.volume(box)) == u"m^3"
    end
end

@testset "edges to grid" begin
    elat = -90:90:90
    elats = [-90:90:90,
             -90.0:90.0:90.0,
             [-90,0,90],
             [-90.0,0.0,90.0],
             (-90.0:90.0:90.0) * u"°",
             [-90.0,0.0,90.0] * u"°",
             (-90:90:90) * u"°",
             [-90,0,90] * u"°",]
    elons = [-180:90:180,
             [-180,-90,0,90,180],
             (-180:90:180) * u"°",
             [-180,-90,0,90,180] * u"°"]
    edepths = [[0, 100, 200, 500, 1000, 3000] * u"m",
               [0.0, 100.0, 200.0, 500.0, 1000.0, 3000.0] * u"m",
               [0.0, 100.0, 200.0, 500.0, 1000.0, 3000.0],
               [0, 100, 200, 500, 1000, 3000],
               [0, 0.1, 0.2, 0.5, 1, 3] * u"km"]
    base_grd = OceanGrid(elats[1], elons[1], edepths[1])
    @testset "Base constructors" for elat in elats, elon in elons, edepth in edepths
        grd = OceanGrid(elat,elon,edepth)
        @test grd isa OceanRectilinearGrid
        @test grd isa OceanGrid
        @test grd == base_grd
        show(grd)
        @testset "OceanGridBox" for box in grd
            @test box isa OceanGridBox
        end
    end
    @testset "Functions" begin
        x = ones(count(iswet(base_grd)))
        @testset "base functions" begin
            @test iswet(base_grd) isa BitArray{3}
            @test latvec(base_grd) isa Vector
            @test lonvec(base_grd) isa Vector
            @test depthvec(base_grd) isa Vector
            @test topdepthvec(base_grd) isa Vector
            @test bottomdepthvec(base_grd) isa Vector
            @test volumevec(base_grd) isa Vector
            @test isseafloorvec(base_grd) isa BitArray{1}
            @test latlondepthvecs(base_grd) isa Tuple{T,T,T} where {T<:Vector}
            x3D = rearrange_into_3Darray(x, base_grd)
            @test x3D isa Array{Float64,3}
            x2 = vectorize(x3D, base_grd)
            @test x2 == x
        end
        @testset "interpolations" begin
        end
        @testset "utils" begin
            @test buildIabove(base_grd) isa SparseMatrixCSC
            @test buildIbelow(base_grd) isa SparseMatrixCSC
        end
        @testset "Slices and views" begin
            @test horizontalslice(x, base_grd, depth=100) isa SubArray{Float64,2}
            @test meridionalslice(x, base_grd, lon=330) isa PermutedDimsArray{Float64,2}
            @test zonalslice(x, base_grd, lat=15) isa PermutedDimsArray{Float64,2}
            @test depthprofile(x, base_grd, lonlat=(330,15)) isa SubArray{Float64,1}
        end
        @testset "Integrals, averages, and stds" begin
            @test ∫dxdy(x, base_grd) isa Vector{T} where {T<:Quantity}
            @test ∫dz(x, base_grd) isa Array{T,2} where {T<:Quantity}
            @test ∫dy(x, base_grd) isa Array{T,2} where {T<:Quantity}
            @test ∫dx(x, base_grd) isa Array{T,2} where {T<:Quantity}
            @test ∫dxdydz(x, base_grd) isa Quantity
            @test horizontalmean(x, base_grd) isa Vector{Float64}
            @test verticalmean(x, base_grd) isa Array{Float64,2}
            @test meridionalmean(x, base_grd) isa Array{Float64,2}
            @test zonalmean(x, base_grd) isa Array{Float64,2}
            @test totalmean(x, base_grd) isa Float64
            @test horizontalstd(x, base_grd) isa Vector{Float64}
            @test verticalstd(x, base_grd) isa Array{Float64,2}
            @test meridionalstd(x, base_grd) isa Array{Float64,2}
            @test zonalstd(x, base_grd) isa Array{Float64,2}
            @test totalstd(x, base_grd) isa Float64
            @test zcumsum(x, base_grd) isa Vector{Float64}
        end
        @testset "Regridding" begin
            lat, lon = collect(-80:20:80), collect(15:30:360)
            y2D = cosd.(lat) * sind.(lon)'
            @test regrid(y2D, lat, lon, base_grd) isa Array{Float64,2}
            @test_broken regridandpaintsurface(y2d, lat, lon, base_grd) isa Array{Float64,3}
        end
    end

end
