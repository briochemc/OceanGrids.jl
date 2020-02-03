module OceanGrids

using Unitful, Interpolations
using LinearAlgebra, SparseArrays

abstract type OceanGrid end

"""
    OceanRectilinearGrid

Ocean rectilinear grid.
Here I use the slightly abuse the term rectilinear,
because the grid supposedly represents layers around the (curved) globe.
Here "rectilinear" means that the grid is defined entirely by just its
`lat`, `lon`, `depth`, `δlat`, `δlon`, and `δdepth` vectors.
For example, the Ocean Circulation Inverse Model (OCIM) version 1
grid is rectilinear in that sense.
"""
struct OceanRectilinearGrid <: OceanGrid
    lat          # °
    lon          # °
    depth        # m
    δlat         # °
    δlon         # °
    δdepth       # m
    lat_3D       # °
    lon_3D       # °
    depth_3D     # m
    δy           # m
    δx_3D        # m
    δy_3D        # m
    δz_3D        # m
    volume_3D    # m³
    depth_top    # m
    depth_top_3D # m
    A_2D         # m²
    wet3D::BitArray
    nlon
    nlat
    ndepth
    nboxes
end

"""
    OceanCurvilinearGrid

Ocean curvilinear grid.
The grid is curvilinear if it can be accessed using Cartesian Indices.
I.e., it maps to a rectilinear grid.
This type of grid typically cannot be inferred from only
`lat`, `lon`, `depth`, `δlat`, `δlon`, and `δdepth` vectors.
Instead, the full 3-dimensional arrays
`lat_3D`, `lon_3D`, `D`epth_3D`, `δx_3D`, `δy_3D`, and `δz_3D`
are required.
For exaample, the displaced-pole grid of the POP model is curcvilinear.
"""
struct OceanCurvilinearGrid <: OceanGrid
    lat_3D       # °
    lon_3D       # °
    depth_3D     # m
    δx_3D        # m
    δy_3D        # m
    δz_3D        # m
    volume_3D    # m³
    depth_top    # m
    depth_top_3D # m
    A_2D         # m²
    wet3D::BitArray
    nlon
    nlat
    ndepth
    nboxes
end

const TU = AbstractArray{<:Quantity}
const T1 = AbstractArray{<:Number}

"""
    OceanGrid(elon::T, elat::U, edepth::V; R=6371.0u"km") where {T<:TU, U<:TU, V<:TU}

Returns an `OceanRectilinearGrid` with boxes whose edges are defined by the
`elon`, `elat`, and `edepth` vectors.
The globe radius can be changed with the keyword `R` (default value 6371 km)
The 3D array of wet boxes, `wet3D` is true everywhere by default.
"""
function OceanGrid(elon::TU, elat::TU, edepth::TU; R=6371.0u"km")
    unitful_data, unitless_data = generate_grid_data_no_wet3D(elon, elat, edepth, R)
    #
    nlat, nlon, ndepth = unitless_data[[2,1,3]]
    wet3D = trues(nlat, nlon, ndepth)
    return OceanRectilinearGrid(unitful_data..., wet3D, unitless_data...)
end

"""
    OceanGrid(elon::T, elat::U, edepth::V, wet3D::BitArray; R=6371.0u"km")

Same as `OceanGrid(elon, elat, edepth; R)` but with specified `wet3D`.
"""
function OceanGrid(elon::TU, elat::TU, edepth::TU, wet3D::BitArray; R=6371.0u"km")
    unitful_data, unitless_data = generate_grid_data_no_wet3D(elon, elat, edepth, R)
    return OceanRectilinearGrid(unitful_data..., wet3D, unitless_data...)
end

function generate_grid_data_no_wet3D(elon::TU, elat::TU, edepth::TU, R)
    R = R |> u"m"
    edepth = edepth .|> u"m"
    elat = elat .|> u"°"
    edon = elon .|> u"°"
    # δ objects
    nlon, nlat, ndepth = length(elon) - 1, length(elat) - 1, length(edepth) - 1
    nboxes = nlon * nlat * ndepth
    lon = edges_to_centers(elon)
    lat = edges_to_centers(elat)
    depth = edges_to_centers(edepth)
    δlon, δlat, δdepth = diff(elon), diff(elat), diff(edepth)
    # depth objects
    δz_3D = repeat(reshape(δdepth, (1,1,ndepth)), outer=(nlat,nlon,1))
    depth_3D = repeat(reshape(depth, (1,1,ndepth)), outer=(nlat,nlon,1))
    depth_top = cumsum(δdepth) - δdepth
    depth_top_3D = repeat(reshape(depth_top, (1,1,ndepth)), outer=(nlat,nlon,1))
    # lat objects
    lat_3D = repeat(reshape(lat, (nlat,1,1)), outer=(1,nlon,ndepth))
    δy = R * δlat ./ 360u"°"
    δy_3D = repeat(reshape(δy, (nlat,1,1)), outer=(1,nlon,ndepth))
    # lon objects
    lon_3D = repeat(reshape(lon, (1,nlon,1)), outer=(nlat,1,ndepth))
    # For δx_3D, first calculate the area and then infer the equivalent distance in meters
    A_2D = R^2 * abs.(sin.(elat[1:end-1]) - sin.(elat[2:end])) * ustrip.(δlon .|> u"rad")'
    A_3D = repeat(A_2D, outer=(1, 1, ndepth))
    δx_3D = A_3D ./ δy_3D
    # volume
    volume_3D = δx_3D .* δy_3D .* δz_3D
    unitful_data = (lat, lon, depth, δlat, δlon, δdepth, lat_3D, lon_3D, depth_3D,
                    δy, δx_3D, δy_3D, δz_3D, volume_3D, depth_top, depth_top_3D, A_2D)
    unitless_data = (nlon, nlat, ndepth, nboxes)
    return unitful_data, unitless_data
end

# Convert units in case not provided
function OceanGrid(elon::TU, elat::TU, edepth::T1; R=6371.0u"km")
    @warn "`edepth` provided without units — assuming meters"
    return OceanGrid(elon, elat, edepth * u"m"; R=R)
end
function OceanGrid(elon::T1, elat::TU, edepth::TU; R=6371.0u"km")
    @warn "`elon` provided without units — assuming degrees"
    return OceanGrid(elon * u"°", elat, edepth; R=R)
end
function OceanGrid(elon::TU, elat::T1, edepth::TU; R=6371.0u"km")
    @warn "`elat` provided without units — assuming degrees"
    return OceanGrid(elon, elat * u"°", edepth; R=R)
end
function OceanGrid(elon::T1, elat::T1, edepth::TU; R=6371.0u"km")
    @warn "`elat` and `elon` provided without units — assuming degrees"
    return OceanGrid(elon * u"°", elat * u"°", edepth; R=R)
end
function OceanGrid(elon::T1, elat::T1, edepth::T1; R=6371.0u"km")
    @warn "`elat`, `elon`, and `edepth` provided without units — assuming degrees and meters"
    return OceanGrid(elon * u"°", elat * u"°", edepth * u"m"; R=R)
end
function OceanGrid(elon::TU, elat::T1, edepth::T1; R=6371.0u"km")
    @warn "`elat` and `edepth` provided without units — assuming and meters"
    return OceanGrid(elon, elat * u"°", edepth * u"m"; R=R)
end
function OceanGrid(elon::T1, elat::TU, edepth::T1; R=6371.0u"km")
    @warn "`elon` and `edepth` provided without units — assuming degrees and meters"
    return OceanGrid(elon * u"°", elat, edepth * u"m"; R=R)
end

edges_to_centers(x::Vector) = 0.5 * (x[1:end-1] + x[2:end])
edges_to_centers(x::AbstractRange) = x[1:end-1] .+ 0.5step(x)

"""
    OceanGrid(nlat::Int, nlon::Int, ndepth::Int)

Returns a regularly spaced `OceanRectilinearGrid` with size `nlat`, `nlon`, and `ndepth`.
"""
function OceanGrid(nlat::Int, nlon::Int, ndepth::Int)
    elat = range(-90,90,length=nlat+1) * u"°"
    elon = range(0,360,length=nlon+1) * u"°"
    edepth = range(0,3682,length=ndepth+1) * u"m"
    return OceanGrid(elon, elat, edepth)
end

function Base.show(io::IO, g::OceanGrid)
    println("OceanGrid of size $(g.nlat)×$(g.nlon)×$(g.ndepth) (lat×lon×depth)")
end

"""
    OceanGridBox

Ocean grid box.
Each grid can be looped over, where the `OceanGridBox` are the elements of the grid.
This useful to investigate specific boxes in the grid.
"""
struct OceanGridBox
    I
    lat          # °
    lon          # °
    depth        # m
    δlat         # °
    δlon         # °
    δdepth       # m
    δx           # m
    δy           # m
    δz           # m
    volume       # m³
    depth_top    # m
    A            # m²
    iswet::Bool
end

"""
    box(g::OceanGrid, i, j, k)

Accesses the individual box of `g::OceanGrid` at index `(i,j,k)`.
Each grid can be looped over, where the `OceanGridBox` are the elements of the grid.
This useful to investigate specific boxes in the grid.
"""
function box(g::OceanGrid, i, j, k)
    return OceanGridBox(
                        CartesianIndex(i,j,k),
                        g.lat[i],
                        g.lon[j],
                        g.depth[k],
                        g.δlat[i],
                        g.δlon[j],
                        g.δdepth[k],
                        g.δx_3D[i,j,k],
                        g.δy_3D[i,j,k],
                        g.δz_3D[i,j,k],
                        g.volume_3D[i,j,k],
                        g.depth_top[k],
                        g.A_2D[i,j],
                        g.wet3D[i,j,k]
                       )
end

"""
    box(g::OceanGrid, I)

Accesses the individual box of `g::OceanGrid` at index `I`.
Each grid can be looped over, where the `OceanGridBox` are the elements of the grid.
This useful to investigate specific boxes in the grid.
"""
function box(g::OceanGrid, I)
    i,j,k = CartesianIndices((g.nlat,g.nlon,g.ndepth))[I].I
    return box(g::OceanGrid, i, j, k)
end

Base.getindex(g::OceanGrid, I) = box(g,I)

"""
    Base.size(g::OceanGrid)

Size of the grid.
"""
Base.size(g::OceanGrid) = g.nlat, g.nlon, g.ndepth

"""
    Base.length(g::OceanGrid)

Size of the grid.
"""
Base.length(g::OceanGrid) = g.nlat * g.nlon * g.ndepth

function Base.show(io::IO, b::OceanGridBox)
    println("$(wet_print(b.iswet)) OceanGridBox at $(b.I):")
    println("  location: $(round(b.lat,digits=1))N, $(round(b.lon,digits=1))E")
    println("  depth: $(round(b.depth,digits=1))")
    println("  size: $(round(b.δx |> u"km",digits=1)) × $(round(b.δy |> u"km",digits=1)) × $(round(b.δz,digits=1)) (δx × δy × δz)")
end

wet_print(iswet) = iswet ? "Wet" : "Dry"

Base.round(q::Quantity; digits=0) = round(q |> ustrip, digits=digits) * unit(q)


area(b::OceanGridBox) = b.A |> u"km^2"
volume(b::OceanGridBox) = b.volume

Base.iterate(g::OceanGrid) = box(g::OceanGrid, 1), 1
function Base.iterate(g::OceanGrid, i)
    if i == g.nboxes
        return nothing
    else
        return box(g::OceanGrid, i+1), i+1
    end
end

Base.isequal(g₁::T, g₂::T) where {T<:OceanGrid} = prod([getfield(g₁,f) == getfield(g₂,f) for f in fieldnames(T)])
Base.:(==)(g₁::T, g₂::T) where {T<:OceanGrid} = isequal(g₁, g₂)

export OceanGrid, OceanCurvilinearGrid, OceanRectilinearGrid, box, OceanGridBox

"""
    iswet(grid)

Returns the 3D BitArray of wet boxes of the grid.
"""
iswet(g::T) where {T<:OceanGrid} = g.wet3D
export iswet

"""
    latvec(grid)

Returns the vector (wet points) of latitudes (units stripped).
"""
latvec(g::T) where {T<:OceanGrid} = ustrip.(g.lat_3D[iswet(g)])

"""
    lonvec(grid)

Returns the vector (wet points) of longitudes (units stripped).
"""
lonvec(g::T) where {T<:OceanGrid} = ustrip.(g.lon_3D[iswet(g)])

"""
    depthvec(grid)

Returns the vector (wet points) of depths (units stripped).
"""
depthvec(g::T) where {T<:OceanGrid} = ustrip.(g.depth_3D[iswet(g)])

"""
    latlondepthvecs(g)

Returns `latvec(g), lonvec(g), depthvec(g)`.
"""
latlondepthvecs(g::T) where {T<:OceanGrid} = latvec(g), lonvec(g), depthvec(g)

export latvec, lonvec, depthvec, latlondepthvecs

#===================================
1D Vector <-> 3D array conversions
===================================#
"""
    rearrange_into_3Darray(x, wet3D::BitArray)

Returns a 3D array of `x` rearranged to the `true` entries of `wet3D`.
Entries where `wet3D` is `false` are filled with `NaN`s.
"""
function rearrange_into_3Darray(x, wet3D::BitArray)
    iwet = findall(wet3D)
    x3d = fill(NaN, size(wet3D))
    x3d[iwet] .= x
    return x3d
end
"""
    rearrange_into_3Darray(x, grid)

Returns a 3D array of `x` rearranged to the wet boxes of the grid.
"""
rearrange_into_3Darray(x, grid) = rearrange_into_3Darray(x, grid.wet3D)

export rearrange_into_3Darray

#============
Interpolation
============#
"""
    interpolate(x, g)

Returns an "interpolation object" for fast interpolation.

It must be queried within the range of the grid in longitude...
"""
function Interpolations.interpolate(x, g::T) where {T<:OceanGrid}
    knots = (ustrip.(g.lat), periodic_longitude(ustrip.(g.lon)), ustrip.(g.depth))
    A = periodic_longitude(rearrange_into_3Darray(x, g))
    itp = interpolate(knots, A, Gridded(Constant()))
    return extrapolate(itp, (Flat(), Periodic(), Flat())) # extrapolate periodically
end
function Interpolations.interpolate(x, g::OceanGrid, lat, lon, depth; itp=interpolate(x,g))
    return itp(ustrip(lat), ustrip.(lon), ustrip.(depth))
end
function Interpolations.interpolate(x, g::OceanGrid, lats::Vector, lons::Vector, depths::Vector; itp=interpolate(x,g))
    return [itp(y,x,z) for (y,x,z) in zip(ustrip.(lats), ustrip.(lons), ustrip.(depths))]
end
function Interpolations.interpolate(x, g::OceanGrid, MD::NamedTuple; itp=interpolate(x,g))
    return interpolate(x, g, MD.lat, MD.lon, MD.depth; itp=itp)
end
function Interpolations.interpolate(x, g::OceanGrid, obs; itp=interpolate(x,g))
    return interpolate(x, g, obs.metadata; itp=itp)
end
Interpolations.interpolate(x, g::OceanGrid, ::Missing; itp=nothing) = missing
export interpolate

"""
    iswet(g, lat, lon, depth)
    iswet(g, lats, lons, depths)
    iswet(g, metadata)

Returns the indices of the provided locations that are "in" a wet box of `g`.

Technically, this uses a nearest neighbour interpolation algorithm, so the bounding boxes
will not work perfectly.
If AIBECS will solely rely on this nearest neighbour interpolation,
then it might be a good idea to replace Interpolations.jl with NearestNeighbors.jl.
"""
function iswet(g, args...; itp=interpolate(1:count(iswet(g)),g))
    J = interpolate(1:count(iswet(g)), g, args...; itp=itp)
    return findall((!isnan).(J))
end
iswet(g, ::Missing; itp=nothing) = missing

"""
    interpolationmatrix(g, lats, lons, depths)
    interpolationmatrix(g, metadata)

Returns the sparse matrix `M` such that `M*x` is the model vector `x` interpolated
onto the provided locations/metadata.

Technically, this requires a linear interpolation, in the sense that
`interpolation(x) = M * x` for some `M`, which seems to require a nearest neighbor
interpolation.
If AIBECS will solely rely on this nearest neighbour interpolation,
then it might be a good idea to replace Interpolations.jl with NearestNeighbors.jl.
Also, the few observations that lie in model land (instead of water) are discarded.

In the future, a more generic approach will use an sparse-aware AD for any type of
interpolation, and use an interpolation that somehow does not discard observations
made at locations not inside an ocean grid cell.
"""
function interpolationmatrix(g, args...; itp=interpolate(1:count(iswet(g)),g))
    J = interpolate(1:count(iswet(g)), g, args...; itp=itp)
    iwet = findall((!isnan).(J))
    J = Int.(J[iwet])
    I = collect(1:length(iwet))
    return sparse(I, J, true, length(I), count(iswet(g)))
end
interpolationmatrix(g, ::Missing; itp=nothing) = missing
export interpolationmatrix


#=
    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    lons = range(ustrip(grd.lon[1]), ustrip(grd.lon[1])+360, length=length(grd.lon)+1)
    stp = Interpolations.scale(itp, lats, lons, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude


    xs = range(0, 2π, length=11)
    ys = range(-π/2, π/2, length=6)
    A = f.(xs,ys')
    itp = interpolate(A, BSpline(Linear())) # interpolate linearly between the data points
    stp = scale(itp, xs, ys) # re-scale to the actual domain
    etp = extrapolate(stp, (Periodic(), Throw())) # extrapolate periodically
=#

periodic_longitude(x3D::Array{T,3}) where T = view(x3D,:,[1:size(x3D,2); 1],:)
periodic_longitude(v::Vector{T}) where T = [v; v[1]+360]

#===================================
Off-diaongal matrices
===================================#

"""
    buildIbelow(wet3D, iwet)

Builds the shifted-diagonal sparse matrix of the indices of below neighbours.

Ibelow[i,j] = 1 if the box represented by the linear index i
lies directly below the box represented by the linear index j.
"""
function buildIbelow(wet3D, iwet)
    nlat, nlon, ndepth = size(wet3D)
    n = nlon * nlat * (ndepth + 1)
    In = sparse(I, n, n)
    idx = zeros(Int, nlat, nlon, ndepth + 1)
    idx[:] = 1:n
    idx .= idx[:, :, [2:ndepth + 1; 1]]      # downward shift
    return In[idx[:], :][iwet, iwet]
end
"""
    buildIbelow(grid)

Builds the shifted-diagonal sparse matrix of the indices of below neighbours for `grid`.
See `buildIbelow(wet3D, iwet)`.
"""
buildIbelow(grid) = buildIbelow(grid.wet3D, indices_of_wet_boxes(grid))
export buildIbelow

"""
    buildIabove(wet3D, iwet)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours.

Iabove[i,j] = 1 if the box represented by the linear index i
lies directly above the box represented by the linear index j.
"""
buildIabove(wet3D, iwet) = copy(transpose(buildIbelow(wet3D, iwet)))
"""
    buildIabove(grid)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours for `grid`.
See `buildIabove(wet3D, iwet)`.
"""
buildIabove(grid) = copy(transpose(buildIbelow(grid)))
export buildIabove



#===================================
Functions returning useful arrays,
vectors, and constants
===================================#

"""
    array_of_volumes(grid)

Returns the 3D array of volumes of grid boxes.
"""
array_of_volumes(grid) = grid.volume_3D
export array_of_volumes


"""
    indices_of_wet_boxes(wet3D::BitArray)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(wet3D::BitArray) = findall(vec(wet3D))
"""
    indices_of_wet_boxes(grid)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(grid) = indices_of_wet_boxes(grid.wet3D)
export indices_of_wet_boxes

"""
    number_of_wet_boxes(wet3D::BitArray)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(wet3D::BitArray) = length(indices_of_wet_boxes(wet3D))
"""
    number_of_wet_boxes(grid)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(grid) = number_of_wet_boxes(grid.wet3D)
export number_of_wet_boxes

"""
    vector_of_volumes(grid)

Returns the vector of volumes of wet boxes.
"""
vector_of_volumes(grid) = array_of_volumes(grid)[indices_of_wet_boxes(grid.wet3D)]
export vector_of_volumes

"""
    vector_of_depths(grid)

Returns the vector of depths of the center of wet boxes.
"""
vector_of_depths(grid) = grid.depth_3D[indices_of_wet_boxes(grid.wet3D)]
export vector_of_depths

"""
    vector_of_depths(grid)

Returns the vector of depths of the top of wet boxes.
"""
vector_of_top_depths(grid) = grid.depth_top_3D[indices_of_wet_boxes(grid.wet3D)]
export vector_of_top_depths

#===================================
Functions return weighted norms
===================================#

"""
    weighted_norm²(x, w)

Returns the square of the weighted norm of `x` using weights `w`.
"""
weighted_norm²(x, w) = transpose(x) * Diagonal(w) * x

"""
    weighted_norm(x, v)

Returns the weighted norm of `x` using weights `w`.
"""
weighted_norm(x, w) = sqrt(weighted_norm²(x, w))

"""
    weighted_mean(x, w)

Returns the weighted mean of `x` using weights `w`.
"""
weighted_mean(x, w) = transpose(w) * x / sum(w)

"""
    vmean(x, w, I)

Returns the weighted mean of `x` using weights `w`, but only over indices `I`.
This is useful to get the observed mean.
(Because there are some grid boxes without observations.)
"""
weighted_mean(x, w, I) = transpose(w[I]) * x[I] / sum(w[I])
export weighted_norm²


end # module
