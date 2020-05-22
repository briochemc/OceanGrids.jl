
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
    topdepthvec(grid)

Returns the vector (wet points) of top depths (units stripped).
(top depths = depth of the top of the box.)
"""
topdepthvec(g::T) where {T<:OceanGrid} = ustrip.(g.depth_top_3D[iswet(g)])

"""
    bottomdepthvec(grid)

Returns the vector (wet points) of bottom depths (units stripped).
(bottom depths = depth of the bottom of the box.)
"""
bottomdepthvec(g::T) where {T<:OceanGrid} = 2depthvec(g) - topdepthvec(g)
function isseafloorvec(g::T) where {T<:OceanGrid}
    nlat, nlon, ndepth = size(g)
    iswetext = cat(iswet(g), falses(nlat,nlon,1), dims=3)
    isseafloor3D =  @views iswetext[:,:,1:ndepth] .& .!iswetext[:,:,2:ndepth+1]
    return isseafloor3D[iswet(g)]
end

"""
    latlondepthvecs(g)

Returns `latvec(g), lonvec(g), depthvec(g)`.
"""
latlondepthvecs(g::T) where {T<:OceanGrid} = latvec(g), lonvec(g), depthvec(g)

export latvec, lonvec, depthvec, latlondepthvecs, topdepthvec, bottomdepthvec
export isseafloorvec

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
    x3d = fill(NaN * unit(eltype(x)), size(wet3D))
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

Returns a `BitArray` of the provided locations that are "in" a wet box of `g`.

Technically, this uses a nearest neighbour interpolation algorithm, so the bounding boxes
will not work perfectly.
If AIBECS will solely rely on this nearest neighbour interpolation,
then it might be a good idea to replace Interpolations.jl with NearestNeighbors.jl.
"""
function iswet(g, args...; itp=interpolate(1:count(iswet(g)),g))
    J = interpolate(1:count(iswet(g)), g, args...; itp=itp)
    return (!isnan).(J)
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
buildIabove(grid) = buildIabove(grid.wet3D, indices_of_wet_boxes(grid))
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
    vector_of_top_depths(grid)

Returns the vector of depths of the top of wet boxes.
"""
vector_of_top_depths(grid) = grid.depth_top_3D[indices_of_wet_boxes(grid.wet3D)]
export vector_of_top_depths

"""
    vector_of_bottom_depths(grid)

Returns the vector of depths of the bottom of wet boxes.
"""
vector_of_bottom_depths(grid) = 2vector_of_depths(grid) - vector_of_top_depths(grid)
export vector_of_bottom_depths

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





#===============
Slices and views
===============#
# depth helper functions
finddepthindex(depth::Quantity, grd) = findmin(abs.(grd.depth .- depth))[2]
convertdepth(depth::Real) = depth * u"m"
convertdepth(depth::Quantity{U, Unitful.𝐋, V}) where {U,V} = depth
convertdepth(x) = error("Not a valid depth")
finddepthindex(depths::Tuple, grd) = ((finddepthindex(d, grd) for d in depths)...,)
convertdepth(depths::Tuple) = ((convertdepth(d) for d in depths)...,)
# lon helper functions
findlonindex(lon::Quantity, grd) = findmin(abs.(grd.lon .- lon))[2]
convertlon(x::Real) = mod(x * u"°", 360u"°")
convertlon(x::Quantity) = unit(x) == Unitful.° ? mod(x, 360u"°") : error("Not a valid lon")
# lat helper functions
findlatindex(lat::Quantity, grd) = findmin(abs.(grd.lat .- lat))[2]
convertlat(y::Real) = y * u"°"
convertlat(y::Quantity) = unit(y) == Unitful.° ? y : error("Not a valid lat")
"""
    HorizontalSlice(x, grd; depth)

Returns a horizontal slice of tracer `x` at depth `depth`.
"""
function horizontalslice(x3D, grd; depth=nothing)
    isnothing(depth) && error("you must specify the depth, e.g., `depth=100`")
    depth = convertdepth(depth)
    iz = finddepthindex(depth, grd)
    return view(x3D, :, :, iz)
end
"""
    meridionalslice(x, grd; lon)

Returns a meridional slice of tracer `x` at longitude `lon`.
"""
function meridionalslice(x3D, grd; lon=nothing)
    isnothing(lon) && error("you must specify the longitude, e.g., `lon=100`")
    lon = convertlon(lon)
    ix = findlonindex(lon, grd)
    return permutedims(view(x3D,:,ix,:), [2,1])
end
"""
    zonalslice(x, grd; lat)

Returns a zonal slice of tracer `x` at latitidue `lat`.
"""
function zonalslice(x3D, grd; lat=nothing)
    isnothing(lat) && error("you must specify the latitude, e.g., `lat=-30`")
    lat = convertlat(lat)
    iy = findlatindex(lat, grd)
    return permutedims(view(x3D,iy,:,:), [2,1])
end

"""
    depthprofile(x, grd; lonlat)

Returns the profile of tracer `x` interpolated at `lonlat=(x,y)` coordinates.
"""
function depthprofile(x, grd; lonlat=nothing)
    isnothing(lonlat) && error("you must specify the coordinates, e.g., `lonlat=(-30,30)`")
    lon, lat = lonlat
    lon, lat = convertlon(lon), convertlat(lat)
    x3D = rearrange_into_3Darray(x, grd)
    ndepth = length(grd.depth)
    knots = (grd.lat, grd.lon, 1:ndepth)
    itp = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    return itp(lat, lon, 1:ndepth)
end
export depthprofile


#=====================
Integrals and averages
=====================#
"""
    integral(x, grd; mask=iswet(grd), dims=(1,2,3))

Returns a "generic" integral of tracer `x` along `dims`.
"""
function integral(x3D, grd, mask3D=iswet(grd); dims=(1,2,3))
    w3D = (1 ∈ dims ? grd.δy_3D : 1) .* (2 ∈ dims ? grd.δx_3D : 1) .* (3 ∈ dims ? grd.δz_3D : 1)
    return dropdims(nansum(x3D .* w3D .* mask3D, dims=dims) ./ any(.!isnan.(mask3D), dims=dims), dims=dims)
end
nansum(x; kwargs...) = sum(x .* .!isnan.(x); kwargs...)
∫dxdydz(x3D, grd, mask3D=iswet(grd)) = integral(x3D, grd, mask3D)
∫dx(x3D, grd, mask3D=iswet(grd))     = integral(x3D, grd, mask3D; dims=2)
∫dy(x3D, grd, mask3D=iswet(grd))     = integral(x3D, grd, mask3D; dims=1)
∫dz(x3D, grd, mask3D=iswet(grd))     = integral(x3D, grd, mask3D; dims=3)
∫dxdy(x3D, grd, mask3D=iswet(grd))   = integral(x3D, grd, mask3D; dims=(1,2))
∫dxdz(x3D, grd, mask3D=iswet(grd))   = integral(x3D, grd, mask3D; dims=(2,3))
∫dydz(x3D, grd, mask3D=iswet(grd))   = integral(x3D, grd, mask3D; dims=(1,3))

# integral convenience aliases
totalintegral = ∫dxdydz
∫dV = ∫dxdydz
∫d³r = ∫dxdydz
zonalintegral = ∫dx
meridionalintegral = ∫dy
verticalintegral = ∫dz
horizontalintegral = ∫dxdy
export ∫dxdydz, ∫dV, ∫d³r, ∫dx, ∫dy, ∫dz, ∫dxdy, ∫dydz, ∫dxdz

# Create mean/average functions from integral functions
for f in (:horizontal, :vertical, :meridional, :zonal, :total)
    fint = Symbol(string(f, "integral"))
    fmean = Symbol(string(f, "mean"))
    favg = Symbol(string(f, "average"))
    @eval begin
        $fmean(x, grd, args...; kwargs...) = $fint(x, grd, args...; kwargs...) ./ $fint(1.0, grd, args...; kwargs...)
        $favg = $fmean
        export $favg
    end
end

# Allow all the integral/average functions to work when `x` is supplied as a vector
for f in (:∫dxdy,    :horizontalmean, :horizontalslice,
          :∫dz,      :verticalmean,
          :∫dy,      :meridionalmean, :meridionalslice,
          :∫dx,      :zonalmean,      :zonalslice,
          :∫dxdydz,  :totalmean)
    @eval begin
        $f(x::Vector, grd, mask=1; kwargs...) = $f(rearrange_into_3Darray(x, grd), grd, rearrange_into_3Darray(mask, grd); kwargs...)
        export $f
    end
end

