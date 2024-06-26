

const AVec = AbstractVector
const AMat = AbstractMatrix
const AArr = AbstractArray

"""
    iswet(grd)

Returns the 3D BitArray of wet boxes of the grid.
"""
iswet(g::T) where {T<:OceanGrid} = g.wet3D
export iswet

"""
    latvec(grd)

Returns the vector (wet points) of latitudes (units stripped).
"""
latvec(g::T) where {T<:OceanGrid} = ustrip.(g.lat_3D[iswet(g)])

"""
    lonvec(grd)

Returns the vector (wet points) of longitudes (units stripped).
"""
lonvec(g::T) where {T<:OceanGrid} = ustrip.(g.lon_3D[iswet(g)])

"""
    depthvec(grd)

Returns the vector (wet points) of depths (units stripped).
"""
depthvec(g::T) where {T<:OceanGrid} = ustrip.(g.depth_3D[iswet(g)])

"""
    topdepthvec(grd)

Returns the vector (wet points) of top depths (units stripped).
(top depths = depth of the top of the box.)
"""
topdepthvec(g::T) where {T<:OceanGrid} = ustrip.(g.depth_top_3D[iswet(g)])

"""
    bottomdepthvec(grd)

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
    rearrange_into_3Darray(x, grd)

Returns a 3D array of `x` rearranged to the wet boxes of the grid.
"""
rearrange_into_3Darray(x, grd) = rearrange_into_3Darray(x, grd.wet3D)

"""
    vectorize(x3D, grd)

Returns the 1D vector corresponding to the 3D field of `x3D` at the wet boxes of `grd`.

This function essentially returns `x3D[iswet(grd)]` but with some additional checks.
"""
function vectorize(x3D::AArr{T,3} where T, grd)
    size(x3D) ≠ size(grd) && error("The size of your 3D array does not match the grid size")
    x3D[iswet(grd)]
end

export rearrange_into_3Darray, vectorize

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
function Interpolations.interpolate(x, g::T) where {T<:OceanCurvilinearGrid}
    # This is weird: I am overloading Interpolations.interpolate with a NearestNeighbors function,
    # but it's to keep things working. TODO: Think of a better way to do all this
    brutetree = BruteTree([lonvec(g)'; latvec(g)'; depthvec(g)'])
    return (y,x,z) -> nn(brutetree, [x'; y'; z'])[1]
end
function Interpolations.interpolate(x, g::OceanGrid, lat, lon, depth; itp=interpolate(x,g))
    return itp(ustrip.(lat), ustrip.(lon), ustrip.(depth))
end
function Interpolations.interpolate(x, g::OceanGrid, lats::Vector, lons::Vector, depths::Vector; itp=interpolate(x,g))
    return [itp(y,x,z) for (y,x,z) in zip(ustrip.(lats), ustrip.(lons), ustrip.(depths))]
end
function Interpolations.interpolate(x, g::OceanGrid, obs; itp=interpolate(x,g))
    return interpolate(x, g, obs.lat, obs.lon, obs.depth; itp=itp)
end
Interpolations.interpolate(x, g::OceanGrid, ::Missing; itp=nothing) = missing
export interpolate

"""
    iswet(g, lat, lon, depth)
    iswet(g, lats, lons, depths)
    iswet(g, obs)

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
    interpolationmatrix(g, obs)

Returns the sparse matrix `M` such that `M*x` is the model vector `x` interpolated
onto the provided locations in the obs table.

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

periodic_longitude(x3D::AArr{T,3}) where T = view(x3D,:,[1:size(x3D,2); 1],:)
periodic_longitude(v::AVec{T}) where T = [v; v[1]+360]

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
    buildIbelow(grd)

Builds the shifted-diagonal sparse matrix of the indices of below neighbours for `grid`.
See `buildIbelow(wet3D, iwet)`.
"""
buildIbelow(grd) = buildIbelow(grd.wet3D, indices_of_wet_boxes(grd))
export buildIbelow

"""
    buildIabove(wet3D, iwet)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours.

Iabove[i,j] = 1 if the box represented by the linear index i
lies directly above the box represented by the linear index j.
"""
buildIabove(wet3D, iwet) = copy(transpose(buildIbelow(wet3D, iwet)))
"""
    buildIabove(grd)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours for `grd`.
See `buildIabove(wet3D, iwet)`.
"""
buildIabove(grd) = buildIabove(grd.wet3D, indices_of_wet_boxes(grd))
export buildIabove



#===================================
Functions returning useful arrays,
vectors, and constants
===================================#

"""
    array_of_volumes(grd)

Returns the 3D array of volumes of grd boxes.
"""
array_of_volumes(grd) = grd.volume_3D
export array_of_volumes


"""
    indices_of_wet_boxes(wet3D::BitArray)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(wet3D::BitArray) = findall(vec(wet3D))
"""
    indices_of_wet_boxes(grd)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(grd) = indices_of_wet_boxes(grd.wet3D)
export indices_of_wet_boxes

"""
    number_of_wet_boxes(wet3D::BitArray)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(wet3D::BitArray) = sum(wet3D)
"""
    number_of_wet_boxes(grd)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(grd) = number_of_wet_boxes(grd.wet3D)
export number_of_wet_boxes

"""
    vector_of_volumes(grid)

Returns the vector of volumes of wet boxes.
"""
vector_of_volumes(grd) = array_of_volumes(grd)[indices_of_wet_boxes(grd.wet3D)]
volumevec(grd) = ustrip.(vector_of_volumes(grd))
export vector_of_volumes, volumevec

"""
    vector_of_depths(grd)

Returns the vector of depths of the center of wet boxes.
"""
vector_of_depths(grd) = grd.depth_3D[indices_of_wet_boxes(grd.wet3D)]
export vector_of_depths

"""
    vector_of_top_depths(grd)

Returns the vector of depths of the top of wet boxes.
"""
vector_of_top_depths(grd) = grd.depth_top_3D[indices_of_wet_boxes(grd.wet3D)]
export vector_of_top_depths

"""
    vector_of_bottom_depths(grd)

Returns the vector of depths of the bottom of wet boxes.
"""
vector_of_bottom_depths(grd) = 2vector_of_depths(grd) - vector_of_top_depths(grd)
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
convertlon(x::Real) = mod(x * °, 360°)
convertlon(x::Quantity) = unit(x) == ° ? mod(x, 360°) : error("Not a valid lon")
# lat helper functions
findlatindex(lat::Quantity, grd) = findmin(abs.(grd.lat .- lat))[2]
convertlat(y::Real) = y * °
convertlat(y::Quantity) = unit(y) == ° ? y : error("Not a valid lat")
export convertlat, convertlon
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
    return PermutedDimsArray(view(x3D,:,ix,:), (2,1))
end
"""
    zonalslice(x, grd; lat)

Returns a zonal slice of tracer `x` at latitidue `lat`.
"""
function zonalslice(x3D, grd; lat=nothing)
    isnothing(lat) && error("you must specify the latitude, e.g., `lat=-30`")
    lat = convertlat(lat)
    iy = findlatindex(lat, grd)
    return PermutedDimsArray(view(x3D,iy,:,:), (2,1))
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
    iy, ix = findlatindex(lat, grd), findlonindex(lon, grd)
    return view(x3D,iy,ix,:)
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
∫dxdydz(x3D, grd, mask3D=iswet(grd)) = integral(x3D, grd, mask3D)[1]
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
    fstd = Symbol(string(f, "std"))
    @eval begin
        $fmean(x, grd, args...; kwargs...) = $fint(x, grd, args...; kwargs...) ./ $fint(1.0, grd, args...; kwargs...)
        $favg = $fmean
        $fstd(x, grd, args...; kwargs...) = sqrt.($favg(abs2.(x - paint3D($favg(x, grd, args...; kwargs...), grd)), grd, args...; kwargs...))
        export $favg
    end
end

"""
    paint3D(x, grd)

Returns a 3D array of the size of `grd` that copies the values of `x` along the dimensions it lacks.

This is useful, e.g., when `x` is a 1D profile, a 2D map, or a 2D zonal average and you want to "copy"
it onto the full 3D array.
"""
function paint3D(x, grd)
    length(unique(size(grd))) ≠ length(size(grd)) && @warn "paint3D might not work because of the grid size"
    dimsout = [findfirst(nx .== collect(size(grd))) for nx in size(x)]
    szin, szout = [1,1,1], collect(size(grd))
    szout[dimsout] .= 1
    szin[dimsout] .= size(grd)[dimsout]
    return repeat(reshape(x, Tuple(szin)), outer=szout)
end
paint3D(x::Number, grd) = x * iswet(grd)
paint3D(x::AArr{T,3} where T, grd) = x

# Allow all the integral/average functions to work when `x` is supplied as a vector
for f in (:∫dxdy,    :horizontalmean, :horizontalstd,
          :∫dz,      :verticalmean,   :verticalstd,
          :∫dy,      :meridionalmean, :meridionalstd,
          :∫dx,      :zonalmean,      :zonalstd,
          :∫dxdydz,  :totalmean,      :totalstd)
    @eval begin
        $f(x::AVec, grd, mask=1; kwargs...) = $f(rearrange_into_3Darray(x, grd), grd, rearrange_into_3Darray(mask, grd); kwargs...)
        export $f
    end
end



#================
Vertical sections
================#

# Functions extending 3D arrays and 1D longitudinal vectors
# to be used by interpolation functions.
# These assume that the longitude is in (0,360) and that
# the longitude is stored in the 2nd dimension (lat,lon,depth)
function lonextend(x3D::AArr{T,3}) where T
    x3Dext = Array{T,3}(undef, (size(x3D) .+ (0,2,0))...)
    x3Dext[:,2:end-1,:] .= x3D
    x3Dext[:,1,:] .= x3D[:,end,:]
    x3Dext[:,end,:] .= x3D[:,1,:]
    return x3Dext
end
function cyclicallon(lon::AVec{T}) where T<:Quantity
    lonext = Vector{T}(undef, length(lon) + 2)
    lonext[2:end-1] .= lon
    lonext[1] = lon[end] - 360°
    lonext[end] = lon[1] + 360°
    return lonext
end
cyclicallon(lon::AVec) = cyclicallon(convertlon.(lon))
# TODO make the function below more generic and be capable of taking in other
# representations than just OceanographyCruises.jl `CruiseTrack`
"""
    verticalsection(x3D, grd; ct=nothing)

Returns a tuple of `(d, depth, x2D)` where the 2D array is the vertical section that follows the cruise track `ct` and `d` is the distance.
This is made to work with the OceanographyCruises.jl package.
"""
function verticalsection(x3D, grd; ct=nothing)
    isnothing(ct) && error("you must include the cruise track with `ct=...`")
    # Create the interpolation object
    x3D = lonextend(x3D)
    knots = ustrip.((grd.lat, cyclicallon(grd.lon), grd.depth))
    itp = interpolate(knots, x3D, Gridded(Constant()))
    # convert lat and lon (TODO check if necessary)
    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [convertlon(st.lon) for st in ct.stations]
    # compute distance
    distance = vcat(0.0u"km", cumsum(diff(ct)))
    # remove identical stations (lat,lon)
    idx = findall(diff(distance) .> 0u"km")
    push!(idx, length(ct))
    d = view(distance, idx)
    y = view(ctlats, idx)
    x = view(ctlons, idx)
    return d, grd.depth, [itp(ustrip.((lat, lon, d))...) for d in grd.depth, (lat,lon) in zip(y, x)]
end

function verticalsection2(x3D, grd; ct=nothing)
    isnothing(ct) && error("you must include the cruise track with `ct=...`")
    # convert lat and lon (TODO check if necessary)
    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [convertlon(st.lon) for st in ct.stations]
    bl = baselon(ctlons)
    ctlons = shiftlon(ctlons, bl=bl)
    # compute distance
    distance = vcat(0.0u"km", cumsum(diff(ct)))
    # remove identical stations (lat,lon)
    idx = findall(diff(distance) .> 0u"km")
    push!(idx, length(ct))
    d = view(distance, idx)
    y = view(ctlats, idx)
    x = view(ctlons, idx)
	# intersections
    xe = lonedges(grd)
    ye = latedges(grd)
	dix = intersections(d, x, shiftlon(xe, bl=bl))
	diy = intersections(d, y, ye)
	ext_d = [d[1]; sort(vcat(dix, diy)); d[end]] # extended distance with intersections values
	# mid points (to get z value between intersections)
	mid_d = [0.5(ext_d[i] + ext_d[i+1]) for i in 1:length(ext_d)-1]
	itpx = LinearInterpolation(d, x)
    itpy = LinearInterpolation(d, y)
	x2 = itpx(mid_d)
    y2 = itpy(mid_d)
	ix = [searchsortedfirst(xe, shiftlon(x, bl=0°), lt=≤) - 1 for x in x2]
    iy = [searchsortedfirst(ye, y) - 1 for y in y2]
    return ext_d, depthedges(grd), [x3D[j,i,k] for k in 1:length(grd.depth), (j,i) in zip(iy,ix)]
end

# aux functions for vertical section delimited by edges of boxes
edges(x::AVec, dx::AVec) = [x .- 0.5 .* dx; x[end] + 0.5dx[end]]
edges(grd::OceanGrid) = lonedges(grd), latedges(grd), depthedges(grd)
lonedges(grd::OceanGrid) = edges(grd.lon, grd.δlon)
latedges(grd::OceanGrid) = edges(grd.lat, grd.δlat)
depthedges(grd::OceanGrid) = edges(grd.depth, grd.δdepth)

# Function to make δlon from lon (approximation)
function edges_from_midpoints(midpoints::Vector{T}, lims) where T
	N = length(midpoints)
    I_left = [I zeros(N)]
    I_right = [zeros(N) I]
    A = (I_left + I_right) / 2
    M = [A; 1 zeros(N)'; zeros(N)' 1]
    edges = M \ [midpoints; collect(lims)]
    return edges
end

#======================================================
Functions for plotting scatter transect on top of section heatmap
======================================================#
"""
    intersections(x::Vector{T}, y, ylevs)

Returns an array of values in `x` coordinates of the intersections
of the linearly interpolated plot of `y` as a function of `x` with the `ylevs`.
"""
function intersections(x::AArr{T}, y, ylevs) where T
	out = T[]
	for i in 1:length(x) - 1
		mini_x = view(x, i:i+1)
		mini_y = view(y, i:i+1)
		idx = sortperm(mini_y)
		a, b = view(mini_y, idx)
        if a ≠ b
		    itp = LinearInterpolation(mini_y[idx], mini_x[idx], extrapolation_bc=Line())
		    for x_intersect in itp(ylevs[a .< ylevs .≤ b])
		    	push!(out, x_intersect)
		    end
        end
	end
	return out
end
shiftlon(x; bl=baselon(x)) = @. mod(mod(x - bl, 360°), 360°) + bl # enforces 0 ≤ m < 360 (https://github.com/JuliaLang/julia/issues/36310)
baselon(x) = (any(0° .≤ x .< 90°) && any(270° .≤ x .< 360°)) ? -180.0° : 0.0°




#=======================================================
Allow functions to work when `x` is supplied as a vector
=======================================================#

for f in (:horizontalslice, :meridionalslice, :zonalslice, :verticalsection, :verticalsection2)
    @eval begin
        $f(x::AVec, grd; kwargs...) = $f(rearrange_into_3Darray(x, grd), grd; kwargs...)
        export $f
    end
end



#=======================================================
Functions to grid/interpolate
=======================================================#

function lonextend(x2D::AArr{T,2}) where T
    x2Dext = Array{eltype(x2D),2}(undef, (size(x2D) .+ (0,2))...)
    x2Dext[:,2:end-1] .= x2D
    x2Dext[:,1] .= x2D[:,end]
    x2Dext[:,end] .= x2D[:,1]
    return x2Dext
end

"""
    regrid(x2D, lat, lon, grd)
    regrid(x3D, lat, lon, depth, grd)

Returns `x2D` (or `x3D`) interpolated onto `grd`.
"""
function regrid(x2D::AArr{T,2} where T, lat, lon, grd)
    size(x2D) ≠ (length(lat), length(lon)) && error("Dimensions of input and lat/lon don't match")
    # need to rearrange in case data was in (-180, 180)
    lon = convertlon.(lon)
    ilon = sortperm(lon)
    sort!(lon)
    x2D = view(x2D, :, ilon)
    knots = convertlat.(lat), cyclicallon(lon)
    x2D = lonextend(x2D)
    itp = LinearInterpolation(knots, x2D, extrapolation_bc = Flat())
    return [itp(y, x) for y in grd.lat, x in grd.lon]
end
function regrid(x3D::AArr{T,3} where T, lat, lon, depth, grd; interpolate_nans=false)
    size(x3D) ≠ (length(lat), length(lon), length(depth)) && error("Dimensions of input and lat/lon/depth don't match")
    if interpolate_nans
        any(isnan, x3D) && @warn "Be aware that NaNs will propagate with interpolation!"
    else
        x3D = inpaint(x3D)
        @warn "NaNs have been inpainted!"
    end
    x3D = lonextend(x3D)
    knots = convertlat.(lat), cyclicallon(lon), convertdepth.(depth)
    itp = LinearInterpolation(knots, x3D, extrapolation_bc = Flat())
    return [itp(y, x, z) for y in grd.lat, x in grd.lon, z in grd.depth]
end

#TODO Rename to rebin
function regrid(vs::AVec{T}, lats, lons, depths, grd) where T
    out = zeros(T, count(iswet(grd)))
    for (i, v) in zip(regrid_indices(lats, lons, depths, tree(grd)), vs)
        out[i] += v
    end
    return out
end
function regrid_indices(lats, lons, depths, tree)
    nind = knn(tree, [mod.(lons, 360)'; lats'; depths'], 1)[1]
    return [ind[1] for ind in nind]
end
tree(grd) = KDTree([lonvec(grd)'; latvec(grd)'; depthvec(grd)'])

export regrid


function regridandpaintsurface(x2D::AArr{T,2} where T, lat, lon, grd)
    x2D2 = regrid(x2D, lat, lon, grd)
    x3D = zeros(size(grd))
    x3D[:,:,1] .= x2D2
    return vectorize(x3D, grd)
end
export regridandpaintsurface

"""
    zcumsum(x, grd)

Returns the cumulative sum of `x` along the `z` dimension.
"""
function zcumsum(x, grd)
    x3D = rearrange_into_3Darray(x, grd)
    return vectorize(cumsum(x3D, dims=3), grd)
end
export zcumsum
