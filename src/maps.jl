################################################################################

"""Abstract type representing the ordering of pixels in a Healpix map.
See also `Ring` and `Nest`.
"""
abstract type Order end

"""The `Ring` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `Nest`.)
"""
abstract type Ring <: Order end

"""The `Nest` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `Ring`.)
"""
abstract type Nest <: Order end

"""
    AbstractMap{T, O} <: AbstractArray{T, 1}

An abstract type representing an Healpix map without a specified
ordering. This can be used to implement multiple dispatch when you
don't care about the ordering of a map."""
abstract type AbstractMap{T, O} <: AbstractArray{T, 1} end
abstract type OrderedMap{T, O <: Order} <: AbstractMap{T, O} end

"""
    struct Map{T, O <: Order} <: AbstractMap{T, O}

A Healpix map. The type `T` is used for the value of the pixels in
a map, and it can be anything (even a string!). The type `O` is used
to specify the ordering of the pixels, and it can either be
`Ring` or `Nest`.

A `Map` type contains the following fields:

- `pixels`: array of pixels
- `resolution`: instance of a `Resolution` object

"""
struct Map{T, O} <: OrderedMap{T, O}
    pixels::Array{T}
    resolution::Resolution
end

"""
Create a map with the specified array of pixels, without redundancy of above
syntax.
"""
function Map{O}(healpixels::Vector{T}) where {T, O}
	Map{T, O}(healpixels, Resolution(npix2nside(length(healpixels))))
end

Map(healpixels) = Map{Nest}(healpixels)

function Map{T, O}(nside::Integer) where {T, O}
	Map{O}(Vector{T}(undef, nside2npix(nside)))
end
Map{T}(nside::Integer) where T = Map{T, Nest}(nside)



###########################################################################
# Map & MaskedMap initialization functions

struct MaskedMap{T, O} <: OrderedMap{T, O}
	pixels::Vector{T}
	inds::Vector{Int64}
	resolution::Resolution
end

# Convenience constructor that infers type
function MaskedMap{O}(pixels::Vector{T}, inds, resolution) where {T,O}
    MaskedMap{T,O}(pixels, inds, resolution)
end

# Convenience constructor that defualts to nest
function MaskedMap(pixels::Vector{T}, inds, resolution) where {T,O}
    MaskedMap{Nest}(pixels, inds, resolution)
end

# Unintialized constructor
function MaskedMap{T, O}(inds, resolution) where {T,O}
	MaskedMap{O}(Vector{T}(undef, length(inds)), inds, resolution)
end

# Unintialized constructor that defualts to Nest
function MaskedMap{T}(inds, resolution) where {T}
	MaskedMap(Vector{T}(undef, length(inds)), inds, resolution)
end

# construct from map
function MaskedMap(m::Map{T,O}, inds) where {T,O}
    MaskedMap{O}(m.pixels[inds], inds, m.resolution)
end

function MaskedMap(m::Map{T, O}, θlims, ϕlims) where {T,O}
	inds = Vector{Int64}()
	for i in eachindex(m.pixels)
		θ, ϕ  = get_latlong_deg(m, i)
		if (θlims[1] <= θ <= θlims[2]) & (ϕlims[1] <= ϕ <= ϕlims[2])
			push!(inds, i)
		end
	end
	MaskedMap(m, inds)
end

# MaskedMap to Map
function Map(mm::MaskedMap{T, O}) where {T, O}
	healpixels = fill(NaN, mm.resolution.numOfPixels)
	healpixels[mm.inds] = mm.pixels
	Map{O}(healpixels)
end

################################################################################
# Base definitions
import Base: size, IndexStyle, getindex, setindex!, eachindex, similar

size(m::AbstractMap{T, O}) where {T, O} = (length(m.pixels),)

IndexStyle(::Type{<:AbstractMap{T, O}}) where {T, O} = IndexLinear()

function getindex(m::AbstractMap{T, O}, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i]
end

function setindex!(m::AbstractMap{T, O}, val, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i] = val
end

eachindex(mm::MaskedMap) = mm.inds

similar(m::Map{T,O}) where {T,O} = Map{T,O}(m.resolution.nside)
similar(m::Map{T,O}, nside::Integer) where {T,O} = Map{T,O}(nside)

similar(mm::MaskedMap{T,O}) where {T,O} = MaskedMap{T,O}(mm.inds, mm.resolution)
similar(mm::MaskedMap{T,O}, nside::Integer) where {T,O} = similar(resize(mm,
                                                                  nside))

################################################################################


import Base: +, -, *, /

+(a::AbstractMap{T, O}, b::AbstractMap{T, O}) where {T <: Number, O} = AbstractMap{T, O}(a.pixels .+ b.pixels)
-(a::AbstractMap{T, O}, b::AbstractMap{T, O}) where {T <: Number, O} = AbstractMap{T, O}(a.pixels .- b.pixels)
*(a::AbstractMap{T, O}, b::AbstractMap{T, O}) where {T <: Number, O} = AbstractMap{T, O}(a.pixels .* b.pixels)
/(a::AbstractMap{T, O}, b::AbstractMap{T, O}) where {T <: Number, O} = AbstractMap{T, O}(a.pixels ./ b.pixels)

+(a::AbstractMap{T, O}, b::Number) where {T <: Number, O} = AbstractMap{T, O}(a.pixels .+ b)
-(a::AbstractMap{T, O}, b::Number) where {T <: Number, O} = a + (-b)
*(a::AbstractMap{T, O}, b::Number) where {T <: Number, O} = AbstractMap{T, O}(a.pixels .* b)
/(a::AbstractMap{T, O}, b::Number) where {T <: Number, O} = AbstractMap{T, O}(a.pixels ./ b)

+(a::Number, b::AbstractMap{T, O}) where {T <: Number, O} = b + a
-(a::Number, b::AbstractMap{T, O}) where {T <: Number, O} = b + (-a)
*(a::Number, b::AbstractMap{T, O}) where {T <: Number, O} = b * a
/(a::Number, b::AbstractMap{T, O}) where {T <: Number, O} = AbstractMap{T, O}(a ./ b.pixels)

###########################################################################
# reordering and coordinate querying


function reorder_inds(m::AbstractMap{T,O}) where {T,O}
    reorder_func = (O ≡ Ring ? ring2nest : nest2ring)
    sortperm([reorder_func(m.resolution, i) for i in eachindex(m)])
end

function reorder(m::Map{T,O}) where {T,O}
    Map{T, (O ≡ Ring ? Nest : Ring)}(m.pixels[reorder_inds(m)])
end
function reorder(mm::MaskedMap{T,O}) where {T,O}
    inds = reorder_inds(mm)
    MaskedMap{T, (O ≡ Ring ? Nest : Ring)}(mm.pixels[inds], inds, mm.resolution)
end

function get_latlong_deg(m::OrderedMap, ind)
	colat, lon = pix2ang(m, ind)
	rad2deg.((colat2lat(colat), lon))
end

function getcoords(m::OrderedMap)
	coords = zeros((2, length(m)))
	for (i, healpixel_ind) in enumerate(eachindex(m))
		coords[:,i] .= get_latlong_deg(m, healpixel_ind)
	end
	coords
end
function getcoords(inds, nside)
	coords = zeros((2, length(inds)))
	for (i, healpixel_ind) in enumerate(inds)
		coords[:,i] .= get_latlong_deg(m, healpixel_ind)
	end
	return coords
end

function resize(m::Map{T,O}, nside) where {T, O}
	m[isequal.(m,NaN)] .= UNSEEN
	m_resized = Map{O}(healpy.ud_grade(m.pixels, nside,
									   order_in=(O ≡ Ring ? "RING" : "NEST")))
	m_resized[isequal.(m_resized,UNSEEN)] .= NaN
	return m_resized
end

function resize(mm::MaskedMap{T,O}, nside) where {T, O}
	resized_map = resize(Map(mm), nside)
	resized_inds = findall(pix -> !isequal(pix, NaN), resized_map)
	MaskedMap{T,O}(resized_map[resized_inds], resized_inds, mm.resolution)
end

###########################################################################

"""
    conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},
                                                 map2::Map{S, O2}) -> Bool

Determine if two Healpix maps are "conformables", i.e., if their
shape and ordering are the same.
"""
conformables(map1::OrderedMap{T, Ring}, map2::OrderedMap{S, Ring}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::OrderedMap{T, Nest}, map2::OrderedMap{S, Nest}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::OrderedMap{T, O1},
             map2::OrderedMap{S, O2}) where {T, S, O1 <: Order, O2 <: Order} = false
