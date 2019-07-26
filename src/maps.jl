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
    Map{T, O <: Order}(nside) -> Map{T, O}

Create an empty map with the specified NSIDE.
"""
function Map{T, O}(nside::Number) where {T, O}
	Map{T,O}(zeros(T, nside2npix(nside)), Resolution(nside))
end

function Map{T, O}(healpixels::Vector{T}) where {T, O}
	Map{T, O}(healpixels, Resolution(npix2nside(length(healpixels))))
end

"""
Create a map with the specified array of pixels, without redundancy of above
syntax.
"""
function Map{O}(healpixels::Vector{T}) where {T, O}
	Map{T, O}(healpixels, Resolution(npix2nside(length(healpixels))))
end

"""
Convenience constructor--default Nest
"""
Map(healpixels::Vector{T}) where T = Map{Nest}(healpixels)

###########################################################################
# Map & MaskedMap initialization functions

struct MaskedMap{T, O} <: OrderedMap{T, O}
	pixels::Vector{T}
	inds::Vector{Int64}
	resolution::Resolution
end

# MaskedMap to Map
function Map(mm::MaskedMap{T, O}) where {T, O}
	healpixels = fill(NaN, mm.resolution.numOfPixels)
	healpixels[mm.inds] = mm.pixels
	Map{O}(healpixels)
end

# Map to MaskedMap
function MaskedMap(m::Map{T, O}, θlims, ϕlims) where {T,O}
	inds = Vector{Int64}()
	for i in eachindex(m.pixels)
		θ, ϕ  = get_latlong_deg(m, i)
		if (θlims[1] <= θ <= θlims[2]) & (ϕlims[1] <= ϕ <= ϕlims[2])
			push!(inds, i)
		end
	end
	MaskedMap{T,O}(m.pixels[inds], inds, m.resolution)
end

# Pixel vector to MaskedMap
function MaskedMap{O}(healpixels::Vector{T}, θlims, ϕlims) where {T, O}
	MaskedMap(Map{O}(healpixels), θlims, ϕlims)
end

# Pixel vector to defaults to ringOrdering
function MaskedMap(healpixels::Vector{T}, θlims, ϕlims) where T <: Number
	MaskedMap{Nested}(healpixels, θlims, ϕlims)
end

###########################################################################
# reordering and coordinate querying

Base.eachindex(mm::MaskedMap) = mm.inds
function reorder_inds(m::AbstractMap{T,O}) where {T,O}
    reorder_func = (O ≡ Ring ? ring2nest : nest2ring)
    sortperm([reorder_func(m.resolution, i) for i in eachindex(m)])
end

function reorder(m::Map{T,O}) where {T,O}
    Map{T, (O ≡ Ring ? Nest : Ring)}(m.pixels[reorder_inds(m)])
end
function reorder(mm::MaskedMap{T,O}) where {T,O}
    inds = reorder_inds(mm)
    Map{T, (O ≡ Ring ? Nest : Ring)}(m.pixels[inds], inds, mm.resolution)
end

function get_latlong_deg(m::OrderedMap, ind)
	colat, lon = pix2ang(m, ind)
	rad2deg.((colat2lat(colat), lon))
end

function getcoords(m::OrderedMap)
	coords = zeros((2, length(m)))
	healpixel_inds = get_healpixel_inds(m)
	for i in eachindex(healpixel_inds)
		coords[:,i] .= get_latlong_deg(m, healpixel_inds[i])
	end
	coords
end

function resize(m::Map{T,O}, nside) where {T, O}
    order =
	Map{O}(healpy.ud_grade(m.pixels, nside, order_in=(O ≡ Ring ? "RING" : "NEST")))
end
function resize(mm::MaskedMap{T,O}, nside) where {T, O}
	θlims, ϕlims = extrema(getcoords(mm), dims=2)
	MaskedMap(resize(Map(mm), nside), θlims, ϕlims)
end

###########################################################################

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

################################################################################
# Iterator interface

Base.size(m::AbstractMap{T, O}) where {T, O} = (length(m.pixels),)

Base.IndexStyle(::Type{<:AbstractMap{T, O}}) where {T, O} = IndexLinear()

function getindex(m::AbstractMap{T, O}, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i]
end

function setindex!(m::AbstractMap{T, O}, val, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i] = val
end

################################################################################

"""
    readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column, t::Type{T})
    readMapFromFITS{T <: Number}(fileName::String, column, t::Type{T})

Read a Healpix map from the specified (1-base indexed) column in a
FITS file. The values will be read as numbers of type T. If the code
fails, FITSIO will raise an exception. (Refer to the FITSIO library
for more information.)
"""

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
