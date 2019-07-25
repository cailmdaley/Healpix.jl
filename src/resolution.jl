# Definition of the composite type "Resolution"

using Printf

"""
    struct Resolution

`Resolution` objects are needed to perform a number of pixel-related
functions, e.g., convert a direction into a pixel number and vice
versa.

The fields of a `Resolution` object are the following:

- `nside`: the NSIDE parameter
- `nsideTimesTwo`: 2 * NSIDE
- `nsideTimesFour`: 4 * NSIDE
- `numOfPixels`: number of pixels in the map
- `order`: order of the map
- `pixelsPerFace`: number of pixels in each Healpix face
- `ncap`
- `fact2`
- `fact1`
"""
struct Resolution
    nside::Int
    nsideTimesTwo::Int
    nsideTimesFour::Int
    numOfPixels::Int

    order::Int
    pixelsPerFace::Int
    ncap::Int
    fact2::Float64
    fact1::Float64
end

# Customize printing
Base.show(io::IO, r::Resolution) = @printf(io, "Healpix resolution(NSIDE = %d)", r.nside)

################################################################################

"""
    Resolution(nside) -> Resolution

Create a `Resolution` object, given a value for `NSIDE`.
"""
function Resolution(nside)
    (1 ≤ nside ≤ NSIDE_MAX) || throw(DomainError(nside))
    # The expression (nside & (nside - 1)) == 0 is a quick check for
    # detecting if nside is a power of two or not
    (nside & (nside - 1) == 0) || throw(DomainError())

    order          = ilog2(nside)
    pixelsPerFace  = nside^2
    numOfPixels    = 12pixelsPerFace
    ncap           = 2 * nside * (nside - 1)
    fact2          = 4 / numOfPixels
    fact1          = 2 * nside * fact2

    result = Resolution(nside,
                        2nside,
                        4nside,
                        numOfPixels,
                        order,
                        pixelsPerFace,
                        ncap,
                        fact2,
                        fact1)

end

########################################################################

"""
    nsideok(nside::Integer) -> Bool

Check whether `nside` is a valid `NSIDE` parameter.
"""
nsideok(nside::Integer) = (nside > 0) && ((nside) & (nside - 1) == 0)

########################################################################

"""
    nside2npix(nside::Integer) -> Integer

Return the number of pixels for a Healpix map with the specified
`NSIDE` value. If `NSIDE` is not an integer power of two, the function
throws a `DomainError` exception.
"""
function nside2npix(nside::Integer)
    (nside > 0) || throw(DomainError(nside, "`NSIDE` is not a positive number"))

    nsidelog2 = round(Int, log2(nside))
    (2^nsidelog2 == nside) || throw(DomainError(nside, "`NSIDE` is not an integer power of two"))

    12(nside^2)
end

########################################################################

"""
    npix2nside(npix::Integer) -> Integer

Given the number of pixels in a Healpix map, return the `NSIDE`
resolution parameter. If the number is invalid, throw a `DomainError`
exception.
"""
function npix2nside(npix::Integer)
    (npix % 12 == 0) || throw(DomainError(npix, "Invalid number of pixels"))

    square_root = sqrt(npix / 12)
    (square_root^2 == npix / 12) || throw(DomainError(npix, "Invalid number of pixels"))

    convert(Int, round(square_root))
end

################################################################################

"""
    nside2pixarea(nside::Integer) -> Real

Return the solid angle of a pixel in a map with the specified `NSIDE` parameter.
The result is expressed in steradians.
"""
nside2pixarea(nside::Integer) = 4π / nside2npix(nside)

################################################################################

"""
    nside2resol(nside::Integer) -> Real

Return the approximate resolution of a map with the specified `NSIDE`. The
resolution is expressed in radians, and it is the square root of the pixel
size.
"""
nside2resol(nside::Integer) = sqrt(nside2pixarea(nside))

################################################################################
