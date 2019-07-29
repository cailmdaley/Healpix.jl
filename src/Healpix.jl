module Healpix
###
export nsideok, nside2pixarea, nside2resol
export Resolution, nside2npix, npix2nside
export ang2pixNest, ang2pixRing, pix2angNest, pix2angRing
export vec2pixNest, vec2pixRing, pix2vecNest, pix2vecRing
export pix2ringpos
export Order, Ring, Nest
export AbstractMap, Map, OrderedMap, MaskedMap
export getcoords, resize, reorder
export readMaskedMap, saveMaskedMap, make_postage_stamps
export ang2vec, vec2ang, ang2pix, pix2ang
export readMapFromFITS, savePixelsToFITS, saveToFITS, conformables
export ringWeightPath, readWeightRing
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms, almIndexL0, almIndex, readAlmFromFITS
export getringinfo!, getringinfo, getinterpolRing
export pix2xyfRing, xyf2pixRing, pix2xyfNest, xyf2pixNest
export ring2nest, nest2ring

import FITSIO
import Base: getindex, setindex!
import DelimitedFiles: readdlm, writedlm

const NSIDE_MAX = 8192

using PyCall; const healpy = PyNULL()
function __init__()
    copy!(healpy, pyimport("healpy"))
end
export healpy

include("math.jl")
include("datatables.jl")
include("resolution.jl")
include("interp.jl")
include("xyf.jl")


include("maps.jl")
include("pixels.jl")
include("projections.jl")
include("alm.jl")
include("mapmaking.jl")
include("fits.jl")


end
