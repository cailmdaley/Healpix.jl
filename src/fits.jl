"""
    readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column, t::Type{T})
    readMapFromFITS{T <: Number}(fileName::String, column, t::Type{T})

Read a Healpix map from the specified (1-base indexed) column in a
FITS file. The values will be read as numbers of type T. If the code
fails, FITSIO will raise an exception. (Refer to the FITSIO library
for more information.)
"""

function readMapFromFITS(f::FITSIO.FITSFile, column, t::Type{T}) where
                        {T <: Number}
    value, comment = FITSIO.fits_read_keyword(f, "NSIDE")
    nside = parse(Int, value)

    value, comment = FITSIO.fits_read_keyword(f, "ORDERING")
    ringOrdering = uppercase(strip(value[2:end-1])) == "RING"

    repeat = (FITSIO.fits_get_coltype(f, column))[2]
    nrows = FITSIO.fits_get_num_rows(f)

    if repeat * nrows != nside2npix(nside)
        error("Wrong number of pixels in column $column of FITS file (NSIDE=$nside)")
    end

    if ringOrdering
        result = Map{T, Ring}(Array{T}(undef, nside2npix(nside)))
    else
        result = Map{T, Nest}(Array{T}(undef, nside2npix(nside)))
    end
    FITSIO.fits_read_col(f, column, 1, 1, result.pixels)

    result
end

function readMapFromFITS(fileName::AbstractString,
                         column,
                         t::Type{T}) where {T <: Number}
    f = FITSIO.fits_open_table(fileName)
    result = readMapFromFITS(f, column, t)
    FITSIO.fits_close_file(f)

    result
end

################################################################################

"""
    savePixelsToFITS(map::Map{T}, f::FITSIO.FITSFile, column) where {T <: Number}

Save the pixels of `map` into the column with index/name `column` in the FITS
file, which must have been already opened.
"""
function savePixelsToFITS(map::AbstractMap{T, O},
                          f::FITSIO.FITSFile,
                          column) where {T <: Number, O}

    FITSIO.fits_update_key(f, "PIXTYPE", "HEALPIX",
                           "HEALPIX pixelisation")
    FITSIO.fits_update_key(f, "NSIDE", map.resolution.nside,
                           "Value of NSIDE")
    FITSIO.fits_update_key(f, "FIRSTPIX", 1,
                           "First pixel (1 based)")
    FITSIO.fits_update_key(f, "LASTPIX", map.resolution.numOfPixels,
                           "Last pixel (1 based)")
    FITSIO.fits_update_key(f, "INDXSCHM", "IMPLICIT",
                           "Indexing: IMPLICIT or EXPLICIT")
    FITSIO.fits_write_col(f, column, 1, 1, map.pixels)

end

"""
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        f::FITSIO.FITSFile,
                                        column)
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        fileName::String,
                                        typechar="D",
                                        unit="",
                                        extname="MAP")

Save a Healpix map in the specified (1-based index) column in a FITS
file. If the code fails, FITSIO will raise an exception. (Refer to the
FITSIO library for more information.)
"""
function saveToFITS(map::OrderedMap{T, Ring},
                    f::FITSIO.FITSFile,
                    column) where {T <: Number}

    FITSIO.fits_update_key(f, "ORDERING", "RING")
    savePixelsToFITS(map, f, column)

end

function saveToFITS(map::OrderedMap{T, Nest},
                    f::FITSIO.FITSFile,
                    column) where {T <: Number}

    FITSIO.fits_update_key(f, "ORDERING", "NEST")
    savePixelsToFITS(map, f, column)

end

"""
    saveToFITS(map::Map{T, O}, filename::AbstractString, typechar="D", unit="", extname="MAP") where {T <: Number, O <: Order}

Save a map into a FITS file. The name of the file is specified in
`filename`; if it begins with `!`, existing files will be overwritten
without warning. The parameter `typechar` specifies the data type to
be used in the FITS file: the default (`D`) will save 64-bit
floating-point values. See the CFITSIO documentation for other
values. The keyword `unit` specifies the measure unit used for the
pixels in the map. The keyword `extname` specifies the name of the HDU
where the map pixels will be written.

"""
function saveToFITS(map::OrderedMap{T, O},
                    fileName::AbstractString;
                    typechar="D",
                    unit="",
                    extname="MAP") where {T <: Number, O <: Order}

    f = FITSIO.fits_create_file(fileName)
    try
        FITSIO.fits_create_binary_tbl(f, 0, [("PIXELS", "1$typechar", unit)], extname)
        saveToFITS(map, f, 1)
    finally
        FITSIO.fits_close_file(f)
    end

end

################################################################################
