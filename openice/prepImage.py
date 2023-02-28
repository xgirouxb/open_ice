# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Collection of functions required to prepare imagery for analysis
# 
# Functions designed for Landsat 7 + 8 TOA and Sentinel 2 TOA
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Global water masking function

# Import global water JRC layer
globalWater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')\
                          .select('occurrence').gte(80)

# Import Large Scale International Boundary Polygons for NA, convert to raster
filterNorthAm = ee.Filter.inList('country_co', ['CA', 'US'])
northAm = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017").filter(filterNorthAm)
northAmCoastline = ee.Image(0).byte().paint(northAm, 1)

# Water mask function
def waterMask(img):
    """
    This function updates an image's mask so it only contains freshwater pixels.
    Uses permanent water (occurence 80-100 %) in the JRC global surface water dataset
    and removes water pixels in marine coastal areas (using political boundary coastlines) 
    see manual: https://storage.googleapis.com/global-surface-water/downloads_ancillary/DataUsersGuidev2021.pdf
    """
    return img.updateMask(globalWater.updateMask(northAmCoastline))

# ---------------------------------------------------------------------------- #
# Tile mask function
def tileMask(imgCol, tile):
    """
    This function masks out the areas outside of the provided tile (geometry),
    for each img in an imgCol.
    """
    # Create mask of tile
    mask = ee.Image(0).byte().paint(tile, 1)
    # Function to map over collection
    def mapTileMask(img):
        # update image mask
        return img.updateMask(mask)
    # return imgCol
    return imgCol.map(mapTileMask)

# ---------------------------------------------------------------------------- #
# Add dates (day of year and year) as bands
def addDate(img):
    """
    This function adds both the year and the day of year (doy) as bands to an image.
    """
    # Get the day of year and year, convert to unsigned 16 bit integer
    doy = img.date().getRelative('day', 'year').add(1).toUint16()
    year = img.date().get('year').toUint16()
    # Create single banded image for each
    doyBand = ee.Image(doy).rename(['doy'])
    yearBand = ee.Image(year).rename(['year'])
    # Add bands to image
    return img.addBands([doyBand, yearBand])

# ---------------------------------------------------------------------------- #
# Add day of year as band
def addDoy(img):
    """
    This function adds the day of year (doy) as a bands to an image.
    """
    # Get the day of year, convert to unsigned 16 bit integer
    doy = img.date().getRelative('day', 'year').add(1).toUint16()
    # Create image with doy as band
    doyBand = ee.Image(doy).rename(['doy'])
    # Add band to image
    return img.addBands(doyBand)

# ---------------------------------------------------------------------------- #
# Add dates (day of year and year) as bands
def addYear(img):
    """
    This function adds the year as a band to an image.
    """
    # Get the year, convert to unsigned 16 bit integer
    year = img.date().get('year').toUint16()
    # Create single banded image
    yearBand = ee.Image(year).rename(['year'])
    # Add band to image
    return img.addBands(yearBand)
    
# ---------------------------------------------------------------------------- #
# Add fractional year since epoch as band
def addFracYear(img):
    """
    This function converts seconds since epoch to fractional year
    and adds it as a band to an image.
    """
    # Fractional years since the epoch
    fracYear = img.date().difference(ee.Date('1970-01-01'), 'year').add(1970)
    # Create single band image, rename
    fracYearBand = ee.Image(fracYear).float().rename(['fracYear'])
    # Add band
    return img.addBands(fracYearBand)

# ---------------------------------------------------------------------------- #
# Add date (YYYY-MM-dd) as property
def addDateProp(img):
    """
    This function adds the date (format YYYY-MM-dd) as a property to an image.
    """
    # Extract date from 'system:time_start'
    date = img.date().format('YYYY-MM-dd')
    # Add property
    return img.set('date', date)