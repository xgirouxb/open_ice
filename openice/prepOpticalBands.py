# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Collection of functions required to prepare optical bands for analyses,
# and add bands or indices specific to optical imagery
# 
# Functions designed for Landsat 8 and Sentinel 2 imagery
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Setup required parameters

# Bands for this script
l7Bands  =   [  'B1',    'B2',  'B3',  'B4',    'B5',    'B7'] # Landsat 7
l8Bands  =   [  'B2',    'B3',  'B4',  'B5',    'B6',    'B7'] # Landsat 8
s2Bands  =   [  'B2',    'B3',  'B4',  'B8',   'B11',   'B12'] # Sentinel 2
stdBands =   ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'] # Common band names

# ---------------------------------------------------------------------------- #
# Rename bands, mask out bad values, and gap fill Landsat 7 TOA
def prepL7TOA(img):
    '''
    This function prepares Landsat 7 Top of Atmosphere images
    for further analysis by masking edge pixels with incomplete
    data, gap filling 60m of L7 SLC off data, and renaming bands.
    '''
    # Rename bands
    img = img.select(l7Bands, stdBands)
  
    # Mask out bad reflectance values (< 0 or >1)
    toa = img.updateMask(img.gte(0).And(img.lte(1)))
  
    # Mask out incomplete pixels
    toa = toa.updateMask(toa.mask().reduce(ee.Reducer.min()))
  
    # Fill L7 SLC off gaps using focal mean
    toa = ee.Image(toa.focalMean(2, 'circle', 'pixels', 8)).blend(toa)

    # Return clipped and masked, set product property
    return img.addBands(toa, stdBands, True)\
              .clip(toa.geometry().buffer(-2000))\
              .set('product', 'l7toa')

# ---------------------------------------------------------------------------- #
# Scale, rename, and mask out bad values in Landsat 8 TOA
def prepL8TOA(img):
    """
    This function prepares Landsat 8 Top of Atmosphere images
    for further analysis by scaling bands B2 to B7 to
    a 0 - 1 range and masking out bad pixels.
    """
    # Select/rename bands
    img = img.select(l8Bands, stdBands);
  
    # Mask out bad reflectance values (< 0 or >1)
    toa = img.updateMask(img.gte(0).And(img.lte(1)))
  
    # Mask out incomplete pixels
    toa = toa.updateMask(toa.mask().reduce(ee.Reducer.min()))
  
    # Return, copy properties, set product property
    return img.addBands(toa, stdBands, True)\
              .set('product', 'l8toa')
        
# ---------------------------------------------------------------------------- #
# Scale, rename, and mask out bad values in Sentinel 2 TOA
def prepS2TOA(img):
    """
    This function prepares Sentinel 2 Top of Atmosphere Reflectance images
    for further analysis by scaling bands B2, B3, B4, B8, B11, and B12 to
    a 0 - 1 range and masking out bad pixels.
    """
    # Select/rename bands
    img = img.select(s2Bands, stdBands)
  
    # Mask out bad reflectance values (< 0 or >10000)
    toa = img.updateMask(img.gte(0).And(img.lte(10000)))
  
    # Mask out incomplete pixels
    toa = toa.updateMask(toa.mask().reduce(ee.Reducer.min()))
  
    # Scale to range 0 - 1
    toa = toa.unitScale(0, 10000)
                
    # Return, copy properties, set product property
    return img.addBands(toa, stdBands, True)\
              .set('product', 's2toa')

# ---------------------------------------------------------------------------- #
# Add Normalized Difference Snow Index
def addNDSI(img):
    """
    This function adds the Normalized Difference Snow Index 
    (aka modified Normalized Difference Water Index, Xu 2006)
    as a band to optical images with standard band names (green and swir1).
    The values are scaled to [0, 1] range for use in classifiers.   
    """
    # Compute NDSI, scale to [0-1] range for classifiers
    ndsi = img.normalizedDifference(['green','swir1'])\
              .unitScale(-1, 1)\
              .rename('ndsi')
    return img.addBands(ndsi)