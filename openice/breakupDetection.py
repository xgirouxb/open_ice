# ---------------------------------------------------------------------------- #
# Process estimates of pixel breakup date on a grid based system
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Import required packages
import prepImage as prepImg
import prepOpticalBands as prepOpt
import classifyIce as classIce
import logisticFilter as logisticFilter
import sequenceDetection as seqDetect

# ---------------------------------------------------------------------------- #
# Define breakup detection routine
def breakupDetection(tile, year, expDirectory, expFilename, cloudThresh = 90, globalWater = True, logFilter = True):
    
    """
    Main OPEN-ICE function, launches a Google Earth Engine task that exports
    spring breakup image.
    
    INPUT
    tile: ee.Geometry() of study area
    year: integer [2013, 2021], year of interest
    expDirectory: string, directory in which to save asset (e.g., 'myFolder')
    expFilename: string, name of created asset (e.g., 'myBreakupImg')
    cloudThresh: integer [0, 100], filter images from L7, L8, and S2 collection based on
                 cloud cover (default is exclude images with > 90% clouds)
    globalWater: logical, limit analysis to pixels with >80% persistence in JRC global
                 water layer (default is True, but should be manually set to False
                 for tiles with latitudes > 70N)
    logFilter: logical, apply logistic regression to each pixel and remove potential
               misclassifications, i.e., residuals > 0.85 (default is True)
               
    OUTPUT
    ee.Image() of spring breakup with following bands at 30 metre resolution
    - breakupDate: day of year pixel transitioned form ice to water
    - R2: r-squared of logistic temporal filter
    - nPixels: number of total ice and water obs between Feb 15th and Oct 1st of 'year'
    - breakupGap: days elapsed between last observed ice and first observed water
    - year: year of spring breakup
    """

    # ---------------------------------------------------------------------------- #
    # STEP 1: Setup parameters for script

    # Time parameters
    poiStart = ee.Date.fromYMD(year, 2, 15)
    poiEnd = ee.Date.fromYMD(year, 10, 1)

    # ---------------------------------------------------------------------------- #
    # STEP 2: Filter image collections, apply masks, prep bands
    
    # LANDSAT 7 TOP OF ATMOSPHERE REFLECTANCE ------------------------------------ #

    # Build a Landsat 7 TOA Collection for roi and poi
    l7toa = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterDate(poiStart, poiEnd)\
                                                         .filterBounds(tile)\
                                                         .filter(ee.Filter.lt('CLOUD_COVER', cloudThresh))
    
    # Mask for current tile
    l7toa = prepImg.tileMask(l7toa, tile)
    
    # If globalWater is True, apply global water mask using JRC layer 80% occurence
    if globalWater:
        l7toa = l7toa.map(prepImg.waterMask)
    
    # Prep Landsat 7 TOA bands
    l7toa = l7toa.map(prepOpt.prepL7TOA)\
                 .map(prepOpt.addNDSI)
                 
    # Repeat masking steps to remove data interpolated outside areas on interest 
    l7toa = prepImg.tileMask(l7toa, tile)
    if globalWater:
        l7toa = l7toa.map(prepImg.waterMask)

    # LANDSAT 8 SURFACE REFLECTANCE ---------------------------------------------- #

    # Build a Landsat 8 TOA collection for roi and poi
    l8toa = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterDate(poiStart, poiEnd)\
                                                         .filterBounds(tile)\
                                                         .filter(ee.Filter.lt('CLOUD_COVER', cloudThresh))

    # Mask for current tile
    l8toa = prepImg.tileMask(l8toa, tile)
    
    # If globalWater is True, apply global water mask using JRC layer 80% occurence
    if globalWater:
        l8toa = l8toa.map(prepImg.waterMask)
    
    # Prep Landsat 8 TOA bands
    l8toa = l8toa.map(prepOpt.prepL8TOA)\
                 .map(prepOpt.addNDSI)

    # SENTINEL 2 TOP OF ATMOSPHERE REFLECTANCE ----------------------------------- #

    # Build a Sentinel 2 TOA collection for roi and poi    
    s2toa = ee.ImageCollection('COPERNICUS/S2').filterDate(poiStart, poiEnd)\
                                               .filterBounds(tile)\
                                               .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudThresh))
    
    # Mask for current tile
    s2toa = prepImg.tileMask(s2toa, tile)
    
    # If globalWater is True, apply global water mask using JRC layer 80% occurence
    if globalWater:
        s2toa = s2toa.map(prepImg.waterMask)
    
    # Prep Sentinel 2 TOA 
    s2toa = s2toa.map(prepOpt.prepS2TOA)\
                 .map(prepOpt.addNDSI)
    
    # Remove duplicates (select most recently generated image)
    s2toa = ee.ImageCollection(s2toa.sort('GENERATION_TIME', False)\
                                    .distinct(['system:time_start', 'MGRS_TILE']))

    # ---------------------------------------------------------------------------- #
    # STEP 3: Add bands, classify water (0)/ ice (1), merge collections

    # Classify
    l7toa = l7toa.map(classIce.classIceL7TOA)
    l8toa = l8toa.map(classIce.classIceL8TOA)
    s2toa = s2toa.map(classIce.classIceS2TOA)

    # Merge both image collections
    imgCol = ee.ImageCollection(l7toa.merge(l8toa).merge(s2toa))\
                .sort('system:time_start', False)
    
    # ---------------------------------------------------------------------------- #
    # STEP 4: Temporal filter using logistic regression of ice ~ time

    # PREP FOR TIME SERIES ANALYSIS ---------------------------------------------- #
    
    # Apply logistic filter if logFilter == True
    
    if logFilter:
    
        # Sort in reverse chrono order
        imgCol = imgCol.sort('system:time_start', False)

        # Add fractional year for input as time variable in logistic regression
        imgCol = imgCol.map(prepImg.addFracYear)

        # FIT LOGISTIC REGRESSION, APPLY FILTER -------------------------------------- #
        # The logistic regression informs us on pixels that are outliers during
        # non-transition periods of ice phenology. Pixels with high residuals can
        # be masked out, as they constitute misclassifications, specifically 
        # water detected during winter, or ice detected during summer. The cutoff
        # value is high (0.85) so that possible misclassifications occuring during
        # breakup period are kept in case they represent a true succession of
        # changed states (e.g. ice breaks, then wind shifts ice back to that
        # pixel form anther area in the lake).

        # Fit logistic regression of ice ~ time
        imgCol = logisticFilter.fitLogisticRegression(imgCol, 'classIce', 'fracYear')

        # Mask out winter water and summer ice (high residuals)
        imgCol = imgCol.map(logisticFilter.temporalFilter)

        # Sort chronologically
        imgCol = imgCol.sort('system:time_start', True)

        # Count non-masked pixels in stack
        emptyImg = ee.Image(0).clip(tile).rename('classIce')
        nPixels = imgCol.select('classIce')\
                        .merge(emptyImg)\
                        .reduce(ee.Reducer.count())\
                        .subtract(1)\
                        .toUint16()\
                        .rename(['nPixels'])

        # Get R^2 of logistic filter
        rSquared = imgCol.select('R2')\
                         .reduce(ee.Reducer.firstNonNull())\
                         .multiply(100)\
                         .round()\
                         .toUint16()\
                         .rename(['R2'])
    
    else:
        
        # Count non-masked pixels in stack
        emptyImg = ee.Image(0).clip(tile).rename('classIce')
        nPixels = imgCol.select('classIce')\
                        .merge(emptyImg)\
                        .reduce(ee.Reducer.count())\
                        .subtract(1)\
                        .toUint16()\
                        .rename(['nPixels'])

        # Create empty image to add as R2 band
        rSquared = ee.Image().clip(tile)\
                             .toUint16()\
                             .rename(['R2'])

    # ---------------------------------------------------------------------------- #
    # STEP 5: Change detection of pixel transition

    # New sequence detection functions in module "sequenceDetection.py"
    iceBreakupImg = seqDetect.detectSpringBreakup(imgCol, poiStart, tile)

    # Compute date of pixel transition (sequence of ice-water-water)
    # average of t_1 and t_2, when:
    # classIce = 0
    # classIce_1 = 0  ---> pixel transition first obs date of open water
    # classIce_2 = 1
    # breakup date:
    breakupDate = ee.Image(iceBreakupImg)\
                    .select('t_1')\
                    .rename(['breakupDate'])\
                    .updateMask(ee.Image(iceBreakupImg).select('seqDetected').eq(1))
    
    breakupGap = ee.Image(iceBreakupImg).select('t_1')\
                    .subtract(ee.Image(iceBreakupImg).select('t_2'))\
                    .rename(['breakupGap'])\
                    .toUint16()\
                    .updateMask(ee.Image(iceBreakupImg).select('seqDetected').eq(1))
    
    # Add breakup date as band
    breakupDate = breakupDate.addBands(rSquared)\
                             .addBands(nPixels)\
                             .addBands(breakupGap)\
                             .set('year', year)

    # ---------------------------------------------------------------------------- #
    # STEP 6: Launch export task
    task = ee.batch.Export.image.toAsset(
        image = breakupDate, 
        description = expFilename,
        assetId = expDirectory + '/' + expFilename,
        maxPixels = 1e13,
        scale = 30
    )
    task.start()