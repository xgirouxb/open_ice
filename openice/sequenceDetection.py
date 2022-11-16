# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Collection of functions to detect a sequence of transition between
# ice presence (1)/absence (0)

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Import required packages
import datetime

# ---------------------------------------------------------------------------- #
# Spring breakup detection function
# Function to iterate: for each 'image' in a time series, update 'iceBreakupImg'
# until an ice-water-water sequence is flagged
def breakupDate(image, iceBreakupImg):
    """
    This function is passed to iterate(), iceBreakupImg is updated by each
    image in a time series until an ice-water-water sequence is flagged
    INPUTS
    image: image in an image collection of ice-water classifications
           sorted in chronological order, requires bands 'classIce'
           (ice(1)/water(0) classification) and 't' (time as doy)
    iceBreakupImg: initial state defined above
    OUTPUT
    iceBreakupImg updated at each iteration, logs updated variables from 
    current image ('classIce', 't'), the 2 previous images at t-1 
    ('classIce_1', 't_1') and t-2 ('classIce_2', 't_2'), and the flag 
    indicating if an ice-water-water sequence was detected ('seqDetected')   
    """
    # IF classIce band in a given pixel of 'image' HAS data
    pixelHasData = image.select('classIce').mask()  
    # AND IF the iceBreakupImg HAS NOT detected an ice-water-water sequence in that pixel
    noBreakupDetected = ee.Image(iceBreakupImg).select('seqDetected').eq(1).Not()
    # then flag the pixel as okay to be updated
    pixels2Update = pixelHasData.add(noBreakupDetected).eq(2) 

    # update mask over all bands in current image
    currentImg = image.select(['classIce', 't']).updateMask(pixels2Update)

    # AT TIME 0 ------------------------------------------------- #
    # where the currentImg has pixels to contribute, update iceLag0
    iceLag0 = ee.Image(iceBreakupImg).select(['classIce', 't'])
    updatedIceLag0 = iceLag0.where(pixels2Update, currentImg)\
                            .rename(['classIce', 't'])

    # AT TIME -1 ------------------------------------------------ #
    # where the currentImg has pixels to contribute, update iceLag1
    iceLag1 = ee.Image(iceBreakupImg).select(['classIce_1', 't_1'],
                                             ['classIce', 't'])
    updatedIceLag1 = iceLag1.where(pixels2Update, iceLag0)\
                            .rename(['classIce_1', 't_1'])

    # AT TIME -2 ------------------------------------------------ #
    # where the currentImg has pixels to contribute, update iceLag2
    iceLag2 = ee.Image(iceBreakupImg).select(['classIce_2', 't_2'],
                                             ['classIce', 't'])
    updatedIceLag2 = iceLag2.where(pixels2Update, iceLag1)\
                            .rename(['classIce_2', 't_2'])

    # last step --> has an ice-water-water sequence been identified
    cond1 = updatedIceLag0.select('classIce').eq(0) # water at t0
    cond2 = updatedIceLag1.select('classIce_1').eq(0) # water at t_1
    cond3 = updatedIceLag2.select('classIce_2').eq(1) # ice at t_2
    # is each condition in sequence met ?
    seqDetected = cond1.add(cond2).add(cond3).eq(3)\
                       .toUint16()\
                       .rename(['seqDetected'])

    # reassemble updated bands to form iceBreakupImg
    updatedImg = updatedIceLag0.addBands(updatedIceLag1)\
                               .addBands(updatedIceLag2)\
                               .addBands(seqDetected)

    # return iceBreakupImg with updated pixels
    return ee.Image(iceBreakupImg).addBands(updatedImg,
                                            [ 't', 'classIce',
                                            't_1', 'classIce_1',
                                            't_2', 'classIce_2',
                                            'seqDetected'],
                                            True) # overwrite

# ---------------------------------------------------------------------------- #
# Function to prep image collection for change detection

def prep4ChangeDetection(img):
    
    # Get the day of year
    doy = img.date().getRelative('day', 'year').add(1)
    # Create image with doy as band, convert to unsigned 16 bit integer
    t = ee.Image(doy).toUint16().rename(['t'])
    
    # Keep ice, convert to integer
    return img.select(['classIce']).toUint16().rename(['classIce'])\
              .addBands(t)

# ---------------------------------------------------------------------------- #
# Function to prep image collection for change detection
def detectSpringBreakup(imgCol, poiStart, tile):
    """
    imgCol: image collection with ice-water classifications
    poiStart: date (datetime object 'YYYY-MM-DD') at which to start 
    tile: region of interest (ee.Geometry)
    """
    # Setup Image Collection for iterated function            
    imgCol = imgCol.map(prep4ChangeDetection).sort('system:time_start')

    # Setup the initial all ice image 'iceBreakupImg', start at poiStart
    date = poiStart.getRelative('day', 'year').add(1)
    firstImgDate = ee.Image(date).toUint16().rename(['t'])
    firstImgSeqDetect = ee.Image(0).toUint16().rename(['seqDetected'])                           
    allIceImg = ee.Image(1).toUint16().rename(['classIce']).addBands(firstImgDate)
    iceBreakupImg = allIceImg.addBands(allIceImg)\
                             .addBands(allIceImg)\
                             .addBands(firstImgSeqDetect)\
                             .clip(tile)
    
    # Iterate through imgCol, update iceBreakupImg using breakupDate()
    # Cast to image (return type of iterate is unknown)
    iceBreakupImg = ee.Image(imgCol.iterate(breakupDate, iceBreakupImg))
    
    # First observation after poiStart needs to be ice, if it is water mask it out
    firstObsIsIce = imgCol.select(['classIce'])\
                          .reduce(ee.Reducer.firstNonNull())\
                          .eq(1)
    iceBreakupImg = iceBreakupImg.updateMask(firstObsIsIce)
    
    # Return iceBreakupImg
    return iceBreakupImg