# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Collection of functions to create dummy water to ensure logistic regression
# descends to 0 even pixel is dominated by bad observations

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Import required packages
import prepImage as prepImg

# ---------------------------------------------------------------------------- #
# Add 'dummy' property to image 
def addPropDummyYes(img):
    '''
    This function sets property 'dummy' to 'YES'.
    '''
    return img.set('dummy', 'YES')

def addPropDummyNo(img):
    '''
    This function sets property 'dummy' to 'NO'.
    '''
    return img.set('dummy', 'NO')

# ---------------------------------------------------------------------------- #
# Create dummy ice (1) or water (0) image and set timestamp
def dummyWaterImg(datestamp):
    '''
    This function creates a dummy image of water (0) and sets
    the 'system:time_start' property to the input datestamp.
    '''
    # Create dummy water (0) image
    dummyWater = ee.Image(0).toUint16().rename('classIce')
    
    # Return with timestamp
    return dummyWater.set('system:time_start', datestamp)

# ---------------------------------------------------------------------------- #
# Create sequence of water (0) images
def dummyCollection(year, tile, nDummy):
    '''
    This function creates an ImageCollection composed of nDummy water (0)images
    in the 4 weeks prior to poiEnd (Oct 1st). Images are masked to 'tile' extent,
    and to JRC global water layer (and North America coastline).
    '''
    
    # Create sequence of datestamps
    waterStamps = ee.List.sequence(
        start = ee.Date.fromYMD(year, 9, 1).millis(),
        end = ee.Date.fromYMD(year, 10, 1).millis(),
        count = nDummy
    )
    
    # Create ImageCollection of dummy water images
    waterDummyCol = ee.ImageCollection(waterStamps.map(dummyWaterImg))
  
    # Mask tile and water
    waterDummyCol = waterDummyCol.map(prepImg.waterMask)
    waterDummyCol = prepImg.tileMask(waterDummyCol, tile)
  
    # Return collection
    return waterDummyCol

# ---------------------------------------------------------------------------- #
# Add dummy water to collection of real ice-water observations
def addDummyData(imgCol, year, tile):
    '''
    This function adds dummy water (0) to imgCol to period preceding poiEnd 
    for spring breakup of 'year'. Dummy data is masked to tile extent.
    Real and dummy data flagged using 'dummy' property.
    '''
    # Get number of dummy water imgs to produce (approx 10% of true observations)
    nDummy = imgCol.size().divide(10).round()
    
    # Generate dummy data 
    dummyData = dummyCollection(year, tile, nDummy)
    
    # Add dummy property flag
    dummyData = dummyData.map(addPropDummyYes)
    imgCol = imgCol.map(addPropDummyNo)
  
    # Return merged
    return imgCol.merge(dummyData)