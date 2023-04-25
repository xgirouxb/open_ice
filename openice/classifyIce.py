# ---------------------------------------------------------------------------- #
# Load/initialize Earth Engine API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# LANDSAT 7 TOP OF ATMOSPHERE CLASSIFICATION

# ---------------------------------------------------------------------------- #
# Define classifier parameters for L7 TOA

# Bands required for classification
classBandsL7TOA = ['blue', 'swir2']

# Define Landsat 7 TOA classification tree
treeStringL7TOA = "\n".join([  
"1) root 927984 618656 0 (0.333333333 0.333333333 0.333333333)  ",
"  2) swir2< 0.08979684 615416 307150 0 (0.500906704 0.490882590 0.008210706)  ",
"    4) blue< 0.1603449 313364  11715 0 (0.962615361 0.036312403 0.001072235) *",
"    5) blue>=0.1603449 302052  11334 1 (0.021906824 0.962476660 0.015616516) *",
"  3) swir2>=0.08979684 312568   8293 9 (0.003397661 0.023134166 0.973468173) *",
])

# Define L7 TOA decision tree classifier
classL7TOA = ee.Classifier.decisionTree(treeStringL7TOA)

# ---------------------------------------------------------------------------- #
# Define classification function for L7 TOA
def classIceL7TOA(img):
    """
    This function classifies Landsat 7 TOA images into ice, water, and clouds.
    
    INPUT: Landsat 7 TOA image with 'blue' and 'swir2' standard band names
           required for ice-water-cloud classification
    OUTPUT: Returns image with band of ice presence (1) or absence (0),
            and no data where classifier detected clouds
    """
    
    # Classify into ice-water-clouds
    classified =  img.select(classBandsL7TOA)\
                     .classify(classL7TOA)
    # Create cloud mask with buffer
    cloud = classified.eq(9)
    cloudMask = cloud.fastDistanceTransform().gt(30)
    # Create band for ice, mask clouds detected by classifier
    ice = classified.eq(1).toUint16()\
                    .updateMask(cloudMask)\
                    .rename(['classIce']) 
    # Return ice, save properties
    return img.addBands(ice)\
              .select(['classIce'])

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# LANDSAT 8 TOP OF ATMOSPHERE REFLECTANCE CLASSIFICATION

# ---------------------------------------------------------------------------- #
# Define classifier parameters for L8 TOA

# Bands required for classification
classBandsL8TOA = ['blue', 'ndsi']

# Define Landsat 8 TOA classification tree
treeStringL8TOA = "\n".join([ 
"1) root 1174044 782696 0 (3.333333e-01 3.333333e-01 3.333333e-01)  ",
"  2) blue< 0.1432013 393101   1826 0 (9.953549e-01 4.550993e-03 9.412339e-05) *",
"  3) blue>=0.1432013 780943 389632 9 (9.347673e-05 4.988315e-01 5.010750e-01)  ",
"    6) ndsi>=0.8493385 388665    327 1 (5.917693e-05 9.991587e-01 7.821646e-04) *",
"    7) ndsi< 0.8493385 392278   1271 9 (1.274606e-04 3.112589e-03 9.967600e-01) *",
])

# Define L8 TOA decision tree classifier
classL8TOA = ee.Classifier.decisionTree(treeStringL8TOA)

# ---------------------------------------------------------------------------- #
# Define classification function for L8 TOA
def classIceL8TOA(img):
    """
    This function classifies Landsat 8 TOA images into ice, water, and clouds.
    
    INPUT: Landsat 8 TOA image with 'blue' and 'ndsi' standard band names
           required for ice-water-cloud classification
    OUTPUT: Returns image with band of ice presence (1) or absence (0),
            and no data where classifier detected clouds
    """
    
    # Classify into ice-water-clouds
    classified =  img.select(classBandsL8TOA)\
                     .classify(classL8TOA)
    # Create cloud mask with buffer
    cloud = classified.eq(9)
    cloudMask = cloud.fastDistanceTransform().gt(30)
    # Create band for ice, mask clouds detected by classifier
    ice = classified.eq(1).toUint16()\
                    .updateMask(cloudMask)\
                    .rename(['classIce']) 
    # Return ice, save properties
    return img.addBands(ice)\
              .select(['classIce'])

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# SENTINEL 2 TOP OF ATMOSPHERE REFLECTANCE

# ---------------------------------------------------------------------------- #
# Define classifier parameters for S2 TOA

# Bands required for classification
classBandsS2TOA = ['blue', 'ndsi']

# Define Sentinel 2 TOA classification tree
treeStringS2TOA = "\n".join([  
"1) root 902238 601492 0 (3.333333e-01 3.333333e-01 3.333333e-01)  ",
"  2) blue< 0.1414 307714   7105 0 (9.769104e-01 2.300513e-02 8.449404e-05) *",
"  3) blue>=0.1414 594524 293804 9 (2.304365e-04 4.939531e-01 5.058164e-01)  ",
"    6) ndsi>=0.8172837 299419   7080 1 (6.679603e-05 9.763542e-01 2.357900e-02) *",
"    7) ndsi< 0.8172837 295105   1445 9 (3.964691e-04 4.500093e-03 9.951034e-01) *",
])

# Define S2 TOA decision tree classifier
classS2TOA = ee.Classifier.decisionTree(treeStringS2TOA)

# ---------------------------------------------------------------------------- #
# Define classification function for S2 TOA
def classIceS2TOA(img):
    """
    This function classifies Sentinel 2 TOA images into ice, water, and clouds.
    
    INPUT: Sentinel 2 TOA image with 'blue' and 'ndsi' standard band names
           required for ice-water-cloud classification
    OUTPUT: Returns image with band of ice presence (1) or absence (0),
            and no data where classifier detected clouds
    """
    
    # Classify into ice-water-clouds
    classified = img.select(classBandsS2TOA)\
                    .classify(classS2TOA)
    # Create cloud mask with buffer
    cloud = classified.eq(9)
    cloudMask = cloud.fastDistanceTransform().gt(30)
    # Create band for ice, mask out clouds detected by classifiers
    ice = classified.eq(1).toUint16()\
                    .updateMask(cloudMask)\
                    .rename(['classIce']) 
    # Return ice, save properties
    return img.addBands(ice)\
              .select(['classIce'])