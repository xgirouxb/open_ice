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
"1) root 826599 551066 0 (0.333333333 0.333333333 0.333333333)  ",
"  2) swir2< 0.08972447 549042 274065 0 (0.500830538 0.491060793 0.008108669)  ",
"    4) blue< 0.164955 279825   7972 0 (0.971510766 0.027481462 0.001007773) *",
"    5) blue>=0.164955 269217   7294 1 (0.011604022 0.972906614 0.015489364) *",
"  3) swir2>=0.08972447 277557   6476 9 (0.002003192 0.021328952 0.976667856) *",
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
"1) root 1173627 782418 0 (3.333333e-01 3.333333e-01 3.333333e-01)  ",
"  2) blue< 0.1435298 393042   1892 0 (9.951863e-01 4.750128e-03 6.360643e-05) *",
"  3) blue>=0.1435298 780585 389401 9 (7.558434e-05 4.987823e-01 5.011421e-01)  ",
"    6) ndsi>=0.848531 388513    342 1 (4.633050e-05 9.991197e-01 8.339489e-04) *",
"    7) ndsi< 0.848531 392072   1212 9 (1.045726e-04 2.986696e-03 9.969087e-01) *",
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
"1) root 901695 601130 0 (3.333333e-01 3.333333e-01 3.333333e-01)  ",
"  2) blue< 0.1408 307367   6950 0 (9.773886e-01 2.253332e-02 7.808255e-05) *",
"  3) blue>=0.1408 594328 293787 9 (2.490207e-04 4.940689e-01 5.056820e-01)  ",
"    6) ndsi>=0.8190071 298871   6850 1 (6.022665e-05 9.770804e-01 2.285936e-02) *",
"    7) ndsi< 0.8190071 295457   1748 9 (4.399963e-04 5.476262e-03 9.940837e-01) *",
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