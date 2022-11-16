# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Collection of functions to fit a logistic regression of 
# ice presence(1)/absence(0) ~ time (fraction year) and filter out observations
# with high residuals (likely misclassifications)

# ---------------------------------------------------------------------------- #
# Initialize EE API
import ee
ee.Initialize()

# ---------------------------------------------------------------------------- #
# Function to fit logistic regression
def fitLogisticRegression(imgCol, y, x):
    """
    This function fits a linearized (y = β0 + β1 * x + Ɛ) version
    of a logistic regression:
    
    ln(y/(1-y)) = β0 + β1*x

    and solves for β0 and β1 using linearRegression().
    
    INPUTS 
    imgCol: Image Collection containing dependant and independant variables,
            to which function appends outputs
    y:      Dependant variable in imgCol, binary response band
            (ice presence (1) or absence (0))
    x:      Independant variables in imgCol, explanatory band 
            (fractional year for time series)
    OUTPUT 
    Bands added to images in imgCol: B0, B1, fitted, residuals, R2
    """
    
    # Function to prepare inputs for logistic regression
    def prepLogisticInputs(img):
        # transform to avoid infinity
        yt = img.select(y)\
                .remap([0,1], [0.001,0.999])
        # calculate logit of dependant variable ln(y/(1-y))
        logity = yt.divide(ee.Image(1.0).subtract(yt))\
                   .log()\
                   .rename(['logity'])
        # add constant of 1 to estimate β0
        constant = ee.Image.constant(1)
        # add bands to each image
        return img.addBands(logity)\
                  .addBands(constant, ['constant'], True)
    
    # Apply function to each img in imgCol
    logisticInputs = imgCol.map(prepLogisticInputs)
    
    # Fit logistic regression
    logisticRegression = logisticInputs.select(['constant', x, 'logity'])\
                                       .reduce(ee.Reducer.linearRegression(2, 1))
    
    # Extract coefficients estimated by regression
    coefficients = logisticRegression.select('coefficients')\
                                     .arrayProject([0])\
                                     .arrayFlatten([['constant', x]])\
                                     .rename(['B0', 'B1'])
    
    # Define function that adds β coefficients to imgCol
    def addBetaCoefficients(img):
        return img.addBands(coefficients.select('B0'))\
                  .addBands(coefficients.select('B1'))
    
    # Define function to calculate fitted values based on coefficients
    def addFitted(img):
        # exp(β0 + β1*x)
        top = img.select('B0')\
                 .add(img.select('B1').multiply(img.select([x])))\
                 .exp()
        # (1 + exp(β0 + β1*x))
        bottom = ee.Image(1.0).add(top)
        # ICEt = exp(β0 + β1*x) / (1 + exp(β0 + β1*x))
        fitted = top.divide(bottom).rename(['fitted'])
        # mask values fitted to masked pixels
        mask = img.select(y).lte(1)
        fitted = fitted.updateMask(mask)
        # return fitted values with same mask as y
        return img.addBands(fitted)
    
    # Define function to add residuals 
    def addResiduals(img):
        # calculate residuals (observed - fitted)
        residuals = img.select('fitted')\
                       .subtract(img.select(y))\
                       .rename(['residuals'])
        return img.addBands(residuals)
    
    # Fit model and add coefficents, fitted values, and residuals
    outputImgCol = imgCol.map(addBetaCoefficients)\
                         .map(addFitted)\
                         .map(addResiduals)
    
    # Define function to calculate R-squared
    def addRsquared(imgColLogReg):
        # calculate yflat
        yflat = imgColLogReg.select(y).mean()
        # Sum of Squares Regression (SSR) = ∑n_i (yhat_i−yflat)^2
        def yhat_yflat(img):
            return img.select('fitted')\
                      .subtract(yflat)\
                      .pow(2)
        SSR = imgColLogReg.map(yhat_yflat)\
                          .sum()
        # Total sum of squares (SSTO) = ∑n_i(y_i−yflat)^2
        def yi_yflat(img):
            return img.select(y)\
                      .subtract(yflat)\
                      .pow(2)
        SSTO = imgColLogReg.map(yi_yflat)\
                           .sum()
        # Calculate, add R2
        def addCalcR2(img):
            # Calculate R2
            R2 = SSR.divide(SSTO).rename(['R2'])
            return img.addBands(R2)
            
        return imgColLogReg.map(addCalcR2)
    
    # Add r-squared
    outputImgCol = addRsquared(outputImgCol)
    
    # Return collection with added B0, B1, fitted, residuals, R2
    return outputImgCol


# ---------------------------------------------------------------------------- #
# Function to filter out observations with high residuals
def temporalFilter(img):
    """
    This function masks out the observations with high residuals
    for each img in an imgCol.
    """
    mask = img.select('residuals').abs().lte(0.85)
    
    return img.updateMask(mask)\
              .copyProperties(img)