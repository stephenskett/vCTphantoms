# -*- coding: utf-8 -*-
"""
Functions for applying some simple classes of defomation to input 2-d numpy
arrays. Designed primarily for use in conjunction with vCTphantom objects.

Deformation matrices currently implemented are:
- linear stretching along an arbitrary axis (stretch_matrix);
- shearing by a defined shear factor (shear_matrix).

Created by Stephen Skett (unless otherwise specified) as part of the KCL MSc
Project, "Development of a commissioning and QA framework for deformable image
registration in radiotherapy planning."

"""

import numpy as np
import textwrap as tw
from scipy import ndimage as ndimg
import math
#import pdb

##############################################################################

def apply_deformation(inputarray, transformmatrix, **deformkws):
    """
    Perform arbitrary linear transformation to a 2d greyscale image, as
    specified by the supplied transform matrix.
    
    Return value is a dictionary containing the transformed output array
    (key='output_data') and a Boolean value (key='success_status') indicating
    whether the transform was successful (True) or unsuccessful (False).
    
    """
    transformsuccess = False
    try:
        # Test the transform matrix and image/data grid dimensions.
        matrixdims = transformmatrix.shape
        if len(matrixdims) != 2 or (matrixdims[0] != 3 or
                                    matrixdims[1] != 3):
            raise MatrixDimsError
        griddims = inputarray.shape
        if len(griddims) != 2:
            raise GridDimsError
        # Test whether the transform matrix is singular (i.e. whether an
        # inverse is defined):
        # - If matrix is not singular, can apply the inverse matrix to the
        #   specified output coordinates and interpolate over the input space
        #   to find the pixel/grid-data values.
        # - If matrix is singular, cannot calculate inverse; transform cannot
        #   be applied.
        matrixinverse = transformmatrix.I
    except(MatrixDimsError):
        errtext = ("The matrix is not a 3x3 square (this is required for " +
                   "2-d transformations in homogeneous coordinates); " +
                   "cannot proceed. [Hint: if no translation is required, " +
                   "set the bottom row & far-right column to (0 0 1).]")
        print tw.fill(errtext)
        outputarray = inputarray
    except(GridDimsError):
        errtext = ("The dimensions of the grid or axis data provided are " +
                   "invalid: the image- or data-grid must be 2-dimensional.")
        print tw.fill(errtext)
        outputarray = inputarray
    except(np.linalg.LinAlgError):
        # Singular transform matrix - inverse is undefined.
        print "The transformation matrix has no inverse; cannot proceed."
        outputarray = inputarray
    except:
        print "Transformation matrix invalid."
        outputarray = inputarray
    else:
        inputxaxis = range(griddims[1])
        inputyaxis = range(griddims[0])
        # Establish which deformations keywords have been entered; if
        # supplied keyword values are invalid, use defaults.
        try:
            outputdims = deformkws['output_array_dims']
            if len(outputdims) != 2:
                errtext = ("Invalid output array dimensions specified " +
                           "(must be 2-dimensional); will use the input " +
                           "array dimensions instead.")
                print tw.fill(errtext)
                raise GridDimsError
        except:
            outputdims = griddims
        try:
            # Get output axis scaling factors and define output axis ranges.
            outputscalefactors = deformkws['output_axis_scaling']
            if len(outputscalefactors) != 2:
                errtext = ("Invalid output axis scale factors specified " +
                           "(must be 2-dimensional); will use the unscaled" +
                           "input axes.")
                print tw.fill(errtext)
                raise GridDimsError
        except:
            outputxaxis = range(outputdims[1])
            outputyaxis = range(outputdims[0])
        else:
            outputxaxis = [outputscalefactors[1]*i
                               for i in range(outputdims[1])]
            outputyaxis = [outputscalefactors[0]*j
                               for j in range(outputdims[0])]
        try:
            outputdtype = deformkws['output_dtype']
            if isinstance(outputdtype, basestring):
                if not 'uint' in outputdtype:
                    print "Supplied dtype must be an unsigned integer format."
                    raise OutputTypeError
                outdtbits = outputdtype.split('uint')[1]
                if not outdtbits in ['8','16','32','64']:
                    print "Supplied dtype must have a valid number of bits."
                    raise OutputTypeError
            else:
                if not 'uint' in outputdtype.str:
                    print "Supplied dtype must be an unsigned integer format."
                    raise OutputTypeError
        except:
            outputdtype = 'uint64'
        # Initialise arrays of x & y coordinates required for interpolation.
        outputinposxcoords = np.zeros(outputdims, dtype='float')
        outputinposycoords = np.zeros(outputdims, dtype='float')
        # Apply inverse matrix to each output pixel coordinate, to get arrays
        # of output pixel coordinates in the input array frame of reference.
        inposvector = outposvector = np.mat(np.ones((3,1)))
        for i in range(outputdims[1]):
            for j in range(outputdims[0]):
                outposvector[0] = outputxaxis[i]
                outposvector[1] = outputyaxis[j]
                inposvector = matrixinverse*outposvector
                if inposvector[0] < inputxaxis[0]:
                    inposvector[0] = inputxaxis[0]
                if inposvector[1] < inputyaxis[0]:
                    inposvector[1] = inputyaxis[0]
                outputinposxcoords[j, i] = inposvector[0]
                outputinposycoords[j, i] = inposvector[1]
        # Calculate array of required output pixel values by interpolation.
        outputarray = ndimg.map_coordinates(inputarray,
                                            [outputinposycoords,
                                             outputinposxcoords], order=1)
        # Report that transform was applied successfully.
        transformsuccess = True
    finally:
        return outputarray, transformsuccess

##############################################################################

def stretch_matrix(k, offset, axisangle, angle_type = 'degrees'):
    """
    Generates a transformation matrix for linear stretching by a factor of k
    along an arbitrary axis angle (relative to the horizontal).
    
    """
    # Get stretching axis angle in radians.
    if angle_type == {"radians", "rads", "rad", "r"}:
        angle_rads = axisangle
    else:
        angle_rads = axisangle*math.pi/180
    
    # Check that stretching parameter are in an appropriate format
    if (not isinstance(k, (int, float)) or 
        not isinstance(offset, (list, tuple)) or
        not isinstance(axisangle, (int, float))):
            errtext = ("The stretch factor and stretching axis angle must " +
                       "both be numeric, and the positional offset must be " +
                       "a (2-element) list or tuple.")
            print tw.fill(errtext)
            stretchmat = np.mat(0)
    elif (not isinstance(offset[0], (int, float)) or
          not isinstance(offset[1], (int, float))):
            errtext = ("The positional offset elements in x & y must be " +
                       "numeric.")
            print tw.fill(errtext)
            stretchmat = np.mat(0)
    else:
        # Set up stretching transformation matrix
        stretchmat = np.mat(np.zeros((3, 3)))
        stretchmat[0,0] = (k-1)*(math.cos(angle_rads))**2 + 1
        stretchmat[0,1] = stretchmat[1,0] = (k-1)*(math.sin(angle_rads)*
                                                   math.cos(angle_rads))
        stretchmat[0,2] = (1-k)*(offset[1]*(math.cos(angle_rads))**2 +
                                 offset[0]*math.sin(angle_rads)*
                                           math.cos(angle_rads))
        stretchmat[1,1] = (k-1)*(math.sin(angle_rads))**2 + 1
        stretchmat[1,2] = (1-k)*(offset[1]*math.sin(angle_rads)*
                                           math.cos(angle_rads) +
                                 offset[0]*(math.sin(angle_rads))**2)
        stretchmat[2,2] = 1
    
    return stretchmat

##############################################################################

def shear_matrix(k, axisangle, angle_type = 'degrees'):
    """
    Performs linear shearing by a factor of k along an arbitrary axis angle
    (relative to the horizontal).
    
    Return value is a dictionary containing the sheared output array
    (key='output_data') and a Boolean value (key='success_status') indicating
    whether the transform was successful (True) or unsuccessful (False).
    
    """
    # Get stretching axis angle in radians.
    if angle_type == {"radians", "rads", "rad", "r"}:
        angle_rads = axisangle
    else:
        angle_rads = axisangle*math.pi/180
    
    # Check that stretching parameter are in an appropriate format
    if (not isinstance(k, (int, float)) or 
        not isinstance(axisangle, (int, float))):
            errtext = ("The shear factor and shearing axis angle must both " +
                       "be numeric.")
            print tw.fill(errtext)
            shearmat = np.mat(0)
    else:
        # Set up stretching transformation matrix
        shearmat = np.mat(np.zeros((3, 3)))
        shearmat[0,0] = 1 + k*math.sin(angle_rads)*math.cos(angle_rads)
        shearmat[0,1] = k*(math.cos(angle_rads))**2
        shearmat[1,0] = -k*(math.sin(angle_rads))**2
        shearmat[1,1] = 1 - k*math.sin(angle_rads)*math.cos(angle_rads)
        shearmat[2,2] = 1
    
    return shearmat

##############################################################################

# Other transformations go here...

##############################################################################

class MatrixDimsError(Exception):
    def __init__(self):
        pass

class GridDimsError(Exception):
    def __init__(self):
        pass

class OutputDimsError(Exception):
    def __init__(self):
        pass

class OutputTypeError(Exception):
    def __init__(self):
        pass

##############################################################################