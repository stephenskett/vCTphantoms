# -*- coding: utf-8 -*-
"""
Class definitions for 'vCTphantom', and for 'vCTshape' and its subclasses:
  - 'vCTcuboid';
  - 'vCTsphere';
  - 'vCTcylinder.
N.B. Currently the shapes are defined in terms of pixel dimensions rather than
physical distance in the image. Since it is intended to use equivalent pixel
size and slice spacing, this is fine. However, if either the pixel size or
slice spacing are changed in the intended output PyDICOM dataset, the shapes
will not be represented correctly in the images.

Version 2.0 info:
- added shape-specifications defaults & valid yes-responses dictionaries to
  vCTphantom init method;
- added 'vCTstudy' & 'vCTseries' classes to hold data pertaining to a
  particular DICOM study or series involving vCTphantom objects (e.g. instance
  UIDs, times, dates and information pertaining to any structures in the
  phantom or deformations applied to it).
- added methods for performing slice-wise deformations on the pixel data of a
  given vCTphantom object (using the new 'image_deform' module).

Created by Stephen Skett as part of the KCL MSc Project, "Development of a
commissioning and QA framework for deformable image registration in
radiotherapy planning."

"""

from __future__ import print_function
import numpy as np
import textwrap as tw
import random
import dcmfunctions as df
import vCTfunctions as vCTf
from datetime import datetime as dt
import image_deform as imdef
#import pdb

##############################################################################

class vCTphantom(object):
    # Initialisation method: store input dimensions and create 3D zero matrix. 
    def __init__(self, rows, columns, numslices, seriesnum=1,
                 patientID="", patientname="", bits=16):
        # Input phantom basic geometry & ID attributes.
        self.rows = vCTf.nint(rows)
        self.columns = vCTf.nint(columns)
        self.slices = vCTf.nint(numslices)
        self.seriesnum = vCTf.nint(seriesnum)
        self.patientID = str(patientID)
        self.patientname = str(patientname)
        
        # Set up phantom reference data attributes.
        self.ref_data = vCTrefdata()
        
        # Determine required phantom dtype & initialise pixel data.
        self.bits = vCTf.nint(bits)
        if self.ref_data.dtypedict.has_key(self.bits):
            self.dtype = self.ref_data.dtypedict[self.bits]
        else:
            warntxt = ("Invalid phantom bit depth entered; 16-bit unsigned "+
                       "integer data type will be used.")
            print(tw.fill(warntxt))
            self.dtype = self.ref_data.dtypedict[16]
        self.reset_pixels()
        self.phantom_shapes = None
        self.phantom_deformations = None
        self.phantom_noise_blur = None
    
    def __str__(self):
        """
        String output method: return pixel data array.
        [N.B. don't need array output method, as self.pixel_data is an array.]
        
        """
        return str(self.pixel_data)
    
    def __additemspec(self, specdict, spectype):
        """
        Adds the specification dictionary for a given vCTphantom-associated
        item (i.e. a vCTshape, a phantom deformation or a noise/blurring
        characteristic) to the appropriate phantom dictionary attribute.
        
        """
        # Check that item specification is valid (i.e. specdict is a Python
        # dictionary & spectype is interpretable as one of the 3 item types).
        if not isinstance(specdict, dict):
            print("Invalid item specification: exiting...")
            return
        if not spectype in self.ref_data.itemstypelist:
            if spectype in self.ref_data.shapetypelist:
                spectype = 'shape'
            elif spectype in self.ref_data.deformtypelist:
                spectype = 'deformation'
            elif spectype in self.ref_data.nbtypelist:
                spectype = 'noiseblur'
            else:
                print("Invalid item specification: exiting...")
                return
        # Loop through the item types & update the relevant item dictionary.
        itemspecsdict = dict(zip(self.ref_data.itemstypelist,
                                 [self.phantom_shapes,
                                  self.phantom_deformations,
                                  self.phantom_noise_blur]))
        for k in itemspecsdict.keys():
            if spectype==k:
                if not isinstance(itemspecsdict[k], dict):
                    itemspecsdict[k] = {1: specdict}
                else:
                    itemnum = 2
                    itemspecadded = False
                    while not itemspecadded:
                        try:
                            if itemnum not in itemspecsdict[k]:
                                raise KeyError
                        except:
                            itemspecsdict[k][itemnum] = specdict
                            itemspecadded = True
                        else:
                            itemnum += 1
                break
        # Store the updated item dictionaries to the relevant vCTphantom
        # object attributes.
        self.phantom_shapes = itemspecsdict['shape']
        self.phantom_deformations = itemspecsdict['deformation']
        self.phantom_noise_blur = itemspecsdict['noiseblur']
    
    def _noise_blur(self, noiseorblur, nbparameter):
        """
        Semi-private method which performs noise or blurring on vCTphantom
        object's pixel data.
        
        The user should not typically use this method directly, but instead
        make calls to either of the 'apply_blurring' or 'add_random_noise'
        methods, both of which are wrappers for this function.
        
        """
        if not noiseorblur in self.ref_data.nbtypelist:
            print("Invalid noise/blurring specification: exiting...")
            return
        if not isinstance(nbparameter, int):
            if isinstance(nbparameter, float):
                nbparameter = vCTf.nint(nbparameter)
            else:
                if noiseorblur=='noise':
                    print("Specified amplitude must be numeric.")
                else:
                    print("Specified kernel sigma value must be numeric.")
                return
        if nbparameter == 0:
            return
        if noiseorblur=='noise':
            self.pixel_data=vCTf.random_additive_noise(self.pixel_data,
                                                       nbparameter,
                                                       phantombits=self.bits)
            nbspec = {'class':'noise', 'amplitude':nbparameter}
        else:
            self.pixel_data=vCTf.gaussian_blur_pixels(self.pixel_data,
                                                      nbparameter)
            nbspec = {'class':'blur', 'sigma':nbparameter}
        # Add noise/blur specification dict to the phantom_noise_blur dict.
        self.__additemspec(nbspec, 'noiseblur')
    
    def reset_pixels(self, newdtype=None):
        """
        Method for setting/resetting pixel data to default (all zeros).
        This method also allows the user to control the dtype of the
        initialized phantom pixel data array.
        
        N.B. if no new pixel data dtype is specified, or the specified value
        is invalid, the original dtype will be used by default.
        
        """
        if newdtype == None or not newdtype in self.ref_data.dtypedict:
            newdtype = self.dtype
        else:
            self.dtype = newdtype
        self.pixel_data = np.zeros((self.slices,
                                    self.rows,
                                    self.columns), dtype=self.dtype)
    
    def get_slice(self, slicenum):
        """
        Method for retrieving a given slice of pixel data.
        
        N.B. Python slicing convention has been used here, i.e. the slices in
        the pixel data are referenced by indices [0] to [self.slices - 1]. 
        
        """
        if not isinstance(slicenum, (int, float)):
            # If slicenum is not numeric, return an empty array.
            print("Specified slice number must be numeric.")
            return np.empty(0)
        if slicenum < self.slices - 0.5:
            # Select the appropriate slice within the range of the number of
            # slices present. If a non-integer value of slicenum is supplied,
            # use nearest neighbour interpolation.
            slicenum = vCTf.nint(slicenum)
            return self.pixel_data[slicenum,:,:]
        else:
            # If slicenum is out of range, return an empty array.
            errtext = ("Specified slice number is out of range. (N.B. the " +
                       "standard Python indexing convention is used here, " +
                       "i.e. slices indexed from [0] to [numslices-1].)")
            print(tw.fill(errtext))
            return np.empty(0)
    
    def write_slice(self, slicenum, slicepixeldata):
        """
        Method for overwriting the pixel data for a given image slice. Returns
        a Boolean value, indicating whether or not the slice was written
        successfully (True) or not (False).
        
        N.B. Python slicing convention has been used here, i.e. the slices in
        the pixel data are referenced by indices [0] to [self.slices - 1]. 
        
        """
        if not isinstance(slicenum, (int, float)):
            print("Specified slice number must be numeric.")
            return False
        elif not slicenum < self.slices - 0.5:
            errtext = ("Specified slice number is out of range. (N.B. the " +
                       "standard Python indexing convention is used here, " +
                       "i.e. slices indexed from [0] to [numslices-1].)")
            print(tw.fill(errtext))
            return False
        elif not isinstance(slicepixeldata, np.ndarray):
            print("Specified slice pixel data must be a numpy array.")
            return False
        elif slicepixeldata.dtype.kind in 'OSU':
            print("Specified slice pixel array must contain numeric data.")
            return False
        elif slicepixeldata.shape != (self.rows, self.columns):
            print("Specified slice pixel array has invalid dimensions.")
            return False
        else:
            # Select the appropriate slice within the range of the number of
            # slices present. If a non-integer value of slicenum is supplied,
            # use nearest neighbour interpolation.
            slicenum = vCTf.nint(slicenum)
            self.pixel_data[slicenum,:,:] = slicepixeldata.astype(self.dtype)
            return True
    
    def apply_blurring(self, kernel_sigma=1):
        """
        Performs Gaussian blurring, using a kernel with size specified by the
        Gaussian sigma value (in pixels), on the phantom pixel data.
        
        """
        self._noise_blur('blur', kernel_sigma)
    
    def add_random_noise(self, amplitude=100):
        """
        Introduces random additive noise to the phantom pixel data (i.e. adds
        a random positive integer to each pixel). The random numbers will be
        uniformly-distributed between zero and the maximum value (specified by
        the 'amplitude' keyword).
        
        The bit-depth of the pixel data is not changed by this process: where
        the value of pixel + noise exceeds the maximum allowed, the new pixel
        value will be truncated.
        
        """
        self._noise_blur('noise', amplitude)
    
    def add_shape(self, shapetype, shapevars, pixelval=1024):
        """
        Method for adding shapes to the 3-D pixel data.
        
        """
        shapetypelist = self.ref_data.shapetypelist
        shapevarslens = self.ref_data.shapetypevarslens
        shapevarsdict = dict(zip(shapetypelist,shapevarslens))
        try:
            # Check the first element of shapevars list is a valid (z, y, x)
            # pixel position (i.e. list/tuple with exactly 3 elements). Check
            # that the specified shape type is valid (i.e. in shapetypelist),
            # and that the right number of shape parameters were supplied.
            if len(shapevars[0]) != 3:
                raise ShapePosnError
            if shapetype not in shapetypelist:
                raise ShapeTypeError(shapetypelist)
            shapevarslen = len(shapevars)
            if shapevarslen not in shapevarsdict[shapetype]:
                raise ShapeVarsError
        except ShapePosnError:
            print('Invalid shape position coordinate.')
            return
        except ShapeTypeError as STerr:
            errtext = ('Invalid shape type specification. You specified a ' +
                       'shape of type \'' + shapetype + '\'. This must be ' +
                       'one of the following: ' + str(STerr) + '.')
            print(tw.fill(errtext))
            return
        except:
            print('Invalid shape definition variables.')
            return
        
        # Determine the appropriate shape class, then create the relevant
        # shape object and define the required shape parameters. Populate the
        # shape's specification dictionary.
        if shapetype in shapetypelist[0:4]:
            shapeclass = 'cuboid'
            if shapetype == 'square':
                ftlcorner = shapevars[0]
                width = depth = shapevars[1]
                length = 1
            elif shapetype == 'rectangle':
                ftlcorner = shapevars[0]
                width = shapevars[1]
                depth = shapevars[2]
                length = 1
            elif shapetype == 'cube':
                ftlcorner = shapevars[0]
                width = depth = length = shapevars[1]
            else:
                ftlcorner = shapevars[0]
                width = shapevars[1]
                depth = shapevars[2]
                length = shapevars[3]
            shapespec = {'class': 'cuboid', 'ftlcorner': ftlcorner,
                         'width': width, 'depth': depth, 'length': length,
                         'pixelval': pixelval}
        elif shapetype in shapetypelist[4:7]:
            shapeclass = 'cylinder'
            if shapetype == 'circle':
                centroid = shapevars[0]
                xradius = yradius = shapevars[1]
                length = 1
            elif shapetype == 'ellipse':
                centroid = shapevars[0]
                xradius = shapevars[1]
                yradius = shapevars[2]
                length = 1
            elif (shapetype == 'cylinder' and shapevarslen == 3):
                centroid = shapevars[0]
                xradius = yradius = shapevars[1]
                length = shapevars[2]
            elif (shapetype == 'cylinder' and shapevarslen == 4):
                centroid = shapevars[0]
                xradius = shapevars[1]
                yradius = shapevars[2]
                length = shapevars[3]
            shapespec = {'class': 'cylinder', 'centroid': centroid,
                         'xradius': xradius, 'yradius': yradius,
                         'length': length, 'pixelval': pixelval}
        else:
            shapeclass = 'sphere'
            centroid = shapevars[0]
            radius = shapevars[1]
            shapespec = {'class': 'sphere', 'centroid': centroid,
                         'radius': radius, 'pixelval': pixelval}
        
        # Define relevant shape and update pixel data in self.
        if shapeclass == 'cuboid':
           shapeobj=vCTcuboid(self,ftlcorner,width,depth,length,pixelval)
        elif shapeclass == 'cylinder':
           shapeobj=vCTcylinder(self,centroid,xradius,yradius,length,pixelval)
        elif shapeclass == 'sphere':
           shapeobj=vCTsphere(self,centroid,radius,pixelval)
        progbarnum = 0
        numpixels = len(shapeobj.pixels_list)
        print("Updating phantom pixel data...")
        for i in range(numpixels):
            coord = shapeobj.pixels_list[i]
            self.pixel_data[coord] = pixelval
            progbarnum = vCTf.increment_progbar(i, numpixels, progbarnum)
        print("")
        
        # Add the shape specification dictionary to the phantom_shapes dict.
        self.__additemspec(shapespec, 'shape')
    
    def deform_phantom(self, deformtype, deformvars):
        """
        Method for performing a given type of slicewise deformation on the
        phantom image pixel data.
        
        """
        deftypelist = self.ref_data.deformtypelist
        defvarslens = self.ref_data.deformtypevarslens
        defvarsdict = dict(zip(deftypelist,defvarslens))
        try:
            # As for adding shapes: first check that the specified
            # deformation type is valid (i.e. in deftypelist), and that
            # the right number of deformation parameters were supplied.
            if deformtype not in deftypelist:
                raise DefTypeError(deftypelist)
            defvarslen = len(deformvars)
            if defvarslen not in defvarsdict[deformtype]:
                raise DefVarsError
            # Check that supplied deformation variables are numeric.
            for i in range(defvarslen):
                if not isinstance(deformvars[i], (int, float)):
                    raise DefVarsError
        except DefTypeError as DTerr:
            errtext = ('Invalid deformation type specification. You ' +
                       'specified a deformation of type \'' + deformtype +
                       '\'. This must be one of the following: ' +
                       str(DTerr) + '.')
            print(tw.fill(errtext))
            return False
        except:
            print('Invalid deformation variables specified.')
            return False
        
        # For each deformation type, assign the relevant variables and define
        # the transformation matrix.
        if deformtype == 'stretch':
            k = deformvars[0]
            xoffset = deformvars[1]
            yoffset = deformvars[2]
            axisangle = deformvars[3]
            transformmatrix = imdef.stretch_matrix(k, [yoffset, xoffset],
                                                   axisangle)
            deformspec = {'class': 'stretch', 'kfactor': k,
                          'xoffset': xoffset, 'yoffset': yoffset,
                          'axisangle': axisangle}
        elif deformtype == 'shear':
            k = deformvars[0]
            axisangle = deformvars[1]
            transformmatrix = imdef.shear_matrix(k, axisangle)
            deformspec = {'class': 'shear', 'kfactor': k,
                          'axisangle': axisangle}
        else:
            # Add code here if/when other deformation classes are implemented.
            pass
        
        # Perform the deformation operation for each slice in the phantom.
        if not len(transformmatrix) == 1:
            print("Deforming slices...")
            progbarnum = 0
            for i in range(self.slices):
                # Print progress bar to command line.
                progbarnum = vCTf.increment_progbar(i, self.slices, progbarnum)
                # Get pixel data for image slice i.
                slicei = self.get_slice(i)
                # Efficiency improvement: if all slice pixels are same, slice
                # data will not change, so needn't bother doing deformation.
                if slicei.max == slicei.min:
                    continue
                # Do deformation...
                slicei, deformstatus = imdef.apply_deformation(
                                                 slicei, transformmatrix,
                                                 output_dtype = self.dtype)
                # If deformation worked, write new pixel data for slice i;
                # if not, stop deforming slices.
                if deformstatus:
                    self.write_slice(i, slicei)
                else:
                    break
            print("")
        else:
            print("Invalid transformation matrix.")
            return False
        
        # Add the deformation specification dictionary to the
        # phantom_deformations dict.
        if deformstatus:
            self.__additemspec(deformspec, 'deformation')
            return True
        else:
            print("There was a problem with the deformation.")
            return False

##############################################################################

class vCTshape(object):
    def __init__(self, phantomobj, shapeclass=None, pixelval=1024):
        """
        Initialisation method: create dummy shape-pixels list, and get data
        about the grid & bit depth from the phantom object.
        
        """
        self.pixels_list = []
        self.phantom_rows = phantomobj.rows
        self.phantom_columns = phantomobj.columns
        self.phantom_slices = phantomobj.slices
        self.phantom_bits = phantomobj.bits
        self.phantom_dtype = phantomobj.dtype
        if shapeclass not in phantomobj.ref_data.shapeclasslist:
            self.shape_class = None
        else:
            self.shape_class = shapeclass
        if not isinstance(pixelval, (int, float)):
            self.pixelval = None
        elif pixelval > 2**self.phantom_bits - 1:
            self.pixelval = 2**self.phantom_bits - 1
        else:
            self.pixelval = vCTf.nint(pixelval)
    
    def __repr__(self):
        """
        String representation method: return shape-pixels list.
        [N.B. don't need an array output function, as self.pixels_list is
        array-like.]
        
        """
        return str(self.pixels_list)
    
    def __str__(self):
        """
        'Pretty' string representation method: return string of shape-pixels 
        list, in a more readable format.
        
        """
        pixels_list = self.pixels_list
        n = len(pixels_list)
        if n > 8:
            outputstr = ('[' + str(pixels_list[0]) + ',\n' +
                         ' ' + str(pixels_list[1]) + ',\n' +
                         ' ' + str(pixels_list[2]) + ',\n' +
                         ' ...,\n' +
                         ' ' + str(pixels_list[n-3]) + ',\n' +
                         ' ' + str(pixels_list[n-2]) + ',\n' +
                         ' ' + str(pixels_list[n-1]) + ']')
        else:
            outputstr = '['
            for i in range(n):
                if i != 0:
                    outputstr = outputstr + ' '
                outputstr = outputstr + str(pixels_list[i])
                if i != n-1:
                    outputstr = outputstr + ',\n'
                else:
                    outputstr = outputstr + ']'
        return outputstr

##############################################################################

class vCTcuboid(vCTshape):
    def __init__(self,phantomobj,ftlcorner,width,depth,length,pixelval=1024):
        """
        Initialisation method for vCTcuboid objects; perform the standard
        vCTshape initialisation, then add cuboid-specific attributes and
        create the pixel-coordinates-and-values list for the desired cuboid
        using the define_cuboid method.
        
        """
        super(vCTcuboid, self).__init__(phantomobj, shapeclass='cuboid',
                                        pixelval=pixelval)
        self.ftlcorner = (vCTf.nint(ftlcorner[0]),
                          vCTf.nint(ftlcorner[1]),
                          vCTf.nint(ftlcorner[2]))
        self.width = vCTf.nint(width)
        self.depth = vCTf.nint(depth)
        self.length = vCTf.nint(length)
        self.define_cuboid()
        
    def define_cuboid(self,ftlcorner=None,width=None,depth=None,length=None):
        """
        Create pixels list for simple square, rectangle, cube or cuboid
        phantoms from user-specified dimensions.
     
        """
        self.pixels_list = []
        if ftlcorner != None:
            try:
                self.ftlcorner = (vCTf.nint(ftlcorner[0]),
                                  vCTf.nint(ftlcorner[1]),
                                  vCTf.nint(ftlcorner[2]))
            except:
                print('Invalid shape position coordinate.')
                return
        if width != None: self.width = vCTf.nint(width)
        if depth != None: self.depth = vCTf.nint(depth)
        if length != None: self.length = vCTf.nint(length)
        print("Generating list of shape pixels...")
        progbarnum = 0
        for c in range(self.length):
            for a in range(self.width):
                for b in range(self.depth):
                    x = a + self.ftlcorner[2]
                    y = b + self.ftlcorner[1]
                    z = c + self.ftlcorner[0]
                    if (x >= 0 and x < self.phantom_columns and
                        y >= 0 and y < self.phantom_rows and
                        z >= 0 and z < self.phantom_slices):
                            self.pixels_list.append((z,y,x))
            progbarnum = vCTf.increment_progbar(c, self.length, progbarnum)
        print("")

##############################################################################

class vCTcylinder(vCTshape):
    def __init__(self,phantomobj,centroid,xradius,yradius,length,
                 pixelval=1024):
        """
        Initialisation method for vCTcylinder objects; perform the standard
        vCTshape initialisation, then add cylinder-specific attributes and
        create the pixel-coordinates-and-values list for the desired cylinder
        using the define_cylinder method.
        
        """
        super(vCTcylinder, self).__init__(phantomobj, shapeclass='cylinder',
                                          pixelval=pixelval)
        self.centroid = (vCTf.nint(centroid[0]),
                         vCTf.nint(centroid[1]),
                         vCTf.nint(centroid[2]))
        self.rx = vCTf.nint(xradius)
        self.ry = vCTf.nint(yradius)
        self.length = vCTf.nint(length)
        self.define_cylinder()
    
    def define_cylinder(self,centroid=None,xradius=None,yradius=None,
                        length=None):
        """
        Create pixels list for circle, ellipse or cylinder phantoms from
        user-specified dimensions.
        
        """
        self.pixels_list = []
        if centroid != None:
            try:
                self.centroid = (vCTf.nint(centroid[0]),
                                 vCTf.nint(centroid[1]),
                                 vCTf.nint(centroid[2]))
            except:
                print('Invalid shape position coordinate.')
                return
        if xradius != None: self.rx = vCTf.nint(xradius)
        if yradius != None: self.ry = vCTf.nint(yradius)
        if length != None:
            halflength = vCTf.nint(length/2)
            self.length = length
        else:
            halflength = vCTf.nint(self.length/2)
        print("Generating list of shape pixels...")
        progbarnum = 0
        for c in range(self.length):
            for a in range(2*(self.rx+1)):
                for b in range(2*(self.ry+1)):
                    xrel = a - self.rx
                    yrel = b - self.ry
                    zrel = c - halflength
                    x = self.centroid[2] + xrel
                    y = self.centroid[1] + yrel
                    z = self.centroid[0] + zrel
                    if (x >= 0 and x < self.phantom_columns and
                        y >= 0 and y < self.phantom_rows and
                        z >= 0 and z < self.phantom_slices and
                        ((float(xrel)/float(self.rx))**2 +
                         (float(yrel)/float(self.ry))**2 < 1.0)):
                            self.pixels_list.append((z,y,x))
            progbarnum = vCTf.increment_progbar(c, self.length, progbarnum)
        print("")

##############################################################################

class vCTsphere(vCTshape):
    def __init__(self,phantomobj,centroid,radius,pixelval=1024):
        """
        Initialisation method for vCTsphere objects; perform the standard
        vCTshape initialisation, then add sphere-specific attributes and
        create the pixel-coordinates-and-values list for the desired sphere
        using the define_sphere method.
        
        """
        super(vCTsphere, self).__init__(phantomobj, shapeclass='sphere',
                                        pixelval=pixelval)
        self.centroid = (vCTf.nint(centroid[0]),
                         vCTf.nint(centroid[1]),
                         vCTf.nint(centroid[2]))
        self.r = int(round(radius))
        self.define_sphere()
        
    def define_sphere(self,centroid=None,radius=None):
        """
        Create pixels list for simple sphere structure from user-specified
        dimensions.
        
        """
        self.pixels_list = []
        if centroid != None:
            try:
                self.centroid = (vCTf.nint(centroid[0]),
                                 vCTf.nint(centroid[1]),
                                 vCTf.nint(centroid[2]))
            except:
                print('Invalid shape position coordinate.')
                return
        if radius != None: self.r = vCTf.nint(radius)
        print("Generating list of shape pixels...")
        progbarnum = 0
        for c in range(2*(self.r+1)):
            for a in range(2*(self.r+1)):
                for b in range(2*(self.r+1)):
                    xrel = a - self.r
                    yrel = b - self.r
                    zrel = c - self.r
                    x = self.centroid[2] + xrel
                    y = self.centroid[1] + yrel
                    z = self.centroid[0] + zrel
                    if (x >= 0 and x < self.phantom_columns and
                        y >= 0 and y < self.phantom_rows and
                        z >= 0 and z < self.phantom_slices and
                        xrel**2 + yrel**2 + zrel**2 < self.r**2):
                            self.pixels_list.append((z,y,x))
            progbarnum = vCTf.increment_progbar(c, 2*(self.r+1), progbarnum)
        print("")

##############################################################################

class vCTseries(object):
    def __init__(self, phantomobj, studyrandoms=None, seriesdate=None,
                 seriestime=None, seriesFoRuid=None):
        """
        Initialisation method for vCTseries objects:
        - store information about the structures in the supplied vCTphantom
          object (i.e. the list of phantom shapes, including data about
          their shape class, geometric parameters & pixel value) - it is
          possible this may change for different series within a vCT study;
        - if other series-specific information is provided, store it;
          otherwise, generate new study-specific UID randoms, Series Instance
          UID, Series Date & Series Time.
        
        N.B. this does not check whether the supplied keyword values are valid
        (i.e. 'studyrandoms' should be a string of characters containing only
        digits & dots, and starting & ending in a digit; dates & times
        should be numberic strings of the form YYYYMMDD & HHmmss.ffffff
        respectively). Whilst the make_new_uid function can accommodate
        invalid randoms-strings to a certain extent, but users should take
        care when specifying study-specific UID randoms.
        
        """
        try:
            phantomshapes = phantomobj.phantom_shapes
            phantomdefs = phantomobj.phantom_deformations
            phantomnoiseblur = phantomobj.phantom_noise_blur
            seriesnum = phantomobj.seriesnum
        except:
            # Couldn't read phantom object shape and/or deformation dicts
            #  => phantom object invalid.
            # Return a 'blank' series object, i.e. all attributes set to None.
            self.seriesnum = self.seriesuid = self.seriesFoRuid = None
            self.studyrandoms = self.seriesdate = self.seriestime = None
            self.phantom_shapes = self.phantom_deformations = None
            self.phantom_noise_blur = None
        else:
            # Populate series object attributes.
            if studyrandoms == None:
                studyrandoms = (str(random.randrange(1, 999)) + '.' +
                                str(random.randrange(1, 999)))
            if seriesFoRuid == None: seriesFoRuid = df.make_new_uid()
            if seriesdate == None: seriesdate = dt.now().strftime("%Y%m%d")
            if seriestime == None: seriestime = dt.now().strftime("%H%M%S.%f")
            self.seriesnum = seriesnum
            self.studyrandoms = studyrandoms
            self.seriesuid = df.make_new_uid(specifiedrandoms = studyrandoms)
            self.seriesFoRuid = seriesFoRuid
            self.seriesdate = seriesdate
            self.seriestime = seriestime
            self.phantom_shapes = phantomshapes
            self.phantom_deformations = phantomdefs
            self.phantom_noise_blur = phantomnoiseblur
    
    def __str__(self):
        """
        String output method: return Study Instance UID as string.
        
        """
        return str(self.seriesuid)

##############################################################################

class vCTstudy(object):
    def __init__(self, studyrandoms=None, studydate=None, studytime=None):
        """
        Initialisation method for vCTstudy objects: if specific study
        information is provided, store it; otherwise, generate new study-
        specific UID randoms, Study Instance UID, Study Date & Study Time.
        
        N.B. this does not check whether the supplied keyword values are valid
        (i.e. 'studyrandoms' should be a string of characters containing only
        digits & dots, and starting & ending in a digit; dates & times
        should be numberic strings of the form YYYYMMDD & HHmmss.ffffff
        respectively). Whilst the make_new_uid function can accommodate
        invalid randoms-strings to a certain extent, but users should take
        care when specifying study-specific UID randoms.
        
        """
        if studyrandoms == None:
            studyrandoms = (str(random.randrange(1, 999)) + '.' +
                            str(random.randrange(1, 999)))
        if studydate == None: studydate = dt.now().strftime("%Y%m%d")
        if studytime == None: studytime = dt.now().strftime("%H%M%S.%f")
        
        self.studyrandoms = studyrandoms
        self.studyuid = df.make_new_uid(specifiedrandoms = studyrandoms,
                                        addtimestamp = False)
        self.studydate = studydate
        self.studytime = studytime
        self.numseries = 0
        self.seriesobjects = []
    
    def add_series(self, seriesobj):
        """
        Add a new series object to the list of series included in the study.
        
        """
        try:
            if seriesobj.studyrandoms != self.studyrandoms:
                raise TypeError
            for obj in self.seriesobjects:
                if (seriesobj.seriesuid == obj.seriesuid or
                    seriesobj.seriesnum == obj.seriesuid):
                        raise SeriesIDError
        except(TypeError):
            errtext = ("Series not added: the series object provided " +
                       "appears to be from a different (i.e. the study-" +
                       "specific UID randoms do not match this study).")
            print(tw.fill(errtext))
        except(SeriesIDError):
            errtext = ("Series not added: there is a clash with an " +
                       "existing series in the study (i.e. either the " +
                       "series UID or series number match).")
            print(tw.fill(errtext))
        except:
            print("Series not added: the specified series object is invalid.")
        else:
            self.numseries += 1
            self.seriesobjects.append(seriesobj)
    
    def __str__(self):
        """
        String output method: return Study Instance UID as string.
        
        """
        return str(self.studyuid)
    
    def __array__(self):
        """
        Array output method: return nested list with 3 elements, containing:
        1) the Study Instance UID;
        2) the number of series in the study;
        3) a list containing each of the Series Instance UIDs.
        
        """
        output = [[self.studyuid]]
        output.append([self.numseries])
        output.append([])
        for i in range(self.numseries):
            output[2].append(self.seriesuids[i])
        return output

##############################################################################

class vCTrefdata(object):
    def __init__(self):
        """
        Compile all the necessary reference data needed for the vCTphantom
        creation process.
        
        """
        self.shapeclasslist = ['cuboid','cylinder','sphere']
        self.shapetypelist = ['square','rectangle','cube','cuboid',
                              'circle','ellipse','cylinder','sphere']
        self.shapetypevarslens = [[2],[3],[2],[4],[2],[3],[3, 4],[2]]
        self.shapetypespecs = [
            ['Front-Top-Left Corner','Side Length','Pixel Value'],
            ['Front-Top-Left Corner','Width','Depth','Pixel Value'],
            ['Front-Top-Left Corner','Side Length','Pixel Value'],
            ['Front-Top-Left Corner','Width','Depth','Length','Pixel Value'],
            ['Centroid Position','Radius','Pixel Value'],
            ['Centroid Position','Horizontal Radius','Vertical Radius',
             'Pixel Value'],
            ['Centroid Position','Horizontal Radius',
             'Vertical Radius (Optional)','Length','Pixel Value'],
            ['Centroid Position','Radius','Pixel Value']]
        self.shapespecinputdefaults = {'Front-Top-Left Corner': '',
                                       'Centroid Position':     '',
                                       'Side Length':           '',
                                       'Width':                 '',
                                       'Depth':                 '',
                                       'Radius':                '',
                                       'Horizontal Radius':     '',
                                       'Vertical Radius':       '',
                                       'Length':                '',
                                       'Pixel Value':           '1024'}
        self.yesdict = {'y', 'yes', 'yeh', 'yeah', 'yep', 'yup', 'yus', 'ok',
                        'sure', 'aye', 'oui', 'ja', 'alrighty', 'make it so'}
        self.alphabetdict = {1:'a', 2:'b', 3:'c', 4:'d', 5:'e', 6:'f',
                             7:'g', 8:'h', 9:'i', 10:'j', 11:'k', 12:'l',
                             13:'m', 14:'n', 15:'o', 16:'p', 17:'q', 18:'r',
                             19:'s', 20:'t', 21:'u', 22:'v', 23:'w', 24:'x',
                             25:'y', 26:'z'}
        self.deformtypelist = ['stretch', 'shear']
        self.deformtypevarslens = [[4],[2]]
        self.nbtypelist = ['noise','blur','blurring']
        self.dtypedict = {8:'uint8', 16:'uint16', 32:'uint32'}
        self.itemstypelist = ['shape','deformation','noiseblur']

##############################################################################

class ShapeTypeError(Exception):
    def __init__(self, shapetypelist):
        self.typelist = ', '.join(shapetypelist)
    def __str__(self):
        return self.typelist

class ShapePosnError(Exception):
    def __init__(self):
        pass

class ShapeVarsError(Exception):
    def __init__(self):
        pass

class DefTypeError(Exception):
    def __init__(self, deftypelist):
        self.typelist = ', '.join(deftypelist)
    def __str__(self):
        return self.typelist

class DefVarsError(Exception):
    def __init__(self):
        pass

class SeriesIDError(Exception):
    def __init__(self):
        pass

##############################################################################