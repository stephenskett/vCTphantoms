# -*- coding: utf-8 -*-
"""
Set of very simple functions, to be used by vCTphantom objects.

Version 2.1 info:
- relocated a number of functions from this module into a new module
  (vCT_interactive_creation) to ensure clarity of module scope.

Created by Stephen Skett as part of the KCL MSc Project, "Development of a
commissioning and QA framework for deformable image registration in
radiotherapy planning."

"""
from __future__ import print_function
import numpy as np
from scipy import ndimage as ndimg

##############################################################################

def nint(num):
    """
    Rounds the specified numeric scalar value (num) to the nearest integer,
    and returns the result in integer format.
    
    """
    if isinstance(num, (int, float)):
        return int(round(num, 0))
    else:
        return

##############################################################################

def increment_progbar(i, n, progbarnum):
    if not (isinstance(i, int) and isinstance(n, int) and
            isinstance(progbarnum, int)):
        return
    else:
        while True:
            if 70.0*i/n >= progbarnum:
                print("|", end="")
                progbarnum += 1
            else:
                break
        return progbarnum

##############################################################################

def random_additive_noise(phantompixels, amplitude, phantombits=16):
    """
    Adds random noise to an image slice (i.e. adds a small random number to
    each pixel). The 'amplitude' keyword specifies the maximum amplitude of
    these uniformly-distributed additive noise values.
    
    """
    if not (isinstance(phantompixels, np.ndarray) or
            isinstance(amplitude, int) or isinstance(phantombits, int)):
        print("Invalid noise specification...")
        return
    imdtype = phantompixels.dtype
    imshape = phantompixels.shape
    phantompixels = phantompixels + np.random.randint(0, amplitude, imshape)
    np.clip(phantompixels, 0, 2**phantombits-1, out=phantompixels)
    return phantompixels.astype(imdtype)

##############################################################################

def gaussian_blur_pixels(phantompixels, kernel_sigma):
    if not (isinstance(phantompixels, np.ndarray) or
            isinstance(kernel_sigma, int)):
        print("Invalid Gaussian blurring specification...")
        return
    return ndimg.gaussian_filter(phantompixels, kernel_sigma)

##############################################################################