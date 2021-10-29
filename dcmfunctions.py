# -*- coding: utf-8 -*-
"""
Some general functions for manipulating DICOM files in Python.

Created by Stephen Skett (unless otherwise specified) as part of the KCL MSc
Project, "Development of a commissioning and QA framework for deformable image
registration in radiotherapy planning."

"""

from __future__ import print_function
import dicom
import random
import numpy as np
import textwrap as tw
from fnmatch import fnmatch
from datetime import datetime as dt
from dicom.tag import Tag
from dicom.UID import UID
from dicom.datadict import dictionaryVR as VR
from os import listdir as ldir
import os

##############################################################################

def setup_DICOM_defaults(ds,imtype='vCT'):
    """
    Initialise imported DICOM file to project defaults.
    
    """
    try:
        fulltagslist = list(ds.keys())
        metatagslist = list(ds.file_meta.keys())
    except:
        print('The dataset provided must be a valid PyDICOM object.')
    
    # Set up initialisation variables.
    inittagslist = []
    initvalslist = []
    initVRslist = []
    initdate = dt.now().strftime("%Y%m%d")
    inittime = dt.now().strftime("%H%M%S.%f")
    instanceUID = make_new_uid(UIDformat='evalstr')
    initpixmat = np.zeros((1,1),dtype='uint16')
    
    # Read in initialisation file (location depending on computer in use).
    computername = os.environ['COMPUTERNAME']
    if "PCID" in computername:
        #PCH file location:
        inittagspath = r'J:/USERS/sketts/KCL MSc Project/Python Code/'
    elif computername=="STEPHEN-PC":
        #Home-PC file location:
        inittagspath = (r'E:/My Documents/Medical Physics/STP/KCL MSc ' +
                        r'Project/Python Code/')
    try:
        inittagsfname = inittagspath + "InitPyDICOM" + imtype + ".txt"
        inittagsfile = open(inittagsfname)
    except:
        print('Invalid initialisation image specification.')
    else:
        # Read tags from initialisation file.
        for line in inittagsfile:
            line = line.strip()
            if line[0] == '(':
                linetag = line[0:16]
                linetag = Tag(eval(linetag))
                lineVR = VR(linetag)
                lineval = line[19:]
                inittagslist.append(linetag)
                initvalslist.append(lineval)
                initVRslist.append(lineVR)
        
        # Delete any unwanted tags in PyDICOM header.
        for elemtag in metatagslist:
            if elemtag not in inittagslist:
                del ds.file_meta[elemtag]
        for elemtag in fulltagslist:
            if elemtag not in inittagslist:
                del ds[elemtag]
        
        # Change PyDICOM element values/add elements as required.
        for elemtag, elemval, elemVR in zip(inittagslist,
                                            initvalslist,
                                            initVRslist):
            if str(elemtag) == '(7fe0, 0010)':
                ds = replace_pixel_data(ds, initpixmat)
            else:
                if elemval == '#':
                    if (str(elemtag) == '(0002, 0003)' or
                        str(elemtag) == '(0008, 0018)'):
                        elemval = instanceUID
                    elif elemVR == 'DA':
                        elemval = initdate
                    elif elemVR == 'TM':
                        elemval = inittime
                # Set DICOM file meta-information...
                if str(elemtag).split(',')[0] == '(0002':
                    ds.file_meta.add_new(elemtag,
                                         elemVR,
                                         element_value_format(elemval,elemVR))
                # ...then set DICOM header information
                else:
                    # Skip the rows & columns entries -- these will be edited
                    # along with the pixel data itself.
                    if (str(elemtag) != '(0028, 0010)' and
                        str(elemtag) != '(0028, 0011)' and
                        str(elemtag) != '(fffc, fffc)'):
                            ds.add_new(elemtag,
                                       elemVR,
                                       element_value_format(elemval, elemVR))
        
        # Check Transfer Syntax UID and adjust Boolean test attributes.
        newTS = str(ds.file_meta.TransferSyntaxUID)
        allowedTStypes = {'Implicit VR Little Endian',
                          'Explicit VR Little Endian',
                          'Deflated Explicit VR Little Endian'}
        if newTS in allowedTStypes:
            ds.is_little_endian = True
            if newTS.split(' VR ')[0] == 'Implicit':
                ds.is_implicit_VR = True
            else:
                ds.is_implicit_VR = False
        else:
            errtext = ('Specified transfer syntax is invalid or not ' +
                       'supported; DICOM default (Implicit VR Little' +
                       ' Endian) used.')
            print(tw.fill(errtext))
            newTSUID = UID('1.2.840.10008.1.2')
            ds.file_meta.TransferSyntaxUID = newTSUID
            ds.is_little_endian = True
            ds.is_implicit_VR = True
    finally:
        return ds

##############################################################################

def element_value_format(elemval, elemVR):
    """
    Get element value string from initialisation text file in right format to
    set new tag value.
    
    N.B. this implementation does not support DICOM 'Sequence' (SQ) elements.
    
    """
    elemlen = len(elemval)
    numericVRs = {'DS','FD','FL','IS','OF','SL','SS','UI','UL','US'}
    if elemlen > 1:
        # Remove leading and trailing apostrophes...
        if elemval[0] == "'" and elemval[elemlen-1] == "'":
            newelemval = elemval[1:elemlen-1]
            # ...and evaluate for any of the VRs requiring numerical input.
            if elemVR in numericVRs:
                newelemval = eval(newelemval)
        # If data element is a list, evaluate the value string...
        elif elemval[0] == '[' and elemval[elemlen-1] == ']':
            newelemval = eval(elemval)
        # Also evaluate for any of the VRs requiring numerical input...
        elif elemVR in numericVRs:
            newelemval = eval(elemval)
        # ...otherwise, use simple input string format.
        else:
            newelemval = elemval
    elif elemVR in numericVRs:
        newelemval = eval(elemval)
    else:
        newelemval = elemval
    return newelemval

##############################################################################

def make_new_uid(UIDformat = 'UIDobj',
                 specifiedrandoms = None,
                 addtimestamp = True):
    """
    Adapted from a routine created by Jamie Fairfoul as part of the
    'cloneDCM' program.    
    
    General purpose function to generate a unique DICOM UID, using the PyDICOM
    root and appending two random numbers and a timestamp.
    
    'specifiedrandoms':
        If this is not None, it should be either a string of characters,
        containing only digits and dots, or a single floating point number.
        This will be used instead of generating new random numbers.
    'addtimestamp':
        If True, a timestamp will be added at the end of the UID.

    N.B. Not guaranteed unique, but has a very, very good chance.
    
    """
    root = dicom.UID.pydicom_root_UID
    if specifiedrandoms == None:
        randomsstr = (str(random.randrange(1, 999)) + '.' +
                      str(random.randrange(1, 999)))
    else:
        randomslist = list(str(specifiedrandoms))
        dotsdigits = [str(x) for x in range(10)]
        dotsdigits.append('.')
        nondotdiginds = []
        for elemnum in range(len(randomslist)):
            strelem = randomslist[elemnum]
            if strelem not in dotsdigits:
                nondotdiginds.append(elemnum)
        for elemnum in nondotdiginds:
            randomslist.pop(elemnum)
        endnums = [0, len(randomslist)-1]
        for elemnum in endnums:
            if randomslist[elemnum] == '.':
                randomslist.pop(elemnum)
        randomsstr = ''.join(randomslist)
    if addtimestamp:
        timestamp = dt.now().strftime("%Y%m%d%H%M%S%f")
        uidstring = root + randomsstr + '.' + timestamp
    else:
        uidstring = root + randomsstr
    
    if UIDformat in {'evalstring','evalstr','UIDstring','UIstring','UIstr'}:
        uidstring = 'UID(\'' + uidstring + '\')'
        return uidstring
    else:
        new_uid = UID(uidstring)
        return new_uid

##############################################################################

def replace_pixel_data(ds,newpixelmatrix):
    """
    Replace pixel data from DICOM dataset 'ds' with 'newpixelmatrix'.
    Return adjusted DICOM dataset.
    
    N.B. This doesn't include generating a new UID for the image - this
    should be done separately as required.
    
    """
    try:
        matrixdims = newpixelmatrix.ndim
        if matrixdims != 2:
            raise MatrixDimsError(matrixdims)
        newdtypename = newpixelmatrix.dtype.name
        accepteddtypes = {'int8','uint8','int16','uint16'}
        if newdtypename not in accepteddtypes:
            raise MatrixTypeError(newdtypename)
    except (NameError, AttributeError):
        print('The input pixel matrix should be a NumPy ndarray object.')
    except MatrixDimsError as m:
        if int(m) == 1:
            errtext = ('The input pixel matrix had 1 dimension. ' +
                       'Please provide a 2-dimensional array instead.')
            print(tw.fill(errtext))
        else:
            errtext = ('The input pixel matrix had ' + str(m) + ' dimensions.'
                       + 'Please provide a 2-dimensional array instead.')
            print(tw.fill(errtext))
    except MatrixTypeError as t:
        errtext = ('The data type of the input pixels is ' + str(t) +
                   ', which is not valid. (Must be 8-bit or 16-bit integer.)')
        print(tw.fill(errtext))
    else:
        try:
            dcmcols = ds.Columns
            dcmrows = ds.Rows
            dcmdtypename = ds.pixel_array.dtype.name
        except:
            print('The dataset input is not a valid PyDicom object.')
        else:
            newcols = newpixelmatrix.shape[1]
            newrows = newpixelmatrix.shape[0]
            if dcmcols != newcols:
                ds.Columns = newcols
            if dcmrows != newrows:
                ds.Rows = newrows
            if dcmdtypename != newdtypename:
                ds.pixel_array.dtype = newdtypename
                newdtypebits = 8*newpixelmatrix.dtype.itemsize
                ds.BitsAllocated = newdtypebits
                ds.BitsStored = newdtypebits
                ds.HighBit = newdtypebits - 1
            ds.PixelData = newpixelmatrix.tostring()
    finally:
        return ds

##############################################################################

def change_data_element(ds, tagstring, newelemval):
    """
    Changes the value of the data element specified by 'tagstring' to
    'newelemval' for the dataset, 'ds'. Returns the adjusted dataset.
    
    """
    try:
        tagparts = tagstring[1:len(tagstring)-1].split(',')
        tagparts = [part.strip() for part in tagparts]
        if len(tagstring.split('x')) == 1:
            tagstring = ('(0x' + tagparts[0] + 
                         ',0x' + tagparts[1] + ')')
        elemtag = Tag(eval(tagstring))
        elemVR = VR(elemtag)
    except:
        print('The specified tag is not valid.')
    else:
        if tagparts[0] == '0002':
            try:
                ds.file_meta[elemtag].value = element_value_format(newelemval,
                                                                   elemVR)
            except:
                ds.file_meta.add_new(elemtag, elemVR,
                                     element_value_format(newelemval, elemVR))
        else:
            try:
                ds[elemtag].value = element_value_format(newelemval, elemVR)
            except:
                ds.add_new(elemtag, elemVR,
                           element_value_format(newelemval, elemVR))
    finally:
        return ds

##############################################################################

def delete_data_element(ds, tagstring):
    """
    Deletes the data element specified by 'tagstring' for the dataset, 'ds'.
    Returns the adjusted dataset.
    
    """
    try:
        tagparts = tagstring[1:len(tagstring)-1].split(',')
        tagparts = [part.strip() for part in tagparts]
        if len(tagstring.split('x')) == 1:
            tagstring = ('(0x' + tagparts[0] + 
                         ',0x' + tagparts[1] + ')')
        elemtag = Tag(eval(tagstring))
    except:
        print('The specified tag is not valid.')
    else:
        if tagparts[0] == '0002':
            try:
                del ds.file_meta[elemtag]
            except:
                pass
        else:
            try:
                del ds[elemtag]
            except:
                pass
    finally:
        return ds

##############################################################################

def change_instance_UID(ds, newinstUID, UIDformat = 'UIDobj'):
    """
    Changes the SOP Instance UID of the provided DICOM dataset (ds), in both
    the file meta and the header, to the specified value (newinstUID).
    
    UIDformat:
        may be used to specify the format of the input UID; there are 3 main
        options.
        Simple string, selected using one of the following: 'string', 'str',
        'rawstring', 'rawstr':
            the provided newinstUID is a string containing simply the text of
            the UID, e.g. '1.2.826.0.1.3680043.8.498.1'.
        Evaluation string, selected using one of the following: 'evalstring',
        'evalstr', 'UIDstring', 'UIstring', 'UIstr':
            the provided newinstUID is a string of the form UID('...') - i.e.
            one which can be evaluated to a UID using the eval command.
        UID object (default), selected using any other value:
            the provided newinstUID is a PyDICOM UID object.
    
    """
    if UIDformat in {'evalstring','evalstr','UIDstring','UIstring','UIstr'}:
        UIDstring = newinstUID
        ds = change_data_element(ds, '(0002,0003)', UIDstring)
        ds = change_data_element(ds, '(0008,0018)', UIDstring)
    elif UIDformat in {'string','str','rawstring','rawstr'}:
        UIDstring = eval('UID(\'' + newinstUID + '\')')
        ds = change_data_element(ds, '(0002,0003)', UIDstring)
        ds = change_data_element(ds, '(0008,0018)', UIDstring)
    else:
        ds.file_meta.MediaStorageSOPInstanceUID = newinstUID
        ds.SOPInstanceUID = newinstUID
    return ds

##############################################################################

def get_dcm_files_list(filespath):
    """
    From a given input directory, return a list of filenames for all DICOM
    files present; if no DICOM files present, return None.
    
    """
    if not os.path.exists(filespath):
        print('The specified DICOM input file path is not valid.')
        return
    else:
        dcmfileslist = fileslist = ldir(filespath)
        for fname in fileslist:
            if not fnmatch(fname, '*.dcm'):
                dcmfileslist.pop(fileslist.index(fname))
        if dcmfileslist == []:
            dcmfileslist = None
        else:
            if not (filespath[len(filespath)-1] == "/" or
                    filespath[len(filespath)-1] == "\\"):
                filespath = filespath + r'/'
            dcmfileslist = [filespath + fname for fname in dcmfileslist]
        return dcmfileslist

##############################################################################

class MatrixDimsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)
    def __int__(self):
        return self.value

class MatrixTypeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

##############################################################################