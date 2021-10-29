# -*- coding: utf-8 -*-
"""
Set of functions which enable the user to create vCTphantoms, containing user-
specified sets of vCTshapes and phantom deformations easily and efficiently,
and/or save them to/load them from DICOM image file series.

The phantoms are created and manipulated interactively using the simple dialog
boxes generated by the fedit function of the formlayout package.

Version 2.1 info:
- the functions in this new module have been moved out of the vCTfunctions
  module, to negate any confusion about the differing roles of the two sets
  of functions (i.e. functions in this module make direct calls to the
  vCTphantoms module, whereas those in vCTfunctions do not).

Created by Stephen Skett as part of the KCL MSc Project, "Development of a
commissioning and QA framework for deformable image registration in
radiotherapy planning."

"""
from __future__ import print_function
import dicom
import numpy as np
import dcmfunctions as df
from datetime import datetime as dt
from os import makedirs as mdir
import os
import vCTphantoms as vCTp
import vCTfunctions as vCTf
from scipy import ndimage as ndimg
from formlayout import fedit
import textwrap as tw

#import pdb

##############################################################################

def save_vCTphantom_series(vCTphantomobj, phantomsavepath, vCTstudyobj=None,
                           pixelscalefactor=1):
    """
    For a given vCTphantom object, save all image slices to DICOM files.
    The vCT phantom may represent one series which is part of a larger study
    set: in such case, include specified study- & series-specific elements.
    
    N.B. currently only phantoms using square pixel matrices are allowed.
    
    """
    # Read in test DICOM file to PyDICOM ds (location depending on computer
    # in use).
    computername = os.environ['COMPUTERNAME']
    if "PCID" in computername:
        #PCH file location:
        dcmfname = r'J:/USERS/sketts/KCL MSc Project/Python Code/inittest.dcm'
    elif computername=="STEPHEN-PC":
        #Home-PC file location:
        dcmfname = (r'E:/My Documents/Medical Physics/STP/KCL MSc Project/' +
                    r'Python Code/inittest.dcm')
    else:
        print("Unrecognised calling computer: exiting...")
        return False, vCTstudyobj
    ds = dicom.read_file(dcmfname)
    
    # Test whether phantom is square.
    if vCTphantomobj.rows == vCTphantomobj.columns:
        matrixsize = vCTphantomobj.rows
    else:
        errtext = ("Cannot save series: phantoms with non-square slices " +
                   "not currently handled.")
        print(tw.fill(errtext))
        return False, vCTstudyobj
    
    # Generate initialised PyDICOM ds
    ds = df.setup_DICOM_defaults(ds)
    del ds.AdditionalPatientHistory
    
    # PyDICOM ds header: define study-specific variables.
    if not isinstance(vCTstudyobj, vCTp.vCTstudy):
        vCTstudyobj = vCTp.vCTstudy()
    nslices = vCTphantomobj.slices
    ds.PatientName = vCTphantomobj.patientname
    ds.PatientID = vCTphantomobj.patientID
    ds.PatientBirthDate = vCTstudyobj.studydate
    ds.StudyDate = vCTstudyobj.studydate
    ds.StudyTime = vCTstudyobj.studytime
    ds.StudyInstanceUID = vCTstudyobj.studyuid
    ds.DataCollectionDiameter = round(matrixsize*0.71, -1)
    ds.ReconstructionDiameter = ds.DataCollectionDiameter
    if pixelscalefactor == 1:
        pixelsize = ds.SpacingBetweenSlices
    else:
        pixelsize = ds.SpacingBetweenSlices*pixelscalefactor
        ygrid, xgrid = (np.mgrid[0:vCTf.nint(matrixsize/pixelscalefactor),
                                 0:vCTf.nint(matrixsize/pixelscalefactor)]
                                                       *pixelscalefactor)
    ds.PixelSpacing[0] = ds.PixelSpacing[1] = pixelsize
    ds.ImagePositionPatient[0] = -0.5*pixelsize*matrixsize
    ds.ImagePositionPatient[1] = -0.5*pixelsize*matrixsize
    
    # PyDICOM ds header: define series-specific variables.
    studyrandoms = vCTstudyobj.studyrandoms
    vCTseriesobj = vCTp.vCTseries(vCTphantomobj, studyrandoms=studyrandoms)
    seriesnum = vCTseriesobj.seriesnum
    if seriesnum==None:
        errtext = ("Cannot save series: check that the supplied " +
                   "vCTphantom object is valid.")
        print(tw.fill(errtext))
        return False, vCTstudyobj
    ds.FrameOfReferenceUID = vCTseriesobj.seriesFoRuid
    ds.SeriesNumber = ds.AcquisitionNumber = seriesnum
    ds.SeriesInstanceUID = vCTseriesobj.seriesuid
    ds.SeriesTime = dt.now().strftime("%H%M%S.%f")
    
    # Set up save directory.
    dcmsavepath = phantomsavepath + "Series" + str(seriesnum) + r'/'
    if not os.path.exists(dcmsavepath):
        mdir(dcmsavepath)
    
    # Loop through slices in phantom object, saving each to DICOM file.
    phantnameparts = ds.PatientName.split('^')
    slicezeroloc = ds.SliceLocation
    progbarnum = 0
    savesuccess = []
    for i in range(nslices):
        slicei = vCTphantomobj.get_slice(i)
        if pixelscalefactor != 1:
            slicei = ndimg.map_coordinates(slicei, [ygrid, xgrid], order=1)
        ds = df.replace_pixel_data(ds, slicei)
        ds.ContentTime = dt.now().strftime("%H%M%S.%f")
        ds.InstanceNumber = i + 1
        ds.SliceLocation = slicezeroloc+(i-nslices/2)*ds.SpacingBetweenSlices
        ds.ImagePositionPatient[2] = ds.SliceLocation
        ds = df.change_instance_UID(ds, df.make_new_uid())
        try:
            ds.save_as(dcmsavepath + 
                       phantnameparts[0] + "_" + 
                       phantnameparts[1] + "-vCT_" + 
                       str(seriesnum) + "." + str(i) + ".dcm",
                       WriteLikeOriginal = False)
        except:
            errtext = ("Save problem: check that specified phantom Patient " +
                       "Name is valid (must contain a '^' character)")
            print("")
            print(tw.fill(errtext))
            savesuccess.append(False)
            break
        else:
            progbarnum = vCTf.increment_progbar(i, nslices, progbarnum)
            savesuccess.append(True)
    print("")
    
    # Return list: first element indicates whether all images in the series
    # were saved correctly, and if not, at what point in the series the DICOM
    # save failed; second element is the up-to-date study object.
    if savesuccess[len(savesuccess)-1]:
        vCTstudyobj.add_series(vCTseriesobj)
    return savesuccess, vCTstudyobj

##############################################################################

def dcmread_vCTphantom_series(dcmseriespath, seriesnum=1, phantombits=None):
    """
    Reads data into a vCTphantom from a DICOM images series.
    
    N.B.1. This assumes that all DICOM images saved in the specified directory
    are from the same DICOM series.
    N.B.2. Currently only phantoms using square pixel matrices are allowed.
    
    """
    dcmfileslist = df.get_dcm_files_list(dcmseriespath)
    try:
        if dcmfileslist == None: raise IOError
        numslices = len(dcmfileslist)
        # Set up vCTphantom object based on 1st file in dcmfileslist. Check
        # that the pixel data is valid for inclusion in a vCTphantom object.
        sliceds = dicom.read_file(dcmfileslist[0])
        seriesUID = sliceds.SeriesInstanceUID
        vCTphantomobj = get_vCTphantom_slice_from_ds(sliceds,
                                                     seriesnum=seriesnum,
                                                     numslices=numslices,
                                                     phantombits=phantombits)
        # Loop through all remaining files in dcmfileslist, adding to the
        # vCTphantom object's pixel data each time. For each file, check
        # that the DICOM Series Instance UIDs are the same.
        for i in range(numslices):
            sliceds = dicom.read_file(dcmfileslist[i])
            sliceiseriesUID = sliceds.SeriesInstanceUID
            if sliceiseriesUID != seriesUID: raise SeriesUIDError
            vCTphantomobj = get_vCTphantom_slice_from_ds(sliceds,
                                                         vCTphantomobj
                                                             =vCTphantomobj)
    except IOError:
        print("")
        errtext = ("No DICOM files found, or number of slices not " +
                   "specified when required; cannot add data to vCTphantom.")
        print(tw.fill(errtext))
        vCTphantomobj = None
        readsuccess = False
    except SliceDimsError:
        print("")
        print("Slice not added: only square slice pixel matrices allowed.")
        vCTphantomobj = None
        readsuccess = False
    except PixelDimsError:
        print("")
        print("Slice not added: the image pixel scaling does not match.")
    except SeriesUIDError:
        print("")
        print("Slice not added: must be from the same series.")
        vCTphantomobj = None
        readsuccess = False
    except:
        print("")
        errtext = ("Slice not added: check that slice dataset and/or " +
                   "vCTphantom object are valid.")
        print(tw.fill(errtext))
        vCTphantomobj = None
        readsuccess = False
    else:
        readsuccess = True
    finally:
        return vCTphantomobj, readsuccess

##############################################################################

def get_vCTphantom_slice_from_ds(ds, vCTphantomobj=None, seriesnum=1,
                                 numslices=None, phantombits=None):
    """
    Get slice data from a PyDICOM dataset and add to a vCTphantom object.
    
    """
    slicepixdims = ds.pixel_array.shape
    if ds.PixelSpacing[0] == ds.PixelSpacing[1]:
        pixelscaling = ds.PixelSpacing[0]/ds.SpacingBetweenSlices
    else:
        raise SliceDimsError
    matrixsize = slicepixdims[0]*pixelscaling
    if vCTphantomobj==None:
        if phantombits==None: phantombits = 16
        if numslices==None: raise IOError
        vCTphantomobj = vCTp.vCTphantom(matrixsize, matrixsize, numslices,
                                        seriesnum=seriesnum,
                                        patientID=ds.PatientID,
                                        patientname=ds.PatientName,
                                        bits=phantombits)
        vCTphantomobj.dcmpixelscaling = pixelscaling
    else:
        if hasattr(vCTphantomobj,'dcmpixelscaling'):
            if pixelscaling != vCTphantomobj.dcmpixelscaling:
                raise PixelDimsError
        else:
            if pixelscaling != 1:
                raise PixelDimsError
            else:
                vCTphantomobj.dcmpixelscaling = 1
    if pixelscaling == 1:
        slicepixels = ds.pixel_array
    else:
        ygrid, xgrid = (np.mgrid[0:vCTf.nint(matrixsize),
                                 0:vCTf.nint(matrixsize)]/pixelscaling)
        slicepixels = ndimg.map_coordinates(ds.pixel_array, [ygrid, xgrid],
                                            order=1)
    slicenumber = ds.InstanceNumber - 1
    writesuccess = vCTphantomobj.write_slice(slicenumber, slicepixels)
    if not writesuccess: raise Exception
    return vCTphantomobj

##############################################################################

def subphase_announcer(subphasenum, subphasetext, phantomobj=None):
    """
    Returns a string containing an appropriately named & labelled sub-phase
    announcement.
    
    """
    if phantomobj==None:
        alphabetdict = vCTp.vCTrefdata().alphabetdict
    else:
        alphabetdict = phantomobj.ref_data.alphabetdict
    if subphasenum <= 26:
        subphasestr = ("---Subphase " + alphabetdict[subphasenum] + ": " +
                       subphasetext)
    else:
        subphasestr = ("---Subphase " + alphabetdict[subphasenum/26] +
                       alphabetdict[subphasenum%26] + ": " + subphasetext)
    return subphasestr     

##############################################################################

def interactive_noise_blur(phantomobj, subphasenum=1):
    """
    Function which allows the user to add random, uniformly-distributed noise
    and/or perform Gaussian blurring to the pixel data of a vCTphantom,
    interactively specifying the parameters using fedit dialog boxes.
    
    N.B. the try...except...else structures in this function are included to
    catch errors resulting from the user cancelling out of the fedit window.
    If the fedit is cancelled, nothing happens; otherwise, the specified
    niose/blurring is added.
    
    """
    
    # Perform Gaussian blurring to phantom pixels, with user-specified sigma
    blurcomment = ('N.B. if blurring has already been applied earlier, ' +
                   'this is not recommended.')
    try:
        sigma = fedit([('Gaussian Blur Kernel Sigma (0.0 = no blurring)',
                        0.0)], title = 'Apply Gaussian blur?',
                      comment = tw.fill(blurcomment, width=60))[0]
    except:
        pass
    else:
        print(subphase_announcer(subphasenum, "add blurring..."))
        subphasenum += 1
        phantomobj.apply_blurring(sigma)
    
    # Add random noise to the image
    noisecomment = ('N.B. if random noise has already been added ' +
                    'earlier, this is not recommended.')
    try:
        amplitude = fedit([('Additive noise max. amplitude (0 = no noise)',
                            0)], title = 'Add random noise?',
                          comment = tw.fill(noisecomment, width=60))[0]
    except:
        pass
    else:
        print(subphase_announcer(subphasenum, "add noise..."))
        subphasenum += 1
        phantomobj.add_random_noise(amplitude)
    
    return phantomobj, subphasenum

##############################################################################

def interactive_add_shapes(phantomobj, studyobj=None, **savekws):
    """
    Function which allows the user to add shapes to a vCTphantom object,
    interactively specifying the shape type(s) and parameters using fedit
    dialog boxes. 
    
    """
    addmoreshapes = True
    firstshapeloop = True
    subphasenum = 1
    allsavestatus = True
    shapetypelist = phantomobj.ref_data.shapetypelist
    shapetypespecs = phantomobj.ref_data.shapetypespecs
    shapespecsdict = dict(zip(shapetypelist, shapetypespecs))
    shapetypefield = [('Shape Type', [6] + shapetypelist)]
    savedialogtitle = ('Save the new phantom image set?')
    savedialogcomment = ('N.B. This can also be done later if desired, ' +
                         'after further shapes have been added.')
    
    #pdb.set_trace()
    
    while addmoreshapes:
        # Ask user whether (more) shapes are required...
        if firstshapeloop:
            feditdata = fedit([('Add a shape? (Yes or no)', 'Yes')],
                              title = 'Add shapes?')
            firstshapeloop = False
        else:
            feditdata = fedit([('Add another shape? (Yes or no)', 'Yes')],
                              title = 'Add more shapes?')
        if feditdata == None:
            print('Operation cancelled: exiting.')
            addmoreshapes = False
            break
        addshapetest = feditdata[0]
        if addshapetest.lower() not in phantomobj.ref_data.yesdict:
            addmoreshapes = False
            break
        
        # Subphase announcement! (Check subphasenum & do arithmetic if needed)
        print(subphase_announcer(subphasenum, "select & add shapes..."))
        subphasenum += 1
        
        # ...If so, ask user what shape type.
        feditdata = fedit(shapetypefield, title='Enter desired shape type')
        if feditdata == None:
            print('Operation cancelled: exiting.')
            addmoreshapes = False
            break
        shapenum = feditdata[0]
        shapetype = shapetypelist[shapenum]
        shapespecslist = shapespecsdict[shapetype]
        
        # Set up fedit default strings (from the vCTphantom object's
        # 'shapespecinputdefaults' dictionary) & get user input shape params.
        shapespecdefs = []
        for param in shapespecslist:
            if param == 'Vertical Radius (Optional)':
                param = 'Vertical Radius'
            shapespecdefs.append(
                    phantomobj.ref_data.shapespecinputdefaults[param])
        shapespecfields = zip(shapespecslist, shapespecdefs)
        shapespecdata = fedit(shapespecfields,
                              title = 'Enter shape specification parameters')
        if shapespecdata == None:
            print('Operation cancelled: exiting.')
            addmoreshapes = False
            break
        
        # Sort out shape params and set new shapespecinputdefaults strings.
        # (N.B. this doesn't check whether the right number of shape params
        # have been supplied; the vCTphantoms 'add_shape' method will throw up
        # errors if the wrong number of params is entered -- but then, the
        # user is expected to be at least a little bit intelligent.)
        for entry in shapespecdata:
            entryind = shapespecdata.index(entry)
            if entry != '':
                try:
                    entry = eval(entry)
                    if entryind == 0:
                        if not isinstance(entry, (list, tuple)):
                            raise vCTp.ShapePosnError
                        for coord in entry:
                            if not isinstance(coord, (int, float)):
                                raise vCTp.ShapeVarsError
                    else:
                        if not isinstance(entry, (int, float)):
                            raise vCTp.ShapeVarsError
                except vCTp.ShapePosnError:
                    errtext = ("Shape addition failed: the shape position " +
                               "variable (either 'Front-Top-Left Corner' " +
                               "or 'Centroid Position') must be a Python " +
                               "sequence (either a list or a tuple).")
                    print(tw.fill(errtext))
                    break
                except vCTp.ShapeVarsError:
                    errtext = ("Shape addition failed: shape specification " +
                               "entries must be numeric.")
                    print(tw.fill(errtext))
                    break
            param = shapespecslist[entryind]
            if param == 'Vertical Radius (Optional)':
                phantomobj.ref_data.shapespecinputdefaults[
                    'Vertical Radius'] = shapespecdata[entryind]
            else:
                phantomobj.ref_data.shapespecinputdefaults[
                    param] = shapespecdata[entryind]
            shapespecdata[entryind] = entry
        if shapetype == 'cylinder' and shapespecdata[2] == '':
            shapespecdata.pop(2)
        
        # Add shape to phantom and report success/failure to user
        phantomobj.add_shape(shapetype,
                             shapespecdata[0:(len(shapespecdata)-1)],
                             pixelval = shapespecdata[len(shapespecdata)-1])
        print("Success! Your specified " + shapetype + " was added.")
        
        # Ask user if they wish to save the new image set
        if not savekws.has_key('phantomsavepath'):
            continue
        else:
            phantomsavepath = savekws['phantomsavepath']
        try:
            saveimagestest = fedit([('Save image set? (Yes or no)','Yes')],
                                   title = savedialogtitle,
                                   comment = savedialogcomment)[0]
            if saveimagestest.lower() not in phantomobj.ref_data.yesdict:
                raise IOError
            pixelscalefactor = savekws['pixelscalefactor']
        except KeyError:
            pixelscalefactor = 1
        except (IOError, TypeError):
            continue
        except:
            break
        
        # Subphase announcement! (Check subphasenum & do arithmetic if needed)
        print(subphase_announcer(subphasenum, "save images if required..."))
        subphasenum += 1
        
        # Save images and get updated vCTstudy object
        savestatus, studyobj = save_vCTphantom_series(phantomobj,
                                                      phantomsavepath,
                                                      vCTstudyobj=studyobj,
                                                      pixelscalefactor
                                                        =pixelscalefactor)
        numsaves = len(savestatus)
        if not savestatus[numsaves-1]:
            print("DICOM save failed at series " + 
                  str(phantomobj.seriesnum) + ", file " +
                  str(numsaves) + "; exiting.")
            allsavestatus = False
            addmoreshapes = False
            break
        else:
            phantomobj.seriesnum += 1
        
    # Add random noise and/or blurring as required
    phantomobj, subphasenum = interactive_noise_blur(phantomobj,
                                                     subphasenum=subphasenum)
    
    # Ask user if they wish to save the new image set
    try:
        phantomsavepath = savekws['phantomsavepath']
        saveimagestest = fedit([('Save image set? (Yes or no)','Yes')],
                               title = savedialogtitle,
                               comment = savedialogcomment)[0]
        if saveimagestest.lower() not in phantomobj.ref_data.yesdict:
            raise IOError 
    except:
        pass
    else:
        try:
            pixelscalefactor = savekws['pixelscalefactor']
        except:
            pixelscalefactor = 1
        # Save images and get updated vCTstudy object
        print(subphase_announcer(subphasenum, "save images if required..."))
        savestatus, studyobj = save_vCTphantom_series(phantomobj,
                                                      phantomsavepath,
                                                      vCTstudyobj=studyobj,
                                                      pixelscalefactor
                                                        =pixelscalefactor)
        numsaves = len(savestatus)
        if not savestatus[numsaves-1]:
            print("DICOM save failed at series " + 
                  str(phantomobj.seriesnum) + ", file " +
                  str(numsaves) + "; exiting.")
            allsavestatus = False
        else:
            phantomobj.seriesnum += 1
    finally:
        return phantomobj, studyobj, allsavestatus

##############################################################################

def interactive_deform_phantom(phantomobj, studyobj, **savekws):
    """
    Function which allows the user to perform deformations to the pixel data
    of a vCTphantom object, interactively specifying the deformation type(s)
    and parameters using fedit dialog boxes. 
    
    """
    domoredeformations = True
    firstdeformationloop = True
    deformtypelist = phantomobj.ref_data.deformtypelist
    subphasenum = 1
    allsavestatus = True
    
    while domoredeformations:
        # Ask user whether (more) defomations are required...
        if firstdeformationloop:
            feditdata = fedit([('Apply deformation to the image? ' +
                                '(Yes or no)', 'Yes')],
                              title = 'Apply deformation?')
            firstdeformationloop = False
        else:
            feditdata = fedit([('Apply another deformation to the image? ' +
                                '(Yes or no)', 'No')],
                               title = 'Apply another deformation?')
        if feditdata == None:
            print('Operation cancelled: exiting.')
            domoredeformations = False
            break
        deformtest = feditdata[0]
        if deformtest.lower() not in phantomobj.ref_data.yesdict:
            domoredeformations = False
            break
        
        # Subphase announcement! (Check subphasenum & do arithmetic if needed)
        print(subphase_announcer(subphasenum,
                                 "select & apply image deformation..."))
        subphasenum += 1
        
        # Ask user what deformation type they require
        deformtypefield = [('Deformation Type', [0] + deformtypelist)]
        deformtitle = 'Enter desired deformation type'
        feditdata = fedit(deformtypefield, title=deformtitle)
        if feditdata == None:
            print('Operation cancelled: exiting.')
            domoredeformations = False
            break
        deformnum = feditdata[0]
        deformtype = deformtypelist[deformnum]
        if deformtype == 'stretch':
            deformtitle = 'Please specify the stretching parameters'
            deformfields = [('Stretch factor, k              ', 1.0),
                            ('Horizontal offset (pixels)     ', 0.0),
                            ('Vertical offset (pixels)       ', 0.0),
                            ('Stretching axis angle (degrees)', 0.0)]
            deformvars = fedit(deformfields, title = deformtitle)
        elif deformtype == 'shear':
            deformtitle = 'Please specify the shearing parameters'
            deformfields = [('Shear factor, k              ', 1.0),
                            ('Shearing axis angle (degrees)', 0.0)]
            deformvars = fedit(deformfields, title = deformtitle)
        else:
            # Add code here if/when other deformations are implemented.
            pass
        if deformvars == None:
            print('Operation cancelled: exiting.')
            domoredeformations = False
            break
        
        # Call routines for adding deformations to phantom...
        print("Performing deformation...")
        deformsuccess = phantomobj.deform_phantom(deformtype, deformvars)
        if deformsuccess:
            print("Success! Your specified " + deformtype + " was applied.")
        else:
            continue
        
        # Add noise and/or blurring...
        phantomobj, subphasenum = interactive_noise_blur(phantomobj,
                                                         subphasenum
                                                           =subphasenum)
        
        # Ask user if they wish to save the deformed image set
        if not savekws.has_key('phantomsavepath'):
            continue
        else:
            phantomsavepath = savekws['phantomsavepath']
        savedialogtitle = ('Save the deformed phantom image set?')
        savedialogcomment = ('N.B. This can also be done later if desired, ' +
                             'after further image deformations have been ' +
                             'applied.')
        try:
            saveimagestest = fedit([('Save image set? (Yes or no)','Yes')],
                                   title = savedialogtitle,
                                   comment = savedialogcomment)[0]
            if saveimagestest.lower() not in phantomobj.ref_data.yesdict:
                raise IOError
            pixelscalefactor = savekws['pixelscalefactor']
        except KeyError:
            pixelscalefactor = 1
        except:
            break
        
        # Subphase announcement! (Check subphasenum & do arithmetic if needed)
        print(subphase_announcer(subphasenum, "save images if required..."))
        subphasenum += 1
        
        # Save images and get updated vCTstudy object
        savestatus, studyobj = save_vCTphantom_series(phantomobj,
                                   phantomsavepath, vCTstudyobj=studyobj,
                                   pixelscalefactor=pixelscalefactor)
        numsaves = len(savestatus)
        if not savestatus[numsaves-1]:
            print("DICOM save failed at series " + 
                  str(phantomobj.seriesnum) + ", file " +
                  str(numsaves) + "; exiting.")
            allsavestatus = False
            domoredeformations = False
            break
        else:
            phantomobj.seriesnum += 1
    
    # Return tuple: first element indicates whether all images in the study
    # were saved correctly; second element is the updated study object; third
    # element is a boolean value indicating whether all files saved correctly
    # [if an error occurred during saving, returns False; if saving was not
    # requested, returns True].
    return phantomobj, studyobj, allsavestatus

##############################################################################

class SliceDimsError(Exception):
    def __init__(self):
        pass

class PixelDimsError(Exception):
    def __init__(self):
        pass

class SeriesUIDError(Exception):
    def __init__(self):
        pass

##############################################################################