# -*- coding: utf-8 -*-
"""
Script for automating the process of interactive vCTphantom creation, saving
the image sets to DICOM file series and producing a text file output which
provides useful information about the image series in the vCTphantom study.

Version 2.1 info:
- create_vCT_phantom now allows the user to apply deformations to the initial
  generated phantom and saves the resulting data to a new DICOM image set;
  this process may be repeated as many times as required.
- certain portions of the script code have been replaced by functions
  (interactive_add_shapes & interactive_deform_phantom) located in the new
  'vCT_interactive_creation' module; this should make the code more readable,
  and the script's workflow more understandable.

N.B. currently phantoms may only be created using square pixel matrices for
each slice.

Created by Stephen Skett as part of the KCL MSc Project, "Development of a
commissioning and QA framework for deformable image registration in
radiotherapy planning."

"""
from __future__ import print_function
import os
import sys
import vCTphantoms as vCTp
import vCT_interactive_creation as vCTi
from formlayout import fedit
import textwrap as tw

# New phase announcement!
print("")
print("Phase 1 of 3: creating vCTphantom object...")

# Set up fedit defaults for phantom data, depending on computer in use
computername = os.environ['COMPUTERNAME']
if "PCID" in computername:
    #PCH default location:
    defrootpath = r'N:/Prosoma33 Testing/TestPatients/vCTphantoms/'
elif computername=="STEPHEN-PC":
    #Home-PC file location:
    defrootpath = (r'E:/My Documents/Medical Physics/STP/KCL MSc Project/' +
                   r'vCTphantoms/')
else:
    print("Unrecognised calling computer: exiting...")
    sys.exit()

# Choose whether to load from DICOM file series or create new vCTphantom
print("---Subphase a: load DICOM series or make new phantom?")
newphantomfield = [('Load existing DICOM series or create new phantom?',
                    [1,'load existing series','create new phantom'])]
try:
    makenewphantom = fedit(newphantomfield,
                           title='Load from DICOM or create new phantom?')[0]
except:
    sys.exit()

if makenewphantom:
    # New subphase announcement!
    print("---Subphase b: get initial phantom specifications",end=" ")
    
    # Determine appropriate phantom naming defaults
    pidusedtest = True
    i = 0
    while pidusedtest:
        i += 1
        # Consider adjusting...
        defpatid = 'DIRphantom' + '%02d' % i
        if not os.path.exists(defrootpath + defpatid):
            pidusedtest = False
    defpatname = 'ZZZvCT^DIRphantom' + '%02d' % i
    
    # Get user-input phantom data using fedit
    phantomfields = [('Phantoms Root Path:             ', defrootpath),
                     ('Patient Name:                   ', defpatname),
                     ('Patient ID:                     ', defpatid),
                     ('Matrix Size:                    ', 128),
                     ('Number of Slices:               ', 100),
                     ('DICOM output pixel-size scaling:', 0.25)]
    phantomdata = fedit(phantomfields, title = 'Get phantom data',
                        comment = 'Specify phantom definition variables')
    if phantomdata == None: sys.exit()
    
    # Make vCTphantom object
    try:
        matrixsize = int(phantomdata[3])
        nslices = int(phantomdata[4])
        dcmpixelscaling = phantomdata[5]
        if matrixsize < 1 or nslices < 1 or dcmpixelscaling <= 0:
            raise ValueError
    except:
        print("")
        errtext = ('Matrix size and number of slices must be positive ' +
                   'integers, and DICOM output pixel-size scale factor ' +
                   'must be a real number greater than zero.')
        print(tw.fill(errtext))
        sys.exit()
    phantomsrootpath = str(phantomdata[0])
    phantomname = str(phantomdata[1])
    phantompatID = str(phantomdata[2])
    myphantom = vCTp.vCTphantom(matrixsize, matrixsize, nslices,
                                patientID=phantompatID,
                                patientname=phantomname)
    
    # Get DICOM file save details
    phantnameparts = phantomname.split('^')
    studypath = phantomsrootpath + phantnameparts[1] + r'/'
else:
    # New subphase announcement!
    print("---Subphase b: select location of DICOM series")
    
    # Get DICOM load directory data
    phantomfields = [('Phantoms Root Path:  ', defrootpath),
                     ('Study Directory:     ', 'DIRphantom01'),
                     ('Series Sub-Directory:', 'Series1')]
    phantomdata = fedit(phantomfields, title='Get phantom data',
                        comment='Specify phantom DICOM series load directory')
    if phantomdata == None: sys.exit()
    
    # New subphase announcement!
    print("---Subphase c: load DICOM series",end=" ")
    
    # Get phantom load data, and loop to find next available series number
    studypath = str(phantomdata[0]) + str(phantomdata[1]) + r'/'
    phantomloadpath = studypath + str(phantomdata[2]) + r'/'
    seriesusedtest = True
    i = 0
    while seriesusedtest:
        i+=1
        seriespath = "Series" + str(i)
        if not os.path.exists(studypath + seriespath):
            seriesusedtest = False
    # Load data from DICOM series into vCTphantom object
    myphantom, loadsuccess = vCTi.dcmread_vCTphantom_series(phantomloadpath,
                                                            seriesnum=i)
    if not loadsuccess: sys.exit()
    dcmpixelscaling = myphantom.dcmpixelscaling
print("...Done!")

# New phase announcement!
print("Phase 2 of 3: creating phantom structures...")

# Add shapes to phantom as required...
myphantom, mystudy = vCTi.interactive_add_shapes(myphantom,
                                                 phantomsavepath=studypath,
                                                 pixelscalefactor
                                                   =dcmpixelscaling)[:2]

# New phase announcement!
print("Phase 3 of 3: performing image deformations...")

# Perform deformations on phantom as required...
mystudy, allsavedok = vCTi.interactive_deform_phantom(myphantom, mystudy,
                                                      phantomsavepath
                                                        =studypath,
                                                      pixelscalefactor
                                                        =dcmpixelscaling)[1:]

# Finalise the vCT phantom study:
# - print completion confirmation message for the user to the command line;
# - create/amend study info text file, containing study & series UIDs, and
#   lists of the phantom shapes, deformations and any noise or blurring in
#   each series.
if allsavedok and isinstance(mystudy, vCTp.vCTstudy):
    studyinfotxtfpath = studypath + 'studyinfo.txt'
    numseries = mystudy.numseries
    print('All done! '+str(numseries)+' series saved; the series UID', end="")
    if makenewphantom:
        f = open(studyinfotxtfpath, 'w')
        f.write('Study summary: \n')
        f.write(' - study instance UID = '+str(mystudy)+'\n')
        f.write(' - study contains '+str(numseries)+' DICOM series; ')
    else:
        f = open(studyinfotxtfpath, 'a')
        f.write('\n')
        f.write('---Study amended--- \n')
        f.write(' - Phantom object loaded from '+str(phantomdata[2])+'\n')
        f.write(' - '+str(numseries)+' additional DICOM series added; ')
    if numseries == 1:
        print(' is:')
        f.write('its instance UID is: \n')
    else:
        print('s are:')
        f.write('their instance UIDs are: \n')
    for i in range(numseries):
        seriesiobj = mystudy.seriesobjects[i]
        print('    ' + str(seriesiobj))
        f.write('   -> Series ' + str(seriesiobj.seriesnum) + ' -- ' +
                                  str(seriesiobj) + '\n')
    if makenewphantom:
        f.write(' - vCTphantom image dimensions: \n')
        f.write('   -> matrix size = ' + str(myphantom.rows) + '\n')
        f.write('   -> number of slices = ' + str(myphantom.slices) + '\n')
        f.write(' - DICOM output pixel-size scaling factor = ' +
                str(dcmpixelscaling) + '\n')
        f.write('\n')
        f.write('Series details: \n')
    else:
        f.write('\n')
        f.write('Additional series details: \n')
    # Then write lists of shape, deformation and noise/blurring specs for
    # each series.
    firstnoise = True
    firstblur = True
    for i in range(numseries):
        seriesiobj = mystudy.seriesobjects[i]
        f.write(' - Series ' + str(seriesiobj.seriesnum) + ': \n')
        if seriesiobj.phantom_shapes == None:
            if makenewphantom:
                f.write('   -> No phantom shapes. \n')
            else:
                f.write('   -> No additional phantom shapes. \n')
        else:
            for j in range(1,len(seriesiobj.phantom_shapes)+1):
                phantomshapejdict = seriesiobj.phantom_shapes[j]
                if makenewphantom:
                    f.write('   -> Shape ' + str(j) + ': \n')
                else:
                    f.write('   -> Additional shape ' + str(j) + ': \n')
                f.write('      * class: ' + phantomshapejdict['class'] + '\n')
                f.write('      * parameters: \n')
                for k in phantomshapejdict.keys():
                    if k != 'class':
                        f.write('         # ' + str(k) + ' - ' +
                                str(phantomshapejdict[k]) + '\n')
        if seriesiobj.phantom_deformations == None:
            if makenewphantom:
                f.write('   -> No phantom deformations. \n')
            else:
                f.write('   -> No additional phantom deformations. \n')
        else:
            for j in range(1,len(seriesiobj.phantom_deformations)+1):
                phantomdefjdict = seriesiobj.phantom_deformations[j]
                if makenewphantom:
                    f.write('   -> Deformation ' + str(j) + ': \n')
                else:
                    f.write('   -> Additional deformation ' + str(j) + ': \n')
                f.write('      * class: ' + phantomdefjdict['class'] + '\n')
                f.write('      * parameters: \n')
                for k in phantomdefjdict.keys():
                    if k != 'class':
                        f.write('         # ' + str(k) + ' - ' +
                                str(phantomdefjdict[k]) + '\n')
        if seriesiobj.phantom_noise_blur == None:
            if makenewphantom:
                f.write('   -> No noise/blurring added. \n')
            else:
                f.write('   -> No further noise/blurring added. \n')
        else:
            for j in range(1,len(seriesiobj.phantom_noise_blur)+1):
                nbjdict = seriesiobj.phantom_noise_blur[j]
                nbjclass = nbjdict['class']
                if nbjclass=='blur' or nbjclass=='blurring':
                    if firstblur:
                        f.write('   -> Blurring added with sigma (pixels) ' +
                                '= ' + str(nbjdict['sigma']) + '\n')
                        firstblur = False
                    else:
                        f.write('   -> Further blurring added with sigma ' +
                                '(pixels) = ' + str(nbjdict['sigma']) + '\n')
                elif nbjclass=='noise':
                    if firstnoise:
                        f.write('   -> Noise added with amplitude (pixels) ' +
                                '= ' + str(nbjdict['amplitude']) + '\n')
                        firstnoise = False
                    else:
                        f.write('   -> Further noise added with amplitude ' +
                                '(pixels) = ' + str(nbjdict['amplitude']) +
                                '\n')
                else:
                    f.write('   -> Unrecognised noise/blurring added; ' +
                            'specification: \n')
                    f.write('      * class: ' + nbjclass + '\n')
                    f.write('      * parameters: \n')
                    for k in nbjdict.keys():
                        if k != 'class':
                            f.write('         # ' + str(k) + ' = ' +
                                    str(nbjdict[k]) + '\n')
    print('See study information text file for more about phantom contents.')
else:
    print('No new series saved.')
print("")