import casac
from tasks import *
from taskinit import *
import pylab as pl
import os
from math import pi,floor
import pyfits as fits
import glob
import numpy as np

useSMAScanNumbers = False
default(importfitsidi)
default(flagdata)
# Definitions needed to create a SOURCE table in the final MS:
sourceTableDesc = {
'DIRECTION':             {'comment': 'Direction (e.g. RA, DEC).',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': 1,
                          'option': 5,
                          'shape': np.array([2], dtype=np.int32),
                          'valueType': 'double'},
'PROPER_MOTION':         {'comment': 'Proper motion',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': 1,
                          'option': 5,
                          'shape': np.array([2], dtype=np.int32),
                          'valueType': 'double'},
'CALIBRATION_GROUP':     {'comment': 'Number of grouping for calibration purpose.',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'int'},
'CODE':                  {'comment': 'Special characteristics of source, e.g. Bandpass calibrator',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'string'},
'INTERVAL':              {'comment': 'Interval of time for which this set of parameters is accurate',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'double'},
'NAME':                  {'comment': 'Name of source as given during observations',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'string'},
'NUM_LINES':             {'comment': 'Number of spectral lines',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'int'},
'SOURCE_ID':             {'comment': 'Source id',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'int'},
'SPECTRAL_WINDOW_ID':    {'comment': 'ID for this spectral window setup',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'int'},
'TIME':                  {'comment': 'Midpoint of time for which this set of parameters is accurate.',
                          'dataManagerGroup': 'StandardStMan',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'option': 0,
                          'valueType': 'double'},
'POSITION':              {'comment': 'Position (e.g. for solar system objects',
                          'dataManagerGroup': 'source table col st man',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': -1,
                          'option': 0,
                          'valueType': 'double'},
'TRANSITION':            {'comment': 'Line Transition name',
                          'dataManagerGroup': 'source table col st man',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': -1,
                          'option': 0,
                          'valueType': 'string'},
'REST_FREQUENCY':        {'comment': 'Line rest frequency',
                          'dataManagerGroup': 'source table col st man',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': -1,
                          'option': 0,
                          'valueType': 'double'},
'SYSVEL':                {'comment': 'Systemic velocity at reference',
                          'dataManagerGroup': 'source table col st man',
                          'dataManagerType': 'StandardStMan',
                          'maxlen': 0,
                          'ndim': -1,
                          'option': 0,
                          'valueType': 'double'},
 '_define_hypercolumn_': {}}

#sourceTableDesc = {'CALIBRATION_GROUP': {'comment': 'Number of grouping for calibration purpose.',
#                       'dataManagerGroup': 'StandardStMan',
#                       'dataManagerType': 'StandardStMan',
#                       'maxlen': 0,
#                       'option': 0,
#                       'valueType': 'int'},
# 'CODE': {'comment': 'Special characteristics of source, e.g. Bandpass calibrator',
#          'dataManagerGroup': 'StandardStMan',
#          'dataManagerType': 'StandardStMan',
#          'maxlen': 0,
#          'option': 0,
#          'valueType': 'string'},
# 'DIRECTION': {'comment': 'Direction (e.g. RA, DEC).',
#               'dataManagerGroup': 'StandardStMan',
#               'dataManagerType': 'StandardStMan',
#               'maxlen': 0,
#               'ndim': 1,
#               'option': 5,
#               'shape': np.array([2], dtype=np.int32),
#               'valueType': 'double'},
# 'INTERVAL': {'comment': 'Interval of time for which this set of parameters is accurate',
#              'dataManagerGroup': 'StandardStMan',
#              'dataManagerType': 'StandardStMan',
#              'maxlen': 0,
#              'option': 0,
#              'valueType': 'double'},
# 'NAME': {'comment': 'Name of source as given during observations',
#          'dataManagerGroup': 'StandardStMan',
#          'dataManagerType': 'StandardStMan',
#          'maxlen': 0,
#          'option': 0,
#          'valueType': 'string'},
# 'NUM_LINES': {'comment': 'Number of spectral lines',
#               'dataManagerGroup': 'StandardStMan',
#               'dataManagerType': 'StandardStMan',
#               'maxlen': 0,
#               'option': 0,
#               'valueType': 'int'},
# 'POSITION': {'comment': 'Position (e.g. for solar system objects',
#              'dataManagerGroup': 'source table col st man',
#              'dataManagerType': 'StandardStMan',
#              'maxlen': 0,
#              'ndim': -1,
#              'option': 0,
#              'valueType': 'double'},
# 'PROPER_MOTION': {'comment': 'Proper motion',
#                   'dataManagerGroup': 'StandardStMan',
#                   'dataManagerType': 'StandardStMan',
#                   'maxlen': 0,
#                   'ndim': 1,
#                   'option': 5,
#                   'shape': np.array([2], dtype=np.int32),
#                   'valueType': 'double'},
# 'REST_FREQUENCY': {'comment': 'Line rest frequency',
#                    'dataManagerGroup': 'source table col st man',
#                    'dataManagerType': 'StandardStMan',
#                    'maxlen': 0,
#                    'ndim': -1,
#                    'option': 0,
#                    'valueType': 'double'},
# 'SOURCE_ID': {'comment': 'Source id',
#               'dataManagerGroup': 'StandardStMan',
#               'dataManagerType': 'StandardStMan',
#               'maxlen': 0,
#               'option': 0,
#               'valueType': 'int'},
# 'SPECTRAL_WINDOW_ID': {'comment': 'ID for this spectral window setup',
#                        'dataManagerGroup': 'StandardStMan',
#                        'dataManagerType': 'StandardStMan',
#                        'maxlen': 0,
#                        'option': 0,
#                        'valueType': 'int'},
# 'SYSVEL': {'comment': 'Systemic velocity at reference',
#            'dataManagerGroup': 'source table col st man',
#            'dataManagerType': 'StandardStMan',
#            'maxlen': 0,
#            'ndim': -1,
#            'option': 0,
#            'valueType': 'double'},
# 'TIME': {'comment': 'Midpoint of time for which this set of parameters is accurate.',
#          'dataManagerGroup': 'StandardStMan',
#          'dataManagerType': 'StandardStMan',
#          'maxlen': 0,
#          'option': 0,
#          'valueType': 'double'},
# 'TRANSITION': {'comment': 'Line Transition name',
#                'dataManagerGroup': 'source table col st man',
#                'dataManagerType': 'StandardStMan',
#                'maxlen': 0,
#                'ndim': -1,
#                'option': 0,
#                'valueType': 'string'},
# '_define_hypercolumn_': {}}



sourceDMInfo = {'*1': {'COLUMNS': np.array(['DIRECTION', 'PROPER_MOTION', 'CALIBRATION_GROUP', 'CODE', 'INTERVAL',
       'NAME', 'NUM_LINES', 'SOURCE_ID', 'SPECTRAL_WINDOW_ID',
       'TIME'], 
      dtype='|S19'),
        'NAME': 'StandardStMan',
        'SEQNR': 0,
        'SPEC': {'ActualCacheSize': 2,
                 'BUCKETSIZE': 2816,
                 'IndexLength': 166,
                 'PERSCACHESIZE': 2},
        'TYPE': 'StandardStMan'},
 '*2': {'COLUMNS': np.array(['POSITION', 'TRANSITION', 'REST_FREQUENCY', 'SYSVEL'], 
      dtype='|S15'),
        'NAME': 'source table col st man',
        'SEQNR': 1,
        'SPEC': {'ActualCacheSize': 2,
                 'BUCKETSIZE': 1152,
                 'IndexLength': 166,
                 'PERSCACHESIZE': 2},
        'TYPE': 'StandardStMan'}}

#sourceDMInfo = {'*1': {'COLUMNS': np.array(['CALIBRATION_GROUP', 'CODE', 'DIRECTION', 'INTERVAL', 'NAME',
#       'NUM_LINES', 'PROPER_MOTION', 'SOURCE_ID', 'SPECTRAL_WINDOW_ID',
#       'TIME'], 
#      dtype='|S19'),
#        'NAME': 'StandardStMan',
#        'SEQNR': 0,
#        'SPEC': {'ActualCacheSize': 2,
#                 'BUCKETSIZE': 2816,
#                 'IndexLength': 166,
#                 'PERSCACHESIZE': 2},
#        'TYPE': 'StandardStMan'},
# '*2': {'COLUMNS': np.array(['POSITION', 'REST_FREQUENCY', 'SYSVEL', 'TRANSITION'], 
#      dtype='|S15'),
#        'NAME': 'source table col st man',
#        'SEQNR': 1,
#        'SPEC': {'ActualCacheSize': 2,
#                 'BUCKETSIZE': 1152,
#                 'IndexLength': 166,
#                 'PERSCACHESIZE': 2},
#        'TYPE': 'StandardStMan'}}

sidebandsToProcess = []
files = glob.glob('tempFITS-IDI*')
bandList=[]
for b in files:
    band=int(b[b.find('.band')+5:])
    if band not in bandList:
        bandList.append(band)
lower = glob.glob('*Lower*')
upper = glob.glob('*Upper*')
if len(lower) > 0:
    sidebandsToProcess.append('Lower')
if len(upper) > 0:
    sidebandsToProcess.append('Upper')

sourceList = []
uniqueSourceNumbers = []
uniqueSourceNames = []
uniqueSourceRAs = []
uniqueSourceDecs = []
try:
    for line in open('sourceTable'):
        tok = line.split()
        sourceList.append(tok[1])
        if tok[0] not in uniqueSourceNumbers:
            uniqueSourceNumbers.append(tok[0])
            uniqueSourceNames.append(tok[1])
            uniqueSourceRAs.append(tok[2])
            uniqueSourceDecs.append(tok[3])
except IOError:
    print 'No source table found - the SrcId field will not be set properly'

# Parameters for flagdata:

for sideband in sidebandsToProcess:
    # First figure out which chunks need to be processed
    chunksToProcess = []
    for i in range(max(bandList)+1):
        try:
            testFile = open('./tempFITS-IDI%s.band%d' % (sideband, i))
            print 'Chunk s%02d FITS-IDI file is present, will process' % (i)
            chunksToProcess.append(i)
            testFile.close()
        except IOError:
            print 'Chunk s%02d FITS-IDI file is missing, will skip' % (i)
    # Read the individual SMA chunks in
    for i in chunksToProcess:
        fitsidifile = './tempFITS-IDI%s.band%d' % (sideband, i)
        mSFileName = '%s_s%02d' % (sideband, i)
        vis = mSFileName
        print 'Reading in chunk s%02d (%s)' % (i, sideband)
        importfitsidi()
        # Flag the bad data points
        flagbackup = False
        mode       ='clip'
        clipzeros  = True
        flagdata()
    # Get the info about with pad each antenna is on
        tb.open(vis+'/ANTENNA', nomodify=False)
        stations = tb.getcol('STATION')
        hduList = fits.open(fitsidifile)
        arhd=hduList[1]
        cards = arhd.header
        for ii in range(len(cards)):
            try:
                if (' was on pad ' in cards[ii]) and not ('Not used' in cards[ii]) and not ('CSO' in cards[ii]) and not ('JCMT' in cards[ii]):
                    tok = cards[ii].split()
                    ant = int(tok[1])
                    pad = tok[5]+' '+tok[6]
                    stations[ant-1] = pad
            except TypeError:
                continue
        tb.putcol('STATION', stations)
        tb.unlock()
        tb.close()

    # Make SYSCAL table with Tsys
        print 'Making a SYSCAL table for '+mSFileName
        tb.fromfits(tablename=mSFileName+'/SYSCAL',fitsfile=fitsidifile,whichhdu=5,nomodify=False)
        tb.clearlocks()
        tb.close()
        os.system('rm '+mSFileName+'/SYSCAL/table.lock')

    #Delete the FITS-IDI files
    os.system('rm ./tempFITS-IDI%s.band*' % sideband)

    if useSMAScanNumbers:
        
        # make scan numbers from unique times
        for i in chunksToProcess:
            tb.open('%s_s%02d' % (sideband, i), nomodify=False)
            times = tb.getcol('TIME')
            timeSet = set(times)
            timeList = list(timeSet)
            timeList.sort()
            print 'Generating scan numbers from times'
            oldScanNo=tb.getcol('SCAN_NUMBER')
            newScanNo = oldScanNo.copy()
            for i in range(len(times)):
                newScanNo[i] = timeList.index(times[i]) + 1
            tb.putcol('SCAN_NUMBER' ,newScanNo)
            tb.unlock()
            tb.close()
    else:
        
        # make scan numbers from source changes
        for i in chunksToProcess:
            scanNumber = 0
            lastField = -1
            tb.open('%s_s%02d' % (sideband, i), nomodify=False)
            fields = tb.getcol('FIELD_ID')
            print 'Generating scan numbers from source changes'
            oldScanNo=tb.getcol('SCAN_NUMBER')
            newScanNo = oldScanNo.copy()
            for i in range(len(newScanNo)):
                if fields[i] != lastField:
                    scanNumber += 1
                    lastField = fields[i]
                newScanNo[i] = scanNumber
            tb.putcol('SCAN_NUMBER' ,newScanNo)
            tb.unlock()
            tb.close()
            print 'Last scan number = ', scanNumber

    # Fix up the Weights (channel 0 doesn't need fixing)
    for i in chunksToProcess:
        if i != 0:
            print 'Fixing weights for chunk s%02d' % (i)
            tb.open('%s_s%02d' % (sideband, i), nomodify=False)
            weights=tb.getcol('WEIGHT_SPECTRUM')
            oneWeightSet=weights[0][0]
            oldWeights=tb.getcol('WEIGHT')
            newWeights=oldWeights.copy()
            for j in range(len(oldWeights[0])):
                newWeights[0][j]= oneWeightSet[j]
            tb.putcol('WEIGHT',newWeights)
            tb.unlock()
            tb.close()

    vis = []
    for i in chunksToProcess:
        vis.append('%s_s%02d' % (sideband, i))
    concatvis = 'MyData%s' % (sideband)
    #Concatonate all the spectral chunks into one MS
    print 'Concatenating the individual chunk MSs into a single MS'
    concat()

    if len(sourceList) > 0:
        print 'Fixing some FIELD columns'
        tb.open(concatvis+'/FIELD',nomodify=False)
        names=tb.getcol('NAME')
        sourceId=tb.getcol('SOURCE_ID')
        newSourceIds=sourceId.copy()
        for i in range(len(names)):
            newSourceIds[i] = sourceList.index(names[i])
        tb.putcol('SOURCE_ID',newSourceIds)
        codes=tb.getcol('CODE')
        newCodes = []
        for i in range(len(codes)):
            newCodes.append('none')
        insertArray = np.array(newCodes, dtype='|S19')
        tb.putcol('CODE', insertArray)
        tb.unlock()
        tb.close()
        print 'Building the SOURCE table'
        tb.create(concatvis+'/SOURCE',tabledesc=sourceTableDesc,dminfo=sourceDMInfo,nrow=len(sourceList)*len(bandList))
        iArray = []
        for i in uniqueSourceNumbers:
            for j in xrange(len(bandList)):
                iArray.append(i)
        insertArray = np.array(iArray,dtype=np.int32)
        tb.putcol('SOURCE_ID', insertArray)
        iArray = []
        for i in uniqueSourceNames:
            for j in xrange(len(bandList)):
                iArray.append(i)
        insertArray = np.array(iArray,dtype='|S19')
        tb.putcol('NAME', insertArray)
        iArray = []
        for i in uniqueSourceNames:
            for j in xrange(len(bandList)):
                iArray.append('none')
        insertArray = np.array(iArray,dtype='|S19')
        tb.putcol('CODE', insertArray)
        iArray = []
        for i in uniqueSourceNumbers:
            for j in xrange(len(bandList)):
                iArray.append(4.38209253e+09)
        insertArray = np.array(iArray, dtype='d')
        tb.putcol('INTERVAL', insertArray)
        iArray = []
        for i in uniqueSourceNumbers:
            for j in xrange(len(bandList)):
                iArray.append(7.03232577e+09)
        insertArray = np.array(iArray, dtype='d')
        tb.putcol('TIME', insertArray)
        iArray = []
        for i in uniqueSourceNumbers:
            for j in bandList:
                iArray.append(j)
        insertArray = np.array(iArray,dtype=np.int32)
        tb.putcol('SPECTRAL_WINDOW_ID', insertArray)
        rAs=[]
        for i in uniqueSourceRAs:
            for j in bandList:
                rAs.append(float(i))
        decs=[]
        for i in uniqueSourceDecs:
            for j in bandList:
                decs.append(float(i))
        insertArray=np.array([rAs,decs],dtype='d')
        tb.putcol('DIRECTION',insertArray)
        tb.unlock()
        tb.close()

    # Get rid of the single-chunk MSs
    for i in chunksToProcess:
        os.system('rm -r -f %s_s%02d' % (sideband, i))
