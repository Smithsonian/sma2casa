import casac
from tasks import *
from taskinit import *
import pylab as pl
import os
from math import pi,floor
import pyfits as fits

useSMAScanNumbers = False

default(importfitsidi)
default(tflagdata)
sidebandsToProcess = ['Lower', 'Upper']
# Parameters for flagdata:

for sideband in sidebandsToProcess:
    # First figure out which chunks need to be processed
    chunksToProcess = []
    for i in range(51):
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
        vis = '%s_s%02d' % (sideband, i)
        print 'Reading in chunk s%02d (%s)' % (i, sideband)
        importfitsidi()
        # Flag the bad data points
        flagbackup = False
        mode       ='clip'
        clipzeros  = True
        tflagdata()
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

    # Get rid of the single-chunk MSs
    for i in chunksToProcess:
        os.system('rm -r -f %s_s%02d' % (sideband, i))
