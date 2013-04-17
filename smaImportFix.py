import casac
from tasks import *
from taskinit import *
import pylab as pl
import os
from math import pi,floor

# First figure out which chunks need to be processed
chunksToProcess = []
for i in range(49):
    try:
        testFile = open('./tempFITS-IDIUpper.band%d' % (i))
        print 'Chunk s%02d FITS-IDI file is present, will process' % (i)
        chunksToProcess.append(i)
    except IOError:
        print 'Chunk s%02d FITS-IDI file is missing, will skip' % (i)

# Parameters for flagdata:
flagbackup = False
mode       ='clip'
clipzeros  = True

# Read the individual SMA chunks in
for i in chunksToProcess:
    fitsidifile = './tempFITS-IDIUpper.band%d' % (i)
    vis = 'upper_s%02d' % (i)
    print 'Reading in chunk s%02d USB' % (i)
    importfitsidi()
    # Flag the bad data points
    flagdata()

#Delete the FITS-IDI files
os.system('rm ./tempFITS-IDIUpper.band*')

# Fix up the Weights (channel 0 doesn't need fixing)
for i in chunksToProcess:
    if i != 0:
        print 'Fixing weights for chunk s%02d' % (i)
        tb.open('upper_s%02d' % (i), nomodify=False)
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
    vis.append('upper_s%02d' % (i))
concatvis = 'MyDataUSB'
#Concatonate all the spectral chunks into one MS
print 'Concatenating the individual chunk MSs into a single MS'
concat()

# Get rid of the single-chunk MSs
for i in chunksToProcess:
    os.system('rm -r -f upper_s%02d' % (i))
