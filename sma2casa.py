#!/sma/SMAusers/taco/myPythonEnv/bin/python
"""
sma2casa.py
===========

This file defines functions which will import SMA data directly into a CASA Measurement Set.
"""

import numpy as np
import makevis
import os, sys, struct, mmap, pyfits, getopt
from astropy.time import Time
from math import sin, cos, pi, sqrt
#import psyco
#psyco.full()

neededFiles = ('antennas', 'bl_read', 'codes_read', 'in_read', 'sch_read', 'sp_read', 'tsys_read', 'projectInfo')

nBands = 0
bandList = []
spSmallDictL = {}
spSmallDictU = {}
spBigDict = {}
antennas = {}
codesDict = {}
inDict = {}
blDictL = {}
blTsysDictL = {}
blDictU = {}
blTsysDictU = {}
sourceDict = {}
maxScan = 10000000
maxWeight = 0.01
numberOfBaselines = 0
antennaList = []
swappingTsys = False
tsysMapping = range(0,11)
trimEdges = False
edgeTrimFraction = 0.1 # Fraction on each edge of a spectral chunk to flag bad
chunkList = range(49)
sidebandList = [0, 1]
pseudoContinuumFrequency = {}
verbose = True
NaN = float('nan')
newFormat = True

def usage():
    print 'Usage: ', sys.argv[0], ' path-to-SMA-data [options]'
    print 'Options: (which must come *after* the path argument)'
    print '\t -c m,n,o [--chunks=m,n,o]\tProcess only chunks m, n and o'
    print '\t -h [--help]\t\t\tPrint this message and exit'
    print '\t -l [--lower]\t\t\tProcess lower sideband only'
    print '\t -p [--percent]\t\t\t%% to trim on band edge (default = %2.0f)' % (edgeTrimFraction*100.0)
    print '\t -s [--silent]\t\t\tRun silently unless an error occurs'
    print '\t -t [--trim]\t\t\tSet the amplitude at chunk edges to 0.0)'
    print "\t -T\t\t\t\t-T n=m means use ant m's Tsys for ant n"
    print '\t -u [--upper]\t\t\tProcess upper sideband only'

def normalize0to360(angle):
    while angle >= 360.0:
        angle -= 360.0
    while angle < 0.0:
        angle += 360.0
    return angle

def toGeocentric(X, Y, Z):
    """
    Transform antenna coordinates to FITS-IDI style GEOCENTRIC coordinates
    """
    SMALong = 2.713594675620429
    SMALat  = 0.3459976585365961
    SMARad  = 6382.248*1000.0
    Rx = SMARad*cos(SMALong)*cos(SMALat)
    Ry = SMARad*sin(SMALong)*cos(SMALat)
    Rz = SMARad*sin(SMALat)
    x = X*cos(SMALong) + Y*sin(SMALong) + Rx
    y = Y*cos(SMALong) - X*sin(SMALong) + Ry
    z = Z + Rz
    return (x, y, z)

def makeInt(data, size):
    tInt = 0
    if newFormat:
        for i in range(size):
            tInt += ord(data[i])<<(i<<3)
    else:
        for i in range(size):
            tInt += ord(data[size-i-1])<<(i<<3)
    return tInt

def makeFloat(data):
    if not newFormat:
        hostByteOrderData = ''
        for i in range(4):
            hostByteOrderData += chr(ord(data[3-i]))
        return (struct.unpack('f', hostByteOrderData[:4]))[0]
    else:
        return (struct.unpack('f', data[:4]))[0]

def makeDouble(data):
    if not newFormat:
        hostByteOrderData = ''
        for i in range(8):
            hostByteOrderData += chr(ord(data[7-i]))
        return (struct.unpack('d', hostByteOrderData[:8]))[0]
    else:
        return (struct.unpack('d', data[:8]))[0]

def read(dataDir):
    global nBands, bandList, antennas, codesDict, inDict, blDictL, blDictU
    global spSmallDictL, spSmallDictU, spBigDict, sourceDict, maxWeight, numberOfBaselines
    global antennaList, blTsysDictL, blTsysDictU, newFormat, pseudoContinuumFrequency

    # Check that the directory contains all the required files
    if verbose:
        print 'checking that the needed files exist in ', dataDir
    dirContents = os.listdir(dataDir)
    for needed in neededFiles:
        if not needed in dirContents:
            print "The directory %s does not have a %s file - aborting" % (dataDir, needed)
            sys.exit(-1)

    ###
    ### Determine if this is an old-format or new-format data file
    ###
    f = open(dataDir+'/in_read', 'rb')
    data = f.read()
    firstInt = makeInt(data[0:], 4)
    if firstInt != 0:
        newFormat = True
    else:
        newFormat = False
    f.close()
    if verbose:
        if newFormat:
            print 'This is a new format data set.'
        else:
            print 'This is an old format data set'

    ###
    ### Read in the antennas file - make a dictionary antennas[ant #] = (X, Y, Z)
    ###
    if verbose:
        print 'Reading antennas file'
    for line in open(dataDir+'/antennas'):
        tok = line.split()
        X = float(tok[1])
        Y = float(tok[2])
        Z = float(tok[3])
        antennas[int(tok[0])] = toGeocentric(X, Y, Z)

    ###
    ###  Read in the codes file (codes_read)
    ###  Build a disctionary, codesDict, which has an entry for each unique type of code
    ###  The entries in that disctionary are themselves dictionaries, with integer keys
    ###  pointing to the various entries for that type of code
    ###
    if verbose:
        print 'Reading codes_read'
    f = open(dataDir+'/codes_read', 'rb')
    data = f.read()
    codesRecLen = 42
    nCodesRecords = len(data)/codesRecLen
    # First make a list with all unique code strings
    codeStrings = []
    for rec in range(nCodesRecords):
        codeString = ''
        for i in range(12):
            if data[codesRecLen*rec + i] < ' ':
                break
            if (data[codesRecLen*rec + i] >= ' ') and (data[codesRecLen*rec + i] <= 'z'):
                codeString += data[codesRecLen*rec + i]
        icode = makeInt(data[codesRecLen*rec + 12:], 2)
        payload = ''
        for i in range(14, 14+26):
            if data[codesRecLen*rec + i] < ' ':
                break
            if (data[codesRecLen*rec + i] >= ' ') and (data[codesRecLen*rec + i] <= 'z'):
                payload += data[codesRecLen*rec + i]
        if codeString not in codeStrings:
            codeStrings.append(codeString)
            codesDict[codeString] = {}
        codesDict[codeString][icode] = payload

    ###
    ### Read in the integration header information (in_read)
    ### Produces a disctionary with integration numbers for the
    ### keys, and tuples holding the entry values
    ###
    if verbose:
        print 'Reading in_read'
    f = open(dataDir+'/in_read', 'rb')
    data = f.read()
    savedSourceInfo = True
    if newFormat:
        inRecLen = 188
    else:
        inRecLen = 132
    nInRecords = len(data)/inRecLen
    for rec in range(nInRecords):
        if newFormat:
            traid     =    makeInt(data[rec*inRecLen +   0:], 4)
            inhid     =    makeInt(data[rec*inRecLen +   4:], 4)
            az        =  makeFloat(data[rec*inRecLen +  12:])
            el        =  makeFloat(data[rec*inRecLen +  16:])
            hA        =  makeFloat(data[rec*inRecLen +  20:])
            iut       =    makeInt(data[rec*inRecLen +  24:], 2)
            iref_time =    makeInt(data[rec*inRecLen +  26:], 2)
            dhrs      = makeDouble(data[rec*inRecLen +  28:])
            vc        =  makeFloat(data[rec*inRecLen +  36:])
            sx        = makeDouble(data[rec*inRecLen +  40:])
            sy        = makeDouble(data[rec*inRecLen +  48:])
            sz        = makeDouble(data[rec*inRecLen +  56:])
            rinteg    =  makeFloat(data[rec*inRecLen +  64:])
            proid     =    makeInt(data[rec*inRecLen +  68:], 4)
            souid     =    makeInt(data[rec*inRecLen +  72:], 4)
            isource   =    makeInt(data[rec*inRecLen +  76:], 2)
            ivrad     =    makeInt(data[rec*inRecLen +  78:], 2)
            offx      =  makeFloat(data[rec*inRecLen +  80:])
            offy      =  makeFloat(data[rec*inRecLen +  84:])
            ira       =    makeInt(data[rec*inRecLen +  88:], 2)
            idec      =    makeInt(data[rec*inRecLen +  90:], 2)
            rar       = makeDouble(data[rec*inRecLen +  92:])
            decr      = makeDouble(data[rec*inRecLen + 100:])
            epoch     =  makeFloat(data[rec*inRecLen + 108:])
            size      =  makeFloat(data[rec*inRecLen + 112:])
        else: # Old format file
            traid     =    makeInt(data[rec*inRecLen +   6:], 4)
            inhid     =    makeInt(data[rec*inRecLen +  10:], 4)
            az        =  makeFloat(data[rec*inRecLen +  20:])
            el        =  makeFloat(data[rec*inRecLen +  24:])
            hA        =  makeFloat(data[rec*inRecLen +  28:])
            iut       =    makeInt(data[rec*inRecLen +  32:], 2)
            iref_time =    makeInt(data[rec*inRecLen +  34:], 2)
            dhrs      = makeDouble(data[rec*inRecLen +  36:])
            vc        =  makeFloat(data[rec*inRecLen +  44:])
            sx        = makeDouble(data[rec*inRecLen +  50:])
            sy        = makeDouble(data[rec*inRecLen +  58:])
            sz        = makeDouble(data[rec*inRecLen +  66:])
            rinteg    =  makeFloat(data[rec*inRecLen +  74:])
            proid     =    makeInt(data[rec*inRecLen +  78:], 4)
            souid     =    makeInt(data[rec*inRecLen +  82:], 4)
            isource   =    makeInt(data[rec*inRecLen +  86:], 2)
            ivrad     =    makeInt(data[rec*inRecLen +  88:], 2)
            offx      =  makeFloat(data[rec*inRecLen +  90:])
            offy      =  makeFloat(data[rec*inRecLen +  94:])
            ira       =    makeInt(data[rec*inRecLen + 100:], 2)
            idec      =    makeInt(data[rec*inRecLen + 102:], 2)
            rar       = makeDouble(data[rec*inRecLen + 104:])
            decr      = makeDouble(data[rec*inRecLen + 112:])
            epoch     =  makeFloat(data[rec*inRecLen + 120:])
            size      =  makeFloat(data[rec*inRecLen + 128:])
        inDict[inhid] = (traid, inhid, az, el, hA, iut, iref_time, dhrs, vc, sx, sy, sz,
                         rinteg, proid, souid, isource, ivrad, offx, offy, ira, idec, rar, decr, epoch, size)
        if souid not in sourceDict:
            # This funny business with savedSourceInfo is a kludge - it keeps the source info from being
            # taken from the first integration header, which usually still has the RA and Dec from the
            # previous scan (and is therefore the header for a bad scan.   This kludge makes the code store
            # the source info from the 2nd integration header, which should have the correct coordinates.
            # but this is not guaranteed, so it should be fixed up later
            if savedSourceInfo:
                savedSourceInfo = False
            else:
                sourceDict[souid] = (codesDict['source'][isource], rar, decr)
                savedSourceInfo = True
        else:
            savedSourceInfo = True

    ###
    ### Read in bl_read
    ###
    if verbose:
        print 'Reading bl_read'
    f = open(dataDir+'/bl_read', 'rb')
    fSize = os.path.getsize(dataDir+'/bl_read')
    if newFormat:
        blRecLen = 158
    else:
        blRecLen = 118
    nBlRecords = fSize/blRecLen
    baselineList = []
    blSidebandDict = {}
    for rec in range(nBlRecords):
        if ((rec % 10000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nBlRecords, 100.0*float(rec)/float(nBlRecords))
            sys.stdout.flush()
        data = f.read(blRecLen)
        if len(data) == blRecLen:
            if newFormat:
                blhid     =    makeInt(data[  0:], 4)
                inhid     =    makeInt(data[  4:], 4)
                isb       =    makeInt(data[  8:], 2)
                ipol      =    makeInt(data[ 10:], 2)
                ant1rx    =    makeInt(data[ 12:], 2)
                ant2rx    =    makeInt(data[ 14:], 2)
                pointing  =    makeInt(data[ 16:], 2)
                irec      =    makeInt(data[ 18:], 2)
                u         =  makeFloat(data[ 20:])
                v         =  makeFloat(data[ 24:])
                w         =  makeFloat(data[ 28:])
                prbl      =  makeFloat(data[ 32:])
                coh       =  makeFloat(data[ 36:])
                avedhrs   = makeDouble(data[ 40:])
                ampave    =  makeFloat(data[ 48:])
                phaave    =  makeFloat(data[ 52:])
                blsid     =    makeInt(data[ 56:], 4)
                ant1      =    makeInt(data[ 60:], 2)
                ant2      =    makeInt(data[ 62:], 2)
                ant1TsysOff =  makeInt(data[ 64:], 4)
                ant2TsysOff =  makeInt(data[ 68:], 4)
            else: # old format
                blhid     =    makeInt(data[  0:], 4)
                inhid     =    makeInt(data[  4:], 4)
                isb       =    makeInt(data[  8:], 2)
                ipol      =    makeInt(data[ 10:], 2)
                ant1rx    =    makeInt(data[ 16:], 2)
                ant2rx    =    makeInt(data[ 18:], 2)
                pointing  =    makeInt(data[ 22:], 2)
                irec      =    makeInt(data[ 24:], 2)
                u         =  makeFloat(data[ 28:])
                v         =  makeFloat(data[ 32:])
                w         =  makeFloat(data[ 36:])
                prbl      =  makeFloat(data[ 40:])
                coh       =  makeFloat(data[ 52:])
                avedhrs   = makeDouble(data[ 72:])
                ampave    =  makeFloat(data[ 80:])
                phaave    =  makeFloat(data[ 84:])
                blsid     =    makeInt(data[ 92:], 4)
                ant1      =    makeInt(data[ 96:], 2)
                ant2      =    makeInt(data[ 98:], 2)
                ant1TsysOff =  0
                ant2TsysOff =  0
            if (ant1, ant2) not in baselineList:
                baselineList.append((ant1, ant2))
                numberOfBaselines +=1
            if ant1 not in antennaList:
                antennaList.append(ant1)
            if ant2 not in antennaList:
                antennaList.append(ant2)
            blSidebandDict[blhid] = isb
            if isb == 0:
                blDictL[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2, ant1TsysOff, ant2TsysOff)
            else:
                blDictU[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2, ant1TsysOff, ant2TsysOff)
    nAnts = len(antennaList)
    if verbose:
        print '%d antennas and %d baselines seen in this data set.' % (nAnts, numberOfBaselines)
    if numberOfBaselines != (nAnts*(nAnts-1)/2):
        print 'Number of baselines seen inconsistant with the number of antennas - aborting'
        sys.exit(-1)

    ###
    ### Pull the Tsys info out of tsys_read
    ###
    if newFormat:
        tsysFile = os.open(dataDir+'/tsys_read', os.O_RDONLY)
        tsysMap = mmap.mmap(tsysFile, 0, prot=mmap.PROT_READ);
        for bl in blDictL:
            ant1Tsys4to6 = makeFloat(tsysMap[blDictL[bl][17]+12:blDictL[bl][17]+16])
            ant1Tsys6to8 = makeFloat(tsysMap[blDictL[bl][17]+28:blDictL[bl][17]+32])
            ant2Tsys4to6 = makeFloat(tsysMap[blDictL[bl][18]+12:blDictL[bl][18]+16])
            ant2Tsys6to8 = makeFloat(tsysMap[blDictL[bl][18]+28:blDictL[bl][18]+32])
            blTsysDictL[bl] = (ant1Tsys4to6, ant1Tsys6to8, ant2Tsys4to6, ant2Tsys6to8)
        for bl in blDictU:
            ant1Tsys4to6 = makeFloat(tsysMap[blDictU[bl][17]+16:blDictU[bl][17]+20])
            ant1Tsys6to8 = makeFloat(tsysMap[blDictU[bl][17]+32:blDictU[bl][17]+36])
            ant2Tsys4to6 = makeFloat(tsysMap[blDictU[bl][18]+16:blDictU[bl][18]+20])
            ant2Tsys6to8 = makeFloat(tsysMap[blDictU[bl][18]+32:blDictU[bl][18]+36])
            blTsysDictU[bl] = (ant1Tsys4to6, ant1Tsys6to8, ant2Tsys4to6, ant2Tsys6to8)
    else: # Old format file
        if verbose:
            print 'Processing old-format Tsys data'
        antLowDict  = {}
        antHighDict = {}
        f = open(dataDir+'/tsys_read', 'rb')
        fSize = os.path.getsize(dataDir+'/tsys_read')
        tsysRecLen = 400
        nTsysRecords = fSize/tsysRecLen
        for rec in range(nTsysRecords):
            if ((rec % 10000) == 0) and verbose:
                print '\t processing record %d of %d (%2.0f%% done)' % (rec, nBlRecords, 100.0*float(rec)/float(nBlRecords))
                sys.stdout.flush()
            data = f.read(tsysRecLen)
            if len(data) == tsysRecLen:
                inhid    =   makeInt(data[  0:], 4)
                iAnt     =   makeInt(data[ 14:], 2)
                tsysLow  = makeFloat(data[ 16:])
                tsysHigh = makeFloat(data[208:])
                if (inhid, iAnt) not in antLowDict:
                    antLowDict[(inhid, iAnt)]  = tsysLow
                    antHighDict[(inhid, iAnt)] = tsysHigh
        for bl in blDictL:
            inhid = blDictL[bl][0]
            ant1  = blDictL[bl][15]
            ant2  = blDictL[bl][16]
            blTsysDictL[bl] = (antLowDict[(inhid, ant1)], antHighDict[(inhid, ant1)], antLowDict[(inhid, ant2)], antHighDict[(inhid, ant2)])
        for bl in blDictU:        
            inhid = blDictU[bl][0]
            ant1  = blDictU[bl][15]
            ant2  = blDictU[bl][16]
            blTsysDictU[bl] = (antLowDict[(inhid, ant1)], antHighDict[(inhid, ant1)], antLowDict[(inhid, ant2)], antHighDict[(inhid, ant2)])
    if swappingTsys:
        print 'Swapping Tsys values'
        # Make a mapping dictionary that maps (intergration, antenna) to Tsys
        mapping = {}
        for bl in blDictL:
            mapping[(blDictL[bl][0], blDictL[bl][15])] = (blTsysDictL[bl][0], blTsysDictL[bl][1])
            mapping[(blDictL[bl][0], blDictL[bl][16])] = (blTsysDictL[bl][2], blTsysDictL[bl][3])
        # Now use that mapping table to perform the Tsys substitution
        for bl in blDictL:
            blTsysDictL[bl] = (mapping[(blDictL[bl][0], tsysMapping[blDictL[bl][15]])][0],mapping[(blDictL[bl][0], tsysMapping[blDictL[bl][15]])][1], mapping[(blDictL[bl][0], tsysMapping[blDictL[bl][16]])][0], mapping[(blDictL[bl][0], tsysMapping[blDictL[bl][16]])][1])
        mapping = {}
        for bl in blDictU:
            mapping[(blDictU[bl][0], blDictU[bl][15])] = (blTsysDictU[bl][0], blTsysDictU[bl][1])
            mapping[(blDictU[bl][0], blDictU[bl][16])] = (blTsysDictU[bl][2], blTsysDictU[bl][3])
        # Now use that mapping table to perform the Tsys substitution
        for bl in blDictU:
            blTsysDictU[bl] = (mapping[(blDictU[bl][0], tsysMapping[blDictU[bl][15]])][0],mapping[(blDictU[bl][0], tsysMapping[blDictU[bl][15]])][1], mapping[(blDictU[bl][0], tsysMapping[blDictU[bl][16]])][0], mapping[(blDictU[bl][0], tsysMapping[blDictU[bl][16]])][1])

    ###
    ### Count the number of spectral bands
    ###
    if verbose:
        print 'Counting the number of bands'
    f = open(dataDir+'/sp_read', 'rb')
    fSize = os.path.getsize(dataDir+'/sp_read')
    if newFormat:
        spRecLen = 188
    else:
        spRecLen = 100
    data = f.read(spRecLen)
    nSpRecords = fSize/spRecLen
    for rec in range(nSpRecords):
        if ((rec % 100000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nSpRecords, 100.0*float(rec)/float(nSpRecords))
            sys.stdout.flush()
        data = f.read(spRecLen)
        if len(data) == spRecLen:
            if newFormat:
                sphid     =    makeInt(data[  0:], 4)
                blhid     =    makeInt(data[  4:], 4)
                inhid     =    makeInt(data[  8:], 4)
                igq       =    makeInt(data[ 12:], 2)
                ipq       =    makeInt(data[ 14:], 2)
                iband     =    makeInt(data[ 16:], 2)
                fsky      = makeDouble(data[ 36:])
                fres      =  makeFloat(data[ 44:])
                wt        =  makeFloat(data[ 84:])
                flags     =    makeInt(data[ 88:], 4)
                nch       =    makeInt(data[ 96:], 2)
                dataoff   =    makeInt(data[100:], 4)
                rfreq     = makeDouble(data[104:])
            else:
                sphid     =    makeInt(data[  0:], 4)
                blhid     =    makeInt(data[  4:], 4)
                inhid     =    makeInt(data[  8:], 4)
                igq       =    makeInt(data[ 12:], 2)
                ipq       =    makeInt(data[ 14:], 2)
                iband     =    makeInt(data[ 16:], 2)
                fsky      = makeDouble(data[ 38:])
                fres      =  makeFloat(data[ 46:])
                wt        =  makeFloat(data[ 58:])
                if wt > 0:
                    flags = 0
                else:
                    flags = -1
                nch       =    makeInt(data[ 68:], 2)
                dataoff   =    makeInt(data[ 72:], 4)
                rfreq     = makeDouble(data[ 82:])
            if (flags != 0) and (wt > 0):
                wt *= -1.0
            spBigDict[(iband, blhid)] = (dataoff, wt)
            if iband not in bandList:
                bandList.append(iband)
            try:
                if (nch == 1) and (blSidebandDict[blhid] not in pseudoContinuumFrequency):
                    pseudoContinuumFrequency[blSidebandDict[blhid]] = fsky*1.0e9
                if blSidebandDict[blhid] == 0:
                    wt = wt/(max(blTsysDictL[blhid][0], blTsysDictL[blhid][1])*max(blTsysDictL[blhid][2], blTsysDictL[blhid][3]))
                    if iband not in spSmallDictL:
                        spSmallDictL[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
                else:
                    wt = wt/(max(blTsysDictU[blhid][0], blTsysDictU[blhid][1])*max(blTsysDictU[blhid][2], blTsysDictU[blhid][3]))
                    if iband not in spSmallDictU:
                        spSmallDictU[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
            except KeyError:
                # This exception can occur if the last record to bl_read was not written, because dataCatcher was interrupted
                wt = -1.0
            if abs(wt) > maxWeight:
                maxWeight = abs(wt)
        if inhid > maxScan:
            break
    nBands = len(bandList)
    bandList.sort()

if len(sys.argv) < 2:
    usage()
    exit(-1)
dataSet = sys.argv[1]
try:
    opts, args = getopt.getopt(sys.argv[2:], "c:hlp:stT:u", ['chunks=', 'help', 'lower', 'percent=', 'silent', 'trim', 'Tsys=', 'upper'])
except getopt.GetoptError as err:
    usage()
    sys.exit(-1)
for o, a in opts:
    if o in ('-c', '--chunks'):
        sChunkList = a.split(',')
        chunkList = []
        for chunk in sChunkList:
            chunkList.append(int(chunk))
    elif o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-l', '--lower'):
        sidebandList = [0]
    elif o in ('-p', '--percent'):
        edgeTrimFraction = float(a)/100.0
        trimEdges = True
    elif o in ('-s', '--silent'):
        verbose = False
    elif o in ('-t', '--trim'):
        trimEdges = True
    elif o in ('-T', '--Tsys'):
        donee = int(a[0])
        doner = int(a[2])
        tsysMapping[donee] = doner
        swappingTsys = True
    elif o in ('-u', '--upper'):
        sidebandList = [1]
    else:
        print o, a
        assert False, 'unhandled option'
if verbose:
    for a in range(1,9):
        if a != tsysMapping[a]:
            print 'Antenna %d will use the Tsys from antenna %d' % (a, tsysMapping[a])

read(dataSet)
for line in open(dataSet+'/projectInfo'):
    projectPI = line.strip()
visFileLen = makevis.open(dataSet+'/sch_read')
for band in bandList:
    for sb in range(2):
        if (band in chunkList) and (sb in sidebandList):
            if sb == 0:
                sideBand = 'Lower'
                blDict = blDictL
                blTsysDict = blTsysDictL
                spSmallDict = spSmallDictL
            else:
                sideBand = 'Upper'
                blDict = blDictU
                blTsysDict = blTsysDictU
                spSmallDict = spSmallDictU
            lowestFSky = spSmallDict[band][2] - 52.0e6
        
            ###
            ### Make the Primary HDU
            ###
            hdu = pyfits.PrimaryHDU()
            hdulist=pyfits.HDUList([hdu])
            header = hdulist[0].header
            header.update('groups', True)
            header.update('gcount', 0)
            header.update('pcount', 0)
            if (band == 0) and ((49 in bandList) or (50 in bandList)):
                header.update('correlat', 'SMA-Hybrid')
            elif band < 49:
                header.update('correlat', 'SMA-Legacy')
            else:
                header.update('correlat', 'SMA-SWARM')
            header.update('fxcorver', 1)
            header.add_history('Translated to FITS-IDI by sma2casa.py')
            header.add_history('Original SMA dataset: '+dataSet[-15:])
#            if newFormat:
#                header.add_history('The SMA data set used the new SWARM format')
#            else:
#                header.add_history('The SMA data set used the old format')
            header.add_history('Project PI: '+projectPI)

            ###
            ### Make the ARRAY_GEOMETRY table
            ###
            refTimeString = codesDict['ref_time'][0][-4:]
            mon = codesDict['ref_time'][0][:3]
            if mon == 'Jan':
                refTimeString += '-01-'
            elif mon == 'Feb':
                refTimeString += '-02-'
            elif mon == 'Mar':
                refTimeString += '-03-'
            elif mon == 'Apr':
                refTimeString += '-04-'
            elif mon == 'May':
                refTimeString += '-05-'
            elif mon == 'Jun':
                refTimeString += '-06-'
            elif mon == 'Jul':
                refTimeString += '-07-'
            elif mon == 'Aug':
                refTimeString += '-08-'
            elif mon == 'Sep':
                refTimeString += '-09-'
            elif mon == 'Oct':
                refTimeString += '-10-'
            elif mon == 'Nov':
                refTimeString += '-11-'
            else:
                refTimeString += '-12-'
            day = int(codesDict['ref_time'][0][5:6])
            refTimeString = '%s%02d' % (refTimeString, day)
            jD = Time(refTimeString, scale='utc').jd
            T = (jD - 2451545.0)/36525.0
            gST = 100.46061837 + 36000.770053608*T + 0.000387933*T**2 - (T**3)/38710000.0
            gST = normalize0to360(gST)
            antNames = np.array(['SMA1    ','SMA2    ','SMA3    ','SMA4    ','SMA5    ','SMA6    ','SMA7    ','SMA8    '])
            antList = []
            velList = []
            antNumList = []
            antStyleList = []
            antOffList = []
            antDiameterList = []
            for ant in range(1, 9):
                antList.append(antennas[ant])
                velList.append((0.0, 0.0, 0.0))
                antNumList.append(ant)
                antStyleList.append(0)
                antOffList.append((0.0, 0.0, 0.0))
                antDiameterList.append(6.0)
            c1 = pyfits.Column(name='ANNAME',   format='8A', array=antNames                        )
            c2 = pyfits.Column(name='STABXYZ',  format='3D', array=antList,         unit='METERS'  )
            c3 = pyfits.Column(name='DERXYZ',   format='3E', array=velList,         unit='METERS/S')
            c4 = pyfits.Column(name='NOSTA',    format='1J', array=antNumList                      )
            c5 = pyfits.Column(name='MNTSTA',   format='1J', array=antStyleList                    )
            c6 = pyfits.Column(name='STAXOF',   format='3E', array=antOffList                      )
            c7 = pyfits.Column(name='DIAMETER', format='1E', array=antDiameterList, unit='METERS'  )
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
            arrayGeometryHDU = pyfits.new_table(coldefs)
            header = arrayGeometryHDU.header
            header.update('extname',  'ARRAY_GEOMETRY')
            header.update('frame',    'GEOCENTRIC')
            header.update('arrnam',   'SMA')
            header.update('tabrev',   1)
            header.update('obscode',  ' ');
            header.update('extver',   1,                         'Array number')
            header.update('arrayx',   0.0,                       'Array origin x coordinate')
            header.update('arrayy',   0.0,                       'Array origin y coordinate')
            header.update('arrayz',   0.0,                       'Array origin z coordinate')
            header.update('numorb',   0,                         'Number of orbital parameters')
            header.update('polarx',   0.0,                       'x coordinate of north pole (arcsec)')
            header.update('polary',   0.0,                       'y coordinate of north pole (arcsec)')
            header.update('freq',     lowestFSky,                'Reference frequency')
            header.update('timesys',  'UTC',                     'Time system')
            header.update('degpdy' ,  0.0,                       'Earth rotation rate')
            header.update('ut1utc' ,  0.0,                       'UT1 - UTC (seconds) - PLACEHOLDER')
            header.update('iatutc' ,  0.0,                       'IAT - UTC (seconds) - PLACEHOLDER')
            header.update('rdate',    refTimeString,             'Reference date')
            header.update('gstia0',   gST,                       'GST at 0h on reference date (degrees)')
            header.update('no_stkd',  1,                         'Number of Stokes parameters')
            header.update('stk_1',    -2,                        'First Stokes parameter')
            header.update('no_band',  1,                         'Number of bands')
            header.update('ref_freq', lowestFSky,                'File reference frequency')
            header.update('ref_pixl', 1,                         'Reference pixel')
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels')
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth')

            ###
            ### Make the SOURCE table
            ###
            iDList = []
            nameList = []
            qualList = []
            calList = []
            freqList = []
            fluxList = []
            rAList = []
            decList = []
            eqList = []
            velList = []
            vtList = []
            vdList = []
            restList = []
            for source in range(1,len(sourceDict)+1):
                iDList.append(source)
                nameList.append(sourceDict[source][0])
                qualList.append(0)
                calList.append('    ')
                freqList.append(1)
                fluxList.append(0.0)
                rAList.append(sourceDict[source][1]*180.0/pi)
                decList.append(sourceDict[source][2]*180.0/pi)
                eqList.append('J2000')
                velList.append(inDict[0][8]*1000.0)
                vtList.append('LSR')
                vdList.append('RADIO')
                restList.append(lowestFSky)
            c1  = pyfits.Column(name='SOURCE_ID',   format='1J',  array=iDList  )
            c2  = pyfits.Column(name='SOURCE',      format='16A', array=nameList)
            c3  = pyfits.Column(name='QUAL',        format='1J',  array=qualList)
            c4  = pyfits.Column(name='CALCODE',     format='4A',  array=calList )
            c5  = pyfits.Column(name='FREQID',      format='1J',  array=freqList)
            c6  = pyfits.Column(name='IFLUX',       format='E1',  array=fluxList)
            c7  = pyfits.Column(name='QFLUX',       format='E1',  array=fluxList)
            c8  = pyfits.Column(name='UFLUX',       format='E1',  array=fluxList)
            c9  = pyfits.Column(name='VFLUX',       format='E1',  array=fluxList)
            c10 = pyfits.Column(name='ALPHA',       format='E1',  array=fluxList)
            c11 = pyfits.Column(name='FREQOFF',     format='E1',  array=fluxList)
            c12 = pyfits.Column(name='RAEPO',       format='1D',  array=rAList  )
            c13 = pyfits.Column(name='DECEPO',      format='1D',  array=decList )
            c14 = pyfits.Column(name='EQUINOX',     format='8A',  array=eqList  )
            c15 = pyfits.Column(name='RAAPP',       format='1D',  array=rAList  )
            c16 = pyfits.Column(name='DECAPP',      format='1D',  array=decList )
            c17 = pyfits.Column(name='SYSVEL',      format='D' ,  array=velList )
            c18 = pyfits.Column(name='VELTYP',      format='8A',  array=vtList  )
            c19 = pyfits.Column(name='VELDEF',      format='8A',  array=vdList  )
            c20 = pyfits.Column(name='RESTFREQ',    format='D' ,  array=restList)
            c21 = pyfits.Column(name='PMRA',        format='D' ,  array=fluxList)
            c22 = pyfits.Column(name='PMDEC',       format='D' ,  array=fluxList)
            c23 = pyfits.Column(name='PARALLAX',    format='E' ,  array=fluxList)
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16,
                                      c17, c18, c19, c20, c21, c22, c23])
            sourceHDU = pyfits.new_table(coldefs)
            header = sourceHDU.header
            header.update('extname',  'SOURCE')
            header.update('tabrev',   1)
            header.update('obscode',  ' ');
            header.update('no_stkd',  1,                         'Number of Stokes parameters')
            header.update('stk_1',    -2,                        'First Stokes parameter')
            header.update('no_band',  1,                         'Number of bands')
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels')
            header.update('ref_freq', lowestFSky,                'File reference frequency')
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth')
            header.update('ref_pixl', 1,                         'Reference pixel')

            ###
            ### Create the FREQUENCY table
            ###
            c1  = pyfits.Column(name='FREQID',         format='1J',  array=[1]  )
            c2  = pyfits.Column(name='BANDFREQ',        format='D',   array=[0.0])
            c3  = pyfits.Column(name='CH_WIDTH',        format='1E' ,  array=[abs(spSmallDict[band][1])])
            c4  = pyfits.Column(name='TOTAL_BANDWIDTH', format='1E',  array=[abs(spSmallDict[band][1])*float(spSmallDict[band][0])])
            if sideBand == 'Upper':
                c5  = pyfits.Column(name='SIDEBAND',        format='1J',  array=[1] )
            else:
                c5  = pyfits.Column(name='SIDEBAND',        format='1J',  array=[-1])
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
            frequencyHDU = pyfits.new_table(coldefs)
            header = frequencyHDU.header
            header.update('extname',  'FREQUENCY')
            header.update('tabrev',   1)
            header.update('obscode',  ' ');
            header.update('no_stkd',  1,                         'Number of Stokes parameters')
            header.update('stk_1',    -2,                        'First Stokes parameter')
            header.update('no_band',  1,                         'Number of bands')
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels')
            header.update('ref_freq', lowestFSky,                'File reference frequency')
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth')
            header.update('ref_pixl', 1,                         'Reference pixel')

            ###
            ### Create the ANTENNA table
            ###
            timeList = []
            timeIntList = []
            arrayList = []
            levelsList = []
            polAList = []
            polBList = []
            polCalAList = []
            for ant in range(1, 9):
                timeList.append(0.0)
                timeIntList.append(1.0)
                arrayList.append(1)
                levelsList.append(2)
                polAList.append('L')
                polBList.append('R')
                polCalAList.append([0.0, 0.0])
            c1  = pyfits.Column(name='TIME',          format='D',  array=timeList   )
            c2  = pyfits.Column(name='TIME_INTERVAL', format='E',  array=timeIntList)
            c3  = pyfits.Column(name='ANNAME',        format='8A', array=antNames   )
            c4  = pyfits.Column(name='ANTENNA_NO',    format='J', array=antNumList )
            c5  = pyfits.Column(name='ARRAY',         format='J',  array=arrayList  )
            c6  = pyfits.Column(name='FREQID',        format='J',  array=arrayList  )
            c7  = pyfits.Column(name='NO_LEVELS',     format='J',  array=levelsList )
            c8  = pyfits.Column(name='POLTYA',        format='1A', array=polAList   )
            c9  = pyfits.Column(name='POLAA',         format='E',  array=timeList   )
            c10 = pyfits.Column(name='POLCALA',       format='2E', array=polCalAList)
            c11 = pyfits.Column(name='POLTYB',        format='1A', array=polBList   )
            c12 = pyfits.Column(name='POLAB',         format='E',  array=timeList   )
            c13 = pyfits.Column(name='POLCALB',       format='2E', array=polCalAList)
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13])
            antennaHDU = pyfits.new_table(coldefs)
            header = antennaHDU.header
            header.update('extname',  'ANTENNA')
            header.update('nopcal',   2)
            header.update('poltype', 'APPROX')

            header.update('tabrev',   1)
            header.update('obscode',  ' ');
            header.update('no_stkd',  1,                         'Number of Stokes parameters')
            header.update('stk_1',    -2,                        'First Stokes parameter')
            header.update('no_band',  1,                         'Number of bands')
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels')
            header.update('ref_freq', lowestFSky,                'File reference frequency')
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth')
            header.update('ref_pixl', 1,                         'Reference pixel')

            ###
            ### Create the SYSTEM_TEMPERATURE table
            ###
            blKeysSorted = sorted(blDict)
#            for bl in blKeysSorted:
#                print bl, blTsysDict[bl][0], blTsysDict[bl][1], blTsysDict[bl][2], blTsysDict[bl][3]
            blPos = 0
            scanDict = {}
            antTsysDict = {}
            inKeysSorted = sorted(inDict)
            antennaListSorted = sorted(antennaList)
            for inh in inKeysSorted:
                scanDict[inh] = (inDict[inh][7]/24.0, inDict[inh][12]/86400.0, inDict[inh][14])
                antCount = 0
                antList = []
                for bl in blKeysSorted[blPos:]:
                    if blDict[bl][0] == inh:
                        ant1 = blDict[bl][15]
                        ant2 = blDict[bl][16]
                        if ant1 not in antList:
                            antList.append(ant1)
                            antCount += 1
                            antTsysDict[(inh, ant1, 1)] = blTsysDict[bl][0]
                            antTsysDict[(inh, ant1, 2)] = blTsysDict[bl][1]
                        if ant2 not in antList:
                            antList.append(ant2)
                            antCount += 1
                            antTsysDict[(inh, ant2, 1)] = blTsysDict[bl][2]
                            antTsysDict[(inh, ant2, 2)] = blTsysDict[bl][3]
                        blPos += 1
                        if antCount == len(antennaList):
                            break
                if antCount != len(antennaList):
                    print 'Only %d antenna entries found for scan %d - aborting' % (antCount, inh)
                    sys.exit(-1)
            timeList = []
            intervalList = []
            sourceList = []
            antList = []
            allOnesList = []
            tsys1List = []
            tsys2List = []
            nanList = []
            for inh in inKeysSorted:
                for ant in antennaListSorted:
                    timeList.append(scanDict[inh][0])
                    intervalList.append(scanDict[inh][1])
                    sourceList.append(scanDict[inh][2])
                    antList.append(ant)
                    allOnesList.append(1)
                    tsys1List.append(antTsysDict[(inh, ant, 1)])
                    tsys2List.append(antTsysDict[(inh, ant, 2)])
                    nanList.append(NaN)
            c1  = pyfits.Column(name='TIME',          format='1D',  array=timeList    )
            c2  = pyfits.Column(name='TIME_INTERVAL', format='1E',  array=intervalList)
            c3  = pyfits.Column(name='SOURCE_ID',     format='1J',  array=sourceList  )
            c4  = pyfits.Column(name='ANTENNA_NO',    format='1J',  array=antList     )
            c5  = pyfits.Column(name='ARRAY',         format='1J',  array=allOnesList )
            c6  = pyfits.Column(name='FREQID',        format='1J',  array=allOnesList )
            c7  = pyfits.Column(name='TSYS_1',        format='E',   array=tsys1List   )
            c8  = pyfits.Column(name='TANT_1',        format='E',   array=nanList     )
            c9  = pyfits.Column(name='TSYS_2',        format='E',   array=tsys2List   )
            c10 = pyfits.Column(name='TANT_2',        format='E',   array=nanList     )
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10])
            tsysHDU = pyfits.new_table(coldefs)
            header = tsysHDU.header
            header.update('extname',  'SYSTEM_TEMPERATURE')
            header.update('tabrev',   1)
            header.update('no_pol', 2)
            header.update('obscode',  ' ');
            header.update('no_stkd',  1,                         'Number of Stokes parameters')
            header.update('stk_1',    -2,                        'First Stokes parameter')
            header.update('no_band',  1,                         'Number of bands')
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels')
            header.update('ref_freq', lowestFSky,                'File reference frequency')
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth')
            header.update('ref_pixl', 1,                         'Reference pixel')

            ###
            ### Create the UV_DATA table
            ###
            if spSmallDict[band][1] < 0:
                reverseChannelOrder = True
            else:
                reverseChannelOrder = False
            scanOffset = 0
            scanNo = 0
            blPos = 0
            if verbose:
                print 'Creating UV_DATA table for band ', band, sideBand, ' with ', spSmallDict[band][0], ' channels'
            uList        = []
            vList        = []
            wList        = []
            jDList       = []
            timeList     = []
            baselineList = []
            arrayList    = []
            sourceList   = []
            freqList     = []
            intList      = []
            matrixList   = []
            nChannels = spSmallDict[band][0]
            if band == 0:
                firstGoodChannel = 0
                lastGoodChannel = 1000000
            else:
                firstGoodChannel = int(float(nChannels)*edgeTrimFraction)
                lastGoodChannel  = int(float(nChannels)*(1.0-edgeTrimFraction) + 0.5)
            while (scanOffset < visFileLen) and (scanNo < maxScan):
                foundBlEntry = False
                foundSpEntry = False
                nBlFound = 0
                scanNo = makevis.scanno(scanOffset, newFormat);
                recSize = makevis.recsize(scanOffset, newFormat)
                if scanNo > 0:
                    for bl in blKeysSorted[blPos:]:
                        if blDict[bl][0] == scanNo:
                            foundBlEntry = True
                            nBlFound += 1
                            if (band, bl) in spBigDict:
                                foundSpEntry = True
                                matrixEntry  = []
                                uList.append(blDict[bl][6]*1000.0/pseudoContinuumFrequency[sb])
                                vList.append(blDict[bl][7]*1000.0/pseudoContinuumFrequency[sb])
                                wList.append(blDict[bl][8]*1000.0/pseudoContinuumFrequency[sb])
                                jDList.append(jD)
                                timeList.append(inDict[scanNo][7]/24.0)
                                baselineList.append(256*blDict[bl][15] + blDict[bl][16])
                                arrayList.append(1)
                                sourceList.append(inDict[scanNo][14])
                                freqList.append(1)
                                intList.append(inDict[scanNo][12])
                                if newFormat:
                                    dataoff = scanOffset+spBigDict[(band, bl)][0] + 8
                                else:
                                    dataoff = scanOffset+spBigDict[(band, bl)][0] + 16
                                scaleExp = makevis.scaleexp(dataoff, newFormat)
                                scale = (2.0**scaleExp) * sqrt(2.0)*130.0*2.0
                                if bl == 0:
                                    ant1Tsys = sqrt(abs(blTsysDict[bl][0]*blTsysDict[bl][1]))
                                    ant2Tsys = sqrt(abs(blTsysDict[bl][2]*blTsysDict[bl][3]))
                                if 24 < bl < 49:
                                    ant1Tsys = blTsysDict[bl][1]
                                    ant2Tsys = blTsysDict[bl][3]
                                else:
                                    ant1Tsys = blTsysDict[bl][0]
                                    ant2Tsys = blTsysDict[bl][2]
                                scale *= sqrt(abs(ant1Tsys*ant2Tsys))
                                if newFormat:
                                    weight = spBigDict[(band, bl)][1]/(maxWeight*ant1Tsys*ant2Tsys)
                                else:
                                    weight = spBigDict[(band, bl)][1]
                                matrixEntry = makevis.convert(nChannels, scale, dataoff, newFormat, weight, trimEdges, firstGoodChannel, lastGoodChannel, reverseChannelOrder)
                                matrixList.append(matrixEntry)
                                if nBlFound == numberOfBaselines:
                                    break
                        blPos += 1
                    if (not foundBlEntry) or (not foundSpEntry):
                        print 'Something not found scanNo = %d, band = %d, foundBl = %d, foundSp = %d' % (scanNo, band, foundBlEntry, foundSpEntry)
                        sys.exit(-1)
                if newFormat:
                    scanOffset += recSize+8
                else:
                    scanOffset += recSize+16

            c10Format = '%dE' % (len(matrixEntry))
            c1  = pyfits.Column(name='UU',          format='1D',       array=uList       , unit='SECONDS')
            c2  = pyfits.Column(name='VV',          format='1D',       array=vList       , unit='SECONDS')
            c3  = pyfits.Column(name='WW',          format='1D',       array=wList       , unit='SECONDS')
            c4  = pyfits.Column(name='DATE',        format='1D',       array=jDList      , unit='DAYS'   )
            c5  = pyfits.Column(name='TIME',        format='1D',       array=timeList    , unit='DAYS'   )
            c6  = pyfits.Column(name='BASELINE',    format='1J',       array=baselineList                )
            c7  = pyfits.Column(name='SOURCE_ID',   format='1J',       array=sourceList                  )
            c8  = pyfits.Column(name='FREQID',      format='1J',       array=freqList                    )
            c9  = pyfits.Column(name='INTTIM',      format='1E',       array=intList     , unit='SECONDS')
            c10 = pyfits.Column(name='FLUX',        format=c10Format,  array=matrixList  , unit='UNCALIB')
            coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10])
            uvDataHDU = pyfits.new_table(coldefs)
            header = uvDataHDU.header
            header.update('extname',  'UV_DATA'                                                          )
            header.update('nmatrix',  1,                          'Number of matrices'                   )
            header.update('tmatx10',  'T'                                                                )
            header.update('maxis',    6,                          'Number of matix axes'                 )
            header.update('maxis1',   3,                          'Number of data points on complex axis')
            header.update('ctype1',   'COMPLEX',                  'Type of axis 1'                       )
            header.update('cdelt1',   1.0                                                                )
            header.update('crval1',   1.0                                                                )
            header.update('crpix1',   1.0                                                                )
            header.update('maxis2',   1,                          'Number of data points on Stokes axis' )
            header.update('ctype2',   'STOKES',                   'Type of axis 2'                       )
            header.update('cdelt2',   0.0                                                                )
            header.update('crval2',   -2.0                                                               )
            header.update('crpix2',   1.0                                                                )
            header.update('maxis3',   spSmallDict[band][0],       'Number of data points on Freq axis'   )
            header.update('ctype3',   'FREQ',                     'Type of axis 3'                       )
            header.update('cdelt3',   abs(spSmallDict[band][1]),  'Channel spacing in Hz'                )
            header.update('crval3',   lowestFSky,                 'Frequency of reference pixel'         )
            header.update('crpix3',   1.0                                                                )
            header.update('maxis4',   1,                          'Number of data points on Band axis  ' )
            header.update('ctype4',   'BAND',                     'Type of axis 4'                       )
            header.update('cdelt4',   1.0                                                                )
            header.update('crval4',   1.0                                                                )
            header.update('crpix4',   1.0                                                                )
            header.update('maxis5',   1,                          'Number of data points on RA axis'     )
            header.update('ctype5',   'RA',                       'Type of axis 5'                       )
            header.update('cdelt5',   1.0                                                                )
            header.update('crval5',   0.0,                        'Must be 0 for a multisource file'     )
            header.update('crpix5',   1.0                                                                )
            header.update('maxis6',   1,                          'Number of data points on Dec axis'    )
            header.update('ctype6',   'DEC',                      'Type of axis 6'                       )
            header.update('cdelt6',   1.0                                                                )
            header.update('crval6',   0.0,                        'Must be 0 for a multisource file'     )
            header.update('crpix6',   1.0                                                                )
            
            header.update('tabrev',   2                                                                  )
            header.update('obscode',  ' '                                                                )
            header.update('no_stkd',  1,                         'Number of Stokes parameters'           )
            header.update('stk_1',    -2,                        'First Stokes parameter'                )
            header.update('no_band',  1,                         'Number of bands'                       )
            header.update('no_chan',  spSmallDict[band][0],      'The number of spectral channels'       )
            header.update('ref_freq', lowestFSky,                'File reference frequency'              )
            header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth'                 )
            header.update('ref_pixl', 1,                         'Reference pixel'                       )
            
#            header.update('weightyp', 'NORMAL',                  'Normal 1/(uncertainty**2) weights'     )
            header.update('telescop', 'SMA',                     'Submillimeter Array, Hawaii'           )
            header.update('observer', 'Taco'                                                             )

            ###
            ### Create the file and write all tables
            ###
            fileName = 'tempFITS-IDI%s.band%d' % (sideBand, band)
            if verbose:
                print 'Writing primary header for file ', fileName
            hdulist.writeto(fileName)
            if verbose:
                print 'Writing ARRAY_GEOMETRY table'
            pyfits.append(fileName, arrayGeometryHDU.data, header=arrayGeometryHDU.header)
            if verbose:
                print 'Writing SOURCE table'
            pyfits.append(fileName, sourceHDU.data,        header=sourceHDU.header       )
            if verbose:
                print 'Writing FREQUENCY table'
            pyfits.append(fileName, frequencyHDU.data,     header=frequencyHDU.header    )
            if verbose:
                print 'Writing ANTENNA table'
            pyfits.append(fileName, antennaHDU.data,       header=antennaHDU.header      )
            if verbose:
                print 'Writing SYSTEM_TEMPERATURE table'
            pyfits.append(fileName, tsysHDU.data,          header=tsysHDU.header         )
            if verbose:
                print 'Writing UV_DATA table'
            pyfits.append(fileName, uvDataHDU.data,        header=uvDataHDU.header       )
