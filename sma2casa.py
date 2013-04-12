#!/sma/SMAusers/taco/myPythonEnv/bin/python
"""
sma2casa.py
===========

This file defines functions which will import SMA data directly into a CASA Measurement Set.
"""

import numpy as np
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
speedOfLight = 2.997925e8
maxScan = 10000000
maxWeight = 0.01
numberOfBaselines = 0
edgeTrimFraction = 0.1 # Fraction on each edge of a spectral chunk to flag bad
chunkList = range(1,49)
sidebandList = [0, 1]
verbose = True

def usage():
    print 'Usage: ', sys.argv[0], ' path-to-SMA-data [options]'
    print 'Options: (which must come *after* the path argument)'
    print '\t -c m,n,o [--chunks=m,n,o]\tProcess only chunks m, n and o'
    print '\t -h [--help]\t\t\tPrint this message and exit'
    print '\t -l [--lower]\t\t\tProcess lower sideband only (default is both)'
    print '\t -s [--silent]\t\t\tRun silently unless an error occurs'
    print '\t -u [--upper]\t\t\tProcess upper sideband only (default is both)'

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
    for i in range(size):
        tInt += ord(data[i])<<(i<<3)
    return tInt

def makeFloat(data):
    return (struct.unpack('f', data[:4]))[0]

def makeDouble(data):
    return (struct.unpack('d', data[:8]))[0]

def read(dataDir):
    global nBands, bandList, antennas, codesDict, inDict, blDictL, blDictU
    global spSmallDictL, spSmallDictU, spBigDict, sourceDict, maxWeight, numberOfBaselines
    global blTsysDictL, blTsysDictU

    # Check that the directory contains all the required files
    if verbose:
        print 'checking that the needed files exist in ', dataDir
    dirContents = os.listdir(dataDir)
    for needed in neededFiles:
        if not needed in dirContents:
            print "The directory %s does not have a %s file - aborting" % (dataDir, needed)
            sys.exit(-1)

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
    inRecLen = 188
    nInRecords = len(data)/inRecLen
    savedSourceInfo = True
    for rec in range(nInRecords):
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
    blRecLen = 158
    nBlRecords = fSize/blRecLen
    baselineList = []
    blSidebandDict = {}
    for rec in range(nBlRecords):
        if ((rec % 10000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nBlRecords, 100.0*float(rec)/float(nBlRecords))
            sys.stdout.flush()
        data = f.read(blRecLen)
        if len(data) == blRecLen:
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
            if (ant1, ant2) not in baselineList:
                baselineList.append((ant1, ant2))
                numberOfBaselines +=1
            ant1TsysOff =  makeInt(data[ 64:], 4)
            ant2TsysOff =  makeInt(data[ 68:], 4)
            blSidebandDict[blhid] = isb
            if isb == 0:
                blDictL[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2, ant1TsysOff, ant2TsysOff)
            else:
                blDictU[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2, ant1TsysOff, ant2TsysOff)
    if verbose:
        print numberOfBaselines, ' baselines seen in this data set.'

    ###
    ### Pull the Tsys info out of tsys_read
    ###
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

    ###
    ### Count the number of spectral bands
    ###
    if verbose:
        print 'Counting the number of bands'
    f = open(dataDir+'/sp_read', 'rb')
    fSize = os.path.getsize(dataDir+'/sp_read')
    spRecLen = 188
    data = f.read(spRecLen)
    nSpRecords = fSize/spRecLen
    for rec in range(nSpRecords):
        if ((rec % 100000) == 0) and verbose:
            print '\t processing record %d of %d (%2.0f%% done)' % (rec, nSpRecords, 100.0*float(rec)/float(nSpRecords))
            sys.stdout.flush()
        data = f.read(spRecLen)
        if len(data) == spRecLen:
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
            if flags != 0:
                wt = 0.0
            if wt > maxWeight:
                maxWeight = wt
            spBigDict[(iband, blhid)] = (dataoff, wt)
            if iband not in bandList:
                bandList.append(iband)
            if blSidebandDict[blhid] == 0:
                if iband not in spSmallDictL:
                    spSmallDictL[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
            else:
                if iband not in spSmallDictU:
                    spSmallDictU[iband] = (nch, fres*1.0e6, fsky*1.0e9, rfreq*1.0e9)
        if inhid > maxScan:
            break
    nBands = len(bandList)
    bandList.sort()

if len(sys.argv) < 2:
    usage()
    exit(-1)
dataSet = sys.argv[1]
try:
    opts, args = getopt.getopt(sys.argv[2:], "c:lhsu", ['chunks=', 'help', 'silent'])
except getopt.GetoptError as err:
    usage()
    sys.exit(-1)
for o, a in opts:
    if o in ('-s', '--silent'):
        verbose = False
    elif o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-l', '--lower'):
        sidebandList = [0]
    elif o in ('-u', '--upper'):
        sidebandList = [1]
    elif o in ('-c', '--chunks'):
        sChunkList = a.split(',')
        chunkList = []
        for chunk in sChunkList:
            chunkList.append(int(chunk))
    else:
        print o, a
        assert False, 'unhandled option'

read(dataSet)
for line in open(dataSet+'/projectInfo'):
    projectPI = line.strip()
visFile = os.open(dataSet+'/sch_read', os.O_RDONLY)
visMap = mmap.mmap(visFile, 0, prot=mmap.PROT_READ);
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
            if (band == 0) and ((49 in bandList) or (50 in bandlist)):
                header.update('correlat', 'SMA-Hybrid')
            elif band < 49:
                header.update('correlat', 'SMA-Legacy')
            else:
                header.update('correlat', 'SMA-SWARM')
            header.update('fxcorver', 1)
            header.add_history('Translated to FITS-IDI by sma2casa.py')
            header.add_history('Original SMA dataset: '+dataSet[-15:])
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
            for ant in range(1,9):
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
            for ant in range(1,9):
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
            blPos = 0
            scanDict = {}
            antTsysDict = {}
            inKeysSorted = sorted(inDict)
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
                        if antCount == 8:
                            break
                if antCount != 8:
                    print 'Only %d antenna entries found for scan %d - aborting' % (antCount, inh)
            timeList = []
            intervalList = []
            sourceList = []
            antList = []
            allOnesList = []
            tsys1List = []
            tsys2List = []
            nanList = []
            NaN = float('nan')
            for inh in inKeysSorted:
                for ant in range(1, 9):
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
                firstGoodChannel = int(float(nChannels)*edgeTrimFraction) - 1
                lastGoodChannel  = int(float(nChannels)*(1.0-edgeTrimFraction) + 0.5) - 1
            while (scanOffset < len(visMap)) and (scanNo < maxScan):
                foundBlEntry = False
                foundSpEntry = False
                nBlFound = 0
                scanNo  = makeInt(visMap[scanOffset:scanOffset+4],   4)
                recSize = makeInt(visMap[scanOffset+4:scanOffset+8], 4)
                if scanNo > 0:
                    for bl in blKeysSorted[blPos:]:
                        if blDict[bl][0] == scanNo:
                            foundBlEntry = True
                            nBlFound += 1
                            if (band, bl) in spBigDict:
                                foundSpEntry = True
                                matrixEntry  = []
                                uList.append(blDict[bl][6]/speedOfLight)
                                vList.append(blDict[bl][7]/speedOfLight)
                                wList.append(blDict[bl][8]/speedOfLight)
                                jDList.append(jD)
                                timeList.append(inDict[scanNo][7]/24.0)
                                baselineList.append(256*blDict[bl][15] + blDict[bl][16])
                                arrayList.append(1)
                                sourceList.append(inDict[scanNo][14])
                                freqList.append(1)
                                intList.append(inDict[scanNo][12])
                                dataoff = scanOffset+spBigDict[(band, bl)][0] + 8
                                scaleExp = makeInt(visMap[dataoff:dataoff+2], 2)
                                if scaleExp > (2**15-1):
                                    scaleExp -= 2**16
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
                                weight = spBigDict[(band, bl)][1]/(maxWeight*ant1Tsys*ant2Tsys)
                                dataoff += 2
                                for i in range(0, nChannels):
                                    real = ord(visMap[dataoff  ])+(ord(visMap[dataoff+1])<<8)
                                    if real > (2**15-1):
                                        real -= 2**16
                                    imag = ord(visMap[dataoff+2])+(ord(visMap[dataoff+3])<<8)
                                    if imag > (2**15-1):
                                        imag -= 2**16
                                    matrixEntry.append(float(real*scale))
                                    matrixEntry.append(float(-imag*scale))
                                    if firstGoodChannel <= i <= lastGoodChannel:
                                        matrixEntry.append(weight)
                                    else:
                                        matrixEntry.append(0.0)
                                    dataoff += 4
                                if reverseChannelOrder and (band != 0):
                                    # This is a chunk with an odd number of LSB downconversions, so channels must
                                    # be reordered, because the FITS-IDI standard demands positive increments in
                                    # frequency between channels.
                                    swappedEntry = []
                                    nSets = len(matrixEntry)/3
                                    for i in range(nSets):
                                        swappedEntry.append(matrixEntry[3*nSets - (3 + i*3)])
                                        swappedEntry.append(matrixEntry[3*nSets - (2 + i*3)])
                                        swappedEntry.append(matrixEntry[3*nSets - (1 + i*3)])
                                    matrixList.append(swappedEntry)
                                else:
                                    matrixList.append(matrixEntry)
                                if nBlFound == numberOfBaselines:
                                    break
                        blPos += 1
                    if (not foundBlEntry) or (not foundSpEntry):
                        print 'Something not found scanNo = %d, band = %d, foundBl = %d, foundSp = %d' % (scanNo, band, foundBlEntry, foundSpEntry)
                        sys.exit(0)
                scanOffset += recSize+8

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
