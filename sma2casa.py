#!/home/rtm/myPythonEnv/bin/python
"""
sma2casa.py
===========

This file defines functions which will import SMA data directly into a CASA Measurement Set.
"""

import numpy as np
import os, sys, struct, mmap, pyfits
from astropy.time import Time
from math import sin, cos, pi
#import psyco
#psyco.full()

neededFiles = ('antennas', 'bl_read', 'codes_read', 'in_read', 'sch_read', 'sp_read', 'tsys_read')

nBands = 0
bandList = []
spSmallDict = {}
spBigDict = {}
antennas = {}
codesDict = {}
inDict = {}
blDictL = {}
blDictU = {}
sourceDict = {}
lowestFSky = 1.0e30
speedOfLight = 2.997925e8
maxScan = 100

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
    global nBands, bandList, antennas, codesDict, inDict, blDictL, blDictU, lowestFSky
    global spSmallDict, spBigDict, sourceDict

    # Check that the directory contains all the required files
    dirContents = os.listdir(dataDir)
    for needed in neededFiles:
        if not needed in dirContents:
            print "The directory %s does not have a %s file - aborting" % (dataDir, needed)
            sys.exit(-1)

    ###
    ### Read in the antennas file - make a dictionary antennas[ant #] = (X, Y, Z)
    ###
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
    print 'Reading in_read'
    f = open(dataDir+'/in_read', 'rb')
    data = f.read()
    inRecLen = 188
    nInRecords = len(data)/inRecLen
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
        rinteg    =  makeFloat(data[rec*inRecLen +  60:])
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
            sourceDict[souid] = (codesDict['source'][isource], rar, decr)

    ###
    ### Read in bl_read
    ###
    print 'Reading bl_read'
    f = open(dataDir+'/bl_read', 'rb')
    fSize = os.path.getsize(dataDir+'/bl_read')
    blRecLen = 158-8
    nBlRecords = fSize/blRecLen
    for rec in range(nBlRecords):
        if (rec % 10000) == 0:
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
            # Put Tsys offsets here!
            if isb == 0:
                blDictL[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2)
            else:
                blDictU[blhid] = (inhid, ipol, ant1rx, ant2rx, pointing, irec, u, v, w, prbl, coh, avedhrs, ampave, phaave,
                                  blsid, ant1, ant2)

    ###
    ### Count the number of spectral bands
    ###
    print 'Counting the number of bands'
    f = open(dataDir+'/sp_read', 'rb')
    fSize = os.path.getsize(dataDir+'/sp_read')
    spRecLen = 188
    data = f.read(spRecLen)
    nSpRecords = fSize/spRecLen
    for rec in range(nSpRecords):
        if (rec % 100000) == 0:
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
            wt        =  makeFloat(data[ 80:])
            flags     =    makeInt(data[ 84:], 4)
            nch       =    makeInt(data[ 96:], 2)
            dataoff   =    makeInt(data[100:], 4)
            rfreq     = makeDouble(data[104:])
            if flags != 0:
                wt = 0.0
            if fsky < lowestFSky:
                lowestFSky = fsky
            spBigDict[(iband, blhid)] = (dataoff, wt)
            if iband not in bandList:
                bandList.append(iband)
                spSmallDict[iband] = (nch, fres*1.0e6, rfreq*1.0e9)
        if inhid > maxScan:
            break
    nBands = len(bandList)
    bandList.sort()

dataSet = '/sma/SMAusers/taco/SWARMTest/sma/rtdata/newFormat/mir_data/130307_08:57:53'
for line in open(dataSet+'/projectInfo'):
    projectPI = line.strip()
read(dataSet)
visFile = os.open(dataSet+'/sch_read', os.O_RDONLY)
visMap = mmap.mmap(visFile, 0, prot=mmap.PROT_READ);
print 'len(visMap) = ', len(visMap)
#for band in bandList:
for band in range(1):
    for sb in range(1):
        if sb == 0:
            sideBand = 'Lower'
            blDict = blDictL
        else:
            sideBand = 'Upper'
            blDict = blDictU
        
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
        header.update('extver',   1,                  'Array number')
        header.update('arrayx',   0.0,                'Array origin x coordinate')
        header.update('arrayy',   0.0,                'Array origin y coordinate')
        header.update('arrayz',   0.0,                'Array origin z coordinate')
        header.update('numorb',   0,                  'Number of orbital parameters')
        header.update('polarx',   0.0,                'x coordinate of north pole (arcsec)')
        header.update('polary',   0.0,                'y coordinate of north pole (arcsec)')
        header.update('freq',     lowestFSky*1.0e9,   'Reference frequency')
        header.update('timesys',  'UTC',              'Time system')
        header.update('degpdy' ,  0.0,                'Earth rotation rate')
        header.update('ut1utc' ,  0.0,                'UT1 - UTC (seconds) - PLACEHOLDER')
        header.update('iatutc' ,  0.0,                'IAT - UTC (seconds) - PLACEHOLDER')
        header.update('rdate',    refTimeString,      'Reference date')
        header.update('gstia0',   gST,                'GST at 0h on reference date (degrees)')
        header.update('no_stkd',  1,                  'Number of Stokes parameters')
        header.update('stk_1',    -2,                 'First Stokes parameter')
        header.update('no_band',  1,                  'Number of bands')
        header.update('ref_freq', lowestFSky*1.0e9,   'File reference frequency')
        header.update('ref_pixl', 1,                  'Reference pixel')
        header.update('no_chan',  spSmallDict[band][0],  'The number of spectral channels')
        header.update('chan_bw',  abs(spSmallDict[band][1]),  'The channel bandwidth')

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
            restList.append(spSmallDict[band][2])
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
        header.update('no_stkd',  1,                  'Number of Stokes parameters')
        header.update('stk_1',    -2,                 'First Stokes parameter')
        header.update('no_band',  1,                  'Number of bands')
        header.update('no_chan',  spSmallDict[band][0],  'The number of spectral channels')
        header.update('ref_freq', lowestFSky*1.0e9,   'File reference frequency')
        header.update('chan_bw',  abs(spSmallDict[band][1]),  'The channel bandwidth')
        header.update('ref_pixl', 1,                  'Reference pixel')

        ###
        ### Create the FREQUENCY table
        ###
        c1  = pyfits.Column(name='FREQ_ID',         format='1J',  array=[1]  )
        c2  = pyfits.Column(name='BANDFREQ',        format='D',   array=[0.0])
        c3  = pyfits.Column(name='CH_WIDTH',        format='E' ,  array=[abs(spSmallDict[band][1])])
        c4  = pyfits.Column(name='TOTAL_BANDWIDTH', format='1J',  array=[abs(spSmallDict[band][1])*float(spSmallDict[band][0])])
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
        header.update('no_stkd',  1,                  'Number of Stokes parameters')
        header.update('stk_1',    -2,                 'First Stokes parameter')
        header.update('no_band',  1,                  'Number of bands')
        header.update('no_chan',  spSmallDict[band][0],  'The number of spectral channels')
        header.update('ref_freq', lowestFSky*1.0e9,   'File reference frequency')
        header.update('chan_bw',  abs(spSmallDict[band][1]),  'The channel bandwidth')
        header.update('ref_pixl', 1,                  'Reference pixel')

        ###
        ### Create the ANTENNA table
        ###
        timeList = []
        timeIntList = []
        arrayList = []
        levelsList = []
        polAList = []
        polBList = []
        for ant in range(1,9):
            timeList.append(0.0)
            timeIntList.append(1.0)
            arrayList.append(1)
            levelsList.append(2)
            polAList.append('L')
            polBList.append('L')
        c1  = pyfits.Column(name='TIME',          format='D',  array=timeList   )
        c2  = pyfits.Column(name='TIME_INTERVAL', format='E',  array=timeIntList)
        c3  = pyfits.Column(name='ANNAME',        format='8A', array=antNames   )
        c4  = pyfits.Column(name='ANTENNA_NO',    format='8A', array=antNumList )
        c5  = pyfits.Column(name='ARRAY',         format='J',  array=arrayList  )
        c6  = pyfits.Column(name='FREQID',        format='J',  array=arrayList  )
        c7  = pyfits.Column(name='NO_LEVELS',     format='J',  array=levelsList )
        c8  = pyfits.Column(name='POLTYA',        format='1A', array=polAList   )
        c9  = pyfits.Column(name='POLAA',         format='E',  array=timeList   )
        c10 = pyfits.Column(name='POLCALA',       format='E',  array=timeList   )
        c11 = pyfits.Column(name='POLTYB',        format='1A', array=polBList   )
        c12 = pyfits.Column(name='POLAB',         format='E',  array=timeList   )
        c13 = pyfits.Column(name='POLCALB',       format='E',  array=timeList   )
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13])
        antennaHDU = pyfits.new_table(coldefs)
        header = antennaHDU.header
        header.update('extname',  'ANTENNA')
        header.update('nopcal',   0)

        header.update('tabrev',   1)
        header.update('obscode',  ' ');
        header.update('no_stkd',  1,                  'Number of Stokes parameters')
        header.update('stk_1',    -2,                 'First Stokes parameter')
        header.update('no_band',  1,                  'Number of bands')
        header.update('no_chan',  spSmallDict[band][0],  'The number of spectral channels')
        header.update('ref_freq', lowestFSky*1.0e9,   'File reference frequency')
        header.update('chan_bw',  abs(spSmallDict[band][1]),  'The channel bandwidth')
        header.update('ref_pixl', 1,                  'Reference pixel')

        ###
        ### Create the UV_DATA table
        ###
        # Map the entire visibilities file into memory, for easy (and fast) access
        scanOffset = 0
        scanNo = 0
        print 'Creating UV_DATA table for band ', band, ' with ', spSmallDict[band][0], ' channels'
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
        weightList   = []
        matrixList   = []
        while (scanOffset < len(visMap)) and (scanNo < maxScan):
            foundBlEntry = False
            foundSpEntry = False
            nBlFound = 0
            scanNo  = makeInt(visMap[scanOffset:scanOffset+4],   4)
            recSize = makeInt(visMap[scanOffset+4:scanOffset+8], 4)
            if scanNo > 0:
                for bl in blDict:
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
                            weightList.append(spBigDict[(band, bl)][1])
                            dataoff = scanOffset+spBigDict[(band, bl)][0] + 8
                            scaleExp = makeInt(visMap[dataoff:dataoff+2], 2)
                            if scaleExp > 2**15:
                                scaleExp -= 2**16
                            scale = 2.0**scaleExp
                            dataoff += 2
                            for i in range(0, spSmallDict[band][0]):
                                real = float(ord(visMap[dataoff  ])+(ord(visMap[dataoff+1])<<8))*scale
                                imag = float(ord(visMap[dataoff+2])+(ord(visMap[dataoff+3])<<8))*scale
                                matrixEntry.append(real)
                                matrixEntry.append(imag)
                            matrixList.append(matrixEntry)
                        if nBlFound == 28:
                            break
                if (not foundBlEntry) or (not foundSpEntry):
                    print 'Something not found scanNo = %d, band = %d, foundBl = %d, foundSp = %d' % (scanNo, band, foundBlEntry, foundSpEntry)
                    sys.exit(0)
            scanOffset += recSize+8

        # We must loop through all the visibilities in the sch_read file, so first we
        # need to figure out how large each record is.

        c12Format = '%dE' % (len(matrixEntry))
        print 'Matrix format: ', c12Format
        c1  = pyfits.Column(name='UU',          format='1D',       array=uList       , unit='SECONDS')
        c2  = pyfits.Column(name='VV',          format='1D',       array=vList       , unit='SECONDS')
        c3  = pyfits.Column(name='WW',          format='1D',       array=wList       , unit='SECONDS')
        c4  = pyfits.Column(name='DATE',        format='1D',       array=jDList      , unit='DAYS'   )
        c5  = pyfits.Column(name='TIME',        format='1D',       array=timeList    , unit='DAYS'   )
        c6  = pyfits.Column(name='BASELINE',    format='1J',       array=baselineList                )
        c7  = pyfits.Column(name='ARRAY',       format='1J',       array=arrayList                   )
        c8  = pyfits.Column(name='SOURCE_ID',   format='1J',       array=sourceList                  )
        c9  = pyfits.Column(name='FREQID',      format='1J',       array=freqList                    )
        c10 = pyfits.Column(name='INTTIM',      format='1E',       array=intList     , unit='SECONDS')
        c11 = pyfits.Column(name='WEIGHT',      format='1E',       array=weightList                  )
        c12 = pyfits.Column(name='FLUX',        format=c12Format,  array=matrixList  , unit='UNCALIB')
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12])
        uvDataHDU = pyfits.new_table(coldefs)
        header = uvDataHDU.header
        header.update('extname',  'UV_DATA'                                                          )
        header.update('nmatrix',  1,                          'Number of matrices'                   )
        header.update('tmatx12',  'T'                                                                )
        header.update('maxis',    6,                          'Number of matix axes'                 )
        header.update('maxis1',   2,                          'Number of data points on complex axis')
        header.update('ctype1',   'COMPLEX',                  'Type of axis 1'                       )
        header.update('cdelt1',   1                                                                  )
        header.update('crval1',   1                                                                  )
        header.update('crpix1',   1                                                                  )
        header.update('maxis2',   1,                          'Number of data points on Stokes axis' )
        header.update('ctype2',   'STOKES',                   'Type of axis 2'                       )
        header.update('cdelt2',   0.0                                                                )
        header.update('crval2',   -2.0                                                               )
        header.update('crpix2',   1.0                                                                )
        header.update('maxis3',   spSmallDict[band][0],       'Number of data points on Freq axis'   )
        header.update('ctype3',   'FREQ',                     'Type of axis 3'                       )
        header.update('cdelt3',   abs(spSmallDict[band][1]),  'Channel spacing in Hz'                )
        header.update('crval3',   spSmallDict[band][2],       'Frequency of reference pixel'         )
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
        header.update('ref_freq', lowestFSky*1.0e9,          'File reference frequency'              )
        header.update('chan_bw',  abs(spSmallDict[band][1]), 'The channel bandwidth'                 )
        header.update('ref_pixl', 1,                         'Reference pixel'                       )

        header.update('telescop', 'SMA',                     'Submillimeter Array, Hawaii'           )
        header.update('observer', 'Taco'                                                             )

        ###
        ### Create the file and write all tables
        ###
        fileName = 'tempFITS-IDI%s.band%d' % (sideBand, band)
        hdulist.writeto(fileName)
        pyfits.append(fileName, arrayGeometryHDU.data, header=arrayGeometryHDU.header)
        pyfits.append(fileName, sourceHDU.data,        header=sourceHDU.header       )
        pyfits.append(fileName, frequencyHDU.data,     header=frequencyHDU.header    )
        pyfits.append(fileName, antennaHDU.data,       header=antennaHDU.header      )
        pyfits.append(fileName, uvDataHDU.data,        header=uvDataHDU.header       )

