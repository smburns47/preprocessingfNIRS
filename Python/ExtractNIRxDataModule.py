#-----------------------------------------------
#step 1 - subfunctions for extracting data
#-----------------------------------------------

import numpy as np
import glob
import re
import sys

from typing import Dict, Any, Union

ignorePattern = r'[^a-zA-Z].*'

def parse(fname):
    f = open(fname, 'r')
    inValue = False
    key = None
    value = None
    result = {}
    for line in f:
        if not inValue:
            if re.match(ignorePattern, line) is not None:
                continue
            keyValue = line.split('=')
            key = keyValue[0]
            value = keyValue[1]
            if value[0] == '"':
                value = value[1:]
                splitVal = value.split('"')
                if len(splitVal) == 1:
                    inValue = True
                else:
                    value = splitVal[0]
                    result[key] = value
                    key = None
                    value = None
            else:
                result[key] = value
                key = None
                value = None
        else:
            splitLine = line.split('"')
            value = value + splitLine[0]
            if len(splitLine) > 1:
                result[key] = value
                key = None
                value = None
                inValue = False
    f.close()
    return result


def extractNIRxData(subjfolder):
    wl1filename = glob.glob(subjfolder + "/" + "*.wl1")
    wl2filename = glob.glob(subjfolder + "/" + "*.wl2")
    hdrfilename = glob.glob(subjfolder + "/" + "*.hdr")
    wl1 = np.genfromtxt(wl1filename[0], delimiter=" ")
    wl2 = np.genfromtxt(wl2filename[0], delimiter=" ")
    d = np.concatenate((wl1,wl2), axis=1)

    hdrDict = parse(hdrfilename[0])

    samprate = float(hdrDict["SamplingRate"])

    sdmask = hdrDict["S-D-Mask"][1:-1]
    sdmask = np.fromstring(sdmask, sep="	")
    sdmask = np.concatenate((sdmask, sdmask), axis=0)
    sd_ind = np.nonzero(sdmask)
    sd_ind = sd_ind[0]

    wavelength1 = int(hdrDict["Wavelengths"][0:3])
    wavelength2 = int(hdrDict["Wavelengths"][4:])
    wavelengths = np.array([wavelength1, wavelength2])

    d = d[:,[sd_ind]]
    d = np.squeeze(d)

    evts = hdrDict["Events"][1:-1]
    if len(evts)>0:
        evts = np.fromstring(evts, sep="    ")
        evtType = evts[1::3]
        timeFrame = evts[2::3]
        evts = np.vstack((evtType, timeFrame)).T
        evts = np.unique(evts,axis=0)
        evts = evts[timeFrame!=0]
        uniqueEvtTypes = np.unique(evtType).size
        s = np.zeros((d.shape[0],uniqueEvtTypes))
        for i in range(uniqueEvtTypes):
            s[[timeFrame], i-1] = 1
    else:
        s = np.zeros((d.shape[0], 1))

    NIRxData = dict()  # type: Dict[str, Union[float, Any]]
    NIRxData['d'] = d
    NIRxData['sd_ind'] = sd_ind
    NIRxData['samprate'] = samprate
    NIRxData['wavelengths'] = wavelengths
    NIRxData['s'] = s
    NIRxData['nSrcs'] = int(hdrDict["Sources"])
    NIRxData['nDets'] = int(hdrDict["Detectors"])
    NIRxData['numchannels'] = d.shape[1]/2

    return NIRxData
