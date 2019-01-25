#-----------------------------------------------
#step 2 - subfunctions for removing bad channels
#-----------------------------------------------

import numpy as np
from scipy import signal


def removeBadChannels(NIRxData, satlength, QCoDthresh):
    d = NIRxData['d']
    numchannels = NIRxData['numchannels']
    numchannels = d.shape[1]/2
    channelmask = np.ones((1, numchannels))
    nanTooLong = int(round(satlength*samprate)) + 1
    for c in range(numchannels):
        channelwl1=d[:,c]
        channelwl2=d[:,c+numchannels]
        #fix NaN points in time series if there are any
        #if too many in a row, mark as bad channel
        usableChannel = 1
        if np.isnan(channelwl1).sum() > 0:
            for i in range(d.shape[0]):
                subVec = channelwl1[i:i+nanTooLong]
                if np.isnan(subVec).sum() == nanTooLong:
                    usableChannel = 0
                else:
                    ok = ~np.isnan(channelwl1)
                    xp = ok.ravel().nonzero()[0]
                    fp = channelwl1[~np.isnan(channelwl1)]
                    x = np.isnan(channelwl1).ravel().nonzero()[0]
                    channelwl1[np.isnan(channelwl1)] = np.interp(x, xp, fp)
                    d[:,c] = channelwl1

        channelmask[c]=usableChannel
        channelmask[c+numchannels]=usableChannel
        usableChannel = 1
        if np.isnan(channelwl2).sum() > 0:
            for i in range(d.shape[0]):
                subVec = channelwl2[i:i+nanTooLong]
                if np.isnan(subVec).sum() == nanTooLong:
                    usableChannel = 0
                else:
                    ok = ~np.isnan(channelwl2)
                    xp = ok.ravel().nonzero()[0]
                    fp = channelwl2[~np.isnan(channelwl2)]
                    x = np.isnan(channelwl2).ravel().nonzero()[0]
                    channelwl2[np.isnan(channelwl2)] = np.interp(x, xp, fp)
                    d[:, c] = channelwl2

        channelmask[c] = usableChannel
        channelmask[c + numchannels] = usableChannel

        f, psd1 = signal.welch(channelwl1, samprate, nperseg=512)
        psdquarters = int(round(len(psd1)) / 4)
        Q1 = sum(psd1[:psdquarters])
        Q3 = sum(psd1[2 * psdquarters + 1:3 * psdquarters])
        QCoD = (Q1 - Q3) / (Q1 + Q3)
        if QCoD < QCoDthresh:
            channelmask[c]=0

        f, psd2 = signal.welch(channelwl2, samprate, nperseg=512)
        psdquarters = int(round(len(psd2)) / 4)
        Q1 = sum(psd2[:psdquarters])
        Q3 = sum(psd2[2 * psdquarters + 1:3 * psdquarters])
        QCoD = (Q1 - Q3) / (Q1 + Q3)
        if QCoD < QCoDthresh:
            channelmask[c] = 0

    return channelmask