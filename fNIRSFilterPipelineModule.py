#-----------------------------------------------
#step 4 - our filtering pipeline
#-----------------------------------------------

import numpy as np
from scipy.signal import butter, lfilter
from scipy import stats

def hmrMotionArtifact(d, samprate, tMotion, tMask, stdThresh, channelmask):
    channelmask=np.concatenate((channelmask,channelmask))
    channelmask=channelmask.flatten()
    channelmask.astype(bool)
    tInclude = np.ones((d.shape[0],1))
    artifactBuffer = int(round(tMask*samprate))
    stddiff = np.std((d[1:,channelmask] - d[:-1, channelmask]), axis=0)
    maxdiff = np.zeros((d.shape[0]-1, channelmask.shape[0]))
    for ii in range(int(round(tMotion*samprate))):
        maxdiff = np.maximum([abs(d[ii+1:, channelmask] - d[:-1, channelmask])],maxdiff)

    maxdiff=np.squeeze(maxdiff)
    mcThresh = stddiff*stdThresh
    mcThresh = np.array(mcThresh)[np.newaxis]
    badInds = maxdiff > (np.ones(maxdiff.shape[0])*np.transpose(mcThresh)).T
    badInds = badInds.max(axis=1)
    badInds = np.where(badInds)[0]


    if np.any(badInds):
        badInds = np.matlib.repmat(np.reshape(badInds, (badInds.shape[0],1), 1, 2*artifactBuffer+1)) + np.matlib.repmat(np.arange(-artifactBuffer,artifactBuffer+1), len(badInds),1)
        badInds = badInds[badInds[:]>0]
        badInds = badInds[badInds[:]<= d.shape[0]-1]
        tInclude[1+badInds]=0

    return tInclude, channelmask


def hmrMotionCorrectPCA(d, tInclude, nSV, channelmask):
    lstNoInc = np.where(tInclude==0)[0]
    if len(lstNoInc)==0:
        dN = d;
        svs = [];
        nSV = 0;
        return dN, svs, nSV


    y = d[lstNoInc[:,None],np.where(channelmask)[0]]
    yc = y
    yo = y

    c = np.matmul(y,y.T)
    [V, St, foo] = np.linalg.svd(c)
    svs = np.diagflat(St) / np.sum(np.diagflat(St))
    svs = np.diag(svs)

    svsc = svs
    for idx in range(1,svs.shape[0]):
        svsc[idx] = svsc[idx-1] + svs[idx]

    ev = np.zeros((svs.shape[0],1))
    ev[:nSV] = 1
    ev = np.diagflat(ev)

    yc = yo - np.matmul(np.matmul(np.matmul(y,V),ev),V.T)

    lstMs = np.where(np.diff(tInclude, axis=0)==-1)
    lstMf = np.where(np.diff(tInclude, axis=0)==1)
    if len(lstMf[0])==0:
        lstMf = len(tInclude)
    if len(lstMs[0])==0:
        lstMs = 1

    if lstMs[0][0]>lstMf[0][0]:
        lstMs = np.array([[1],[lstMs]])
    if lstMs[0][-1]>lstMf[0][-1]:
        np.append(lstMf,len(tInclude))

    lstMb = lstMf[0]-lstMs[0]
    for ii in range(1,len(lstMb)):
        lstMb[ii] = lstMb[ii-1] + lstMb[ii]

    dN = d

    for ii in range (len(channelmask)):
        jj = np.where(channelmask)[0][ii]
        lst = np.arange(lstMs[0][0],lstMf[0][0])
        if lstMs[0][0]>1:
            dN[lst,jj] = yc[:lstMb[0],ii] - yc[0,ii] + dN[lst[0],jj]
        else:
            dN[lst,jj] = yc[:lstMb[0],ii] - yc[lstMb[0],ii] + dN[lst[-1],jj]

        for kk in range(len(lstMf[0])-1):
            lst = np.arange(lstMf[0][kk],lstMs[0][kk+1])
            dN[lst,jj] = d[lst,jj] - d[lst[0],jj] + dN[lst[0],jj]
            lst = np.arange(lstMs[0][kk+1],lstMf[0][kk+1]-1)
            dN[lst,jj] = yc[lstMb[kk]+1:lstMb[kk+1],ii] - yc[lstMb[kk]+1,ii] + dN[lst[0],jj]

        if lstMf[-1]<len(d):
            lst = np.arange(lstMf[0][-1]-1,d.shape[0])
            dN[lst,jj] = d[lst,jj] - d[lst[0],jj] + dN[lst[0], jj]

    return dN, svs, nSV


def bandpass(lpf, hpf, samprate, order=3):
    nyq = 0.5 * samprate
    low = lpf / nyq
    high = hpf / nyq
    b, a = butter(order, [low,high], btype='band')
    return b,a

def bandpassfilt(d, lpf, hpf, samprate, order=3)
    b, a = bandpass(lpf, hpf, samprate, order=order)
    y = lfilter(b, a, d)
    return y

def hmrOD2Conc(d, channelmask, NIRxData):
    numchannels = NIRxData["numchannels"]
    rho = np.array([2.75294831343702,2.59391909620926,2.89390643929921,3.16382322961505,2.62110739729229,2.96172632777674,4.02461249099262,3.28059438345427,4.31176280152085,2.63190764322845,3.08764730907214,2.95608363013052,3.68502671633002,2.56198644288384,3.87389418414280,2.81333836236465,4.68451444199667,4.03947524760419,3.07902875414238,3.17981701964377,2.66445419093828,3.28545594196557,3.18657921386250,3.13663592212414,3.37508719234535,2.83028658464260,2.88482679589375,2.79984685563647,3.65855234981939,3.14913673024192,3.12014715297748,3.42142126543336,3.45008367302570,3.13518319465310,3.15575545703306,3.03938416302363,3.89959747207689,2.70004596591252,3.16815655225146,3.65038524218728,3.39040723363248,2.97196610140026,3.25809998242152,3.46681906555635,3.56571275183508,3.09993024732817,3.16793438087845,3.44560997911601,3.54847655824251,3.79541579326342,2.95548593251323,3.07691677374201,3.04685735288605,2.92218358424758,3.08326627716818,3.36150615411623,3.09732754125042,4.00067912079507,3.07365664932104,3.65964913671669,3.27312919173483,3.19826407937980,3.22496017356654,2.81842556000572,3.80476039695552,3.22916547969144,3.12223843081221,4.02193571565341,3.42735381304277,3.47063482450417,3.41890782130324,3.37461180799098,3.84184129209686,2.94832335985654,3.24329036146342,3.21906248929818,4.02411531296331,3.82763140632397,2.62146145726230,3.23045561538537,2.87596897982358,3.50492594817468,3.16303581351615,3.20377067400470,3.67247430919941,2.99739473474148,3.18679351803099,3.93055230249248,2.85339588972147,3.26121026985840,3.09512239745092,4.55205294886391,4.17585440283239,3.98761339137248,3.50021427915492,4.09472480965336,4.06180343884522,3.98288373609679,3.36189172936905,3.81180505423535,2.13836386052515,3.50957262355404,4.07305723583142,2.99062912860228,2.87721015402533,3.59465890242625,2.84181253657155,1.97004187259817])

    nTpts = d.shape[0]
    dc= []
    einv = np.array([[-0.000255602763275000,0.000546224086947472],[0.000359022063140430,-0.000211256829313717]])
    for idx in range (len(channelmask)/2)
        idx1 = np.where(channelmask==1)[0][idx]
        x = np.ones((nTpts,1))*rho[idx]*[6,6]
        dc[:,:,idx] = np.divide((np.matmul(d[:,[idx,idx+numchannels]],einv)),x)

    dc[:, 3,:] = dc[:, 1,:] + dc[:, 2,:]
    return dc


def hmrIntensity2Conc(d, samprate, hpf, lpf, channelmask, NIRxData):
    ppf = [6,6]
    dm = np.mean(abs(d), axis=0)
    nTpts = d.shape[0]
    dod = np.divide(abs(d),(np.ones((nTpts,1)))*dm)
    dod=-np.log(dod)

    dod = bandpassfilt(dod,samprate,hpf,lpf)

    dc = hmrOD2Conc(d, channelmask, NIRxData)

    return dc


def fNIRSFilterPipeline(NIRxData, channelmask)
    # -----------------------------------------------
    # step 4 - our filtering pipeline
    # -----------------------------------------------
    samprate = NIRxData['samprate']
    d = NIRxData['d']
    tInclude = hmrMotionArtifact(d,samprate,0.5, 2, 10, 5)
    dfiltered,svs, nSV = hmrMotionCorrectPCA(d, tInclude, 2, channelmask)
    dconverted = hmrIntensity2Conc(dfiltered,samprate,0.005,0.5,channelmask,NIRxData)
    dnormed = stats.zscore(dconverted)
    dnormed[:,np.where(channelmask==1)[0]] = dnormed
