%-----------------------------------------------
%step 4 - our filtering pipeline
%-----------------------------------------------
 
function [oxy, deoxy, totaloxy, z_oxy, z_deoxy, z_totaloxy]= fNIRSFilterPipeline(d, SD, samprate)
    %depends on Homer2 package
    warning off; %sometimes hmrIntensity2Conc gives a warning we don't care about here
    %see hmrMotionArtifact in Homer2 documentation for parameter description
    numchannels = size(d,2)/2;
    tInc = hmrMotionArtifact(d, samprate, SD, ones(length(d)), 0.5, 2, 10, 5);
    %see hmrMotionCorrectPCA in Homer2 documentation for parameter description
    [dfiltered,~,~] = hmrMotionCorrectPCA(SD, d, tInc, 2);
    %see hmrIntensity2Conc in Homer2 documentation for parameter description
    [dconverted, ~] = hmrIntensity2Conc(dfiltered, SD, samprate, 0.005, 0.5, [6, 6]);
    dnormed = zscore(dconverted);
    dnormed(:,:,(SD.MeasListAct')==0)=NaN;
    dconverted(:,:,(SD.MeasListAct')==0)=NaN;
    oxy = zeros(length(dconverted), numchannels);
    deoxy = zeros(length(dconverted), numchannels);
    totaloxy = zeros(length(dconverted), numchannels);
    z_oxy = zeros(length(dnormed), numchannels);
    z_deoxy = zeros(length(dnormed), numchannels);
    z_totaloxy = zeros(length(dnormed), numchannels);
    for c = 1:numchannels
        oxy(:,c) = dconverted(:,1,c);
        deoxy(:,c) = dconverted(:,2,c);
        totaloxy(:,c) = dconverted(:,3,c);
        z_oxy(:,c) = dnormed(:,1,c);
        z_deoxy(:,c) = dnormed(:,2,c);
        z_totaloxy(:,c) = dnormed(:,3,c);
    end
    warning on;
end