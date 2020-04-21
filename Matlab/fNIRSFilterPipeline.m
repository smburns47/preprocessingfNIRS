%-----------------------------------------------
%step 4 - our filtering pipeline
%-----------------------------------------------
 
function [new_d, oxy, deoxy, totaloxy, z_oxy, z_deoxy, z_totaloxy]= fNIRSFilterPipeline(d, SD, samprate)
    %depends on Homer2 package
    warning off; %sometimes hmrIntensity2Conc gives a warning we don't care about here
    %see hmrMotionArtifact in Homer2 documentation for parameter description
    numchannels = size(d,2)/2;
    d(:,SD.MeasListAct==0)=NaN;
    [~,tIncCh] = hmrMotionArtifactByChannel(d, samprate, SD, ones(length(d),1), 1, 1, 5, 2);
    
    %4/21/2020
    %Homer2 motionPCA filtering has not worked very well on our recent
    %data, so reverting to a simpler method of correction for
    %discontinutites and changes in volatility, and then generating report
    %for which channels might still have problematic artifacts after processing
    %while I work on learning Homer3. Blog description coming soon.
    dfiltered = BaselineVolatilityCorrection(d, SD, tIncCh);

    %see hmrIntensity2Conc in Homer2 documentation for parameter description
    [dconverted, ~] = hmrIntensity2Conc(dfiltered, SD, samprate, 0.008, 0.2, [6, 6]);
    dnormed = zscore(dconverted);
    oxy = zeros(size(dconverted,1), numchannels);
    deoxy = zeros(size(dconverted,1), numchannels);
    totaloxy = zeros(size(dconverted,1), numchannels);
    z_oxy = zeros(size(dnormed,1), numchannels);
    z_deoxy = zeros(size(dnormed,1), numchannels);
    z_totaloxy = zeros(size(dnormed,1), numchannels);
    new_d = zeros(size(dconverted,1), numchannels*2);
    for c = 1:numchannels
        oxy(:,c) = dconverted(:,1,c);
        deoxy(:,c) = dconverted(:,2,c);
        totaloxy(:,c) = dconverted(:,3,c);
        z_oxy(:,c) = dnormed(:,1,c);
        z_deoxy(:,c) = dnormed(:,2,c);
        z_totaloxy(:,c) = dnormed(:,3,c);
        new_d(:,(c*2)-1) = oxy(:,c);
        new_d(:,c*2) = deoxy(:,c);
    end
    warning on;
end