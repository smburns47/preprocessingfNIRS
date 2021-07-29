%-----------------------------------------------
%step 4 - our filtering pipeline
%-----------------------------------------------
 
function [dconverted, dnormed]= fNIRSFilterPipeline(d, SD, samprate, coords)
    %depends on Homer2 package
    warning off; %sometimes hmrIntensity2Conc gives a warning we don't care about here
    dod = hmrIntensity2OD( d );
    if ~isempty(coords)
        dod = spatialPCFilter(dod, coords);
    end
    dod = hmrBandpassFilt(dod,samprate,0.008,0.2);
    ppf = [6 6];
    if length(ppf)~=length(SD.Lambda)
        disp( sprintf('Warning: ppf must be a vector with the same number of elements as SD.Lambda\n    Using default values of [6 6]') );
        ppf = [6 6];
    end
    dconverted = hmrOD2Conc( dod, SD, ppf );
    
    %see hmrMotionArtifact in Homer2 documentation for parameter description
    %[~,tIncCh] = hmrMotionArtifactByChannel(dod, samprate, SD, ones(length(d),1), 1, 1, 5, 2);
    
    %4/21/2020
    %Homer2 motionPCA filtering has not worked very well on our recent
    %data, so reverting to a simpler method of correction for
    %discontinutites and changes in volatility, and then generating report
    %for which channels might still have problematic artifacts after processing
    %while I work on learning Homer3. Blog description coming soon.
    %dfiltered = BaselineVolatilityCorrection(dod, samprate, SD, tIncCh);

    %see in Homer2 documentation for parameter description
    
    %[dconverted, ~] = hmrIntensity2Conc(dfiltered, SD, samprate, 0.008, 0.2, [6, 6]);
    dnormed = zscore(dconverted);
    warning on;
end