%-----------------------------------------------
%step 2 - subfunctions for removing bad channels
%-----------------------------------------------
%depends on inpaint_nans function, written by John D'Errico 2012

%bad channel defined as any where detector saturation occurs for >2sec, 
%or where power spectrum resembles white noise. 
%Feel free to change these parameters if you have a good reason to do so
%
%reasoning for default choices:
%- if saturation occurs, data will be 'NaN'. But if this only lasts a
%short amount of time (e.g. <8 points=<2 seconds at 4Hz), we can fill in what 
%those data points would have likely been with reasonable confidence.
%
%- power spectrum of the signal shows how many sine waves at each freq.
%make up the raw signal. Good signal should have a large peak at lower
%frequencies. Pure noise will have random numbers of all freqencies. 
%We will use a modified version of the quartile coefficient of
%dispersion
%(https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion)
%to automatically decide which channels have good or bad
%signal. Essentially, it sums the frequency amplitudes in the
%first and fourth quartiles of the frequency range, and then
%compares them via (Q1-Q4)/(Q1+Q4). Larger QCoD is cleaner
%signal. Default threshold is set to 0.6-0.03*samprate based on 
%simulations of how much signal can be recovered at r=0.85 from imposed 
%noise (preprint in prep). Change first parameter to <0.6 to allow more 
%noise in the signal, or change to >0.6 for more stringency.  

function [d, channelmask] = removeBadChannels(d, samprate, satlength, QCoDthresh)
    numchannels = size(d,2)/2;    
    channelmask = ones(1,numchannels);
    for c=1:size(d,2)
        %fix NaN points in time series if there are any
        %if too many in a row, mark as bad channel
        if sum(isnan(d(:,c)))>0
             if ~isempty(strfind(isnan(d(:,c))', true(1,round(satlength*samprate)))) 
                  channelmask(1,c)=0;
             end
             newch = inpaint_nans(d(:,c),4);
             d(:,c)=newch;
        end
    end
    
    for c=1:numchannels
        %check QCoD of each channel
        meand=d(:,c)-mean(d(:,c));
        [psd1,~]=pwelch(meand,[],[],[],samprate);
        psdquarters=round(length(psd1)/4);
        Q1 = sum(psd1(1:psdquarters));
        Q4 = sum(psd1(3*psdquarters+1:end));
        QCoD = (Q1-Q4)/(Q1+Q4);
         if QCoD<QCoDthresh
             channelmask(1,c)=0;
         end
        meand=d(:,c+numchannels)-mean(d(:,c+numchannels));
        [psd2,~]=pwelch(meand,[],[],[],samprate);
        Q1 = sum(psd2(1:psdquarters));
        Q4 = sum(psd2(3*psdquarters+1:end));
        QCoD = (Q1-Q4)/(Q1+Q4);
        if QCoD<QCoDthresh
            channelmask(1,c)=0;
        end
    end 
end
