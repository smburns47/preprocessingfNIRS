function [channelmask, QCoDvector] = testQCoD(testsubjectpath, QCoDthresh, suppressPlot)   
%       plots the power spectral density graphs of the given test subject, 
%       as well as the automatic bad channel detection to make sure your 
%       chosen CoV threshold is performing as you want it to. Does NOT do any 
%       other preprocessing, so channels might be saturated. It's
%       recommended that you do this with several subject folders first
%       before the full preprocessing job.

%INPUTS: testsubjpath: relative path to one scan folder with all raw nirx files
%        QCoDthresh: QCoD threshold to test out
%        suppressPlot: 0 or 1. 0 to let plots happen, 1 to keep them from
%        coming up.
%
%OUTPUTS: the channelmask and vector of QCoD values

    [d, ~, samprate] = extractNIRxData(testsubjectpath);
    numchannels = size(d,2)/2;
    channelmask = ones(1,numchannels);
    QCoDvector = nan(2,numchannels);
    
    %check CoV of each channel
    for c=1:numchannels
        meand=d(:,c)-mean(d(:,c));
        [psd1,~]=pwelch(meand,[],[],[],samprate);
        psdquarters=round(length(psd1)/4);
        Q1 = sum(psd1(1:psdquarters));
        Q3 = sum(psd1(2*psdquarters+1:3*psdquarters));
        QCoD = (Q1-Q3)/(Q1+Q3);
        QCoDvector(1,c)=QCoD;
         if QCoD<QCoDthresh
             channelmask(1,c)=0;
         end
        meand=d(:,c+numchannels)-mean(d(:,c+numchannels));
        [psd2,fs]=pwelch(meand,[],[],[],samprate);
        Q1 = sum(psd2(1:psdquarters));
        Q3 = sum(psd2(2*psdquarters+1:3*psdquarters));
        QCoD = (Q1-Q3)/(Q1+Q3);
        QCoDvector(2,c)=QCoD;
        if QCoD<QCoDthresh
            channelmask(1,c)=0;
        end

        if ~suppressPlot
            figure()
            plot(fs,psd1)
            title(strcat('Channel- ', num2str(c)))
            hold on
            plot(fs,psd2, 'r')
        end
    end
end
