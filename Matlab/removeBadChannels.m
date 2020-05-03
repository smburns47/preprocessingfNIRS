%-----------------------------------------------
%step 2 - subfunctions for removing bad channels
%-----------------------------------------------
%depends on inpaint_nans function, written by John D'Errico 2012

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
