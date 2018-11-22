%-----------------------------------------------
%step 2 - subfunctions for removing bad channels
%-----------------------------------------------
%depends on inpaint_nans function, written by John D'Errico 2012

function channelmask = removeBadChannels(d, samprate, satlength, QCoDthresh)
    numchannels = size(d,2)/2;
    channelmask = ones(1,numchannels);
    for c=1:numchannels
        %fix NaN points in time series if there are any
        %if too many in a row, mark as bad channel
        if sum(isnan(d(:,c)))>0
             if ~isempty(strfind(isnan(d(:,c)), true(1,round(satlength*samprate)))) 
                  channelmask(1,c)=0;
                  channelmask(1,c+numchannels)=0;
             end
             newch = inpaint_nans(d(:,c),4);
             d(:,c)=newch;
        end
        if sum(isnan(d(:,c+numchannels)))>0
             if ~isempty(strfind(isnan(d(:,c+numchannels)), true(1,round(satlength*samprate)))) 
                  channelmask(1,c)=0;
                  channelmask(1,c+numchannels)=0;
             end
             newch = inpaint_nans(d(:,c+numchannels),4);
             d(:,c+numchannels)=newch;
        end

        %check QCoD of each channel
        meand=d(:,c)-mean(d(:,c));
        [psd1,~]=pwelch(meand,[],[],[],samprate);
        psdquarters=round(length(psd1)/4);
        Q1 = sum(psd1(1:psdquarters));
        Q3 = sum(psd1(2*psdquarters+1:3*psdquarters));
        QCoD = (Q1-Q3)/(Q1+Q3);
         if QCoD<QCoDthresh
             channelmask(1,c)=0;
         end
        meand=d(:,c+numchannels)-mean(d(:,c+numchannels));
        [psd2,~]=pwelch(meand,[],[],[],samprate);
        Q1 = sum(psd2(1:psdquarters));
        Q3 = sum(psd2(2*psdquarters+1:3*psdquarters));
        QCoD = (Q1-Q3)/(Q1+Q3);
        if QCoD<QCoDthresh
            channelmask(1,c)=0;
        end
    end 
end
