function channelmask = testCoV(testsubjectpath, QCoDthresh, suppressPlot)   
%       plots the power spectral density graphs of the given test subject, 
%       as well as the automatic bad channel detection to make sure your 
%       chosen CoV threshold is performing as you want it to. Does NOT do any 
%       other preprocessing, so channels might be saturated. It's
%       recommended that you do this with several subject folders first
%       before the full preprocessing job.

%INPUTS: testsubjpath: relative path to the folder with all raw nirx files
%        CoVthresh: CoV threshold to test out (0.5 is default)
%        suppressPlot: 0 or 1. 0 to let plots happen, 1 to keep them from
%        coming up.
%
%OUTPUTS: the channelmask

    [d, ~, samprate] = extractData(testsubjectpath);
    numchannels = size(d,2)/2;
    channelmask = ones(1,numchannels); 
    
    %check CoV of each channel
    for c=1:numchannels
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
        [psd2,fs]=pwelch(meand,[],[],[],samprate);
        Q1 = sum(psd2(1:psdquarters));
        Q3 = sum(psd2(2*psdquarters+1:3*psdquarters));
        QCoD = (Q1-Q3)/(Q1+Q3);
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

%-----------------------------------------------
%step 1 - subfunctions for extracting data
%-----------------------------------------------

function [d, sd_ind, samprate, wavelengths, s] = extractData(path)
    wl1file = dir(strcat(path,'/','*.wl1'));
    wl2file = dir(strcat(path,'/','*.wl2'));
    hdrfile = dir(strcat(path,'/','*.hdr'));
    wl1filename = wl1file(1).name;
    wl2filename = wl2file(1).name;
    hdrfilename = hdrfile(1).name;
    wl1=load(strcat(path,'/',wl1filename));
    wl2=load(strcat(path,'/',wl2filename));
    d=[wl1 wl2];

    fid = fopen(strcat(path,'/',hdrfilename));
    tmp = textscan(fid,'%s','delimiter','\n');%This just reads every line
    hdrString = tmp{1};
    fclose(fid);

    samprate = readSamplingRate(hdrString);
    sd_ind = getSDMask(hdrString);
    wavelengths = readWavelengths(hdrString);
    s = readEvents(hdrString, d);
    d = d(:,sd_ind);        
end

function samprate = readSamplingRate(hdrString)
    keyword = 'SamplingRate=';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
    tmp = hdrString{ind};
    samprate = str2num(tmp(length(keyword)+1:end));
end

function sd_ind = getSDMask(hdrString)
    keyword = 'S-D-Mask="#';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)) + 1; %This gives cell of hdr_str with keyword
    tmp = strfind(hdrString(ind+1:end),'#');
    ind2 = find(~cellfun(@isempty,tmp)) - 1;
    ind2 = ind + ind2(1);
    sd_ind = cell2mat(cellfun(@str2num,hdrString(ind:ind2),'UniformOutput',0));
    sd_ind = sd_ind';
    sd_ind = find([sd_ind(:);sd_ind(:)]);
end

function wavelengths = readWavelengths(hdrString)
    keyword = 'Wavelengths=';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
    tmp = hdrString{ind};
    Wavelength1 = str2double(tmp(length(keyword)+2:length(keyword)+5));
    Wavelength2 = str2double(tmp(end-4:end-1));
    wavelengths = [Wavelength1 Wavelength2];
end

function s = readEvents(hdrString, d)
    keyword = 'Events="#';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)) + 1; %This gives cell of hdr_str with keyword
    tmp = strfind(hdrString(ind+1:end),'#');
    ind2 = find(~cellfun(@isempty,tmp)) - 1;
    ind2 = ind + ind2(1);
    events = cell2mat(cellfun(@str2num,hdrString(ind:ind2),'UniformOutput',0));
    if ~isempty(events)
        events = events(:,2:3);
        markertypes = unique(events(:,1));
        s = zeros(length(d),length(markertypes));
    else
        s = zeros(length(d),1);
    end
end
