%-----------------------------------------------
%step 1 - subfunctions for extracting data
%-----------------------------------------------

function [d, sd_ind, samprate, wavelengths, s] = extractNIRxData(subjfolder)
    wl1file = dir(strcat(subjfolder,'/','*.wl1'));
    wl2file = dir(strcat(subjfolder,'/','*.wl2'));
    hdrfile = dir(strcat(subjfolder,'/','*.hdr'));
    wl1filename = wl1file(1).name;
    wl2filename = wl2file(1).name;
    hdrfilename = hdrfile(1).name;
    wl1=load(strcat(subjfolder,'/',wl1filename));
    wl2=load(strcat(subjfolder,'/',wl2filename));
    d=[wl1 wl2];

    fid = fopen(strcat(subjfolder,'/',hdrfilename));
    tmp = textscan(fid,'%s','delimiter','\n');%This just reads every line
    hdrString = tmp{1};
    fclose(fid);

    samprate = getSamplingRate(hdrString);
    sd_ind = getSDMask(hdrString);
    wavelengths = getWavelengths(hdrString);
    d = d(:,sd_ind);
    s = getEvents(hdrString, d);
end

function samprate = getSamplingRate(hdrString)
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

function wavelengths = getWavelengths(hdrString)
    keyword = 'Wavelengths=';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
    tmp = hdrString{ind};
    Wavelength1 = str2double(tmp(length(keyword)+2:length(keyword)+5));
    Wavelength2 = str2double(tmp(end-4:end-1));
    wavelengths = [Wavelength1 Wavelength2];
end

function s = getEvents(hdrString, d)
    keyword = 'Events="#';
    tmp = strfind(hdrString,keyword);
    ind = find(~cellfun(@isempty,tmp)) + 1; %This gives cell of hdr_str with keyword
    tmp = strfind(hdrString(ind+1:end),'#');
    ind2 = find(~cellfun(@isempty,tmp)) - 1;
    ind2 = ind + ind2(1);
    events = cell2mat(cellfun(@str2num,hdrString(ind:ind2),'UniformOutput',0));
    if ~isempty(events)
        events = events(:,2:3);
        events = unique(events,'rows');
        events = events(events(:,2)~=0,:);
        markertypes = unique(events(:,1));
        s = zeros(size(d,1),length(markertypes));
        for j=1:length(markertypes)
            s(events(:,2),j)=1;
        end
    else
        s = zeros(size(d,1),1);
    end
end