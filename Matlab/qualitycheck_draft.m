dataprefix = 'dyad';
rawdir=uigetdir('','Choose Data Directory');

dyads = 1;

currdir=dir(strcat(rawdir,filesep,dataprefix,'*'));
if length(currdir)<1
    error(['ERROR: No data files found with ',dataprefix,' prefix']);
end

[probefile,probepath] = uigetfile('*_probeInfo.mat','Choose probeInfo File');
load(fullfile(probepath,probefile));
if ~exist('probeInfo','var')
    error('ERROR: Invalid probeInfo file (does not contain a probeInfo object');
end

QCoDthresh = 0.1; %Change to <0.1 to allow more noise, >0.1 for more stringency
subjectmap = zeros(length(currdir),1);
probenumchannels = probeInfo.probes.nChannel0;
channelmap = zeros(probenumchannels,1);

fprintf('\n\t Checking data quality ...\n')
reverseStr = '';
Elapsedtime = tic;
if dyads
    firstdyad=currdir(1).name;
    if isdir(strcat(rawdir,filesep,firstdyad,filesep,'Subject1'))
        for i=1:length(currdir)
            msg = sprintf('\n\t dyad number %d/%d ...',i,length(currdir));
            fprintf([reverseStr,msg]);
            reverseStr = repmat(sprintf('\b'),1,length(msg));      
            
            dyad=currdir(i).name;
            subj1folder = strcat(rawdir,filesep,dyad,filesep,'Subject1');
            subj2folder = strcat(rawdir,filesep,dyad,filesep,'Subject2');
            [d1, ~, samprate1] = extractNIRxData(subj1folder);
            [d2, ~, samprate2] = extractNIRxData(subj2folder);
            numchannels = size(d1,2)/2;
            channelmask = ones(1,numchannels);
            QCoDvector = nan(2,numchannels);
            
            for c=1:numchannels
                meand=d(:,c)-mean(d(:,c));
                [psd1,~]=pwelch(meand,[],[],[],samprate1);
                psdquarters=round(length(psd1)/4);
                Q1 = sum(psd1(1:psdquarters));
                Q4 = sum(psd1(3*psdquarters+1:end));
                QCoD = (Q1-Q4)/(Q1+Q4);
                QCoDvector(1,c)=QCoD;
                if QCoD<QCoDthresh
                    channelmask(1,c)=0;
                end
                meand=d(:,c+numchannels)-mean(d(:,c+numchannels));
                [psd2,fs]=pwelch(meand,[],[],[],samprate1);
                Q1 = sum(psd2(1:psdquarters));
                Q4 = sum(psd2(3*psdquarters+1:end));
                QCoD = (Q1-Q4)/(Q1+Q4);
                QCoDvector(2,c)=QCoD;
                if QCoD<QCoDthresh
                    channelmask(1,c)=0;
                end
            end
            channelmap = sum([channelmap channelmask'],2); 
            subjectmap(i,1) = sum(channelmask);
        end
    else
        firstdyad=currdir(1).name;
        firstdyaddir=dir(strcat(rawdir,filesep,firstdyad,filesep,dataprefix,'*'));
        firstdyadnamelength = length(firstdyad);
        scannames = cell(1,length(firstdyaddir));
        for j=1:length(firstdyaddir)
            scanname = firstdyaddir(j).name(firstdyadnamelength+1:end);
            scannames{j} = scanname;
        end
        channelmaptable = array2table(zeros(probenumchannels,length(firstdyaddir)));
        channelmaptable.Properties.VariableNames= scannames;
        subjectmaptable = array2table(zeros(length(currdir),length(firstdyaddir)));
        subjectmaptable.Properties.VariableNames= scannames;
        for i=1:length(currdir)
            msg = sprintf('\n\t dyad number %d/%d ...',i,length(currdir));
            fprintf([reverseStr,msg]);
            reverseStr = repmat(sprintf('\b'),1,length(msg));      
        
            dyad=currdir(i).name;
            dyadnamelength = length(dyad);
            dyaddir=dir(strcat(rawdir,filesep,dyad,filesep,dataprefix,'*'));
            for j=1:length(dyaddir)
                scanname = dyaddir(j).name;
                subscanname = dyaddir(j).name(dyadnamelength+1:end);
                subj1folder = strcat(rawdir,filesep,dyad,filesep,scanname,filesep,'Subject1');
                subj2folder = strcat(rawdir,filesep,dyad,filesep,scanname,filesep,'Subject2');
                [d1, ~, samprate1] = extractNIRxData(subj1folder);
                [d2, ~, samprate2] = extractNIRxData(subj1folder);
                numchannels = size(d1,2)/2;
                channelmask = ones(1,numchannels);
                QCoDvector = nan(2,numchannels);
            
                for c=1:numchannels
                    meand=d1(:,c)-mean(d1(:,c));
                    [psd1,~]=pwelch(meand,[],[],[],samprate1);
                    psdquarters=round(length(psd1)/4);
                    Q1 = sum(psd1(1:psdquarters));
                    Q4 = sum(psd1(3*psdquarters+1:end));
                    QCoD = (Q1-Q4)/(Q1+Q4);
                    QCoDvector(1,c)=QCoD;
                    if QCoD<QCoDthresh
                        channelmask(1,c)=0;
                    end
                    meand=d1(:,c+numchannels)-mean(d1(:,c+numchannels));
                    [psd2,~]=pwelch(meand,[],[],[],samprate1);
                    Q1 = sum(psd2(1:psdquarters));
                    Q4 = sum(psd2(3*psdquarters+1:end));
                    QCoD = (Q1-Q4)/(Q1+Q4);
                    QCoDvector(2,c)=QCoD;
                    if QCoD<QCoDthresh
                        channelmask(1,c)=0;
                    end
                    meand=d2(:,c)-mean(d2(:,c));
                    [psd1,~]=pwelch(meand,[],[],[],samprate1);
                    psdquarters=round(length(psd1)/4);
                    Q1 = sum(psd1(1:psdquarters));
                    Q4 = sum(psd1(3*psdquarters+1:end));
                    QCoD = (Q1-Q4)/(Q1+Q4);
                    QCoDvector(1,c)=QCoD;
                    if QCoD<QCoDthresh
                        channelmask(1,c)=0;
                    end
                    meand=d2(:,c+numchannels)-mean(d2(:,c+numchannels));
                    [psd2,~]=pwelch(meand,[],[],[],samprate1);
                    Q1 = sum(psd2(1:psdquarters));
                    Q4 = sum(psd2(3*psdquarters+1:end));
                    QCoD = (Q1-Q4)/(Q1+Q4);
                    QCoDvector(2,c)=QCoD;
                    if QCoD<QCoDthresh
                        channelmask(1,c)=0;
                    end
                end
                channelmaptable(:,{subscanname}) = array2table(sum([table2array(channelmaptable(:,{subscanname})) channelmask'],2)); 
                subjectmaptable(i,{subscanname}) = array2table(sum(channelmask));
            end
        end
    end
end
figure();
imagesc(table2array(channelmaptable))
title('Number of dyads with good data in each channel (y) for each scan (x)')
colorbar
set(gca,'XTick',1:5,'XTickLabel',scannames)
figure();
imagesc(table2array(subjectmaptable))
title('Number of good channels in each dyad (y) for each scan (x)')
colorbar
set(gca,'XTick',1:5,'XTickLabel',scannames)
            

 
        