
%path must be in brain_data folder (where your raw data folders are)
%required inpaint_nans, homer2 scripts 

function preprocessingfNIRS(dataprefix, dyads)
%inputs: 
%       dataprefix: string. Prefix of every folder name that should be considered a
%       data folder. E.g., ST for ST_101, ST_102, etc.  
%       dyads: 0 or 1. 1 if hyperscanning, 0 if single subject. 
%
%outputs: preprocessed and .nirs files in a new folder in rawdir called
%           'PreProcessedFiles', sorted by subject

%if you get 'WARNING: Some data points in d are zero...' this is ok.
%this would normally indicate noise in our data, but since we're doing
%motion correction before filtering, our motion correction algorithm might
%force the data into negative numbers while still being good data. you can
%ignore this.

%DEBUGGING TIPS:
%Note that this function is, well, funcitonal, but not 100% optimized for
%everything that could go wrong. A couple problems I have noticed so far,
%and how to fix them if it comes up for you:
%   - If you get an "Index exceeds matrix dimensions" error in
%   hmrMotionArtifact for a subject that's not the first file:
%       Check the SD Mask structure in the .hdr of that subject to see if 
%       it matches the channel structure of the selected probeInfo file. If
%       the wrong probeInfo file was chosen, this will throw the error.
%       also happens if the wrong montage was selected in recording, Simply copy-paste
%       the correct SD Mask and ChannelDistance list into the .hdr file from a
%       subject's .hdr file that had the correct montage.
%   - Delete all false start files from the data directory, or will cause
%   script to error out. 

% added_path = [pwd,'/utils'];
% addpath(added_path);

rawdir=uigetdir('','Choose Data Directory');

currdir=dir(strcat(rawdir,filesep,dataprefix,'*'));
if length(currdir)<1
    error(['ERROR: No data files found with ',dataprefix,' prefix']);
end

supported_devices = {'NIRx-NirScout','NIRx-NirSport1','NIRx-NirSport2','TechEn'};
[device,~] = listdlg('PromptString', 'Select acquisition device:',...
    'SelectionMode', 'single', 'ListString', supported_devices);
if device <= 2
    device=1;
elseif device >= 3
    device=2;
end

if device==1
    [probefile,probepath] = uigetfile('*_probeInfo.mat','Choose probeInfo File');
    load(fullfile(probepath,probefile));
    if ~exist('probeInfo','var')
        error('ERROR: Invalid probeInfo file (does not contain a probeInfo object');
    end
end
multiscan=0;
anytofilter=0;


fprintf('\n\t Preprocessing ...\n')
reverseStr = '';
Elapsedtime = tic;

if dyads
    for i=1:length(currdir)
        dyad=currdir(i).name;
        msg = sprintf('\n\t dyad number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));      
        if exist(strcat(rawdir,filesep,dyad,filesep,'Subject1'),'dir')
            subj1folder = strcat(rawdir,filesep,dyad,filesep,'Subject1');
            subj2folder = strcat(rawdir,filesep,dyad,filesep,'Subject2');
            
            outpath = strcat(rawdir,filesep,'PreProcessedFiles',filesep,dyad);
            if device==2
                outpath = outpath(1:end-5);
            end
            if ~exist(outpath,'dir')
            
                %1) extract data values
                anytofilter=1;
                if device==1
                    [d1, sd_ind1, ~, ~, s1] = extractNIRxData(subj1folder);
                    [d2, sd_ind2, samprate, wavelengths, s2] = extractNIRxData(subj2folder);
                    probenumchannels = probeInfo.probes.nChannel0;
                    datanumchannels = size(d1,2)/2;
                    if probenumchannels~=datanumchannels
                        error('ERROR: number of data channels in subj1 hdr file does not match number of channels in probeInfo file.');
                    end
                    probenumchannels = probeInfo.probes.nChannel0;
                    datanumchannels = size(d2,2)/2;
                    if probenumchannels~=datanumchannels
                        error('ERROR: number of data channels in subj2 hdr file does not match number of channels in probeInfo file.');
                    end
                elseif device==2
                    [d1, samprate, s1, SD1, aux, t] = extractTechEnData(subj1folder);
                    [d2, ~, s2, SD2, ~, ~] = extractTechEnData(subj2folder);
                end
                
                %2) Trim beginning of data to 10s before onset, if there is
                %a lot of dead time before that 
                ssum = sum(s1,2);
                stimmarks = find(ssum);
                if length(stimmarks)>=1
                    begintime = stimmarks(1) - round(samprate*10);
                    if begintime>0
                        d1 = d1(begintime:end,:);
                        d2 = d2(begintime:end,:);
                        s1 = s1(begintime:end,:);
                        s2 = s2(begintime:end,:);
                        stimmarks = stimmarks-begintime;
                    end
                end
                
                %trim off last ten seconds, to remove edge artifacts that
                %might mess up hemodynamics calculation
                d1 = d1((1:end-round(10*samprate)),:);
                d2 = d2((1:end-round(10*samprate)),:);
                s1 = s1((1:end-round(10*samprate)),:);
                s2 = s2((1:end-round(10*samprate)),:);

                %3) identify and remove bad channels
                %bad channel defined as any where detector saturation occurs for >2sec, 
                %or where power spectrum variation is too high. 
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
                %signal. Default threshold is set to 0.1. Change this to <0.1 
                %to allow for greater noise in the signal, or change to >0.1 
                %for more stringency.  
    
                satlength = 2; %in seconds
                QCoDthresh = 0.1;
                channelmask1 = removeBadChannels(d1, samprate, satlength, QCoDthresh);
                channelmask2 = removeBadChannels(d2, samprate, satlength, QCoDthresh);
    
                %4) convert to .nirs format
                if device==1
                    [SD1, ~, ~] = getMiscNirsVars(d1, sd_ind1, samprate, wavelengths, probeInfo, channelmask1);
                    [SD2, aux, t] = getMiscNirsVars(d2, sd_ind2, samprate, wavelengths, probeInfo, channelmask2);
                elseif device==2
                    SD1.MeasListAct = [channelmask1'; channelmask1'];
                    SD1.MeasListVis = SD1.MeasListAct;
                    SD2.MeasListAct = [channelmask2'; channelmask2'];
                    SD2.MeasListVis = SD2.MeasListAct;
                end
                if length(stimmarks)>=1
                    if begintime>0
                        aux = aux(begintime:end,:);
                        t = t(begintime:end);
                    end
                end
                aux = aux((1:end-round(10*samprate)),:);
                t = t((1:end-round(10*samprate)),:);
            
                %5) motion filter, convert to hemodynamic changes
                [new_d1, oxy1, deoxy1, totaloxy1, z_oxy1, z_deoxy1, z_totaloxy1] = fNIRSFilterPipeline(d1, SD1, samprate);
                [new_d2, oxy2, deoxy2, totaloxy2, z_oxy2, z_deoxy2, z_totaloxy2] = fNIRSFilterPipeline(d2, SD2, samprate);

                mkdir(outpath)
                save(strcat(outpath,filesep,dyad,'_subj1_preprocessed.mat'),'oxy1', 'deoxy1', 'totaloxy1','z_oxy1', 'z_deoxy1', 'z_totaloxy1','stimmarks');
                save(strcat(outpath,filesep,dyad,'_subj2_preprocessed.mat'),'oxy2', 'deoxy2', 'totaloxy2','z_oxy2', 'z_deoxy2', 'z_totaloxy2','stimmarks');
                SD=SD1;
                d=new_d1;
                s=s1;
                save(strcat(outpath,'_subj1.nirs'),'aux','d','s','SD','t');
                SD=SD2;
                d=new_d2;
                s=s2;
                save(strcat(outpath,'_subj2.nirs'),'aux','d','s','SD','t');
            end
        else
            multiscan=1;
            dyaddir=dir(strcat(rawdir,filesep,dyad,filesep,dataprefix,'*'));
            for j=1:length(dyaddir)
                scanname = dyaddir(j).name;
                subj1folder = strcat(rawdir,filesep,dyad,filesep,scanname,filesep,'Subject1');
                subj2folder = strcat(rawdir,filesep,dyad,filesep,scanname,filesep,'Subject2');
 
                outpath = strcat(rawdir,filesep,'PreProcessedFiles',filesep,dyad,filesep,scanname);
                if device==2
                    outpath = outpath(1:end-5);
                end
                if ~exist(outpath,'dir')
            
                %1) extract data values
                    anytofilter=1;
                    if device==1
                        [d1, sd_ind1, ~, ~, s1] = extractNIRxData(subj1folder);
                        [d2, sd_ind2, samprate, wavelengths, s2] = extractNIRxData(subj2folder);
                        probenumchannels = probeInfo.probes.nChannel0;
                        datanumchannels = size(d1,2)/2;
                        if probenumchannels~=datanumchannels
                            error('ERROR: number of data channels in subj1 hdr file does not match number of channels in probeInfo file.');
                        end
                        probenumchannels = probeInfo.probes.nChannel0;
                        datanumchannels = size(d2,2)/2;
                        if probenumchannels~=datanumchannels
                            error('ERROR: number of data channels in subj2 hdr file does not match number of channels in probeInfo file.');
                        end
                    elseif device==2
                        [d1, samprate, s1, SD1, aux, t] = extractTechEnData(subj1folder);
                        [d2, ~, s2, SD2, ~, ~] = extractTechEnData(subj2folder);
                    end
    
                    %2) Trim beginning of data to 10s before onset, if there is
                    %a lot of dead time before that 
                    ssum = sum(s1,2);
                    stimmarks = find(ssum);
                    if length(stimmarks)>=1
                        begintime = stimmarks(1) - round(samprate*10);
                        if begintime>0
                            d1 = d1(begintime:end,:);
                            d2 = d2(begintime:end,:);
                            s1 = s1(begintime:end,:);
                            s2 = s2(begintime:end,:);
                            stimmarks = stimmarks-begintime;
                        end
                    end
                    
                    %trim off last ten seconds, to remove edge artifacts that
                    %might mess up hemodynamics calculation
                    d1 = d1((1:end-round(10*samprate)),:);
                    d2 = d2((1:end-round(10*samprate)),:);
                    s1 = s1((1:end-round(10*samprate)),:);
                    s2 = s2((1:end-round(10*samprate)),:);
                    
                    %3) identify and remove bad channels
                    satlength = 2; %in seconds
                    QCoDthresh = 0.1;
                    channelmask1 = removeBadChannels(d1, samprate, satlength, QCoDthresh);
                    channelmask2 = removeBadChannels(d2, samprate, satlength, QCoDthresh);
                    
                    %4)Prep rest of .nirs format
                    if device==1
                        [SD1, ~, ~] = getMiscNirsVars(d1, sd_ind1, samprate, wavelengths, probeInfo, channelmask1);
                        [SD2, aux, t] = getMiscNirsVars(d2, sd_ind2, samprate, wavelengths, probeInfo, channelmask2);
                    elseif device==2
                        SD1.MeasListAct = [channelmask1'; channelmask1'];
                        SD1.MeasListVis = SD1.MeasListAct;
                        SD2.MeasListAct = [channelmask2'; channelmask2'];
                        SD2.MeasListVis = SD2.MeasListAct;
                    end
                    if length(stimmarks)>=1
                        if begintime>0
                            aux = aux(begintime:end,:);
                            t = t(begintime:end);
                        end
                    end
                    aux = aux((1:end-round(10*samprate)),:);
                    t = t((1:end-round(10*samprate)),:);
                    
                    %5) motion filter, convert to hemodynamic changes
                    [new_d1, oxy1, deoxy1, totaloxy1, z_oxy1, z_deoxy1, z_totaloxy1] = fNIRSFilterPipeline(d1, SD1, samprate);
                    [new_d2, oxy2, deoxy2, totaloxy2, z_oxy2, z_deoxy2, z_totaloxy2] = fNIRSFilterPipeline(d2, SD2, samprate);

                    mkdir(outpath)
                    save(strcat(outpath,filesep,scanname,'_subj1_preprocessed.mat'),'oxy1', 'deoxy1', 'totaloxy1','z_oxy1', 'z_deoxy1', 'z_totaloxy1','stimmarks');
                    save(strcat(outpath,filesep,scanname,'_subj2_preprocessed.mat'),'oxy2', 'deoxy2', 'totaloxy2','z_oxy2', 'z_deoxy2', 'z_totaloxy2','stimmarks');
                    SD=SD1;
                    d=new_d1;
                    s=s1;
                    save(strcat(outpath,filesep,scanname,'_subj1.nirs'),'aux','d','s','SD','t');
                    SD=SD2;
                    d=new_d2;
                    s=s2;
                    save(strcat(outpath,filesep,scanname,'_subj2.nirs'),'aux','d','s','SD','t');
                end
            end
        end
    end
    Elapsedtime = toc(Elapsedtime);
    fprintf('\n\t Elapsed time: %g seconds \n', Elapsedtime);
    
else
    %all again but no dyad stuff
    fprintf('\n\t Preprocessing ...\n')
    reverseStr = '';
    Elapsedtime = tic;
    for i=1:length(currdir);
        subj=currdir(i).name;
        subjdir=dir(strcat(rawdir,filesep,subj,filesep,dataprefix,'*'));
        msg = sprintf('\n\t subject number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        %if there is only one scan per participant
        if isempty(subjdir) || ~isdir(strcat(rawdir,filesep,subj,filesep,subjdir(1).name))
            subjfolder = strcat(rawdir,filesep,subj);
            outpath = strcat(rawdir,filesep,'PreProcessedFiles',filesep,subj);
            if device==2
                outpath = outpath(1:end-5);
            end
            if ~exist(outpath,'dir')
                
                %1) extract data values
                anytofilter=1;
                if device==1
                    [d, sd_ind, samprate, wavelengths, s] = extractNIRxData(subjfolder);
                    probenumchannels = probeInfo.probes.nChannel0;
                    datanumchannels = size(d,2)/2;
                    if probenumchannels~=datanumchannels
                        error('ERROR: number of data channels in hdr file does not match number of channels in probeInfo file.');
                    end
                elseif device==2
                    [d, samprate, s, SD, aux, t] = extractTechEnData(subjfolder);
                end
                
                %2) Trim beginning of data to 10s before onset, if there is
                %a lot of dead time before that 
                ssum = sum(s,2);
                stimmarks = find(ssum);
                if length(stimmarks)>=1
                    begintime = stimmarks(1) - round(samprate*10);
                    if begintime>0
                        d = d(begintime:end,:);
                        s = s(begintime:end,:);
                        stimmarks = stimmarks-begintime;
                    end
                end
                
                %trim off last ten seconds, to remove edge artifacts that
                %might mess up hemodynamics calculation
                d = d((1:end-round(10*samprate)),:);
                s = s((1:end-round(10*samprate)),:);
                
                %3) identify and remove bad channels
                satlength = 2; %in seconds
                QCoDthresh = 0.1;
                channelmask = removeBadChannels(d, samprate, satlength, QCoDthresh);
            
                %4)Prep rest of .nirs format
                if device==1
                    [SD, aux, t] = getMiscNirsVars(d, sd_ind, samprate, wavelengths, probeInfo, channelmask);
                elseif device==2
                    SD.MeasListAct = [channelmask'; channelmask'];
                    SD.MeasListVis = SD.MeasListAct;
                end
                if length(stimmarks)>=1
                    if begintime>0
                        aux = aux(begintime:end,:);
                        t = t(begintime:end);
                    end
                end
                aux = aux((1:end-round(10*samprate)),:);
                t = t((1:end-round(10*samprate)),:);
                
                %4) motion filter, convert to hemodynamic changes
                [new_d, oxy, deoxy, totaloxy, z_oxy, z_deoxy, z_totaloxy] = fNIRSFilterPipeline(d, SD, samprate);
            
                mkdir(outpath)
                if device==2
                    subj=subj(1:end-5);
                end
                d=new_d;
                save(strcat(outpath,filesep,subj,'_preprocessed.mat'),'oxy', 'deoxy', 'totaloxy','z_oxy', 'z_deoxy', 'z_totaloxy','stimmarks');
                save(strcat(outpath,filesep,subj,'.nirs'),'aux','d','s','SD','t');
            end
        %if there are more than one scan per participant    
        else
            multiscan=1;
            for j=1:length(subjdir)
                scanname = subjdir(j).name;
                subjfolder = strcat(rawdir,filesep,subj,filesep,scanname);
            
                outpath = strcat(rawdir,filesep,'PreProcessedFiles',filesep,subj,filesep,scanname);
                if device==2
                    outpath = outpath(1:end-5);
                end
                if ~exist(outpath,'dir')
                    
                    %1) extract data values
                    anytofilter=1;
                    if device==1
                        [d, sd_ind, samprate, wavelengths, s] = extractNIRxData(subjfolder);
                        probenumchannels = probeInfo.probes.nChannel0;
                        datanumchannels = size(d,2)/2;
                        if probenumchannels~=datanumchannels
                            error('ERROR: number of data channels in hdr file does not match number of channels in probeInfo file.');
                        end
                    elseif device==2
                        [d, samprate, s, SD, aux, t] = extractTechEnData(subjfolder);
                    end
                    
                    %2) Trim beginning of data to 10s before onset, if there is
                    %a lot of dead time before that 
                    ssum = sum(s,2);
                    stimmarks = find(ssum);
                    if length(stimmarks)>=1
                        begintime = stimmarks(1) - round(samprate*10);
                        if begintime>0
                            d = d(begintime:end,:);
                            s = s(begintime:end,:);
                            stimmarks = stimmarks-begintime;
                        end
                    end
                    %trim off last ten seconds, to remove edge artifacts that
                    %might mess up hemodynamics calculation
                    d = d((1:end-round(10*samprate)),:);
                    s = s((1:end-round(10*samprate)),:);

                    %3) identify and remove bad channels
                    satlength = 2; %in seconds
                    QCoDthresh = 0.1;
                    channelmask = removeBadChannels(d, samprate, satlength, QCoDthresh);
                    
                    %4)Prep rest of .nirs format
                    if device==1
                        [SD, aux, t] = getMiscNirsVars(d, sd_ind, samprate, wavelengths, probeInfo, channelmask);
                    elseif device==2
                        SD.MeasListAct = [channelmask'; channelmask'];
                        SD.MeasListVis = SD.MeasListAct;
                    end
                    if length(stimmarks)>=1
                        if begintime>0
                            aux = aux(begintime:end,:);
                            t = t(begintime:end);
                        end
                    end
                    aux = aux((1:end-round(10*samprate)),:);
                    t = t((1:end-round(10*samprate)),:);

                    %5) motion filter, convert to hemodynamic changes
                    [new_d, oxy, deoxy, totaloxy, z_oxy, z_deoxy, z_totaloxy] = fNIRSFilterPipeline(d, SD, samprate);
            
                    mkdir(outpath)
                    if device==2
                        scanname=scanname(1:end-5);
                    end
                    d=new_d;
                    save(strcat(outpath,filesep,scanname,'_preprocessed.mat'),'oxy', 'deoxy', 'totaloxy','z_oxy', 'z_deoxy', 'z_totaloxy','stimmarks');
                    save(strcat(outpath,filesep,scanname,'.nirs'),'aux','d','s','SD','t');
                end
            end
        end
    end
    Elapsedtime = toc(Elapsedtime);
    fprintf('\n\t Elapsed time: %g seconds\n', Elapsedtime);
end

if anytofilter
    qualityAssessment(dataprefix,dyads,multiscan,size(d,2),samprate,0.1,strcat(rawdir,filesep,'PreProcessedFiles'));
end

end
