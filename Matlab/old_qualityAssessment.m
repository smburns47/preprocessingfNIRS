function qualityAssessment(dataprefix,hyperscan,multiscan,channelnum,samprate,thresh,preprocdir)

%inputs: 
%       dataprefix: string. Prefix of every folder name that should be considered a
%           data folder. E.g., ST for ST_101, ST_102, etc.  
%       dyads: 0 or 1. 1 if hyperscanning, 0 if single subject.
%       multiscan: 0 or 1. If the experiment has multiple scans per person
%       channelnum: integer. Number of channels in the montage
%       samprate: double. Sampling rate in experiment
%       thresh: double, multiple of 0.05. Max amount you want to allow a 
%          synchrony estimate to vary by (measurement error allowance).
%           Default is 0.1 
%       preprocdir: string. Path to PreProcessedFiles directory
%
%outputs: csv data quality reports, per subject and whole dataset, reporting
%           which channels had poor data quality.

if ~exist('preprocdir','var')
    preprocdir=uigetdir('','Choose Preprocessed Data Directory');
    currdir=dir(strcat(preprocdir,filesep,dataprefix,'*'));
else
    currdir=dir(strcat(preprocdir,filesep,dataprefix,'*'));
end
if length(currdir)<1
    error(['ERROR: No data files found with ',dataprefix,' prefix']);
end

thresh_corr_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.03,0.08,0.13,0.22,0.27,0.31,0.34,0.38,0.42,0.47,0.52,0.58,0.66,0.75,0.95,1,1,1,1,1});
thresh_ps_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.08,0.21,0.26,0.30,0.35,0.39,0.43,0.47,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96});
if ~exist('thresh','var')
    thresh=0.1;
end
thresh_corr = thresh_corr_map(thresh);
thresh_ps = thresh_ps_map(thresh);

fprintf('\n\t Generating data quality reports ...\n')
    reverseStr = '';
    
if hyperscan
    qatable_corr = array2table(cell(length(currdir)*2,1));
    qatable_corr.Properties.VariableNames={'subjname'};
    qatable_ps = qatable_corr;
    qatable_noise = qatable_corr;
    chtable_corr = array2table([1:channelnum]');
    chtable_corr.Properties.VariableNames={'channelnum'};
    chtable_ps = chtable_corr;
    chtable_noise = chtable_corr;
    scannames = cell(1,0);
    for i=1:length(currdir)
        msg = sprintf('\n\t dyad number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        dyad=currdir(i).name;
        dyadnamelength = length(dyad);
        dyaddir=dir(strcat(preprocdir,filesep,dyad,filesep,dataprefix,'*'));
        if size(scannames,2)<1
            subjtable1_corr = array2table([1:channelnum]');
            subjtable1_corr.Properties.VariableNames={'channelnum'};
        else
            subjtable1_corr = array2table([[1:channelnum]' zeros(channelnum,size(scannames,2))]);
            varnames = [{'channelnum'} scannames];
            subjtable1_corr.Properties.VariableNames=varnames;
        end
        subjtable1_ps = subjtable1_corr;
        subjtable1_noise = subjtable1_corr;
        subjtable2_corr = subjtable1_corr;
        subjtable2_ps = subjtable1_corr;
        subjtable2_noise = subjtable1_corr;
        qatable_corr.subjname{i*2-1} = strcat(dyad,'_1');
        qatable_corr.subjname{i*2} = strcat(dyad,'_2');
        qatable_ps.subjname{i*2-1} = strcat(dyad,'_1');
        qatable_ps.subjname{i*2} = strcat(dyad,'_2');
        qatable_noise.subjname{i*2-1} = strcat(dyad,'_1');
        qatable_noise.subjname{i*2} = strcat(dyad,'_2');
        if multiscan
            for j=1:length(dyaddir)
                scanname = dyaddir(j).name;
                subscanname = dyaddir(j).name(dyadnamelength+1:end);
                subscanname = regexprep(subscanname,'_','');
                if ~any(strcmp(scannames,subscanname))
                    scannames = [scannames,subscanname];
                    qatable2 = array2table(zeros(length(currdir)*2,1));
                    chtable2 = array2table(zeros(channelnum,1));
                    qatable2.Properties.VariableNames={subscanname};
                    chtable2.Properties.VariableNames={subscanname};
                    qatable_corr = [qatable_corr qatable2];
                    qatable_ps = [qatable_ps qatable2];
                    qatable_noise = [qatable_noise qatable2];
                    chtable_corr = [chtable_corr chtable2];
                    chtable_ps = [chtable_ps chtable2];
                    chtable_noise = [chtable_noise chtable2];
                    subjtable1_corr = [subjtable1_corr chtable2];
                    subjtable1_ps = [subjtable1_ps chtable2];
                    subjtable1_noise = [subjtable1_noise chtable2];
                    subjtable2_corr = [subjtable2_corr chtable2];
                    subjtable2_ps = [subjtable2_ps chtable2];
                    subjtable2_noise = [subjtable2_noise chtable2];
                end
                %subj1
               scandir = dir(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,'*_subj1_preprocessed.mat'));
               if exist(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,scandir(1).name),'file')
                   load(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,scandir(1).name))
                   channelcount_corr = 0;
                   channelcount_ps = 0;
                   noisecount = 0;
                   for k=1:size(z_oxy1,2)
                       trace = z_oxy1(:,k);
                       if ~any(isnan(trace))
                           trace_orig = trace;
                           offset=round(samprate);
                           for datapoint=(offset+6):(length(trace)-6)
                               if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                                   trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                               end
                           end
                           autocorrdiff = corrcoef(trace_orig, trace);
                           autocorrdiff = autocorrdiff(1,2);
                           PS1 = angle(hilbert(trace));
                           PS2 = angle(hilbert(trace_orig));
                           avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                           if autocorrdiff<(1-thresh_corr)
                               subjtable1_corr.(subscanname)(k) = 1;
                               chtable_corr.(subscanname)(k) = chtable_corr.(subscanname)(k) + 1;
                               channelcount_corr=channelcount_corr+1;
                           end
                           if avgPS<(1-thresh_ps)
                                subjtable1_ps.(subscanname)(k) = 1;
                                chtable_ps.(subscanname)(k) = chtable_ps.(subscanname)(k) + 1;
                                channelcount_ps=channelcount_ps+1;
                           end
                       else
                           subjtable1_noise.(subscanname)(k) = 1;
                           chtable_noise.(subscanname)(k) = chtable_noise.(subscanname)(k) + 1;
                           noisecount=noisecount+1;
                       end
                   end
               qatable_corr.(subscanname)(i*2-1) = channelcount_corr;
               qatable_ps.(subscanname)(i*2-1) = channelcount_ps;
               qatable_noise.(subscanname)(i*2-1) = noisecount;

               end
               %subj2
               scandir = dir(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,'*_subj2_preprocessed.mat'));
               if exist(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,scandir(1).name),'file')
                   load(strcat(preprocdir,filesep,dyad,filesep,scanname,filesep,scandir(1).name))
                   channelcount_corr = 0;
                   channelcount_ps = 0;
                   noisecount = 0;
                   for k=1:size(z_oxy2,2)
                       trace = z_oxy2(:,k);
                       if ~any(isnan(trace))
                           trace_orig = trace;
                           offset=round(samprate);
                           for datapoint=(offset+6):(length(trace)-6)
                               if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                                   trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                               end
                           end
                           autocorrdiff = corrcoef(trace_orig, trace);
                           autocorrdiff = autocorrdiff(1,2);
                           PS1 = angle(hilbert(trace));
                           PS2 = angle(hilbert(trace_orig));
                           avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                           if autocorrdiff<(1-thresh_corr)
                               subjtable2_corr.(subscanname)(k) = 1;
                               chtable_corr.(subscanname)(k) = chtable_corr.(subscanname)(k) + 1;
                               channelcount_corr=channelcount_corr+1;
                           end
                           if avgPS<(1-thresh_ps)
                                subjtable2_ps.(subscanname)(k) = 1;
                                chtable_ps.(subscanname)(k) = chtable_ps.(subscanname)(k) + 1;
                                channelcount_ps=channelcount_ps+1;
                           end
                       else
                           subjtable2_noise.(subscanname)(k) = 1;
                           chtable_noise.(subscanname)(k) = chtable_noise.(subscanname)(k) + 1;
                           noisecount=noisecount+1;
                       end
                   end
               qatable_corr.(subscanname)(i*2) = channelcount_corr;
               qatable_ps.(subscanname)(i*2) = channelcount_ps;
               qatable_noise.(subscanname)(i*2) = noisecount;
               end
            end
           outpath_noise=strcat(preprocdir,filesep,dyad,filesep,'preQAreport_subj1_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_flaggedchannels_ps.csv');
           writetable(subjtable1_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable1_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable1_noise,outpath_noise,'Delimiter',',');
           outpath_noise=strcat(preprocdir,filesep,dyad,filesep,'preQAreport_subj2_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_ps.csv');
           writetable(subjtable2_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable2_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable2_noise,outpath_noise,'Delimiter',',');
           
           combined = table2array(subjtable1_corr(:,2:end)) + table2array(subjtable1_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable1_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable1_ps(:,2:end)) + table2array(subjtable1_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable1_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
           
           combined = table2array(subjtable2_corr(:,2:end)) + table2array(subjtable2_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable2_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable2_ps(:,2:end)) + table2array(subjtable2_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable2_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
        else
            subscanname = 'scan';
            if ~any(strcmp(scannames,subscanname))
                scannames = [scannames,subscanname];
                qatable2 = array2table(zeros(length(currdir)*2,1));
                chtable2 = array2table(zeros(channelnum,1));
                qatable2.Properties.VariableNames={subscanname};
                chtable2.Properties.VariableNames={subscanname};
                qatable_corr = [qatable_corr qatable2];
                qatable_ps = [qatable_ps qatable2];
                qatable_noise = [qatable_noise qatable2];
                chtable_corr = [chtable_corr chtable2];
                chtable_ps = [chtable_ps chtable2];
                chtable_noise = [chtable_noise chtable2];
                subjtable1_corr = [subjtable1_corr chtable2];
                subjtable1_ps = [subjtable1_ps chtable2];
                subjtable1_noise = [subjtable1_noise chtable2];
                subjtable2_corr = [subjtable2_corr chtable2];
                subjtable2_ps = [subjtable2_ps chtable2];
                subjtable2_noise = [subjtable2_noise chtable2];
            end
            %subj1
            scandir = dir(strcat(preprocdir,filesep,dyad,filesep,'*_subj1_preprocessed.mat'));
           if exist(strcat(preprocdir,filesep,dyad,filesep,scandir(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,scandir(1).name))
               channelcount_corr = 0;
               channelcount_ps = 0;
               noisecount = 0;
               for k=1:size(z_oxy1,2)
                   trace = z_oxy1(:,k);
                   if ~any(isnan(trace))
                       trace_orig = trace;
                       offset=round(samprate);
                       for datapoint=(offset+6):(length(trace)-6)
                           if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                               trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                           end
                       end
                       autocorrdiff = corrcoef(trace_orig, trace);
                       autocorrdiff = autocorrdiff(1,2);
                       PS1 = angle(hilbert(trace));
                       PS2 = angle(hilbert(trace_orig));
                       avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                       if autocorrdiff<(1-thresh_corr)
                           subjtable1_corr.scan(k) = 1;
                           chtable_corr.scan(k) = chtable_corr.scan(k) + 1;
                           channelcount_corr=channelcount_corr+1;
                       end
                       if avgPS<(1-thresh_ps)
                            subjtable1_ps.scan(k) = 1;
                            chtable_ps.scan(k) = chtable_ps.scan(k) + 1;
                            channelcount_ps=channelcount_ps+1;
                       end
                   else
                      subjtable1_noise.scan(k) = 1;
                      chtable_noise.scan(k) = chtable_noise.scan(k) + 1;
                      noisecount=noisecount+1;
                   end
               end
               qatable_corr(i*2-1,2) = channelcount_corr;
               qatable_ps(i*2-1,2) = channelcount_ps;
               qatable_noise(i*2-1,2) = noisecount;
            end
            %subj2
            scandir = dir(strcat(preprocdir,filesep,dyad,filesep,'*_subj2_preprocessed.mat'));
           if exist(strcat(preprocdir,filesep,dyad,filesep,scandir(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,scandir(1).name))
               channelcount_corr = 0;
               channelcount_ps = 0;
               noisecount = 0;
               for k=1:size(z_oxy2,2)
                   trace = z_oxy2(:,k);
                   if ~any(isnan(trace))
                       trace_orig = trace;
                       offset=round(samprate);
                       for datapoint=(offset+6):(length(trace)-6)
                           if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                               trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                           end
                       end
                       autocorrdiff = corrcoef(trace_orig, trace);
                       autocorrdiff = autocorrdiff(1,2);
                       PS1 = angle(hilbert(trace));
                       PS2 = angle(hilbert(trace_orig));
                       avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                       if autocorrdiff<(1-thresh_corr)
                           subjtable2_corr.scan(k) = 1;
                           chtable_corr.scan(k) = chtable_corr.scan(k) + 1;
                           channelcount_corr=channelcount_corr+1;
                       end
                       if avgPS<(1-thresh_ps)
                            subjtable2_ps.scan(k) = 1;
                            chtable_ps.scan(k) = chtable_ps.scan(k) + 1;
                            channelcount_ps=channelcount_ps+1;
                       end
                   else
                      subjtable2_noise.scan(k) = 1;
                      chtable_noise.scan(k) = chtable_noise.scan(k) + 1;
                      noisecount=noisecount+1;
                   end
               end
               qatable_corr(i*2,2) = channelcount_corr;
               qatable_ps(i*2,2) = channelcount_ps;
               qatable_noise(i*2,2) = noisecount;
           end
           outpath_noise=strcat(preprocdir,filesep,dyad,filesep,'preQAreport_subj2_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_ps.csv');
           writetable(subjtable2_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable2_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable2_noise,outpath_noise,'Delimiter',',');
           outpath_noise=strcat(preprocdir,filesep,dyad,filesep,'preQAreport_subj2_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_flaggedchannels_ps.csv');
           writetable(subjtable2_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable2_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable2_noise,outpath_noise,'Delimiter',',');
           
           combined = table2array(subjtable1_corr(:,2:end)) + table2array(subjtable1_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable1_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable1_ps(:,2:end)) + table2array(subjtable1_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable1_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj1_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
           
           combined = table2array(subjtable2_corr(:,2:end)) + table2array(subjtable2_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable2_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable2_ps(:,2:end)) + table2array(subjtable2_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable2_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,dyad,filesep,'postQAreport_subj2_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
        end
    end

else
    qatable_corr = cell2table(cell(length(currdir),1));
    qatable_corr.Properties.VariableNames={'subjname'};
    qatable_ps = qatable_corr;
    qatable_noise = qatable_corr;
    chtable_corr = array2table([1:channelnum]');
    chtable_corr.Properties.VariableNames={'channelnum'};
    chtable_ps = chtable_corr;
    chtable_noise = chtable_corr;
    scannames = cell(1,0);

    for i=1:length(currdir)
        msg = sprintf('\n\t subj number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        subj=currdir(i).name;
        subjnamelength = length(subj);
        subjdir=dir(strcat(preprocdir,filesep,subj,filesep,dataprefix,'*'));
        if size(scannames,2)<1
            subjtable_corr = array2table([1:channelnum]');
            subjtable_corr.Properties.VariableNames={'channelnum'};
        else
            subjtable_corr = array2table([[1:channelnum]' zeros(channelnum,size(scannames,2))]);
            varnames = [{'channelnum'} scannames];
            subjtable_corr.Properties.VariableNames=varnames;
        end
        subjtable_ps = subjtable_corr;
        subjtable_noise = subjtable_corr;
        qatable_corr.subjname{i} = subj;
        qatable_ps.subjname{i} = subj;
        qatable_noise.subjname{i} = subj;
        if multiscan
            for j=1:length(subjdir)
                scanname = subjdir(j).name;
                subscanname = subjdir(j).name(subjnamelength+1:end);
                subscanname = regexprep(subscanname,'_','');
                if ~any(strcmp(scannames,subscanname))
                    scannames = [scannames,subscanname];
                    qatable2 = array2table(zeros(length(currdir),1));
                    qatable2.Properties.VariableNames={subscanname};
                    chtable2 = array2table(zeros(channelnum,1));
                    chtable2.Properties.VariableNames={subscanname};
                    qatable_corr = [qatable_corr qatable2];
                    qatable_ps = [qatable_ps qatable2];
                    qatable_noise = [qatable_noise qatable2];
                    chtable_corr = [chtable_corr chtable2];
                    chtable_ps = [chtable_ps chtable2];
                    chtable_noise = [chtable_noise chtable2];
                    subjtable_corr = [subjtable_corr chtable2];
                    subjtable_ps = [subjtable_ps chtable2];
                    subjtable_noise = [subjtable_noise chtable2];
                end
                scandir = dir(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,'*_preprocessed.mat'));
                if exist(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,scandir(1).name),'file')
                   load(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,scandir(1).name))
                   channelcount_corr = 0;
                   channelcount_ps = 0;
                   noisecount = 0;
                   for k=1:size(z_oxy,2)
                       trace = z_oxy(:,k);
                       if ~any(isnan(trace))
                           trace_orig = trace;
                           offset=round(samprate);
                           for datapoint=(offset+6):(length(trace)-6)
                               if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                                   trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                               end
                           end
                           autocorrdiff = corrcoef(trace_orig, trace);
                           autocorrdiff = autocorrdiff(1,2);
                           PS1 = angle(hilbert(trace));
                           PS2 = angle(hilbert(trace_orig));
                           avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                           if autocorrdiff<(1-thresh_corr)
                               subjtable_corr.(subscanname)(k) = 1;
                               chtable_corr.(subscanname)(k) = chtable_corr.(subscanname)(k) + 1;
                               channelcount_corr=channelcount_corr+1;
                           end
                           if avgPS<(1-thresh_ps)
                                subjtable_ps.(subscanname)(k) = 1;
                                chtable_ps.(subscanname)(k) = chtable_ps.(subscanname)(k) + 1;
                                channelcount_ps=channelcount_ps+1;
                           end
                       else
                           subjtable_noise.(subscanname)(k) = 1;
                           chtable_noise.(subscanname)(k) = chtable_noise.(subscanname)(k) + 1;
                           noisecount=noisecount+1;
                       end
                   end
                   qatable_corr.(subscanname)(i) = channelcount_corr;
                   qatable_ps.(subscanname)(i) = channelcount_ps;
                   qatable_noise.(subscanname)(i) = noisecount;
                end
            end
           outpath_noise=strcat(preprocdir,filesep,subj,filesep,'preQAreport_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,subj,filesep,'postQAreport_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,subj,filesep,'postQAreport_flaggedchannels_ps.csv');
           writetable(subjtable_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable_noise,outpath_noise,'Delimiter',',');
           combined = table2array(subjtable_corr(:,2:end)) + table2array(subjtable_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,subj,filesep,'postQAreport_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable_ps(:,2:end)) + table2array(subjtable_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,subj,filesep,'postQAreport_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
        else
            subscanname = 'scan';
            if ~any(strcmp(scannames,subscanname))
                scannames = [scannames,subscanname];
                qatable2 = array2table(zeros(length(currdir),1));
                qatable2.Properties.VariableNames={subscanname};
                chtable2 = array2table(zeros(channelnum,1));
                chtable2.Properties.VariableNames={subscanname};
                qatable_corr = [qatable_corr qatable2];
                qatable_ps = [qatable_ps qatable2];
                qatable_noise = [qatable_noise qatable2];
                chtable_corr = [chtable_corr chtable2];
                chtable_ps = [chtable_ps chtable2];
                chtable_noise = [chtable_noise chtable2];
                subjtable_corr = [subjtable_corr chtable2];
                subjtable_ps = [subjtable_ps chtable2];
                subjtable_noise = [subjtable_noise chtable2];
            end
            scandir = dir(strcat(preprocdir,filesep,subj,filesep,'*_preprocessed.mat'));
            if exist(strcat(preprocdir,filesep,subj,filesep,scandir(1).name),'file')
               load(strcat(preprocdir,filesep,subj,filesep,scandir(1).name))
               channelcount_corr = 0;
               channelcount_ps = 0;
               noisecount = 0;
               for k=1:size(z_oxy,2)
                   trace = z_oxy(:,k);
                   if ~any(isnan(trace))
                       trace_orig = trace;
                       offset=round(samprate);
                       for datapoint=(offset+6):(length(trace)-6)
                           if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                               trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+5,1),11);
                           end
                       end
                       autocorrdiff = corrcoef(trace_orig, trace);
                       autocorrdiff = autocorrdiff(1,2);
                       PS1 = angle(hilbert(trace));
                       PS2 = angle(hilbert(trace_orig));
                       avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
                       if autocorrdiff<(1-thresh_corr)
                           subjtable_corr.scan(k) = 1;
                           chtable_corr.scan(k) = chtable_corr.scan(k) + 1;
                           channelcount_corr=channelcount_corr+1;
                       end
                       if avgPS<(1-thresh_ps)
                            subjtable_ps.scan(k) = 1;
                            chtable_ps.scan(k) = chtable_ps.scan(k) + 1;
                            channelcount_ps=channelcount_ps+1;
                       end
                   else
                      subjtable_noise.scan(k) = 1;
                      chtable_noise.scan(k) = chtable_noise.scan(k) + 1;
                      noisecount=noisecount+1;
                   end
                   qatable_corr.scan(i) = channelcount_corr;
                   qatable_ps.scan(i) = channelcount_ps;
                   qatable_noise.scan(i) = noisecount;
               end
            end
           outpath_noise=strcat(preprocdir,filesep,subj,filesep,'preQAreport_noisychannels.csv');
           outpath_corr=strcat(preprocdir,filesep,subj,filesep,'postQAreport_flaggedchannels_corr.csv');
           outpath_ps=strcat(preprocdir,filesep,subj,filesep,'postQAreport_flaggedchannels_ps.csv');
           writetable(subjtable_corr,outpath_corr,'Delimiter',',');
           writetable(subjtable_ps,outpath_ps,'Delimiter',',');
           writetable(subjtable_noise,outpath_noise,'Delimiter',',');
           combined = table2array(subjtable_corr(:,2:end)) + table2array(subjtable_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable_corr;
           combinedtable(:,2:end) = combined;
           outpath_corr=strcat(preprocdir,filesep,subj,filesep,'postQAreport_allbadch_corr.csv');
           writetable(combinedtable,outpath_corr,'Delimiter',',');
           combined = table2array(subjtable_ps(:,2:end)) + table2array(subjtable_noise(:,2:end));
           combined = array2table(combined);
           combinedtable = subjtable_ps;
           combinedtable(:,2:end) = combined;
           outpath_ps=strcat(preprocdir,filesep,subj,filesep,'postQAreport_allbadch_ps.csv');
           writetable(combinedtable,outpath_ps,'Delimiter',',');
        end
    end
end
qaoutpath_noise=strcat(preprocdir,filesep,'preQAreport_allsubj_noisychannels.csv');
qaoutpath_corr=strcat(preprocdir,filesep,'postQAreport_allsubj_flaggedchannels_corr.csv');
qaoutpath_ps=strcat(preprocdir,filesep,'postQAreport_allsubj_flaggedchannels_ps.csv');
writetable(qatable_corr,qaoutpath_corr,'Delimiter',',');
writetable(qatable_ps,qaoutpath_ps,'Delimiter',',');
writetable(qatable_noise,qaoutpath_noise,'Delimiter',',');
combined = table2array(qatable_corr(:,2:end)) + table2array(qatable_noise(:,2:end));
combined = array2table(combined);
combinedtable = qatable_corr;
combinedtable(:,2:end) = combined;
outpath_corr=strcat(preprocdir,filesep,'postQAreport_allsubj_allbadch_corr.csv');
writetable(combinedtable,outpath_corr,'Delimiter',',');
combined = table2array(qatable_ps(:,2:end)) + table2array(qatable_noise(:,2:end));
combined = array2table(combined);
combinedtable = qatable_ps;
combinedtable(:,2:end) = combined;
outpath_ps=strcat(preprocdir,filesep,'postQAreport_allsubj_allbadch_ps.csv');
writetable(combinedtable,outpath_ps,'Delimiter',',');

choutpath_noise=strcat(preprocdir,filesep,'preQAreport_allch_noisychannels.csv');
choutpath_corr=strcat(preprocdir,filesep,'postQAreport_allch_flaggedchannels_corr.csv');
choutpath_ps=strcat(preprocdir,filesep,'postQAreport_allch_flaggedchannels_ps.csv');
writetable(chtable_corr,choutpath_corr,'Delimiter',',');
writetable(chtable_ps,choutpath_ps,'Delimiter',',');
writetable(chtable_noise,choutpath_noise,'Delimiter',',');
combined = table2array(chtable_corr(:,2:end)) + table2array(chtable_noise(:,2:end));
combined = array2table(combined);
combinedtable = chtable_corr;
combinedtable(:,2:end) = combined;
outpath_corr=strcat(preprocdir,filesep,'postQAreport_allch_allbadch_corr.csv');
writetable(combinedtable,outpath_corr,'Delimiter',',');
combined = table2array(chtable_ps(:,2:end)) + table2array(chtable_noise(:,2:end));
combined = array2table(combined);
combinedtable = chtable_ps;
combinedtable(:,2:end) = combined;
outpath_ps=strcat(preprocdir,filesep,'postQAreport_allch_allbadch_ps.csv');
writetable(combinedtable,outpath_ps,'Delimiter',',');

end
