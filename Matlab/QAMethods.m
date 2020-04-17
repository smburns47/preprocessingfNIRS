% dataprefix = 'ST';
% dyads = 0;
% channelnum = 108;
% diff_thresh = 0.85;
% samprate = 1.95;

function QAMethods(dataprefix,dyads,channelnum,samprate,thresh,currdir)

if ~exists('currdir','var')
    preprocdir=uigetdir('','Choose Preprocessed Data Directory');
    currdir=dir(strcat(preprocdir,filesep,dataprefix,'*'));
    if length(currdir)<1
        error(['ERROR: No data files found with ',dataprefix,' prefix']);
    end
end

thresh_corr_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.03,0.08,0.13,0.22,0.27,0.31,0.34,0.38,0.42,0.47,0.52,0.58,0.66,0.75,0.95,1,1,1,1,1});
thresh_ps_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.08,0.21,0.26,0.30,0.35,0.39,0.43,0.47,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96});
if~exists('thresh','var')
    thresh=0.1;
end
thresh_corr = thresh_corr_map(thresh);
thresh_ps = thresh_ps_map(thresh);


if dyads
    %to be written
else
    %to be written - make sure to detect if there are multiple scans or not;
    %maybe pass an argument from main file so you don't have to keep
    %detecting it
    fprintf('\n\t Generating data quality reports ...\n')
    reverseStr = '';
    
%     firstsubj=currdir(1).name;
%     firstsubjdir=dir(strcat(preprocdir,filesep,firstsubj,filesep,dataprefix,'*'));
%     firstsubjnamelength = length(firstsubj);
%     scannames = cell(1,length(firstsubjdir));
%     for j=1:length(firstsubjdir)
%         scanname = firstsubjdir(j).name(firstsubjnamelength+1:end);
%         scannames{j} = scanname;
%     end
%     for x=1:length(scannames)
%         scannames{x} = regexprep(scannames{x},'_','');
%     end
    qatable_corr = array2table(zeros(length(currdir),1));
    qatable_corr.Properties.VariableNames={'subjname'};
    all_subjtables_corr = cell(1,length(currdir));
    qatable_ps = qatable_corr;
    all_subjtables_ps = all_subjtables_corr;
    scannames = cell(1,0);

    for i=1:length(currdir)
        msg = sprintf('\n\t subj number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        subj=currdir(i).name;
        subjnamelength = length(subj);
        subjdir=dir(strcat(preprocdir,filesep,subj,filesep,dataprefix,'*'));
        subjtable_corr = array2table([1:108]');
        subjtable_corr.Properties.VariableNames={'channelnum'};
        subjtable_ps = subjtable_corr;
        qatable_corr.subjname{i} = subj;
        qatable_ps.subjname{i} = subj;
        for j=1:length(subjdir)
            scanname = subjdir(j).name;
            subscanname = subjdir(j).name(subjnamelength+1:end);
            subscanname = regexprep(subscanname,'_','');
            if ~any(strcmp(scannames,subscanname))
                scannames = [scannames,subscanname];
                qatable2 = array2table(zeros(length(currdir),1));
                qatable2.Properties.VariableNames=subscanname;
                qatable_corr = [qatable_corr qatable2];
                qatable_ps = [qatable_ps qatable2];
            end
            scandir = dir(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,dataprefix,'*.mat'));
            if exist(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,scandir(1).name),'file')
               load(strcat(preprocdir,filesep,subj,filesep,scanname,filesep,scandir(1).name))
               channelcount_corr = 0;
               channelcount_ps = 0;
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
                           channelcount_corr=channelcount_corr+1;
                       end
                       if avgPS<(1-thresh_ps)
                            subjtable_ps.(subscanname)(k) = 1;
                            channelcount_ps=channelcount_ps+1;
                       end
                   else
                      channelcount_corr=channelcount_corr+1; 
                      channelcount_ps=channelcount_ps+1;
                   end
               end
               qatable_corr.(subscanname)(i) = channelcount_corr;
               qatable_ps.(subscanname)(i) = channelcount_ps;
               all_subjtables_corr{i} = subjtable_corr;
               all_subjtables_ps{i} = subjtable_ps;
            end
        end
    end

    %--------------------------------------
    dataprefix = 'DYAD';
    dyads = 1;
    channelnum = 23;
    diff_thresh = 0.85;
    samprate = 7.81;

    %function qacheckMotionPerformance(dataprefix,dyads)
    preprocdir=uigetdir('','Choose Data Directory');
    currdir=dir(strcat(preprocdir,filesep,dataprefix,'*'));
    if length(currdir)<1
        error(['ERROR: No data files found with ',dataprefix,' prefix']);
    end

    fprintf('\n\t Checking motion correction performance ...\n')
    reverseStr = '';

    qatable_corr = array2table(zeros(length(currdir)*2,1));
    qatable2 = cell2table(cell(length(currdir)*2,1));
    qatable2.Properties.VariableNames={'subjname'};
    qatable_corr = [qatable2 qatable_corr];
    all_subjtables_corr = cell(1,length(currdir));
    for i=1:length(currdir)
        msg = sprintf('\n\t subj number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        dyad=currdir(i).name;
        dyadnamelength = length(dyad);
        subjtable_corr = array2table(zeros(channelnum,1));
        channelnames = array2table([1:channelnum]');
        channelnames.Properties.VariableNames={'channelnum'};
        subjtable_corr = [channelnames subjtable_corr];
        qatable_corr.subjname{i*2-1} = strcat(dyad,'_1');
        qatable_corr.subjname{i*2} = strcat(dyad,'_2');
        scandir1 = dir(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,dataprefix,'*.mat'));
        scandir2 = dir(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,dataprefix,'*.mat'));
        %subj1
        if exist(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,scandir1(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,scandir1(1).name))
               channelcount = 0;
               for k=1:size(oxy1,2)
                   trace = oxy1(:,k);
                   trace_orig = trace;
                   offset=round(0.5*samprate);
                   for datapoint=(offset+5):(length(trace)-5)
                       if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                           trace(datapoint-4:datapoint+4,1) = linspace(trace(datapoint-5,1),trace(datapoint+5,1),9);
                       end
                   end
                   autocorrdiff = corrcoef(trace_orig, trace);
                   autocorrdiff = autocorrdiff(1,2);
                   if autocorrdiff<diff_thresh
                       subjtable_corr{k,2} = 1;
                       channelcount=channelcount+1;
                   end
               end
               qatable_corr{i*2-1,2} = channelcount;
               all_subjtables_corr{i} = subjtable_corr;
        end
        %subj2
        if exist(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,scandir2(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,scandir2(1).name))
               channelcount = 0;
               for k=1:size(oxy2,2)
                   trace = oxy2(:,k);
                   trace_orig = trace;
                   offset=round(0.5*samprate);
                   for datapoint=(offset+5):(length(trace)-5)
                       if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                           trace(datapoint-4:datapoint+4,1) = linspace(trace(datapoint-5,1),trace(datapoint+5,1),9);
                       end
                   end
                   autocorrdiff = corrcoef(trace_orig, trace);
                   autocorrdiff = autocorrdiff(1,2);
                   if autocorrdiff<diff_thresh
                       subjtable_corr{k,2} = 1;
                       channelcount=channelcount+1;
                   end
               end
               qatable_corr{i*2,2} = channelcount;
               all_subjtables_corr{i} = subjtable_corr;
        end
    end

    %-------------------------------

    dataprefix = 'DYAD';
    dyads = 1;
    channelnum = 23;
    diff_thresh = 0.75;
    samprate = 7.81;

    %function qacheckMotionPerformance(dataprefix,dyads)
    preprocdir=uigetdir('','Choose Data Directory');
    currdir=dir(strcat(preprocdir,filesep,dataprefix,'*'));
    if length(currdir)<1
        error(['ERROR: No data files found with ',dataprefix,' prefix']);
    end

    fprintf('\n\t Checking motion correction performance ...\n')
    reverseStr = '';

    qatable_corr = array2table(zeros(length(currdir)*2,1));
    qatable2 = cell2table(cell(length(currdir)*2,1));
    qatable2.Properties.VariableNames={'subjname'};
    qatable_corr = [qatable2 qatable_corr];
    all_subjtables_corr = cell(1,length(currdir));
    for i=1:length(currdir)
        msg = sprintf('\n\t subj number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        dyad=currdir(i).name;
        dyadnamelength = length(dyad);
        subjtable_corr = array2table(zeros(channelnum,1));
        channelnames = array2table([1:channelnum]');
        channelnames.Properties.VariableNames={'channelnum'};
        subjtable_corr = [channelnames subjtable_corr];
        qatable_corr.subjname{i*2-1} = strcat(dyad,'_1');
        qatable_corr.subjname{i*2} = strcat(dyad,'_2');
        scandir1 = dir(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,dataprefix,'*.mat'));
        scandir2 = dir(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,dataprefix,'*.mat'));
        %subj1
        if exist(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,scandir1(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,'subj1',filesep,scandir1(1).name))
               channelcount = 0;
               for k=1:size(oxy1,2)
                   trace = oxy1(:,k);
                   trace_orig = trace;
                   offset=round(1*samprate);
                   for datapoint=(offset+6):(length(trace)-6)
                       if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                           trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+6,1),11);
                       end
                   end
                   bothfiltered = [trace'; trace_orig'];
                   PS1 = angle(hilbert(bothfiltered(1,:))); % instantaneous phases of channelxtime data
                   PS2 = angle(hilbert(bothfiltered(2,:)));
                   PS_timecourse = 1-sin(abs(PS1-PS2)/2);
                   avg_phasesynch = nanmean(PS_timecourse,2);
                   if avg_phasesynch<diff_thresh
                       subjtable_corr{k,2} = 1;
                       channelcount=channelcount+1;
                   end
               end
               qatable_corr{i*2-1,2} = channelcount;
               all_subjtables_corr{i} = subjtable_corr;
        end
        %subj2
        if exist(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,scandir2(1).name),'file')
               load(strcat(preprocdir,filesep,dyad,filesep,'subj2',filesep,scandir2(1).name))
               channelcount = 0;
               for k=1:size(oxy2,2)
                   trace = oxy2(:,k);
                   trace_orig = trace;
                   offset=round(1*samprate);
                   for datapoint=(offset+6):(length(trace)-6)
                       if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
                           trace(datapoint-5:datapoint+5,1) = linspace(trace(datapoint-6,1),trace(datapoint+6,1),11);
                       end
                   end
                   bothfiltered = [trace'; trace_orig'];
                   PS1 = angle(hilbert(bothfiltered(1,:))); % instantaneous phases of channelxtime data
                   PS2 = angle(hilbert(bothfiltered(2,:)));
                   PS_timecourse = 1-sin(abs(PS1-PS2)/2);
                   avg_phasesynch = nanmean(PS_timecourse,2);
                   if avg_phasesynch<diff_thresh
                       subjtable_corr{k,2} = 1;
                       channelcount=channelcount+1;
                   end
               end
               qatable_corr{i*2,2} = channelcount;
               all_subjtables_corr{i} = subjtable_corr;
        end
    end
       
end        

