function qualityReport(dataprefix,hyperscan,multiscan,scannames,channelnum,preprocdir)

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

fprintf('\n\t Generating data quality reports ...\n')
    reverseStr = '';
  
qatable1a = array2table(cell(1,1));
qatable1a.Properties.VariableNames={'subjname'};
qatable1b = array2table(zeros(1,length(scannames)));
qatable1b.Properties.VariableNames=scannames;
qatable1 = [qatable1a qatable1b];
qatable2 = qatable1;
qatable_copy = qatable1;
chtable1a = array2table([1:channelnum]');
chtable1a.Properties.VariableNames={'channelnum'};
chtable1b = array2table(zeros(channelnum,length(scannames)));
chtable1b.Properties.VariableNames=scannames;
chtable1 = [chtable1a chtable1b];
chtable2 = chtable1;    
    
if hyperscan
    for i=1:length(currdir)
        msg = sprintf('\n\t group number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        group=currdir(i).name;
        groupdir=dir(strcat(preprocdir,filesep,group,filesep,dataprefix,'*'));
        
        for j=1:length(groupdir)
            subjname = groupdir(j).name;
            subjnamelength = length(subjname);
            qatable_copy.subjname{1} = subjname;
            qatable1 = [qatable1; qatable_copy];
            qatable2 = [qatable2; qatable_copy];
            subjdir = dir(strcat(preprocdir,filesep,group,filesep,subjname,filesep,dataprefix,'*'));
            
            if multiscan
                for k=1:length(subjdir)
                    scanname = subjdir(k).name;
                    subscanname = scanname(subjnamelength+1:end);
                    subscanname = regexprep(subscanname,'_','');
                    scandir = dir(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,'*_nonoisych.mat'));
                    if exist(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,scandir(1).name),'file')
                        load(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,scandir(1).name))
                        goodchannels = ~isnan(z_oxy(1,:));
                        sumgoodchannels = sum(goodchannels);
                        qatable1.(subscanname)(end) = sumgoodchannels;
                        chtable1.(subscanname)(:) = chtable1.(subscanname)(:) + goodchannels';
                    end
                    scandir = dir(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,'*_nouncertainch.mat'));
                    if exist(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,scandir(1).name),'file')
                        load(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scanname,filesep,scandir(1).name))
                        goodchannels = ~isnan(z_oxy(1,:));
                        sumgoodchannels = sum(goodchannels);
                        qatable2.(subscanname)(end) = sumgoodchannels;
                        chtable2.(subscanname)(:) = chtable2.(subscanname)(:) + goodchannels';
                    end
                end
                       
            else  
                scandir = dir(strcat(preprocdir,filesep,group,filesep,subjname,filesep,'*_nonoisych.mat'));
                if exist(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scandir(1).name),'file')
                    load(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scandir(1).name)) 
                    goodchannels = ~isnan(z_oxy(1,:));
                    sumgoodchannels = sum(goodchannels);
                    qatable1.scan(end) = sumgoodchannels;
                    chtable1.scan(:) = chtable1.scan(:) + goodchannels';
                end
                scandir = dir(strcat(preprocdir,filesep,group,filesep,subjname,filesep,'*_nouncertainch.mat'));
                if exist(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scandir(1).name),'file')
                    load(strcat(preprocdir,filesep,group,filesep,subjname,filesep,scandir(1).name)) 
                    goodchannels = ~isnan(z_oxy(1,:));
                    sumgoodchannels = sum(goodchannels);
                    qatable2.scan(end) = sumgoodchannels;
                    chtable2.scan(:) = chtable2.scan(:) + goodchannels';
                end
            end
        end
    end

else
    for i=1:length(currdir)
        msg = sprintf('\n\t subject number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        subjname = currdir(i).name;
        subjnamelength = length(subjname);
        qatable_copy.subjname{1} = subjname;
        qatable1 = [qatable1; qatable_copy];
        qatable2 = [qatable2; qatable_copy];
        subjdir = dir(strcat(preprocdir,filesep,subjname,filesep,dataprefix,'*'));

        if multiscan
            for k=1:length(subjdir)
                scanname = subjdir(k).name;
                subscanname = scanname(subjnamelength+1:end);
                subscanname = regexprep(subscanname,'_','');
                scandir = dir(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,'*_nonoisych.mat'));
                if exist(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,scandir(1).name),'file')
                    load(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,scandir(1).name))
                    goodchannels = ~isnan(z_oxy(1,:));
                    sumgoodchannels = sum(goodchannels);
                    qatable1.(subscanname)(end) = sumgoodchannels;
                    chtable1.(subscanname)(:) = chtable1.(subscanname)(:) + goodchannels';
                end
                scandir = dir(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,'*_nouncertainch.mat'));
                if exist(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,scandir(1).name),'file')
                    load(strcat(preprocdir,filesep,subjname,filesep,scanname,filesep,scandir(1).name))
                    goodchannels = ~isnan(z_oxy(1,:));
                    sumgoodchannels = sum(goodchannels);
                    qatable2.(subscanname)(end) = sumgoodchannels;
                    chtable2.(subscanname)(:) = chtable2.(subscanname)(:) + goodchannels';
                end
            end
                       
        else  
            scandir = dir(strcat(preprocdir,filesep,subjname,filesep,'*_nonoisych.mat'));
            if exist(strcat(preprocdir,filesep,subjname,filesep,scandir(1).name),'file')
                load(strcat(preprocdir,filesep,subjname,filesep,scandir(1).name)) 
                goodchannels = ~isnan(z_oxy(1,:));
                sumgoodchannels = sum(goodchannels);
                qatable1.scan(end) = sumgoodchannels;
                chtable1.scan(:) = chtable1.scan(:) + goodchannels';
            end
            scandir = dir(strcat(preprocdir,filesep,subjname,filesep,'*_nouncertainch.mat'));
            if exist(strcat(preprocdir,filesep,subjname,filesep,scandir(1).name),'file')
                load(strcat(preprocdir,filesep,subjname,filesep,scandir(1).name)) 
                goodchannels = ~isnan(z_oxy(1,:));
                sumgoodchannels = sum(goodchannels);
                qatable2.scan(end) = sumgoodchannels;
                chtable2.scan(:) = chtable2.scan(:) + goodchannels';
            end
        end
    end
end

qaoutpath_noisy=strcat(preprocdir,filesep,'QAreport_allsubj_cleanchannels.csv');
qaoutpath_uncertain=strcat(preprocdir,filesep,'QAreport_allsubj_cleanandcertainchannels.csv');
writetable(qatable1,qaoutpath_noisy,'Delimiter',',');
writetable(qatable2,qaoutpath_uncertain,'Delimiter',',');

choutpath_noisy=strcat(preprocdir,filesep,'QAreport_allch_cleanchannels.csv');
choutpath_uncertain=strcat(preprocdir,filesep,'QAreport_allch_cleanandcertainchannels.csv');
writetable(chtable1,choutpath_noisy,'Delimiter',',');
writetable(chtable1,choutpath_uncertain,'Delimiter',',');

end
