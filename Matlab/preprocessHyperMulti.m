function preprocessHyperMulti(dataprefix, currdir, rawdir, device)

fprintf('\n\t Preprocessing ...\n')
reverseStr = '';
Elapsedtime = tic;

scannames = {};

for i=1:length(currdir)
    group=currdir(i).name;  
    groupdir=dir(strcat(rawdir,filesep,group,filesep,dataprefix,'*'));

    for j=1:length(groupdir)
        subjname = groupdir(j).name;
        subjdir = dir(strcat(rawdir,filesep,group,filesep,subjname,filesep,dataprefix,'*'));
            
            for k=1:length(subjdir)
                scanname = subjdir(k).name;
                msg = sprintf('\n\t group %d/%d, subj %d/%d, scan %d/%d ...',i,length(currdir),...
                    j,length(groupdir),k,length(subjdir));
                fprintf([reverseStr,msg]);
                reverseStr = repmat(sprintf('\b'),1,length(msg)); 
                subscanname = scanname(length(subjname)+1:end);
                subscanname = regexprep(subscanname,'_','');
                if ~any(strcmp(scannames,subscanname))
                    scannames = [scannames,subscanname];
                end
                scanfolder = strcat(rawdir,filesep,group,filesep,subjname,filesep,scanname);
                outpath = strcat(rawdir,filesep,'PreProcessedFiles',filesep,group,filesep,subjname,filesep,scanname);

                if ~exist(outpath,'dir')

                %1) extract data values
            if device==1
                [d, sd_ind, samprate, wavelengths, s] = extractNIRxData(scanfolder);
                probenumchannels = probeInfo.probes.nChannel0;
                datanumchannels = size(d,2)/2;
                if probenumchannels~=datanumchannels
                    error('ERROR: number of data channels in hdr file does not match number of channels in probeInfo file.');
                end
            elseif device==2
                [d, samprate, s, SD, aux, t] = extractTechEnData(scanfolder);
            end
            digfile = strcat(scanfolder,filesep,'digpts.txt');
            if device==2 && exist(digfile,'file')
                mni_ch_table = getMNIcoords(digfile, SD);
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

            %3) identify noisy channels
            satlength = 2; %in seconds
            QCoDthresh = 0.1;
            [d, channelmask] = removeBadChannels(d, samprate, satlength, QCoDthresh);
            if device==1
                [SD, aux, t] = getMiscNirsVars(d, sd_ind, samprate, wavelengths, probeInfo, totalmask);
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
            [dconverted, dnormed] = fNIRSFilterPipeline(d, SD, samprate);

            %5) final data quality assessment, remove uncertain channels
            % default is to use phase synchrony to check how impactful
            % remaining spikes are - can change to "corr" as well
            % default QA threshold is 0.1 - amount of measurement error
            % to be allowed in data (out of 1)
            % right now quality assessment is only run on the oxy values
            %TO DO?: in future implementation with GUI, ask usr which signal
            %they plan on analyzing (z scored or no, chromophore)
            qamethod = 'ps';
            thresh = 0.1;
            qamask = qualityAssessment(squeeze(dconverted(:,1,:)),samprate,qamethod,thresh);
            z_qamask = qualityAssessment(squeeze(dnormed(:,1,:)),samprate,qamethod,thresh);
            
            %6) Output results
            mkdir(outpath)
            
            totalmask = channelmask;
            totalmask(~qamask) = 0;
            z_totalmask = channelmask;
            z_totalmask(~z_qamask) = 0;

            numchannels = size(dconverted,3);
            oxy = zeros(size(dconverted,1), numchannels);
            deoxy = zeros(size(dconverted,1), numchannels);
            totaloxy = zeros(size(dconverted,1), numchannels);
            z_oxy = zeros(size(dnormed,1), numchannels);
            z_deoxy = zeros(size(dnormed,1), numchannels);
            z_totaloxy = zeros(size(dnormed,1), numchannels);
            new_d = zeros(size(dconverted,1), numchannels*2);
            for c = 1:numchannels
                oxy(:,c) = dconverted(:,1,c);
                deoxy(:,c) = dconverted(:,2,c);
                totaloxy(:,c) = dconverted(:,3,c);
                z_oxy(:,c) = dnormed(:,1,c);
                z_deoxy(:,c) = dnormed(:,2,c);
                z_totaloxy(:,c) = dnormed(:,3,c);
                new_d(:,(c*2)-1) = oxy(:,c);
                new_d(:,c*2) = deoxy(:,c);
            end
            save(strcat(outpath,filesep,group,'_',subjname,'_preprocessed.mat'),'oxy', 'deoxy', 'totaloxy','z_oxy', 'z_deoxy', 'z_totaloxy','s','samprate','t','SD');
            
            oxy(:,~channelmask) = NaN;
            deoxy(:,~channelmask) = NaN;
            totaloxy(:,~channelmask) = NaN;
            z_oxy(:,~channelmask) = NaN;
            z_deoxy(:,~channelmask) = NaN;
            z_totaloxy(:,~channelmask) = NaN;
            save(strcat(outpath,filesep,group,'_',subjname,'_preprocessed_nonoisych.mat'),'oxy', 'deoxy', 'totaloxy','z_oxy', 'z_deoxy', 'z_totaloxy','s','samprate','t','SD');
            
            oxy(:,~totalmask) = NaN;
            deoxy(:,~totalmask) = NaN;
            totaloxy(:,~totalmask) = NaN;
            z_oxy(:,~z_totalmask) = NaN;
            z_deoxy(:,~z_totalmask) = NaN;
            z_totaloxy(:,~z_totalmask) = NaN;
            save(strcat(outpath,filesep,group,'_',subjname,'_preprocessed_nouncertainch.mat'),'oxy', 'deoxy', 'totaloxy','z_oxy', 'z_deoxy', 'z_totaloxy','s','samprate','t','SD');
            writetable(mni_ch_table,strcat(outpath,filesep,'channel_mnicoords.csv'),'Delimiter',',');
            end
        end
    end
end

preprocdir = strcat(rawdir,filesep,'PreProcessedFiles');
qualityReport(dataprefix,1,1,scannames,numchannels,preprocdir);

    
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds\n', Elapsedtime);

end