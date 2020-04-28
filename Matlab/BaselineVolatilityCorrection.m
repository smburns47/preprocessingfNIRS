function newd = BaselineVolatilityCorrection(d, samprate, SD, tIncCh)

mlAct = SD.MeasListAct; % prune bad channels

newd=d;

lstAct = find(mlAct==1);
fullstd = nanstd(d);

for ii = 1:length(lstAct)
    idx_ch = lstAct(ii);
    lstMA = find(tIncCh(:,idx_ch)==0);   % sublist of motion artifact segments
    
    if ~isempty(lstMA)
        
        % Find indexes of starts and ends of MA segments
        lstMs = find(diff(tIncCh(:,idx_ch))==-1);   % starting indexes of mvt segments
        lstMf = find(diff(tIncCh(:,idx_ch))==1);    % ending indexes of mvt segments
            
        % Case where there's a single MA segment, that either starts at the
        % beginning or ends at the end of the total time duration
        if isempty(lstMf)
            lstMf = size(tIncCh,1);
        end
        if isempty(lstMs)
            lstMs = 1;
        end
        % If any MA segment either starts at the beginning or
        % ends at the end of the total time duration
        if lstMs(1)>lstMf(1)
            lstMs = [1;lstMs];
        end
        if lstMs(end)>lstMf(end)
            lstMf(end+1,1) = size(tIncCh,1);
        end
        
        lstMl = lstMf-lstMs;    % lengths of MA segments
        nbMA = length(lstMl);   % number of MA segments
        
        % Do baseline and volatility correction on each MA segment
        % only include channels in the active meas list
        
        for jj = 1:nbMA
            lst = lstMs(jj):(lstMf(jj));
            dataseg = newd(lst,lstAct(ii));
            badstd = std(dataseg);
            badmean = nanmean(dataseg,1);
            prevmean = nanmean(newd(1:lst(1)-1,lstAct(ii)),1);
            if isnan(prevmean)
                prevmean = nanmean(newd(:,lstAct(ii)),1);
            end
            adjusted = prevmean + ((dataseg - badmean) .* (fullstd(1,lstAct(ii))/badstd));
            newd(lst,lstAct(ii)) = adjusted;
            
            newmean = nanmean(newd(1:lst(end),lstAct(ii)),1);
            if lst(end)~=size(newd,1)
                if length(newd(lst(end)+1:end,lstAct(ii)))<=30*round(samprate)
                    nextmean = nanmean(newd(lst(end)+1:end,lstAct(ii)),1);
                else
                    nextmean = nanmean(newd(lst(end)+1:lst(end)+(30*round(samprate)+1),lstAct(ii)),1);
                end
                meandiff = newmean-nextmean;
                newd(lst(end)+1:end,lstAct(ii)) = newd(lst(end)+1:end,lstAct(ii)) + meandiff;
            end
        end
    end
end