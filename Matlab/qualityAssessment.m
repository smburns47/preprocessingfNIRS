function qamask = qualityAssessment(data,samprate,qamethod,thresh)

thresh_corr_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.03,0.08,0.13,0.22,0.27,0.31,0.34,0.38,0.42,0.47,0.52,0.58,0.66,0.75,0.95,1,1,1,1,1});
thresh_ps_map = containers.Map({0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},...
    {0,0.08,0.21,0.26,0.30,0.35,0.39,0.43,0.47,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96});
if ~exist('thresh','var')
    thresh=0.1;
end
thresh_corr = thresh_corr_map(thresh);
thresh_ps = thresh_ps_map(thresh);

qamask = zeros(1,size(data,2));

for k=1:size(data,2)
    trace = data(:,k);
    if ~any(isnan(trace))
        trace = zscore(trace);
        trace_orig = trace;
        offset=round(samprate);
        bigchanges = zeros(length(trace),1);
        for datapoint=((offset*3)+1):(length(trace)-((offset*3)+1))
           if abs(trace(datapoint-offset,1)-trace(datapoint,1))>3
               bigchanges(datapoint-(offset*3):datapoint+(offset*3),1) = 1;
           end
        end
        lstMs = find(diff(bigchanges(:,1))==1);
        lstMf = find(diff(bigchanges(:,1))==-1);
        if length(lstMs) ~= length(lstMf)
            if length(lstMf)>length(lstMs)
               lstMs = [1;lstMs];
            end
            if length(lstMf)<length(lstMs)
               lstMf(end+1,1) = size(bigchanges,1)-1;
            end
        end
        if ~isempty(lstMs)
            for segment = 1:length(lstMs)
                tracelength = length(trace(lstMs(segment):lstMf(segment),1));
                if lstMs(1) == 1
                    pretrace = nanmean(trace(lstMf(1):lstMf(1)+round(10*samprate)));
                    trace(lstMs(segment):lstMf(segment)) = linspace(pretrace,trace(lstMf(segment)+1),tracelength);
                else
                    trace(lstMs(segment):lstMf(segment)) = linspace(trace(lstMs(segment)-1),trace(lstMf(segment)+1),tracelength);
                end
            end
        end
        
        if strcmp(qamethod,'corr')
            autocorrdiff = corrcoef(trace_orig, trace);
            autocorrdiff = autocorrdiff(1,2);
            if autocorrdiff>=(1-thresh_corr)
               qamask(1,k) = 1;
            end

        elseif strcmp(qamethod,'ps')
            PS1 = angle(hilbert(trace));
            PS2 = angle(hilbert(trace_orig));
            avgPS = nanmean(1-sin(abs(PS1-PS2)/2),1);
            if avgPS>=(1-thresh_ps)
                qamask(1,k) = 1;
            end

        else
            error('ERROR: Invalid quality assessment method argument (in qualityAssessment)');
        end
    end
end