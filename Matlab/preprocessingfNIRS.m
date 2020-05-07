
%path must be in brain_data folder (where your raw data folders are)
%required inpaint_nans, homer2 scripts 

function preprocessingfNIRS(dataprefix, hyperscan, multiscan)
%inputs: 
%       dataprefix: string. Prefix of every folder name that should be considered a
%       data folder. E.g., ST for ST_101, ST_102, etc.  
%       hyperscan: 0 or 1. 1 if hyperscanning, 0 if single subject.
%       multiscan: 0 or 1. 1 if multiple scans per person, 0 if single
%       scan.
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

if hyperscan
    if multiscan
        preprocessHyperMulti(dataprefix, currdir, rawdir);
    else
        preprocessHyperSingle(dataprefix, currdir, rawdir);
    end
else
    if multiscan
        preprocessSoloMulti(dataprefix, currdir, rawdir);
    else
        preprocessSoloSingle(dataprefix, currdir, rawdir);
    end
end

end
