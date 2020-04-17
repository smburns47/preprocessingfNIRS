%--!--  STILL IN DEV DON'T USE --!--%


function [ paths_list,subj ] = findDataPaths(passeddir,dataprefix)
%utility function - find all datafile paths regardless of file structure

dyads=NaN;
multiscan=NaN;
paths_list = cell(0,1);
subj = cell(0,1);

currdir=dir(strcat(passeddir,filesep,dataprefix,'*'));
if length(currdir)<1
    error(['ERROR: No data files found with ',dataprefix,' prefix']);
end


for i=1:length(passeddir)
    currsubj=currdir(i).name;
    subj = {subj; currsubj};
    paths_list = folderDive(strcat(passeddir,filesep,currsubj),paths_list);
    
%     subj1folder = strcat(passeddir,filesep,subj);
%     subj1dir=dir(strcat(passeddir,filesep,subj,filesep,dataprefix,'*'));
%     msg = sprintf('\n\t dyad number %d/%d ...',i,length(currdir));
%     fprintf([reverseStr,msg]);
%     reverseStr = repmat(sprintf('\b'),1,length(msg));      
    
end

end

function [paths] = folderDive(passeddir, priorpaths)
currdir = dir(passeddir);
paths = priorpaths;
if identifyRawNirx(currdir)
    paths = {paths; passeddir};
%elseif identifyRawTechEn(nextlevelfiles)
    %paths = {paths; strcat(passeddir, filesep, nextlevel)};
else
    for i = 1:size(currdir,1)
        paths2 = cell(0,1);
        nextlevelpath = strcat(passeddir,filesep,currdir(i).name);
        if isdir(nextlevelpath)
           paths2 = folderDive(nextlevelpath, paths) ;
           paths = {paths; paths2};
        end
    end
end
end
