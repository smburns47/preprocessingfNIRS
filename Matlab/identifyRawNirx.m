%--!--  STILL IN DEV DONT USE  --!--%

function [ nirxflag ] = identifyRawNirx( passedfiles )
%Identify if passed set of files are raw Nirx data
%   Detailed explanation goes here
nirxflag=0;
wl1=0;
wl2=0;
hdr=0;
probeInfo=0;
fname='';
for x=1:size(passedfiles,1)
    fname = passedfiles(x).name;
    if strfind(fname,'.wl1')
        wl1=1;
    end
    if strfind(fname,'.wl2')
        wl2=1;
    end
    if strfind(fname,'.hdr')
        hdr=1;
    end
    if strfind(fname,'_probeInfo.mat')
        probeInfo=1;
    end
end
if wl1 && wl2 && hdr && probeInfo
    nirxflag=1;
end

end

