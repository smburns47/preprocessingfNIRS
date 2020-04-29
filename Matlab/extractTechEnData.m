%-----------------------------------------------
%step 1 - subfunctions for extracting data
%-----------------------------------------------

function [d, samprate, s, SD, aux, t] = extractTechEnData(subjfolder)
    nirsfile = dir(strcat(subjfolder,'/*.nirs'));
    loadedsubjfile = load(strcat(subjfolder,'/',nirsfile(1).name),'-mat');
    d = loadedsubjfile.d;
    s = loadedsubjfile.s;
    SD = loadedsubjfile.SD;
    aux = loadedsubjfile.aux;
    t = loadedsubjfile.t;
    samprate = 1/mean(diff(t(:)));
end