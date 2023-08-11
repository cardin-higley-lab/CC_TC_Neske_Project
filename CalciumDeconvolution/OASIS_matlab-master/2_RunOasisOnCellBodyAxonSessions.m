%Before running the code below, run "oasis_setup.m"

date = '';
mouse = '';
session = '';

numCBs = size(dir([fullfile('S:/','Imaging','Garrett','FMB208_2PRig',date,mouse,session,'dFFsCellBodies') '/*.mat']),1);

dffCB = cell(numCBs);
for iCB = 1:numCBs
 dffCB{iCB} = cell2mat(struct2cell(load(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','dFFsCellBodies/','dffCB_ROI', num2str(iCB),'.mat'))));
end

numAxons = size(dir([fullfile('S:/','Imaging','Garrett','FMB208_2PRig',date,mouse,session,'dFFsAxons') '/*.mat']),1);

dffAxon = cell(numAxons);
for iAxon = 1:numAxons
 dffAxon{iAxon} = cell2mat(struct2cell(load(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','dFFsAxons/','dffAxon_ROI', num2str(iAxon),'.mat'))));
end

deconvCB = cell(numCBs);
for iCB = 1:numCBs
            deconvCB{iCB} = deconvolveCa(dffCB{iCB}, 'foopsi', 'ar1', 'optimize_pars', true, 'optimize_b', true );  
end

deconvAxon = cell(numAxons);
for iAxon = 1:numAxons
            deconvAxon{iAxon} = deconvolveCa(dffAxon{iAxon}, 'foopsi', 'ar1', 'optimize_pars', true, 'optimize_b', true );  
end

for iCB = 1:numCBs
           [pksCB{iCB},locsCB{iCB}] = findpeaks(deconvCB{iCB});  
end

for iAxon = 1:numAxons
           [pksAxon{iAxon},locsAxon{iAxon}] = findpeaks(deconvAxon{iAxon});  
end

mkdir (strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsCellBodies/'));

for iCB = 1:numCBs
    temp1 = deconvCB{iCB};
    temp2 = locsCB{iCB};
    save(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsCellBodies/','deconvDFFcb',num2str(iCB),'.txt'),'temp1','-ascii');
    save(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsCellBodies/','deconvSpikecb',num2str(iCB),'.txt'),'temp2','-ascii');
end

mkdir (strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsAxons/'));

for iAxon = 1:numAxons
    temp1 = deconvAxon{iAxon};
    temp2 = locsAxon{iAxon};
    save(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsAxons/','deconvDFFaxon',num2str(iAxon),'.txt'),'temp1','-ascii');
    save(strcat('S:/','Imaging','/','Garrett','/','FMB208_2PRig','/',date,'/',mouse,'/',session,'/','deConvdFFsAxons/','deconvSpikeaxon',num2str(iAxon),'.txt'),'temp2','-ascii');
end


