%% automates Ca extraction over all tif files in folder
%% Note: run each segment individually
%% run roigui for first file in list, pick ROIs, define input and output directories 
clear all

%define some parameters

inputFileName='file_00001_00001_moco40_ref1-200.tif'; % name of the file you will pick your ROIs from

mainAnalysisDir='C:/Users/garrett/Desktop/Garrett_Local/071223/Mouse23666/Session2/';
mainInputDir=fullfile(mainAnalysisDir,'moco'); %your input directory
%mainOutputDir=fullfile(mainAnalysisDir,'ROIs');  %your output directory
mainOutputDir=fullfile(mainAnalysisDir,'Ftraces');  %your output directory
    if ~exist(mainOutputDir,'dir'), mkdir(mainOutputDir), end 

fsample=30.03;  %sampling frequency 

%%%%%%%%%%%%
    
indx = strfind(inputFileName,'_');
%fileTag= inputFileName(1:indx(4)-1);  % folder name - can replace with any string you want.  e.g., 'MyName_folder'; 
 
inputDir=fullfile(mainInputDir);
outputDir=fullfile(mainOutputDir);
if (~exist(outputDir,'dir')) mkdir(outputDir); end

D=dir([inputDir,'/*.tif']);
lastFileInd = length(D(not([D.isdir])));  %number of tif files in input_dir

input=fullfile(inputDir,inputFileName);

%% if using new ROI list
Fcorr=32768;
%% 
im=readTifStack(input);

im=im-Fcorr; 
roigui(im); 
%% run Ca extraction for first file in list, and refine ROI choices
f_traces_gui_v3(im,ROI_list,fsample);  

%% using ROI list from above, runs Ca extraction over the rest of the .tif files in input_dir
for fileInd=1:lastFileInd
    inputFileNameBatch=D(fileInd).name; 
    if strcmp(inputFileNameBatch,inputFileName), continue, end   
    input=fullfile(inputDir,inputFileNameBatch);
    
    outputFilename=[inputFileNameBatch(1:end-4)]; 
    im=readTifStack(input);
    im=im-Fcorr; 
    roinogui(im,ROI_list);
    f_traces_nogui_v3(im,ROI_list,fsample,outputFilename,outputDir);
end
save([outputDir,'/ROI_list'],'ROI_list');

%%