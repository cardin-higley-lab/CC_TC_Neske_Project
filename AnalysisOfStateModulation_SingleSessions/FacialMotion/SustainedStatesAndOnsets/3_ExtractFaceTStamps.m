function [tStamps]=extractFaceTStamps(params,wheel,proc,timestamps, timing,indivFigureFolder)
%% This function extract face high arousal on/off transition timestamps 
%params: parameters
%wheel: contains wheel speed
%proc: contains face data
%timestamps/timing: camera/wheel times
%indivFiguresFolder: where plots will be saved 
%Sweyta Lohani 2022 
date = '';
mouse = '';
session = '';
workingDir = '';
%% get face data 
face_Norm=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_faceNorm.mat'))));

% get camera sync times and frame rate
pupil_time = cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_pupilTime.mat'))));
params.fspupilcam = cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_camFrameRate.mat'))));

%% get quiescent bouts
sitOn_final=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_sitOnFinal.mat'))));
sitOff_final=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_sitOffFinal.mat'))));
%% other parameters
%analysis window for events (face and pupil) during quiescence
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset such as airpuff
params.minSitDuration=5;%minimum sit duration during quiescnece state
params.minArousalDuration=0; %minimum face/pupil arousal state (high or low arousal)
params.TimeSinceArousalOn=0;%for face/pupil arousal state, minimum time since  onset
params.TimeBeforeArousalOff=0;%for face/pupil arousal state, minimum time before  offset
params.TimeBeforeArousalOn=4;%for face/pupil arousal state, minimum time of no/low arousal before high arousal state
%% do change point detection on face PC1

%get z thresholds based on face data during quiescence only
b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
zthres_High=quantile(face_Norm,0.60);
zthres_Low=quantile(face_Norm,0.40);

%get on and off timestamps for high and low face movment
[h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
title('FaceHighArousal');
[h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
title('FaceLowArousal');


toDelete=ones(1,length(Face_HighArousal_OnTStamp));
for rj=1:length(Face_HighArousal_OnTStamp)
    tmp = find (Face_HighArousal_OnTStamp(rj)>=sitOn_final & Face_HighArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
Face_HighArousal_On_int=Face_HighArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
Face_HighArousal_Off_int=Face_HighArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

toDelete=ones(1,length(Face_LowArousal_OnTStamp));
for rj=1:length(Face_LowArousal_OnTStamp)
    tmp = find (Face_LowArousal_OnTStamp(rj)>=sitOn_final & Face_LowArousal_OffTStamp(rj)<=sitOff_final);
    toDelete(rj)=isempty(tmp);
end
Face_LowArousal_On_int=Face_LowArousal_OnTStamp(~toDelete)+params.TimeSinceArousalOn;
Face_LowArousal_Off_int=Face_LowArousal_OffTStamp(~toDelete)+params.TimeBeforeArousalOff;

%get high face/pupil states that are preceded by low face state
finalIndex=zeros(1,length(Face_HighArousal_On_int));
for ll=1:length(Face_HighArousal_On_int)
    currOnPoint=Face_HighArousal_On_int(ll);
    tmp=currOnPoint-Face_LowArousal_Off_int;
    idxtmp=find(tmp>=0 &tmp<=0.5); %find states where low arousal ends at less than 0.5s before high arousal onset;
    lowarousallength=Face_LowArousal_Off_int(idxtmp)-Face_LowArousal_On_int(idxtmp);
    if length(lowarousallength)>1 || isempty(lowarousallength), continue, end
    if lowarousallength>=params.TimeBeforeArousalOn
        finalIndex(ll)=1;
        %Face_HighArousal_On_int(ll)=Face_LowArousal_Off_int(idxtmp); %replace aorusal on time with low arousal off time plus two samples
    end
end
Face_HighArousal_On_int1=Face_HighArousal_On_int(logical(finalIndex));
Face_HighArousal_Off_int1=Face_HighArousal_Off_int(logical(finalIndex));


%only select identified high arousal states are at least minimum criterion s long
minimumDuration = params.TimeSinceArousalOn+params.TimeBeforeArousalOff+params.minArousalDuration;
idx3=find((Face_HighArousal_Off_int1-Face_HighArousal_On_int1)>=minimumDuration);
Face_On_final=Face_HighArousal_On_int1(idx3); Face_Off_final=Face_HighArousal_Off_int1(idx3);
%
%% Export results
save(strcat("F:/", date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_whiskOnsetTimes.txt'),'Face_On_final','-ascii');
save(strcat("F:/", date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_whiskOffsetTimes.txt'),'Face_Off_final','-ascii');
%dlmwrite(strcat(date, '_', mouse, "_", session, '_whiskOnsetTimes.txt'), Face_On_final, 'precision', '%.6f')
%dlmwrite(strcat(date, '_', mouse, "_", session, '_whiskOffsetTimes.txt'), Face_Off_final, 'precision', '%.6f')

end