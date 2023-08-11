function[wheelOn_final,wheelOff_final,Face_HighArousal_On_final,Face_HighArousal_Off_final,Face_LowArousal_On_final,Face_LowArousal_Off_final]...
    =getSustainedStates(timing,params,imaging_time,pupil_time,face_Norm,files,indivFigureFolder)
%this function gets sustained locomotion, face high, and face low states 
%written by SL, 2022
date = '';
mouse = '';
session = '';
workingDir = '';
%% get face data 
face_Norm=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_faceNorm.mat'))));
%% 

% get camera sync times and frame rate
pupil_time = cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_pupilTime.mat'))));
params.fspupilcam = cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_camFrameRate.mat'))));
%% get quiescent bouts
sitOn_final=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_sitOnFinal.mat'))));
sitOff_final=cell2mat(struct2cell(load(strcat(workingDir, date, '/', mouse, '/', session, '/', date, '_', mouse, "_", session, '_sitOffFinal.mat'))));
%% Other parameters
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=10;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=10;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=10;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state, could reduce this to 2s if you are doing pearson's correlation (you can extract more states with this duration)
params.minArousalDuration=5; %minimum face/pupil arousal state (high or low arousal),could reduce this to 2s if you are doing pearson's correlation (you can extract more states with this duration)
params.minSitDuration=5;%minimum sit duration during quiescnece state,could reduce this to 2s if you are doing pearson's correlation (you can extract more states with this duration)
params.plotStates=1; %indicate whether a plot of behavioral state timestamps should be generated

%% do change point detection on face to get face high/low movement times during sustained quiescence state
%get Z-thresholds based on face data during quiescence, when mouse isn't moving and when aripuffs are not given
    b1DoPlot=1; blDoPlotDuration=1:15000; smoothWin=1;
   
    zthres_High=quantile(face_Norm,0.60); %upper 60% quantile for face high 
    zthres_Low=quantile(face_Norm,0.40);%lower 40% quantile for face low
    
    %get on and off timestamps for high and low face movment
    [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
    title('FaceHighArousal');
    [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
    title('FaceLowArousal');
    
    %remove outliers (such as caused by grooming) in high face state data by removing states where the state values are greater than 4SD from the whole session average
    znormedface=normalize(face_Norm);   
    maxzdata=nan(1,length(Face_HighArousal_OnTStamp));
    for st=1:length(Face_HighArousal_OnTStamp)
        OnIndx=find(pupil_time==Face_HighArousal_OnTStamp(st));
        OffIndx=find(pupil_time==Face_HighArousal_OffTStamp(st));
        zdata=znormedface(OnIndx:OffIndx);
        maxzdata(st)=nanmax(zdata);
    end
    
    Face_HighArousal_OnTStamp=Face_HighArousal_OnTStamp(maxzdata<4);
    Face_HighArousal_OffTStamp=Face_HighArousal_OffTStamp(maxzdata<4);
    
    % get  face on times if both on and off times occur during
    % sustained queiscence states identified in the previous step. 
    %If on/off times occur beyond the quiescence state, modify the on/off times to be during the quiescence only
    s1=sitOn_final; s2=sitOff_final;
    on_final=cell(1,length(Face_HighArousal_OnTStamp)); off_final=cell(1,length(Face_HighArousal_OnTStamp));
    for rj=1:length(Face_HighArousal_OnTStamp)
        a1=Face_HighArousal_OnTStamp(rj); a2=Face_HighArousal_OffTStamp(rj);
        A=a1-s1; B=a2-s2;
        on_f=nan(1,length(sitOn_final)); off_f=nan(1,length(sitOff_final));
        for rt=1:length(sitOn_final)
            if A(rt)>0,Ac=A(rt); else Ac=0; end
            if B(rt)<0,Bc=abs(B(rt)); else Bc=0; end
            on= s1(rt)+Ac;
            off=s2(rt)-Bc;
            if (off-on)>0,on_f(rt)=on; off_f(rt)=off; end
        end
        on_final{rj}=on_f(find(~isnan(on_f))); off_final{rj}=off_f(find(~isnan(on_f)));
    end
    Face_HighArousalOn_int1=cell2mat(on_final); Face_HighArousalOff_int1=cell2mat(off_final);
    
    on_final=cell(1,length(Face_LowArousal_OnTStamp)); off_final=cell(1,length(Face_LowArousal_OnTStamp));
    for rj=1:length(Face_LowArousal_OnTStamp)
        a1=Face_LowArousal_OnTStamp(rj); a2=Face_LowArousal_OffTStamp(rj);
        A=a1-s1; B=a2-s2;
        on_f=nan(1,length(sitOn_final)); off_f=nan(1,length(sitOff_final));
        for rt=1:length(sitOn_final)
            if A(rt)>0,Ac=A(rt); else Ac=0; end
            if B(rt)<0,Bc=abs(B(rt)); else Bc=0; end
            on= s1(rt)+Ac;
            off=s2(rt)-Bc;
            if (off-on)>0,on_f(rt)=on; off_f(rt)=off; end
        end
        on_final{rj}=on_f(find(~isnan(on_f))); off_final{rj}=off_f(find(~isnan(on_f)));
    end
    Face_LowArousalOn_int1=cell2mat(on_final); Face_LowArousalOff_int1=cell2mat(off_final);
    
    % determine that face high/low arousal/movement time are at least minimum criterion seconds long
    idx3=find((Face_HighArousalOff_int1-Face_HighArousalOn_int1)>=params.minArousalDuration);
    Face_HighArousal_On_final=Face_HighArousalOn_int1(idx3); Face_HighArousal_Off_final=Face_HighArousalOff_int1(idx3);
    
    idx3=find((Face_LowArousalOff_int1-Face_LowArousalOn_int1)>=params.minArousalDuration);
    Face_LowArousal_On_final=Face_LowArousalOn_int1(idx3); Face_LowArousal_Off_final=Face_LowArousalOff_int1(idx3);
    %% Export data
    %dlmwrite(strcat(date, '_', mouse, "_", session, '_whiskStateOnsetTimes.txt'), Face_HighArousal_On_final, 'precision', '%.6f')
    %dlmwrite(strcat(date, '_', mouse, "_", session, '_whiskStateOffsetTimes.txt'), Face_HighArousal_Off_final, 'precision', '%.6f')
    
    %dlmwrite(strcat(date, '_', mouse, "_", session, '_noWhiskStateOnsetTimes.txt'), Face_LowArousal_On_final, 'precision', '%.6f')
    %dlmwrite(strcat(date, '_', mouse, "_", session, '_noWhiskStateOffsetTimes.txt'), Face_LowArousal_Off_final, 'precision', '%.6f')
    Face_LowArousal_On_Trans = transpose(Face_LowArousal_On_final);
    Face_LowArousal_Off_Trans = transpose(Face_LowArousal_Off_final);
    
    Face_HighArousal_On_Trans = transpose(Face_HighArousal_On_final);
    Face_HighArousal_Off_Trans = transpose(Face_HighArousal_Off_final);
    
    save(strcat("F:/", date, '/', mouse, '/', session, '/',date, '_', mouse, "_", session, '_whiskStateOnsetTimes.txt'), 'Face_HighArousal_On_Trans', '-ascii')
    save(strcat("F:/", date, '/', mouse, '/', session, '/',date, '_', mouse, "_", session, '_whiskStateOffsetTimes.txt'), 'Face_HighArousal_Off_Trans', '-ascii')
    
    save(strcat("F:/", date, '/', mouse, '/', session, '/',date, '_', mouse, "_", session, '_noWhiskStateOnsetTimes.txt'), 'Face_LowArousal_On_Trans', '-ascii')
    save(strcat("F:/", date, '/', mouse, '/', session, '/',date, '_', mouse, "_", session, '_noWhiskStateOffsetTimes.txt'), 'Face_LowArousal_Off_Trans', '-ascii')
