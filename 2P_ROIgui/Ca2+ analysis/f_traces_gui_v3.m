function [ROI_list] = f_traces_gui_v3(im,ROI_list,frame_rate)
% This function launches a GUI for processing in vivo two-photon calcium
% imaging data.
%
% Inputs:
% im = an MxNxT matrix containing pixel values for the imaging stack that
% corresponds to the ROIs in ROI_list. This can be the output of the
% function 'readTifStack.m' from the Svoboda lab roigui package
% ROI_list = a structure containing information for ROIs selected from the
% imaging data. This structure is an output of the Svoboda lab GUI
% 'roigui.m'
% frame_rate = frame rate of movie 'im', in hertz
%
% Outputs:
% ROI_list = will automatically update the ROI_list structure with
% information about neuropil location and fluorescence for each ROI and the
% neuropil-subtracted raw fluorescence traces for the first pass of the
% function. ROI_list is updated in the MATLAB workspace everytime the
% "Export Data" button is pressed in the GUI Toolbox.
%
% Notes:
% This GUI requires that the following functions be present in the MATLAB
% file path:
%   neuropil_mask.m
%   dFFfunc.m
%       alphatau.mat
%   foopsi_testfunc.m (TO BE UPDATED when greater foospi functionality is added to the GUI)
%       fast_oopsi.m
%       z1.m
%
%
%
%
%
%
%
%
% Last updated: 30 July 2014, dbarson

tic

%Initialize constant variables and figure handles.
max_im=max(im,[],3);
min_im=min(im,[],3);
DispRange=[prctile(min_im(:),25),prctile(max_im(:),75)]; %SL changed this from 0 and 99 to 25 and 75 
DispRange_local = DispRange;
ROI_window_radius = 25;
hdispfig=-1;
htoolboxfig=-1;
hdispraw=-1;
hdisprawim=-1;
hdispimage=-1;
hdispROI=-1;
hdispCentroid=-1;
reordered = 0;
yfrozen = 0;
isplaying = 0;
ss = size(im);
t = 0:1/frame_rate:(ss(3)-1)/frame_rate;
currROI = 1; %initialization
currFrame = 1;
num_ROIs = length(ROI_list);
pixel_dimension=[1,1];

%set up ROI display structure
d.name='Average';
d.data=mean(im,3);
d.playsize=1;
d.DataRange=[min(d.data(:)),max(d.data(:))];
% d.DataRange=dataRange;
d.DispRange=[prctile(d.data(:),55),prctile(d.data(:),99)]; %SL changed this from 0 and 99 to 25 and 75
d.colormap=1;
dispstruct=d;

%Initialize default settings and do first-pass neuropil calculations.
cell_rad = 10; %in pixels
neuropil_rad = 17; %in pixels
neuropil_alpha = 0.7;
ROI_list = neuropil_mask(im, ROI_list,cell_rad,neuropil_rad,neuropil_alpha);

%Initialize default settings and do first-pass dF/F calculation and spike
%deconvolution.
dFF_filter = [];
tau_0 = 0.2;
tau_1 = .75;
tau_2 = 3;
new_t = [];
new_t_offset = [];
dFF = [];
Nhat_dFF = dFF;


%Initialize GUI Toolbox controls
htoolboxfig=figure('name','Toolbox','position',[500,500,300,200],'Resize','off','MenuBar','none','color',[1,1,1]*0.85,'DeleteFcn',@toolboxfigDelCallback);
uicontrol(htoolboxfig,'Tag','ROIOrder','Style','pushbutton','String','Reorder ROIs','Units','pixel','position',[10,170,130,20],'Callback',@ROIOrderCallback);
uicontrol(htoolboxfig,'Tag','FreezeY','Style','pushbutton','String','Freeze Y-axis','Units','pixel','position',[160,170,130,20],'Callback',@FreezeYCallback);
uicontrol(htoolboxfig,'Tag','text1','units','pixel','position',[10,140,130,20],'Style','text','String','Neuropil Params');
uicontrol(htoolboxfig,'Tag','text2','units','pixel','position',[160,140,130,20],'Style','text','String','Filtering Params');
uicontrol(htoolboxfig,'Tag','text3','units','pixel','position',[5,110,90,20],'Style','text','String','Cell Radius');
uicontrol(htoolboxfig,'Tag','text4','units','pixel','position',[5,80,90,20],'Style','text','String','Neuropil Radius');
uicontrol(htoolboxfig,'Tag','text7','units','pixel','position',[5,50,90,20],'Style','text','String','Neuropil Alpha');
uicontrol(htoolboxfig,'Tag','filter_choice','units','pixel','position',[155,110,90,20],'Style','popupmenu','String',{'exponential','gaussian','no filter'},'Callback',@PopupFilterCallback);
uicontrol(htoolboxfig,'Tag','text5','units','pixel','position',[155,80,90,20],'Style','text','String','tau_1');
uicontrol(htoolboxfig,'Tag','text6','units','pixel','position',[155,50,90,20],'Style','text','String','tau_2');
uicontrol(htoolboxfig,'Tag','EdCellRad','Style','edit','Units','pixel','position',[105,110,40,20],'String',num2str(cell_rad),'Callback',@EdCellRadCallback);
uicontrol(htoolboxfig,'Tag','EdNeuropilRad','Style','edit','Units','pixel','position',[105,80,40,20],'String',num2str(neuropil_rad),'Callback',@EdNeuropilRadCallback);
uicontrol(htoolboxfig,'Tag','EdNeuropilAlpha','Style','edit','Units','pixel','position',[105,50,40,20],'String',num2str(neuropil_alpha),'Callback',@EdNeuropilAlphaCallback);
uicontrol(htoolboxfig,'Tag','EdTau0','Style','edit','Units','pixel','position',[255,110,40,20],'String',num2str(tau_0),'Callback',@EdTau0Callback);
uicontrol(htoolboxfig,'Tag','EdTau1','Style','edit','Units','pixel','position',[255,80,40,20],'String',num2str(tau_1),'Callback',@EdTau1Callback);
uicontrol(htoolboxfig,'Tag','EdTau2','Style','edit','Units','pixel','position',[255,50,40,20],'String',num2str(tau_2),'Callback',@EdTau2Callback);
uicontrol(htoolboxfig,'Tag','BnNeuropil','Style','pushbutton','String','Calculate Neuropil','Units','pixel','position',[10,25,130,15],'Callback',@BnNeuropilCallback);
uicontrol(htoolboxfig,'Tag','BnFiltering','Style','pushbutton','String','Calculate dF/F','Units','pixel','position',[160,25,130,15],'Callback',@BnFilteringCallback);
uicontrol(htoolboxfig,'Tag','BnExportData','Style','pushbutton','String','Export Data','Units','pixel','position',[10,5,130,15],'Callback',@BnExportData);
uicontrol(htoolboxfig,'Tag','BnDeleteROI','Style','pushbutton','String','Delete Selected ROI','Units','pixel','position',[160,5,130,15],'Callback',@BnDeleteROI);
htools = guihandles(htoolboxfig);

%Initialize GUI Display axes and controls
hctrlfig=figure('name','Fluorescence Traces','position',[900,100,600,900],'Resize','off','MenuBar','none','Toolbar','figure','color',[1,1,1]*0.85,'DeleteFcn',@ctrlfigDelCallback);
uicontrol(hctrlfig,'Tag','SlROI','Style','slider','Units','pixel','position',[75,870,400,20],'Max',num_ROIs,'Min',1,'Value',currROI,'SliderStep',[1/(num_ROIs-1),0.1],'Callback',@SlROICallback);
uicontrol(hctrlfig,'Tag','EdROI','Style','edit','Units','pixel','position',[485,870,50,20],'String',num2str(currROI),'Callback',@EdROICallback);
uicontrol(hctrlfig,'Tag','SlFrame','Style','slider','Units','pixel','position',[75,840,400,20],'Max',ss(3),'Min',1,'Value',currFrame,'SliderStep',[1/(ss(3)-1),0.1],'Callback',@SlFrameCallback);
uicontrol(hctrlfig,'Tag','EdFrame','Style','edit','Units','pixel','position',[485,840,50,20],'String',num2str(currFrame),'Callback',@EdFrameCallback);
uicontrol(hctrlfig,'Tag','BnPlay','Style','pushbutton','String','Play','Units','pixel','position',[545,850,45,30],'Callback',@BnPlayCallback);
uicontrol(hctrlfig,'Tag','text1','units','pixel','position',[10,870,55,20],'Style','text','String','ROI');
uicontrol(hctrlfig,'Tag','text2','units','pixel','position',[10,840,55,20],'Style','text','String','Frame');
axes('Tag','ax1','units','pixel','position',[50,690,500,115],'Parent',hctrlfig,'NextPlot','add','Visible','on','XTickLabel',[]);
ylabel('rawF')
axes('Tag','ax2','units','pixel','position',[50,555,500,115],'Parent',hctrlfig,'NextPlot','add','XTickLabel',[]);
ylabel('neuropil')
axes('Tag','ax3','units','pixel','position',[50,420,500,115],'Parent',hctrlfig,'NextPlot','add','XTickLabel',[]);
ylabel('subtractedF')
axes('Tag','ax4','units','pixel','position',[50,285,500,115],'Parent',hctrlfig,'NextPlot','add','XTickLabel',[]);
ylabel('dFF')
%axes('Tag','ax5','units','pixel','position',[50,200,500,70],'Parent',hctrlfig,'NextPlot','add');
%ylabel('foopsi dFF')
%xlabel('time (s)')
axes('Tag','im1','units','pixel','position',[40,10,150,150],'Parent',hctrlfig,'NextPlot','add','Visible','on','XTickLabel',[],'YTickLabel',[],'XLim',[0 2*ROI_window_radius+1],'YLim',[0 2*ROI_window_radius+1]);
ylabel('raw_im')
colormap('gray')
axes('Tag','im2','units','pixel','position',[225,10,150,150],'Parent',hctrlfig,'NextPlot','add','Visible','on','XTickLabel',[],'YTickLabel',[],'XLim',[0 2*ROI_window_radius+1],'YLim',[0 2*ROI_window_radius+1]);
ylabel('ROI');
colormap('gray')
axes('Tag','im3','units','pixel','position',[410,10,150,150],'Parent',hctrlfig,'NextPlot','add','Visible','on','XTickLabel',[],'YTickLabel',[],'XLim',[0 2*ROI_window_radius+1],'YLim',[0 2*ROI_window_radius+1]);
ylabel('neuropil');
colormap('gray')
hobjs=guihandles(hctrlfig);

PopupFilterCallback();
Calculate_dFF_foopsi();

%Initialize trace plots
h_rawF = plot(hobjs.ax1,t,ROI_list(currROI).fmean,t(currFrame),ROI_list(currROI).fmean(currFrame),'o','MarkerFaceColor',[0 1 0],'MarkerSize',5);
h_neuropil = plot(hobjs.ax2,t,ROI_list(currROI).neuropil_fmean,t(currFrame),ROI_list(currROI).neuropil_fmean(currFrame),'o','MarkerFaceColor',[0 1 0],'MarkerSize',5);
h_Fsubtracted = plot(hobjs.ax3,t,ROI_list(currROI).F_neuropilsubtracted,t(currFrame),ROI_list(currROI).F_neuropilsubtracted(currFrame),'o','MarkerFaceColor',[0 1 0],'MarkerSize',5);
h_dFF = plot(hobjs.ax4,new_t,dFF(:,currROI),new_t(currFrame),dFF(currFrame,currROI),'o','MarkerFaceColor',[0 1 0],'MarkerSize',5);
%h_dFF_foopsi = plot(hobjs.ax5,new_t,Nhat_dFF(:,currROI),'r',new_t(currFrame),Nhat_dFF(currFrame,currROI),'o','MarkerFaceColor',[0 1 0],'MarkerSize',5);
linkaxes([hobjs.ax1,hobjs.ax2,hobjs.ax3,hobjs.ax4],'x'); %,hobjs.ax5

%Initialize values of Y axis limits for freezing/unfreezing.
[Ylim_raw] = get(hobjs.ax1,'YLim');
[Ylim_npil] = get(hobjs.ax2,'YLim');
[Ylim_fsub] = get(hobjs.ax3,'YLim');
[Ylim_dFF] = get(hobjs.ax4,'YLim');
%[Ylim_dFF_foopsi] = get(hobjs.ax5,'YLim');

%Initialize display of images in the Traces window
local_FOV_x = round(ROI_list(currROI).centerPos(1))+ROI_window_radius:-1:round(ROI_list(currROI).centerPos(1))-ROI_window_radius;
local_FOV_y = round(ROI_list(currROI).centerPos(2))-ROI_window_radius:round(ROI_list(currROI).centerPos(2))+ROI_window_radius;
local_FOV_x(local_FOV_x < 1 | local_FOV_x > ss(1)) = [];
local_FOV_y(local_FOV_y < 1 | local_FOV_y > ss(2)) = [];
h_raw_im = imagesc(im(local_FOV_x,local_FOV_y,currFrame),'Parent',hobjs.im1);
h_avg_im = imagesc(mean(im(local_FOV_x,local_FOV_y,:),3),'Parent',hobjs.im2);
h_neuropil_im = imagesc(mean(im(local_FOV_x,local_FOV_y,:),3),'Parent',hobjs.im3);hold on;
hdispneuropil=image(zeros([length(local_FOV_x),length(local_FOV_y),3]),'AlphaData',zeros([length(local_FOV_x),length(local_FOV_y)]),'Visible','off');

%Initialize two display windows, the first of which shows average of all movie
%frames with ROIs overlaid, the second of which shows raw data for each
%frame.
CreateNewDispWnd();
drawROI();

if  usejava('awt') % java enabled -> use it to update while dragging
%     l1=handle.listener(hobjs.SlROI,'ActionEvent',@SlROICallback);
%     l2=handle.listener(hobjs.SlFrame,'ActionEvent',@SlFrameCallback);
    l1=addlistener(hobjs.SlROI,'ContinuousValueChange',@SlROICallback);
    l2=addlistener(hobjs.SlFrame,'ContinuousValueChange',@SlFrameCallback);
end

%% The following functions are "Math and Plotting" functions. They do operations on the data matrices and structures and update what is displayed in the GUI windows and axes.
    function Calculate_dFF_foopsi()
        new_t = dFFfunc(t,ROI_list(1).fmean,tau_1,tau_2,tau_0,dFF_filter);
        new_t_offset = length(t)-length(new_t);
        dFF = zeros(length(new_t),num_ROIs);
        Nhat_dFF = dFF;
        for j = 1:num_ROIs
            [new_t dFF(:,j)] = dFFfunc(t,ROI_list(j).F_neuropilsubtracted,tau_1,tau_2,tau_0,dFF_filter);
            % Nhat_dFF(:,j) = foopsi_testfunc(dFF(:,j),frame_rate);
        end
    end

    function UpdateDisplay()
        if ~ishandle(hdispfig)
            CreateNewDispWnd();
        end
        if ~ishandle(hdispraw)
            CreateNewDispWnd();
        end
        
        set(h_rawF,{'YData'},{ROI_list(currROI).fmean;ROI_list(currROI).fmean(currFrame)},{'XData'},{t;t(currFrame)}); %axis('tight');
        set(h_neuropil,{'YData'},{ROI_list(currROI).neuropil_fmean;ROI_list(currROI).neuropil_fmean(currFrame)},{'XData'},{t;t(currFrame)}); %axis('tight');
        set(h_Fsubtracted,{'YData'},{ROI_list(currROI).F_neuropilsubtracted;ROI_list(currROI).F_neuropilsubtracted(currFrame)},{'XData'},{t;t(currFrame)}); %axis('tight');
        
        if currFrame <= new_t_offset
            set(h_dFF,{'YData'},{dFF(:,currROI);dFF(1,currROI)},{'XData'},{new_t;new_t(1)}); %axis('tight');
            %set(h_dFF_foopsi,{'YData'},{Nhat_dFF(:,currROI);Nhat_dFF(1,currROI)},{'XData'},{new_t;new_t(1)}); %axis('tight');
        else
            set(h_dFF,{'YData'},{dFF(:,currROI);dFF(currFrame-new_t_offset,currROI)},{'XData'},{new_t;new_t(currFrame-new_t_offset)}); %axis('tight');
            %set(h_dFF_foopsi,{'YData'},{Nhat_dFF(:,currROI);Nhat_dFF(currFrame-new_t_offset,currROI)},{'XData'},{new_t;new_t(currFrame-new_t_offset)}); %axis('tight');
        end
        
        set(h_raw_im,'CData',im(local_FOV_x,local_FOV_y,currFrame));
        set(get(h_raw_im,'parent'),'CLim',DispRange_local);
        set(hdisprawim,'CData',im(:,:,currFrame));
        set(get(hdisprawim,'parent'),'clim',DispRange);
        drawnow();
    end

    function CreateNewDispWnd()
        hdispfig=figure('name','ROIs','position',[50,100,400,400],'DeleteFcn',@dispfigDelCallback,'colormap',gray);
        hdispimage=imagesc(mean(im,3));
        xlabel('y')
        ylabel('x')
        set(get(hdispimage,'parent'),'clim',dispstruct.DispRange);
        daspect([pixel_dimension,1]);
        hold on;
        hdispROI=image(zeros([ss(1:2),3]),'AlphaData',zeros([ss(1:2)]),'Visible','off');
        hdispCentroid = image(zeros([ss(1:2),3]),'AlphaData',zeros([ss(1:2)]),'Visible','off'); hold off;
        hdispraw = figure('name','Raw Movie','position',[50,500,400,400],'DeleteFcn',@disprawDelCallback,'colormap',gray);
        hdisprawim=imagesc(im(:,:,currFrame));
        set(get(hdisprawim,'parent'),'clim',DispRange);
        %hdispRec=rectangle('Position',[0,0,para.region_size*2+1,para.region_size*2+1],'Visible','off','EdgeColor','r','LineStyle','--','parent',hax);
    end

    function drawROI()
        
        roi_display=zeros([prod(ss(1:2)),3]);
        alpha=zeros(prod(ss(1:2)),1);
        
        sel_ind=currROI;
        for i=1:length(ROI_list)
            rgb=[0.5,1,0.5];
            if sum(sel_ind==i)>0
                roi_display(ROI_list(i).pixel_list,1)=rgb(1);
                roi_display(ROI_list(i).pixel_list,2)=rgb(2);
                roi_display(ROI_list(i).pixel_list,3)=rgb(3);
                alpha(ROI_list(i).pixel_list)=1;
            else
                roi_display(ROI_list(i).pixel_list,1)=1;
                roi_display(ROI_list(i).pixel_list,2)=0;
                roi_display(ROI_list(i).pixel_list,3)=0;
                alpha(ROI_list(i).pixel_list)=1;
            end
        end
        alpha=reshape(alpha,ss(1:2));
        
        roi_display=reshape(roi_display,[ss(1:2),3]);
        set(hdispROI,'CData',roi_display,'AlphaData',alpha);
        set(hdispROI,'Visible','on');
        
        ROI_window_radius = round(1.25*neuropil_rad);
        if ROI_window_radius < 25
            ROI_window_radius = 25;
        end
        local_FOV_x = round(ROI_list(currROI).centerPos(1))+ROI_window_radius:-1:round(ROI_list(currROI).centerPos(1))-ROI_window_radius;
        local_FOV_y = round(ROI_list(currROI).centerPos(2))-ROI_window_radius:round(ROI_list(currROI).centerPos(2))+ROI_window_radius;
        local_FOV_x(local_FOV_x < 1 | local_FOV_x > ss(1)) = [];
        local_FOV_y(local_FOV_y < 1 | local_FOV_y > ss(2)) = [];
        
        neuropil_display=zeros(length(local_FOV_x)*length(local_FOV_y),3);
        alpha_np=zeros(length(local_FOV_x)*length(local_FOV_y),1);
        neuropil_y = ceil(ROI_list(currROI).neuropil_pixel_list/ss(1))-min(local_FOV_y)+1;
        neuropil_x = max(local_FOV_x)-rem(ROI_list(currROI).neuropil_pixel_list,ss(1))+1;
        neuropil_pixel_inds = (neuropil_y-1)*length(local_FOV_x)+neuropil_x;
        neuropil_display(neuropil_pixel_inds,1) = 1;
        alpha_np(neuropil_pixel_inds) = 1;
        alpha_np = reshape(alpha_np,[length(local_FOV_x),length(local_FOV_y)]);
        neuropil_display = reshape(neuropil_display,[length(local_FOV_x),length(local_FOV_y),3]);
        set(h_avg_im,'CData',mean(im(local_FOV_x,local_FOV_y,:),3));
        set(h_neuropil_im,'CData',mean(im(local_FOV_x,local_FOV_y,:),3));
        set(hdispneuropil,'CData',neuropil_display,'AlphaData',alpha_np,'Visible','on');
        max_im_local=max(im(local_FOV_x,local_FOV_y,currFrame),[],3);
        min_im_local=min(im(local_FOV_x,local_FOV_y,currFrame),[],3);
        DispRange_local=[prctile(min_im_local(:),0.0),prctile(max_im_local(:),98)];
        drawnow();
    end

%% All of the below functions are Callbacks for the GUI, which can call one of the "Math and Plotting" functions above. 

    function BnPlayCallback(varargin)
        if isplaying==0
            isplaying=1;
            set(hobjs.BnPlay,'String','Stop');
            while isplaying
                
                currFrame = currFrame+1;
                if currFrame>ss(3)
                    currFrame=1;
                    isplaying=0;
                end
                set(hobjs.EdFrame,'String',num2str(currFrame));
                set(hobjs.SlFrame,'Value',currFrame);
                UpdateDisplay();
            end
            set(hobjs.EdFrame,'String',num2str(currFrame));
            set(hobjs.SlFrame,'Value',currFrame);
            set(hobjs.BnPlay,'String','Play');
            UpdateDisplay();
        else
            isplaying=0;
            set(hobjs.BnPlay,'String','Play');
        end
    end

    function FreezeYCallback(varargin)
        if yfrozen == 0;
            yfrozen = 1;
            set(htools.FreezeY,'String','Unfreeze Y');
            [Ylim_raw] = get(hobjs.ax1,'YLim');
            [Ylim_npil] = get(hobjs.ax2,'YLim');
            [Ylim_fsub] = get(hobjs.ax3,'YLim');
            [Ylim_dFF] = get(hobjs.ax4,'YLim');
            %[Ylim_dFF_foopsi] = get(hobjs.ax5,'YLim');
            set(hobjs.ax1,'YLim',Ylim_raw);
            set(hobjs.ax2,'Ylim',Ylim_npil);
            set(hobjs.ax3,'YLim',Ylim_fsub);
            set(hobjs.ax4,'Ylim',Ylim_dFF);
            %set(hobjs.ax5,'YLim',Ylim_dFF_foopsi);
        else
            yfrozen = 0;
            set(htools.FreezeY,'String','Freeze Y');
            set(hobjs.ax1,'YLimMode','auto');
            set(hobjs.ax2,'YLimMode','auto');
            set(hobjs.ax3,'YLimMode','auto');
            set(hobjs.ax4,'YLimMode','auto');
            %set(hobjs.ax5,'YLimMode','auto');
        end
    end
    function ROIOrderCallback(varargin)
        set(htools.ROIOrder,'String','One moment, please'); drawnow();
        sel_point = inputdlg({'X value','Y value','Radius of display dot (recommended: 3)'},'Enter coordinates (pixels) of point from which ROIs will be ordered',1);
        x_val = round(str2num(sel_point{1}));
        y_val = round(str2num(sel_point{2}));
        dot_rad = str2num(sel_point{3});
        if isempty(sel_point)
            return
        end
        if reordered == 0
            unordered_ROI_list = ROI_list;
            assignin('base','unordered_ROI_list',unordered_ROI_list);
            reordered = 1;
        end
        
        ROI_distances = zeros(num_ROIs,1);
        for i = 1:num_ROIs
            ROI_distances(i) = sqrt((ROI_list(i).centerPos(1)-x_val)^2+(ROI_list(i).centerPos(2)-y_val)^2);
        end
        [~,I]=sort(ROI_distances);
        ordered_ROI_list = ROI_list(I);
        ROI_list = ordered_ROI_list;
        
        centroid_display=zeros(ss(1),ss(2),3);
        alpha_cent=zeros(ss(1),ss(2),1);
        
        for i=1:ss(1)
            for j=1:ss(2)
                if sqrt((i-x_val)^2+(j-y_val)^2) <= dot_rad
                    centroid_display(i,j,1)=1;
                    centroid_display(i,j,2)=1;
                    centroid_display(i,j,3)=0;
                    alpha_cent(i,j)=1;
                end
            end
        end
        set(hdispCentroid,'CData',centroid_display,'AlphaData',alpha_cent)
        set(hdispCentroid,'Visible','on');
        
        Calculate_dFF_foopsi();
        drawROI();
        UpdateDisplay();
        set(htools.ROIOrder,'String','Reorder ROIs');
    end

    function EdCellRadCallback(varargin)
        cell_rad = str2double(get(htools.EdCellRad,'String'));
    end
    function EdNeuropilRadCallback(varargin)
        neuropil_rad = str2double(get(htools.EdNeuropilRad,'String'));
    end
    function EdNeuropilAlphaCallback(varargin)
        neuropil_alpha = str2num(get(htools.EdNeuropilAlpha,'String'));
    end
    function EdTau0Callback(varargin)
        tau_0 = str2num(get(htools.EdTau0,'String'));
    end
    function EdTau1Callback(varargin)
        tau_1 = str2num(get(htools.EdTau1,'String'));
    end
    function EdTau2Callback(varargin)
        tau_2 = str2num(get(htools.EdTau2,'String'));
    end
    function BnNeuropilCallback(varargin)
        set(htools.BnNeuropil,'String','One moment, please'); drawnow();
        ROI_list = neuropil_mask(im, ROI_list,cell_rad,neuropil_rad,neuropil_alpha);
        Calculate_dFF_foopsi();
        drawROI();
        UpdateDisplay();
        set(htools.BnNeuropil,'String','Calculate Neuropil');
    end
    function PopupFilterCallback(varargin)
        dFF_filter_ind = get(htools.filter_choice,'Value');
        dFF_filter_choices = get(htools.filter_choice,'String');
        dFF_filter = dFF_filter_choices{dFF_filter_ind};
    end
    function BnFilteringCallback(varargin)
        Calculate_dFF_foopsi();
        UpdateDisplay();
    end

    function SlROICallback(varargin)
        currROI=round(get(hobjs.SlROI,'Value'));
        set(hobjs.EdROI,'String',num2str(currROI));
        UpdateDisplay();
        drawROI();
    end
    function SlFrameCallback(varargin)
        currFrame=round(get(hobjs.SlFrame,'Value'));
        set(hobjs.EdFrame,'String',num2str(currFrame));
        UpdateDisplay();
    end
    function EdROICallback(varargin)
        currROI = str2double(get(hobjs.EdROI,'String'));
        set(hobjs.SlROI,'Value',currROI);
        UpdateDisplay();
        drawROI();
    end
    function EdFrameCallback(varargin)
        currFrame=str2double(get(hobjs.EdFrame,'String'));
        set(hobjs.SlFrame,'Value',currFrame);
        UpdateDisplay();
    end

    function BnExportData(varargin)
        % Send important data to the workspace
        assignin('base','ROI_list',ROI_list);
        assignin('base','tvec',t);
        assignin('base','tvec_for_dFF',new_t);
        assignin('base','dFF',dFF);
        %assignin('base','Nhat_dFF',Nhat_dFF);
        
        % And save all of this data in a directory
        [fname,pnameOut] = uiputfile('','Choose a path and folder name for the data.');
        rawF=zeros(num_ROIs,ss(3));
        neuropilF = rawF;
        F_subtracted = rawF;
        ROIx = zeros(num_ROIs,1);
        ROIy = ROIx;
        for i=1:num_ROIs
            rawF(i,:)=ROI_list(i).fmean;
            neuropilF(i,:)=ROI_list(i).neuropil_fmean;
            F_subtracted(i,:)=ROI_list(i).F_neuropilsubtracted;
            ROIx(i)=ROI_list(i).centerPos(1);
            ROIy(i)=ROI_list(i).centerPos(2);
        end
        rawF = double(rawF');
        neuropilF = double(neuropilF');
        F_subtracted = double(F_subtracted');
        ROIx = double(ROIx');
        ROIy = double(ROIy');
        tvec = double(t');
        tvec_for_dFF = double(new_t');
        dFF_copy = double(dFF);
        %Nhat_dFF_copy = double(Nhat_dFF);
        
        mkdir(pnameOut,fname);
        dir = cd(strcat(pnameOut,'/',fname));
        save('all_vars.mat','rawF','neuropilF','F_subtracted','ROIx','ROIy','tvec','tvec_for_dFF','dFF_copy');%,'Nhat_dFF_copy');
        save('rawF.txt','rawF','-ascii');
        save('neuropilF.txt','neuropilF','-ascii');
        save('F_neuropilsubtracted.txt','F_subtracted','-ascii');
        save('ROI_xpos.txt','ROIx','-ascii');
        save('ROI_ypos.txt','ROIy','-ascii');
        save('time_vector.txt','tvec','-ascii');
        save('time_vector_for_dFF.txt','tvec_for_dFF','-ascii');
        save('dFF.txt','dFF_copy','-ascii');
        %save('Nhat_dFF.txt','Nhat_dFF_copy','-ascii');
        cd(dir);
    end

    function BnDeleteROI(varargin)
        check = questdlg('Are you sure you want to delete the selected ROI?','Double Check','Cancel');
        switch check
            case 'Yes'
                if currROI == num_ROIs
                    new_ROI_inds = [1:currROI-1];
                    currROI = currROI-1;
                elseif currROI == 1
                    new_ROI_inds = [currROI+1:num_ROIs];
                else
                    new_ROI_inds = [1:currROI-1,currROI+1:num_ROIs];
                end
                new_ROI_list = ROI_list(new_ROI_inds);
                ROI_list = new_ROI_list;
                num_ROIs = length(ROI_list);
            case 'No'
                return
            case 'Cancel'
                return
        end
        set(hobjs.SlROI,'Max',num_ROIs,'Min',1,'Value',currROI,'SliderStep',[1/(num_ROIs-1),0.1]);
        Calculate_dFF_foopsi();
        drawROI();
        UpdateDisplay();
    end

    function dispfigDelCallback(varargin)
        hdispfig=-1;
    end
    function disprawDelCallback(varargin)
        hdispraw=-1;
    end
    function toolboxfigDelCallback(varargin)
        htoolboxfig=-1;
    end
    function ctrlfigDelCallback(varargin)
        if ishandle(hdispfig)
            close(hdispfig);
        end
        if ishandle(htoolboxfig)
            close(htoolboxfig);
        end
        if ishandle(hdispraw)
            close(hdispraw)
        end
    end

toc
end