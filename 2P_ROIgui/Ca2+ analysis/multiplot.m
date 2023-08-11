function multiplot(X,Y,type)

if nargin < 3
    error('Must input two matrices (Xdata and Ydata) and specify if they are row or column vectors.');
end
if size(X,1) ~= size(Y,1) || size(X,2) ~= size(Y,2)
    error('X and Y must be two-dimensional matrices with the same dimensions');
end
if length(size(X)) > 2 || length(size(Y)) > 2
    error('X and Y must be two-dimensional matrices with the same dimensions');
end

switch type
    case 'row'
        num_Data = size(X,1);
    case 'column'
        num_Data = size(X,2);
end

currData = 1;

hctrlfig=figure('name','Data Plotter','position',[900,100,500,550],'MenuBar','none','Toolbar','figure','color',[1,1,1]*0.85);
uicontrol(hctrlfig,'Tag','SlData','Style','slider','Units','pixel','position',[20,500,400,20],'Max',num_Data,'Min',1,'Value',currData,'SliderStep',[1/(num_Data-1),0.1],'Callback',@SlDataCallback);
uicontrol(hctrlfig,'Tag','EdData','Style','edit','Units','pixel','position',[440,500,50,20],'String',num2str(currData),'Callback',@EdDataCallback);
axes('Tag','ax1','units','pixel','position',[50,50,400,400],'Parent',hctrlfig,'NextPlot','add','Visible','on');
hobjs=guihandles(hctrlfig);

if  usejava('awt') % java enabled -> use it to update while dragging
    l1=handle.listener(hobjs.SlData,'ActionEvent',@SlDataCallback);
end

switch type
    case 'row'
        h_Data = plot(X(currData,:),Y(currData,:),'.');
    case 'column'
        h_Data = plot(X(:,currData),Y(:,currData),'.');
end

    function SlDataCallback(varargin)
        currData=round(get(hobjs.SlData,'Value'));
        set(hobjs.EdData,'String',num2str(currData));
        switch type
            case 'row'
                set(h_Data,'XData',X(currData,:),'YData',Y(currData,:));
            case 'column'
                set(h_Data,'XData',X(:,currData),'YData',Y(:,currData));
        end
    end
    function EdDataCallback(varargin)
        currData = str2double(get(hobjs.EdData,'String'));
        set(hobjs.SlData,'Value',currData);
        switch type
            case 'row'
                set(h_Data,'XData',X(currData,:),'YData',Y(currData,:));
            case 'column'
                set(h_Data,'XData',X(:,currData),'YData',Y(:,currData));
        end
    end

end
