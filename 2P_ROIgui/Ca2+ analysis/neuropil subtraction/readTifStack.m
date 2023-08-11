% simple function to read multiframe tif file into data matrix (NxMxT)
% 04.14.09 Tsai-Wen Chen 
% 09.20.10 Improve reading over network by first creating local copy
function data=readTifStack(varargin)

%  readTifStack(filename)
%  readTifStack(filename,index)
%  readTifStack(filename,firstim,lastim)
movelocal=0;
index=[];
if nargin ==0
  [filename, tif_path] = uigetfile('*.tif','choose a file');
  if isequal(filename,0);return;end
  filename = [tif_path filename];  
else
  filename=varargin{1};
end

if nargin == 2 
    index=varargin{2};
end

if nargin ==3
    index=(varargin{2}:varargin{3});
end

if nargin ==4
    index=(varargin{2}:varargin{3});
    movelocal=varargin{4};
end
%%
% local=['C','D','E','F','G'];
if movelocal
    [pathstr, name]=fileparts(filename);
    if isempty(pathstr)
        filename=[pwd,'\',filename];
    end
    dos(['copy "',filename, '"  c:\temp.tif']);
    filename='c:\temp.tif';
    disp('create local copy');
end

info=imfinfo(filename);
nimage=length(info);
if isempty(index)
    index=1:nimage;
end
%%


nread=length(index);
%data=zeros(info(1).Height,info(2).Width,nread);
data=zeros(info(1).Height,info(1).Width,nread,'single');
for i=1:length(index)
    data(:,:,i)=imread(filename,index(i),'Info',info);    
end
