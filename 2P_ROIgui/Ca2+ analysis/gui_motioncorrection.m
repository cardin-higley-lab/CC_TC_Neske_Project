function gui_motioncorrection

global configs
configs.max_shift = 12;
configs.threshold = 0.95;
configs.moving_average = 20;
configs.lims = [40 240];

f = figure('Visible','on', 'Position', [100 100 600 400]);

% Create pop-up menu
txt1 = uicontrol('Style', 'text',...
     'Position', [20 320 80 20],...
     'String', 'max_shift');    

% Create pop-up menu
txt1_A = uicontrol('Style', 'edit',...
     'String', '12',...
     'Position', [20 340 80 50],...
     'Callback', @max_shift);    

       % Create pop-up menu
txt2 = uicontrol('Style', 'text',...
     'Position', [110 320 100 20],...
     'String', 'threshold');    

% Create pop-up menu
txt2_A = uicontrol('Style', 'edit',...
     'String', '0.95',...
     'Position', [110 340 100 50],...
     'Callback', @threshold);  


% Create pop-up menu
txt3 = uicontrol('Style', 'text',...
     'Position', [220 320 100 20],...
     'String', 'mov_av');    

% Create pop-up menu
txt3_A = uicontrol('Style', 'edit',...
     'String', '20',...
     'Position', [220 340 100 50],...
     'Callback', @moving_average);  


% Create pop-up menu
txt4 = uicontrol('Style', 'text',...
     'Position', [330 320 100 20],...
     'String', 'lims');    

% Create pop-up menu
txt4_A = uicontrol('Style', 'edit',...
     'String', '[40 240]',...
     'Position', [330 340 100 50],...
     'Callback', @lims);  

% Create pop-up menu
txt5 = uicontrol('Style', 'pushbutton',...
     'Position', [200 220 80 20],...
     'String', 'Select file', 'Callback', @select_file);    

% Create pop-up menu
txt6_A = uicontrol('Style', 'text',...
     'String', 'GUI for motion correction by Martin Vinck (c) 2014',...
     'Position', [100 100 200 50]);
   
%%%%%%%%% CALLBACKS GUI
function max_shift(source,callbackdata)
global configs
val = str2num(get(source, 'String'));
configs.max_shift = val;

function threshold(source,callbackdata)
global configs
val = str2num(get(source, 'String'));
configs.threshold = val;
       
function moving_average(source,callbackdata)
global configs
val = str2num(get(source, 'String'));
configs.moving_average = val;
       
function lims(source,callbackdata)
global configs
val = str2num(get(source, 'String'));
configs.lims = val;       
       
function select_file(source,callbackdata)
global configs
[referenceMovie, pathname] = uigetfile('All Files');
cfg.filename = 1;
datasets = getdatasets(pathname, {'.tif'}, {''}, cfg);
fprintf('loading in the reference image\n')
dataRef  = readTifStack(fullfile(pathname, referenceMovie));
nMovies  = length(datasets); 
for iMovie = 1:nMovies
  if strcmp(datasets{iMovie}, fullfile(pathname, referenceMovie))
    data = dataRef;
  else
    fprintf('loading in %s\n', datasets{iMovie})    
    data = readTifStack(datasets{iMovie});
  end    
  
  % compute the old average xcorr
  configsnew = configs;
  configsnew.moving_average = 1;
  configsnew.threshold = 0;
  
  [ignore, paramsOrig] = motioncorrection(nanmean(dataRef,3), nanmean(data,3), configsnew);  
 
  figure, 
  subplot(2,1,1)
  imagesc(paramsOrig.cross_corr), colorbar
  title(sprintf('original xcorr for %s', datasets{iMovie}))
  
  % compute the motion correction
  fprintf('running mcorr on %s using \n max_shift = %d, threshold = %2.2f, moving av = %d\n', datasets{iMovie}, configs.max_shift, configs.threshold, configs.moving_average);
  [data, params] = motioncorrection(dataRef, data, configs); % overwrite for memory
  
  % compute the new average xcorr
  [ignore, paramsNew] = motioncorrection(nanmean(dataRef,3), nanmean(data,3), configsnew);  
  
  subplot(2,1,2)
  imagesc(paramsNew.cross_corr), colorbar    
  title(sprintf('new xcorr for %s', datasets{iMovie}))
  
  outputFilename = datasets{iMovie}(1:end-4);
  save(outputFilename, 'data', 'params', '-v7.3');
end
    
%%%% just a subfunction to get the tif files
function [datasets] = getdatasets(inputDir,positiveidentifiers, negativeidentifiers,cfg)
if nargin<4
    cfg.seltype = 'all';
    cfg.filename = 1;
end
if ~isfield(cfg,'seltype'), cfg.seltype = 'all'; end

datasets = {};
dirinfo = dir(inputDir);
for iFile=1:length(dirinfo)
        
    pos = [];
    for k = 1:length(positiveidentifiers)        
        pos(k) = ~isempty(strfind(dirinfo(iFile).name,positiveidentifiers{k}));
    end
    
    if strcmp(cfg.seltype,'all')        
        isFilePos = all(pos);
    elseif strcmp(cfg.seltype,'any')
        isFilePos = any(pos);
    end
    neg = [];
    for k = 1:length(negativeidentifiers)        
        neg(k) = ~isempty(strfind(dirinfo(iFile).name,negativeidentifiers{k}));
    end
    
    if strcmp(cfg.seltype,'all')        
        isFileNeg = all(neg);
    elseif strcmp(cfg.seltype,'any')
        isFileNeg = any(neg);
    end
    
    isValid = isFilePos&&~isFileNeg;
    
    
    if isValid
    	if cfg.filename
	        datasets{end+1} =  fullfile(inputDir,dirinfo(iFile).name);
	    else
            datasets{end+1} =  dirinfo(iFile).name;
        end
    end  
end


% simple function to read multiframe tif file into data matrix (NxMxT)
% 04.14.09 Tsai-Wen Chen 
% 09.20.10 Improve reading over network by first creating local copy
% from mike 
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
%


nread=length(index);
%data=zeros(info(1).Height,info(2).Width,nread);
data=zeros(info(1).Height,info(1).Width,nread,'single');
for i=1:length(index)
    data(:,:,i)=imread(filename,index(i),'Info',info);    
end
 
function [test_images,params] = motioncorrection(ref_images, test_images, cfg)

% MOTIONCORRECTION performs motion correction for a stack of test images based on a stack
% of reference images
% 
% Inputs:
%   ref_images: reference images
%   should be nPixels x mPixels x # refImages
%   test_images: test images, on which we want to do motion correction
%   test_images must be nPixels x mPixels x # testImages
%
%   configurations (cfg):
%     cfg.max_shift    : maximum shift that we want to compute (default = 7)
%
%     cfg.thresholdtype  : 'fraction' or 'abs': 'fraction' uses fraction of maximum 
%     correlation, and absolute simply puts everything below that value to zero.
%     ('fraction' recommended)
%     cfg.threshold below which we set cross-corr values to 0 (default = 0.8)
%
%     cfg.moving_average: number of frames to average over in order to compute
%     shift (default = 1, no moving average); if number is even, we add one to make it odd.
%
% Outputs: 
%
%   test_images: test images, shifted
%   params:      parameters of the cross-correlation for each frame
%     params.x_center_of_mass = center of mass of cross-correlation, x position
%     params.y_center_of_mass = center of mass of cross-correlation, y position
%     params.cross_corr       = cross-correlations between reference and test images
%     params.x_peak           = x position of peak in cross-correlation
%     params.y_peak           = y position of peak in cross-correlation
%
% See MOTIONCORRECTION_TEST for demo

% Copyright (c) Martin Vinck 2014

if ~isfield(cfg, 'max_shift'),      cfg.max_shift  = 12;                 end
if ~isfield(cfg, 'moving_average'), cfg.moving_average = 20;             end
if ~isfield(cfg, 'thresholdtype'),  cfg.thresholdtype = 'fraction';      end
if ~isfield(cfg, 'threshold'),      cfg.threshold = 0.95;                end
if ~isfield(cfg, 'leftlim'),        cfg.leftlim   = 1;                   end
if ~isfield(cfg, 'rightlim'),       cfg.rightlim  = size(ref_images, 2); end

max_shift = cfg.max_shift;
threshold = cfg.threshold;
moving_average = cfg.moving_average;

% ensure that we have an uneven number of frames, just for the convolution
if mod(moving_average,2)==0, moving_average = moving_average + 1;       end

% compute the reference image by averaging
ref_image = nanmean(ref_images,3);

% for the cross-correlation, we ignore the edges
ref_image = ref_image(:,cfg.leftlim:cfg.rightlim,:);
ref_image  = ref_image(max_shift+1:end-max_shift, max_shift+1:end-max_shift);

% preallocate cross-correlation before the loop
n_images = size(test_images,3);
params.cross_corr = zeros(max_shift*2+1, max_shift*2+1,n_images);

% compute the moving average of the images 
test_images_mva = convn(test_images, ones(1,1,moving_average),'same');

% first compute the optimal shift for every so many images
inds = (moving_average+1)/2:moving_average:n_images-(moving_average+1)/2+1;
for i_image = inds
  
  fprintf('Doing image %d out of %d images\n', i_image, n_images);
  % call the main function that computes the cross-correlation and center of mass
  [params_xcorr] = correctimage(max_shift, ref_image, test_images_mva(:,cfg.leftlim:cfg.rightlim,i_image), threshold, cfg);
  
  % store in structure
  params.x_peak(i_image) = params_xcorr.x_peak;
  params.y_peak(i_image) = params_xcorr.y_peak;
  params.x_center_of_mass(i_image) = params_xcorr.x_center_of_mass;
  params.y_center_of_mass(i_image) = params_xcorr.y_center_of_mass;
  params.cross_corr(:,:,i_image)    = params_xcorr.cross_corr;  
end

% now allign frame by frame
for i_image = 1:n_images
  
  % get the index for which we computed the moving average
  indx = inds(nearest(inds,i_image));
  
  % get the parameters for this shift from that moving average; just to avoid confusion: 
  % by convention in a matrix x would go over the columns (dim=2), y over throws (dim=1)
  y_center_of_mass = params.y_center_of_mass(indx);
  x_center_of_mass = params.x_center_of_mass(indx);
  params.x_center_of_mass(i_image) = x_center_of_mass;
  params.y_center_of_mass(i_image) = y_center_of_mass;
  params.x_peak(i_image) = params.x_peak(indx);
  params.y_peak(i_image) = params.y_peak(indx);
  params.cross_corr(:,:,i_image) = params.cross_corr(:,:,indx);
      
  % get the shifted image: no need to cut out the edges completely
  % append with NaNs to the side creating a dummy image that we fill in then
  test_image = test_images(:,:,i_image);
  dum_image  = NaN(size(test_image,1)+max_shift*2, size(test_image,2)+max_shift*2);
  dum_image(max_shift+1:end-max_shift, max_shift+1:end-max_shift) = test_image;
  
  test_image_shifted = dum_image(max_shift + 1 + round(y_center_of_mass): end-max_shift + round(y_center_of_mass),...
  max_shift + 1 + round(x_center_of_mass): end-max_shift + round(x_center_of_mass)); 

  % store it
  test_images(:,:,i_image) = test_image_shifted;  
end

function [params_xcorr] = correctimage(max_shift, ref_image, test_image, threshold, cfg)

% first step: compute the cross-correlations
pearson_r = zeros(max_shift*2+1, max_shift*2+1); % preallocate

x = ref_image(:); % vectorize ref image for correlation

% loop through the various possible shifts
for i_shift = -max_shift : max_shift
  for j_shift = -max_shift : max_shift
    
    % cut out part of the test image
    test_image_new = test_image(max_shift+1+i_shift : end-max_shift+i_shift, ...
      max_shift+j_shift+1 : end-max_shift+j_shift);
    
    y = test_image_new(:);
    
    % ignore NaNs
    sl = ~isnan(x) & ~isnan(y);
        
    % compute correlation
    R = corrcoef(x(sl),y(sl));
    
    % store
    pearson_r(i_shift+max_shift+1, j_shift+max_shift+1) = R(1,2);
  end
end

% second step: threshold depending on input configuration
mx = nanmax(pearson_r(:));
if strcmp(cfg.thresholdtype, 'fraction')
  pearson_r(pearson_r<(mx.*threshold)) = 0;
else
  pearson_r(pearson_r<(threshold)) = 0;
end  
  
% third step: compute the center of mass; first make a grid of x and y
[x,y] = meshgrid(-max_shift:max_shift, -max_shift:max_shift); % deviations from center

% create the centers of mass
x_center_of_mass = sum(x(:).*pearson_r(:))./sum(pearson_r(:));
y_center_of_mass = sum(y(:).*pearson_r(:))./sum(pearson_r(:));

% find the maximum as well
[ignore,ind] = max(pearson_r(:));
[ind_x,ind_y] = ind2sub(size(pearson_r), ind);

% gather the output
params_xcorr.x_center_of_mass = x_center_of_mass;
params_xcorr.y_center_of_mass = y_center_of_mass;
params_xcorr.cross_corr       = pearson_r;
params_xcorr.x_peak           = ind_x;
params_xcorr.y_peak           = ind_y;


%%%%%%%%%% FUNCTION FROM FIELDTRIP TOOLBOX TO COMPUTE NEAREST
function [indx] = nearest(array, val, insideflag, toleranceflag)

% NEAREST return the index of an array nearest to a scalar
%
% Use as
%   [indx] = nearest(array, val, insideflag, toleranceflag)
%
% The second input val can be a scalar, or a [minval maxval] vector for
% limits selection.
%
% If not specified or if left empty, the insideflag and the toleranceflag
% will default to false.
%
% The boolean insideflag can be used to specify whether the value should be
% within the array or not. For example nearest(1:10, -inf) will return 1,
% but nearest(1:10, -inf, true) will return an error because -inf is not
% within the array.
%
% The boolean toleranceflag is used when insideflag is true. It can be used
% to specify whether some tolerance should be allowed for values that are
% just outside the array. For example nearest(1:10, 0.99, true, false) will
% return an error, but nearest(1:10, 0.99, true, true) will return 1. The
% tolerance that is allowed is half the distance between the subsequent
% values in the array.

% Copyright (C) 2002-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see https://urldefense.proofpoint.com/v1/url?u=http://www.ru.nl/neuroimaging/fieldtrip&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=3a40f00280447d942e3f3d180740f485e48ba30a0b1ff344cdb62bbaec2ac99a
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <https://urldefense.proofpoint.com/v1/url?u=http://www.gnu.org/licenses/&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=722a4370dab0cffbdad8b3459629067498773178b3d770ff586609adb5065b01>.
%
% $Id: nearest.m 8534 2013-09-25 08:15:54Z jorhor $

mbreal(array);
mbreal(val);

mbvector(array);
assert(all(~isnan(val)), 'incorrect value (NaN)');

if numel(val)==2
  % interpret this as a range specification like [minval maxval]
  % see also https://urldefense.proofpoint.com/v1/url?u=http://bugzilla.fcdonders.nl/show_bug.cgi?id%3D1431&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=01493f42fd65729f68700058b1015ce76503a9110fd57b4f34dae480bb4adb54
  intervaltol=eps;
  sel = find(array>=val(1) & array<=val(2));
  if isempty(sel)
    error('The limits you selected are outside the range available in the data');
  end
  indx(1) = sel(1);
  indx(2) = sel(end);  
  if indx(1)>1 && abs(array(indx(1)-1)-val(1))<=intervaltol
    indx(1)=indx(1)-1;
  end
  if indx(2)<length(array) && abs(array(indx(2)+1)-val(2))<=intervaltol
    indx(2)=indx(2)+1;
  end
  return
end

mbscalar(val);

if nargin<3 || isempty(insideflag)
  insideflag = false;
end

if nargin<4 || isempty(toleranceflag)
  toleranceflag = false;
end

% ensure that it is a column vector
array = array(:);

% determine the most extreme values in the array
minarray = min(array);
maxarray = max(array);

% do some strict checks whether the value lies within the min-max range
if insideflag
  if ~toleranceflag
    if val<minarray || val>maxarray
      error('the value %g should be within the range of the array from %g to %g', val, minarray, maxarray);
    end
  else
    if ~isequal(array, sort(array))
      error('the input array should be sorted from small to large');
    end
    if numel(array)<2
      error('the input array must have multiple elements to compute the tolerance');
    end
    mintolerance = (array(2)-array(1))/2;
    maxtolerance = (array(end)-array(end-1))/2;
    if val<(minarray-mintolerance) || val>(maxarray+maxtolerance)
      error('the value %g should be within the range of the array from %g to %g with a tolerance of %g and %g on both sides', val, minarray, maxarray, mintolerance, maxtolerance);
    end
  end % toleragceflag
end % insideflag

% FIXME it would be possible to do some soft checks and potentially give a
% warning in case the user did not explicitly specify the inside and
% tolerance flags

% note that [dum, indx] = min([1 1 2]) will return indx=1
% and that  [dum, indx] = max([1 2 2]) will return indx=2
% whereas it is desired to have consistently the match that is most towards the side of the array

if val>maxarray
  % return the last occurence of the largest number
  [dum, indx] = max(flipud(array));
  indx = numel(array) + 1 - indx;
  
elseif val<minarray
  % return the first occurence of the smallest number
  [dum, indx] = min(array);
  
else
  % implements a threshold to correct for errors due to numerical precision
  % see https://urldefense.proofpoint.com/v1/url?u=http://bugzilla.fcdonders.nl/show_bug.cgi?id%3D498&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=27937ff1a911c40a460032a68659053e27ebc6e7232b2d1aadb0b00e0d96ff4c and https://urldefense.proofpoint.com/v1/url?u=http://bugzilla.fcdonders.nl/show_bug.cgi?id%3D1943&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=670ddffb965a6ebb0626519d24ccabd1087d72b130ed2d3ba44c4301720d9434
%   if maxarray==minarray
%     precision = 1;
%   else
%     precision = (maxarray-minarray) / 10^6;
%   end
%   
%   % return the first occurence of the nearest number
%   [dum, indx] = min(round((abs(array(:) - val)./precision)).*precision);

  % use find instead, see https://urldefense.proofpoint.com/v1/url?u=http://bugzilla.fcdonders.nl/show_bug.cgi?id%3D1943&k=dpQisR3avULHgiNaNeY%2Btg%3D%3D%0A&r=Uk5fkA4LxM%2BT8j9lf4736g%3D%3D%0A&m=29oWU0yMceyvecddIUgUQ9x%2BXQPdL0oODKdA%2F5Dzocs%3D%0A&s=670ddffb965a6ebb0626519d24ccabd1087d72b130ed2d3ba44c4301720d9434
  wassorted = true;
  if ~issorted(array)
    wassorted = false;
    [array, xidx] = sort(array);
  end
  
  indx2 = find(array<=val, 1, 'last');
  indx3 = find(array>=val, 1, 'first');
  if abs(array(indx2)-val) <= abs(array(indx3)-val)
    indx = indx2;
  else
    indx = indx3;
  end
  if ~wassorted
    indx = xidx(indx);
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbreal(a)
if ~isreal(a)
  error('Argument to mbreal must be real');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbscalar(a)
if ~all(size(a)==1)
  error('Argument to mbscalar must be scalar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbvector(a)
if ndims(a) > 2 || (size(a, 1) > 1 && size(a, 2) > 1)
  error('Argument to mbvector must be a vector');
end

%%%%
function [zX] = ztransform(X,dim)

sizeX = size(X);
isNotDim = setdiff(1:length(sizeX),dim);
sizeX(isNotDim) = 1;
stdX = nanstd(X,[],dim);
meanX = nanmean(X,dim);
dX = X-repmat(meanX,sizeX);
zX = dX./repmat(stdX,sizeX);









    
    
    




       