function [ROI_list] = neuropil_mask(im, ROI_list,cell_rad,neuropil_rad,neuropil_alpha)

% This script will calculate the mean fluorescence of the neuropil for a
% set diameter around every cell identified using the Svoboda ROIgui. It
% will then subtract the neuropil fluorescence from the corresponding cell
% using a scaling factor defined by "neuropil_alpha." It requires that
% 'ROI_list' (output of ROIgui) and 'im' (imput of ROIgui) be present in
% the MATLAB workspace or the calling function.
% 
% Last updated: 30 July 2014, dbarson

im_dim_m = size(im,1);
im_dim_n = size(im,2);
im_dim_t = size(im,3);

num_pixels = im_dim_m * im_dim_n;
num_ROIs = size(ROI_list,2);

% Find all of the pixels that are between a defined "cell radius" and a
% defined "neuropil radius" from the center of the cell in the selected
% ROI. Exclude all pixels that also belong to any cell ROI in the ROI_list
% structure (i.e. exclude all pixels that define neighboring cells, as well as the selected cell).

cell_pixels = ones(im_dim_m,im_dim_n);
neuropil_pixels = zeros(im_dim_m,im_dim_n,num_ROIs);

for i = 1:num_ROIs
    ROIx = ceil(ROI_list(i).centerPos(1));
    ROIy = ceil(ROI_list(i).centerPos(2));
    cell_pixels(ROI_list(i).pixel_list)= 0;
    for j = 1:im_dim_m
        for k = 1:im_dim_n
            if (ROIx - j)^2 + (ROIy - k)^2 <= cell_rad^2
                cell_pixels(j,k) = 0;
            end
            if (ROIx - j)^2 + (ROIy - k)^2 <= neuropil_rad^2
                neuropil_pixels(j,k,i) = 1;
            end
        end
    end
end

cell_subtracted_neuropil = zeros(im_dim_m,im_dim_n,num_ROIs);
im_2D = reshape(im,num_pixels,im_dim_t);

for i = 1:num_ROIs
   neuropil_temp = neuropil_pixels(:,:,i).*cell_pixels;
   neuropil_pixel_inds = find(neuropil_temp == 1);
   neuropil_traces = im_2D(neuropil_pixel_inds,:);
   neuropil_fmean = mean(neuropil_traces,1);
   ROI_list(i).neuropil_fmean = neuropil_fmean';
   ROI_list(i).neuropil_pixel_list = neuropil_pixel_inds;
   cell_subtracted_neuropil(:,:,i) = neuropil_temp;
   ROI_list(i).F_neuropilsubtracted = ROI_list(i).fmean - neuropil_alpha*ROI_list(i).neuropil_fmean; %neuropil_alpha is an experimentally determined value for the relative contribution of the neuropil to cellular fluorescence
end

% neuropil_flattened = sum(cell_subtracted_neuropil,3);
% neuropil_flattened(neuropil_flattened >= 1) = 1;
% neuropil_flattened_inv = (neuropil_flattened - ones(im_dim_m,im_dim_n));

%% graphical demonstration

% av_im = mean(im,3);
% im_max = max(max(av_im));
% im_min = min(min(av_im));
% im_range = im_max - im_min;
% av_im_normalized = (av_im - im_min)./im_range;
% 
% im_plus_neuropil = av_im_normalized.*neuropil_flattened_inv;
% 
% figure
% h = imagesc(-im_plus_neuropil);
% colormap(gray)


end


