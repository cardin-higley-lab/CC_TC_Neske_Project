% Get and load Ca_soma_dend file
clc; clear all;
addpath 'D:\matlab code'
Dir = uigetdir; %gets directory
files = dir(fullfile(Dir,'Ca*.mat')); % get files
%Files = uigetfile('*.tif', 'Select Multiple Files', 'MultiSelect', 'on' );

load(fullfile(files.folder, files.name));
Ca_dend_soma = Ca_dend_soma_stim3;

soma_dff = Ca_dend_soma{end-1, 4};
plot(soma_dff(1:9000));

dend_avg = Ca_dend_soma{end, 4};
plot(dend_avg(1:9000));

dend_dff = Ca_dend_soma{1, 4};
plot(dend_dff(1:9000));

hold on; plot(dend_dff(1:9000)); plot(soma_dff(1:9000));

% normalize Ca_dff
soma_norm = (soma_dff- ca_min(soma_dff)) ./ (max(soma_dff) - ca_min(soma_dff));
dend_avg_norm = (dend_avg - ca_min(dend_avg)) ./ (max(dend_avg) - ca_min(dend_avg)); 
dend_norm = (dend_dff - ca_min(dend_dff)) ./ (max(dend_dff) - ca_min(dend_dff)); 


hold on; plot(dend_norm(1:9000)); plot(soma_norm(1:9000)); plot(dend_avg_norm(1:9000));

soma_tau_decay = Ca_dend_soma{end-1, 7}* 30;
dend_avg_tau_decay = Ca_dend_soma{end,7} * 30;
dend_tau_decay = Ca_dend_soma{1, 7} *30;


addpath 'D:\matlab code\OASIS_matlab-master'
% Deconvolve soma
[c_soma, s_soma] = deconvolveCa(soma_dff, 'exp2', [80, 40], 'foopsi', 'lambda',0.01, 'shift', 100, 'window', 800);
figure; hold on; plot(s_soma(1:9000)), plot(soma_dff(1:9000)); 

% Deconvolve dend_avg

[c_dend_avg, s_dend_avg] = deconvolveCa(dend_dff_avg, 'exp2', [65, 30], 'foopsi', 'lambda',0.01, 'shift', 100, 'window', 800);
figure; hold on; plot(s_dend_avg(1:9000)), plot(dend_avg(1:9000));


%% Kernel convultion with the soma signal
% 1) reference signal
y = reshape(s_soma, 1, []); T = length(y);
nMax = 200; t = 1:nMax;
%decay = round ((Ca_dend_soma{1, 7}) * 30)/2;
tau_decay = 30+2*(0:5-1);
tau_rise = 10+2*(0:5-1);

y = reshape(s_soma, 1, []);
T = length(y);
plot(y(1:9000));


kernels_soma = {};
n = 0;
for ii = 1:length(tau_decay)
    for jj = 1:length(tau_rise)
        
        kernel = (exp(-t/tau_decay(ii)) .*  (1-exp(-t/tau_rise(jj))));
        
        figure(1); plot(kernel); hold on;
        gt = kernel; ind0 = 1;
        gt = [reshape(gt, 1, nMax), zeros(1, T)];
        gt_max = max(gt);
        gt = gt/gt_max;
        gt1 = [gt(ind0:end), zeros(1, ind0-1)];
        
        y = [y, zeros(1, nMax)];
        conv1_signal = conv(y, gt(1:nMax));
        conv_signal = conv1_signal(1:T) + y(1)*gt1(1:T);
        conv_soma = conv_signal';
        conv_soma_norm = normalize(conv_soma);
        
        b=robustfit(conv_soma_norm, dend_norm);
        conv_signal_sc = conv_soma_norm * b(2);
        dend_reg = (dend_norm) - (conv_soma_norm * b(2)) - b(1);
       
        mean_reg = mean(dend_reg);
        
        n = n +1;
        kernels_soma(n,:) = {sprintf('kernel %d', n) conv_soma_norm  b conv_signal_sc dend_reg};
     end
end

a = 1+5*(0:5-1);
mm = [a; a+1; a+2; a+3; a+4];
for j = 1: length(mm(:,1))
    figure(j);
    t = tiledlayout(5,1);
    for i = mm(j,1):5: mm(j,end)
        nexttile; plot(kernels_soma{i, 5}(1:9000)); hold on; plot(dend_norm(1:9000));% plot(kernels_soma{i, 5}(1:9000))
        
    end
end


%% Kernel convultion with the dend_avg signal
% 1) reference signal
y = reshape(s_dend_avg, 1, []); T = length(y);

nMax = 800; t = 1:nMax;
tau_decay = 40+10*(0:5-1);
tau_rise = 20+5*(0:5-1);
kernels_dend_avg = {};
n = 0;
for ii = 1:length(tau_decay)
    for jj = 1:length(tau_rise)
        
        kernel = (exp(-t/tau_decay(ii)) .*  (1-exp(-t/tau_rise(jj))));
        figure(1); plot(kernel); hold on;
        gt = kernel; ind0 = 1;
        gt = [reshape(gt, 1, nMax), zeros(1, T)];
        gt_max = max(gt);
        gt = gt/gt_max;
        gt1 = [gt(ind0:end), zeros(1, ind0-1)];
       
        y = [y, zeros(1, nMax)];
        conv1_signal = conv(gt(1:nMax), y);
        conv_signal = conv1_signal(1:T) + y(1)*gt1(1:T);
        conv_signal = conv_signal';
        
        conv_signal_norm = (conv_signal- min(conv_signal)) ./ (max(conv_signal) - min(conv_signal));
        
        b=robustfit(conv_signal_norm, dend_norm);
        conv_signal_sc = conv_signal_norm * b(2);
        dend_reg = (dend_norm) - (conv_signal_sc) - b(1);
        mean_reg  = mean(dend_reg);
        
        n = n +1;
        kernels_dend_avg(n,:) = {sprintf('kernel %d', n) kernel' conv_signal conv_signal_norm b conv_signal_sc dend_reg, mean_reg};
        
    end
end

a = 1+5*(0:5-1);
mm = [a; a+1; a+2; a+3; a+4];
for j = 1: length(mm(:,1))
    figure(j);
    t = tiledlayout(5,1);
    for i = mm(j,1):5: mm(j,end)
        nexttile; plot(kernels_dend_avg{i, 4}(1:9000)); hold on; plot(dend_norm(1:9000)); plot(kernels_dend_avg{i, 7}(1:9000))
        
    end
end

%% best fit kernels soma
w = 1;
figure; plot(kernels_soma{w, 4}(1:9000)); hold on; plot(kernels_soma{w, 5}(1:9000))
figure; plot(kernels_soma{w, 4}(1:9000)); hold on; plot(dend_norm(w:9000)); plot(kernels_soma{w, 7}(1:9000))

% best fit kernels dend_avg
ww=11;
figure; plot(kernels_dend_avg{ww, 4}(1:9000)); hold on; plot(kernels_dend_avg{ww, 6}(1:9000))
figure; plot(kernels_dend_avg{ww, 4}(1:9000)); hold on; plot(dend_norm(ww:9000)); plot(kernels_dend_avg{ww, 7}(1:9000))

%% Kernel re convultion with the average dendrtic signal

b=robustfit(dend_avg_norm, dend_norm);
dend_sc = dend_avg_norm * b(2);
dend_regressed = (dend_norm) - (dend_sc);
% avg dendritic signal is removed from dend seg signal by subtracting a scaled version of the dend avg where the scaling factor 
%equals the slope of robust regression 

figure; plot(dend_regressed(1:9000)); hold on; plot(dend_norm(1:9000)); plot(dend_sc(1:9000));


scatter(dend_avg_norm(1:5000),dend_norm(1:5000));
hold on
plot(dend_avg_norm,b(1)+ b(2)*dend_avg_norm,'r-')




%% 


% plotting closed values

a = 1:1:10;
k = [a; a+10; a+20; a+30; a+40; a+50];
for j = 1: length(k(:,1))
    figure(j);
    t = tiledlayout(10,1);
    for i = k(j,1) : 1:k(j,end)
        nexttile; plot(kernels{i, 5}(1:9000)); hold on; plot(dend_norm(1:9000))
        %set(gca,'XTick',[], 'YTick', [])
        %set(gca,'visible','off')
        %title(sprintf('seg %d', k));
    end
end

% plotting distanced values
a = 1+20*(0:5-1);
k = [a; a+6; a+10; a+15; a+19];


for j = 1: length(k(:,1))
    figure(j);
    t = tiledlayout(5,1);
    for i = k(j,1) : 20 :k(j,end)
        nexttile; plot(kernels{i, 5}(1:9000)); hold on; plot(dend_norm(1:9000))
        %set(gca,'XTick',[], 'YTick', [])
        %set(gca,'visible','off')
        %title(sprintf('seg %d', k));
    end
end

%% Kernel convultion with the reference signal
% 1) reference signal

y = reshape(s_soma, 1, []);
T = length(y);
plot(y(1:9000));

% 2) find exponential kernel for each ROIS
nMax = 200; t = 1:nMax;
tau_decay = 45; tau_rise = 25;
kernel = (exp(-t/tau_decay) .*  (1-exp(-t/tau_rise)));
figure; plot(kernel); %hold on; plot(event)

gt = kernel;
ind0 = 1;
gt = [reshape(gt, 1, nMax), zeros(1, T)]; plot(gt)

gt_max = max(gt);
gt = gt/gt_max;
gt1 = [gt(ind0:end), zeros(1, ind0-1)];
plot(gt1); hold on; 

% 3) covolve kernel fit with reference signal
y = [y, zeros(1, nMax)];
conv_signal = conv(y, gt(1:nMax));
conv_soma = conv_signal(1:T) + y(1)*gt1(1:T);
figure;plot(conv_soma(1:9000));
%conv_soma_norm = (conv_soma- mean(conv_soma)) ./ (max(conv_soma) - min(conv_soma));
conv_soma_norm = normalize(conv_soma)

plot(conv_soma_norm(1:9000));
%dend_norm = (dend_dff- mean(dend_dff)) ./ (max(dend_dff) - min(dend_dff));
dend_norm = normalize(dend_dff)

soma_norm = normalize(soma_dff)
figure, plot(conv_soma_norm(1:9000)); hold on; plot(dend_norm(1:9000)), plot(soma_norm(1:9000))


b=robustfit(dend_dff_avg, b);
dend_regressed = (spine) - (dend_dff_avg * b(2)) - b(1);
dend_sc = conv_soma_norm* b(2);
figure(6); plot(dend_regressed(1:9000)); hold on; plot(spine(1:9000));% plot(dend_sc(1:9000));
 box off;
 %set(gca,'XTick',[], 'YTick', [])
  set(gca,'visible','off')
   %title(sprintf('seg %d', k));
  %set(gcf,'Position',[400 400 800 400])
  xlim([700 5760])



conv1_soma_norm = (conv1_signal- min(conv1_signal)) ./ (max(conv1_signal) - min(conv1_signal));

figure, plot(conv1_soma_norm(1:9000)); hold on; plot(dend_norm(1:9000))











%c_norm = normalize(c,'norm',2);
conv1_soma_norm = (conv1_signal- min(conv1_signal)) ./ (max(conv1_signal) - min(conv1_signal));
c1_soma_norm = (c1_soma- min(c1_soma)) ./ (max(c1_soma) - min(c1_soma));
c_soma_norm = (c_soma- min(c1_soma)) ./ (max(c_soma) - min(c_soma));
plot(conv_signal_norm(1:9000));

figure, hold on; plot(dend1_norm(1:9000)); plot(conv1_soma_norm(1:9000)); plot(c1_soma_norm(1:9000))
figure, plot(c1_soma_norm(1:9000)); hold on; plot(dend1_norm(1:9000));


dend1_norm = (dend_dff - ca_min(dend_dff)) ./ (max(dend_dff) - ca_min(dend_dff)); 
soma_norm = (soma_dff- min(soma_dff)) ./ (max(soma_dff) - min(soma_dff));
soma_conv_dend1 = conv_signal_norm;
dend1_F = dend1_norm - conv_signal_norm; 

figure(13)
plot(soma_norm); hold on; plot(conv_signal_norm); plot(dend1_norm)

figure(14)
plot(s_soma); hold on; plot(soma_dff)


%% soma deconvulution 
[c_soma, s_soma, options] = deconvolveCa(soma_dff, 'exp2', 'foopsi', 'lambda', .1, 'shift', 100, 'window', 200, 'smin',1);
[cc, ss] = deconvolveCa(soma_dff, 'exp2', [60, 5], 'thresholded');
[c4, s4] = deconvolveCa(soma_dff, 'exp2', [75, 15], 'foopsi', 'lambda',0.1, 'shift', 100, 'window', 200);

[c6_soma, s6_soma] = deconvolveCa(soma_dff, 'exp2', [50, 8], 'foopsi', 'lambda',0.1, 'shift', 100, 'window', 500);

figure; plot(c1_soma); hold on; plot(c2_soma); plot(c3_soma); plot(c4_soma); plot(c5_soma);  plot(c6_soma)

c4_soma_norm = (c4_soma- min(c4_soma)) ./ (max(c4_soma) - min(c4_soma)); 
c1_soma_norm = (c1_soma- min(c1_soma)) ./ (max(c1_soma) - min(c1_soma));
c2_soma_norm = (c2_soma- min(c2_soma)) ./ (max(c2_soma) - min(c2_soma));
c3_soma_norm = (c3_soma- min(c3_soma)) ./ (max(c3_soma) - min(c3_soma));
figure; plot(c1_soma_norm); hold on; plot(c2_soma_norm); plot(c3_soma_norm); plot(c4_soma_norm)
figure; plot(s4_soma); hold on; 

figure; plot(soma_norm); hold on; plot(c_soma); plot(dend1_norm)
plot(s4); hold on; plot(soma_dff)


% define the parameters options.sn and options.pars for constrained foopsi deconvolution
[c, s,options] = deconvolveCa(soma_dff, 'ar1', 'constrained');

options.sn = GetSn(soma_dff, [.25,.5], 'logmexp');
options.pars = estimate_time_constant(soma_dff,1,options.sn, 5, 0.5);

options.pars = 80 * 30;

[c_soma, s_soma] = deconvolveCa(soma_dff, 'ar1', 'constrained', options.pars);

%[c, s, options] = deconvolveCa(dend_dff, 'exp2', 'foopsi', 'lambda', .4, 'shift', 100, 'window', 200);
[c_roi, s_roi, options] = deconvolveCa(roi_dff, 'exp2', 'foopsi', 'lambda', .4, 'shift', 100, 'window', 200);