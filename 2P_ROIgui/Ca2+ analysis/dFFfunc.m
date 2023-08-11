function [new_t, dFF] = dFFfunc(t,F,tau1,tau2,tau0,filter_type)
% This functions output deltaF/F0 for a calcium imaging ROI. 
% Procedure: 
%   1. Calculate a smoothed F(t). Width of smoothing is defined by tau_1. 
%   2. Calculate the time-dependent baseline F0(t) by taking the minimum value
%   of smoothed F(t) in a time window before t, the width of which is defined
%   by tau_2.
%   3. Calculate the relative change in fluorescence signal R(t) from F(t) and F0(t). R = (F-F0)/F0 
%   4. Apply exponentially weighted moving average filter for R(t) in order
%   to get dF/F(t).
% 
% Inputs:
%   t = vector of time points, in seconds
%   F = fluorescence data from ROI
%   tau1; Width of smoothing window, in seconds, to calculate smoothed F(t)
%   tau2; Size of window, in seconds, before t to select minimum value of smoothed F(t) for F0.
%   tau0; Decay time, in seconds, of expontially weighted moving average filter
%   filter_type = 'no_filter','exponential', or 'gaussian'
% 
% For beginning implementation, try tau1 = 0.75, tau2 = 3, and tau0 = 0.2.
%    
% Code is based off method described in Jia H, et al. "In vivo two-photon
% imaging of sensory evoked dendritic calcium signals in cortical neurons."
% Nature Protocols. 6, 2011, 28 - 35.
% 
% Last updated: 29 July 2014, dbarson
% 
%% Initialization of variables
dt = mean(diff(t));
tau1_imp = round(tau1/dt);
if rem(tau1_imp,2)
    tau1_imp = tau1_imp+1;
end

tau2_imp = round(tau2/dt); 


%% Smoothing of F(t) to calculate baseline F0, and calculation of R(t) = (F(t)-F0(t))/F0(t)

Fx = zeros(length(F),1);
for i=1:length(F)
    
        if i > tau1_imp/2 && i < length(F)-tau1_imp/2
            Fx(i) = (1/tau1_imp)*sum(F(i-tau1_imp/2+1:i+tau1_imp/2));
        elseif i <= tau1_imp/2
            Fx(i) = (2/tau1_imp)*sum(F(i:i+tau1_imp/2-1));
        elseif i >= length(F)-tau1_imp/2
            Fx(i) = (2/tau1_imp)*sum(F(i-tau1_imp/2+1:i));
        end
end

F0 = zeros(length(F)-tau2_imp,1);
R = F0;
for i=1:length(Fx)-tau2_imp;
    F0(i) = min(Fx(i:i+tau2_imp-1));
    R(i) = (F(i+tau2_imp)-F0(i))/F0(i);
end

% figure(6);
% subplot(2,1,1)
% hold on;
% plot(Fx);hold off;
% subplot(2,1,2)
% hold on;
% plot(F0);hold off;

%% Numeric implementation of expontially weighted moving-average filter.
new_t = t(tau2_imp:end-1);

switch filter_type 
    case 'exponential'
        load('alphatau.mat');
        alpha = interp1(num_tau,num_alpha,tau0/dt); 
        dFF = filter(alpha, [1 alpha-1], R);
    case 'gaussian'
        warning('Not yet implemented.');
    case 'no filter'
        dFF = R;
    otherwise 
        warning('Did not recognize filter type. Please select "gaussian", "exponential", or "no filter"');
end
        