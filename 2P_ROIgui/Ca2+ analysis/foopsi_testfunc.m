function [Nhat Phat] = foopsi_testfunc(fluor_trace, im_rate)

F = double(fluor_trace);

% set simulation metadata
T       = length(F); % # of time steps
V.dt    = 1/im_rate;  % time step size

% initialize params
P.a     = 1;    % observation scale
P.b     = 0;    % observation bias
tau     = 1.5;    % decay time constant
P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
P.lam   = 0.1;  % firing rate = lam*dt
P.sig   = 0.1;  % standard deviation of observation noise 

% fast oopsi
[Nhat Phat] = fast_oopsi(F,V,P);

%% plot results
% figure
% tvec=0:V.dt:(T-1)*V.dt;
% h(1)=subplot(211); plot(tvec,F); axis('tight'), ylabel('F (au)')
% h(2)=subplot(212); plot(tvec,Nhat,'r','linewidth',1), axis('tight'), ylabel('fast')
% xlabel('time (sec)') 
% linkaxes(h,'x')


