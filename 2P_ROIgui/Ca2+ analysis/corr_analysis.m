%%
dt = 1/frame_rate;
ss = size(dFF);
num_ROIs = length(ROI_list);
corr_mat = zeros(num_ROIs);
dist_mat = corr_mat;
corr_mat = corr(dFF);

for i = 1:num_ROIs
    for j = 1:num_ROIs
        dist_mat(i,j) = sqrt((ROI_list(i).centerPos(1)-ROI_list(j).centerPos(1))^2 + (ROI_list(i).centerPos(2)-ROI_list(j).centerPos(2))^2);
    end
end

multiplot(dist_mat,corr_mat,'column')

mean_r = mean(corr_mat);
std_r = std(corr_mat);
sig_corrs = zeros(num_ROIs);

for i = 1:num_ROIs
    temp1 = find(corr_mat(i,:) > mean_r(i)+std_r(i));
    sig_corrs(i,temp1) = 1;
end

sig_dists = reshape(sig_corrs.*dist_mat,num_ROIs^2,1);
sig_dists(sig_dists == 0) = [];

figure
hist(sig_dists,[0:50:600]);

%% Make a figure for Mike, 20 September 2014

num_increments = 25;
max_dist = 501;
increment = (max_dist-1)/num_increments;
av_x_all = [];
av_y_all = [];

figure
hold on
for i = 1:num_ROIs
    for j = 1:num_increments
        inds = find(dist_mat(i,:) <= (increment*j + 1) & dist_mat(i,:) > (increment*(j-1) + 1));
        av_x(j) = mean(dist_mat(i,inds));
        av_y(j) = mean(corr_mat(i,inds));
    end
    plot(av_x,av_y,'color',[.8 .8 .8])
    av_x_mat(i,:) = av_x;
    av_y_mat(i,:) = av_y;
end

window = 20;
dist_mat_1d = reshape(dist_mat,num_ROIs^2,1);
corr_mat_1d = reshape(corr_mat,num_ROIs^2,1);
mov_av_corr = zeros(max_dist-window-1,1);
dist_points = window+1:max_dist;

for i = dist_points
    inds = find(dist_mat_1d < i + window & dist_mat_1d >= i-window);
    mov_av_corr(i-window) = mean(corr_mat_1d(inds));
end

plot(dist_points,mov_av_corr,'k')
hold off
    
    

%% Event extraction

dFF_thresh = 0.4;
isi_min = round(0.2/dt);
event_window = 6;
window = round(event_window/dt/2)*2;
num_locs = zeros(num_ROIs,1);
events_array = [];
events_array = cell(num_ROIs,1);
counter = 0;


for i = 1:num_ROIs
    [pks,locs] = findpeaks(dFF(:,i),'MINPEAKHEIGHT',dFF_thresh,'MINPEAKDISTANCE',isi_min);
    num_locs(i) = length(locs);
    temp_event_array = zeros(window,num_locs(i));
    temp_event_code = zeros(1,num_locs(i));
    early_pk_signal = 0;
    late_pk_signal = 0;
    
    for j = 1:num_locs(i)
        if (locs(j)-window/2+1)>0 && (locs(j)+window/2) <= ss(1)
            temp_inds = locs(j)-window/2+1:locs(j)+window/2;
            temp_code = 0;
        elseif (locs(j)-window/2+1)<0
            temp_inds = 1:window;
            temp_code = 1;
            early_pk_signal = 1;
        elseif (locs(j)+window/2)>ss(1)
            temp_inds = ss(1)-window+1:ss(1);
            temp_code = 2;
            late_pk_signal = 1;
        end
%         
%         if temp_inds(end) > ss(1)
%             ind_over = temp_inds(end) - ss(1);
%             temp_inds(end-ind_over:end)= temp_inds(end-ind_over);
%             counter = counter + 1
%         end
%         if temp_inds(1) > 
        temp_event_array(:,j) = dFF(temp_inds,i);
        temp_event_code(j) = temp_code;
    end
    
    if early_pk_signal == 1
        one_inds = find(temp_event_code == 1);
        if length(one_inds) > 1
            temp_event_array(:,one_inds(2:end)) = [];
            temp_event_code(one_inds(2:end)) = [];
            locs(one_inds(2:end)) = [];
            pks(one_inds(2:end)) = [];
        end
    end
    
    if late_pk_signal == 1
        two_inds = find(temp_event_code == 2);
        if length(two_inds) > 1
            temp_event_array(:,two_inds(2:end)) = [];
            temp_event_code(two_inds(2:end)) = [];
            locs(two_inds(2:end)) = [];
            pks(two_inds(2:end)) = [];
        end
    end
    
    events_array{i,1} = temp_event_array; 
    events_array{i,2} = temp_event_code;
    events_array{i,3} = locs;
    events_array{i,4} = pks;
end

step_window_sec = 3;
num_steps = 5;
step_window = round(-step_window_sec/dt:step_window_sec/dt/num_steps:step_window_sec/dt);

for i = 1:num_ROIs
    pk_locs = events_array{i,3};
    num_pks = length(pk_locs);
    r_mat = zeros(num_ROIs,num_pks);
    offset_mat = r_mat;
    for j = 1:num_pks
            r = zeros(1,num_ROIs);
            best_offset = r;
            event_trace = events_array{i,1}(:,j);
            
            for m = 1:length(step_window)
                loc_window = (pk_locs(j)+step_window(m)-window/2+1):(pk_locs(j)+step_window(m)+window/2);
                r_old = r;
                if events_array{i,2}(j) == 0 && loc_window(1)>0 && loc_window(end) <= ss(1)
                    corr_test =  corr(event_trace,dFF(loc_window,:));
                    r = r_old.*(r_old>=corr_test) + corr_test.*(corr_test>r_old);
                    best_offset = best_offset.*(r_old>=corr_test) + step_window(m).*(corr_test>r_old);
                elseif events_array{i,2}(j) == 1 && loc_window(1)>0
                   corr_test =  corr(event_trace,dFF(loc_window,:));
                    r = r_old.*(r_old>=corr_test) + corr_test.*(corr_test>r_old);
                    best_offset = best_offset.*(r_old>=corr_test) + step_window(m).*(corr_test>r_old);
                elseif events_array{i,2}(j) == 2 && loc_window(end) <= ss(1)
                    corr_test =  corr(event_trace,dFF(loc_window,:));
                    r = r_old.*(r_old>=corr_test) + corr_test.*(corr_test>r_old);
                    best_offset = best_offset.*(r_old>=corr_test) + step_window(m).*(corr_test>r_old);
                end
            end
            r_mat(:,j) = r;
            offset_mat(:,j) = best_offset;
    end
    events_array{i,5} = r_mat;
    events_array{i,6} = offset_mat;
end

%%

corr_window = 20; %in seconds
corr_step = 2; %in seconds
corr_substeps = floor(corr_step/dt);
tmax = floor(tvec_for_dFF(end));

for i = 1:num_ROIs
    




end



