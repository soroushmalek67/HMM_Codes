% Select which sleep epoch to analyze (1 = sleep1; 2 = sleep2; 3 = sleep3)
clear all;
Epoch_Select = 3;
Run_Num = 10; %for the save filename
% Select which HMM file to load
ThPer = 60;
DurationLimit = 10;
dataset = '7165_10p';
states = 3;
min = 0;
filename =  sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\New_Result\\%s\\HMM_all_train_DownandUPs_3states_10rep_run10_sleep3Th%d.mat', dataset, ThPer);
% Select which epoch file to load
% epoch_file = 'sleep3_UP_epochs_run7.mat';

load(filename);
load('ts.mat');

output_dir = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\New_Result\\%s' , dataset);



curr_output_dir = sprintf('%s\\NewResult', output_dir);
if ~exist(curr_output_dir, 'dir')
    mkdir(curr_output_dir);
end

curr_fn = [mfilename('fullpath') '.m'];
copyfile(curr_fn, curr_output_dir);


switch Epoch_Select
    case 1
        %         load(epoch_file, 'UP_test')
        savefile = [sprintf('sleep1_Down_UP%d_epochs_run', states-1), num2str(Run_Num),'_Th' num2str(ThPer) '.mat'];
        sleep_start = e.epochs.sleep1(1);
    case 2
        %         load(epoch_file, 'UP_test')
        savefile = [sprintf('sleep2_Down_UP%d_epochs_run', states), num2str(Run_Num),'.mat'];
        sleep_start = e.epochs.sleep2(1);
    case 3
        %         load(epoch_file, 'UP_test')
        load('sleep3_HMM_data.mat')
        savefile = [sprintf('sleep3_Down_UP%d_epochs_runn', states-1), num2str(Run_Num),'_Th' num2str(ThPer) '.mat'];
        
        if exist('e')
            
            sleep_start = e.epochs.sleep3(1);
        else
            sleep_start = epochs.sleep3(1);
        end
        
    otherwise
        error('Wrong Value for Epoch_select')
end

%% Only for the Down satets and 2 UP states %%%

UP1 = [];
UP2 = [];
Down = [];
start_with_UP1 = 0;
start_with_UP2 = 0;
start_with_Down = 0;

UP1_only = 0;
UP2_only = 0;
Down_only = 0;

for i = 1:numEpochs
    
    start_ts(i,1) = HMmodels(i).StartTime ;
    
    state_seq = Vit{i};
    
    if isempty(state_seq)
        continue
    end
    
    state1 = state_seq == 1;
    state2 = state_seq == 2;
    state3 = state_seq == 3;
    
    Down_diff = diff([0 state1 0]);
    UP1_diff = diff([0 state2 0]);
    UP2_diff = diff([0 state3 0]);
    
    Down_start_idx = find(Down_diff==1);
    Down_end_idx = find(Down_diff==-1);
    UP1_start_idx = find(UP1_diff == 1);
    UP1_end_idx = find(UP1_diff == -1) - 1;
    UP2_start_idx = find(UP2_diff == 1);
    UP2_end_idx = find(UP2_diff == -1) - 1;
    
    
    Down_temp = [];
    UP1_temp = [];
    UP2_temp = [];
    
    if ~isempty(UP1_start_idx) && ~isempty(UP2_start_idx) && ~isempty(Down_start_idx)
        if state_seq(6) == 1
            start_with_Down = start_with_Down +1 ;
        elseif state_seq == 2
            start_with_UP1 = start_with_UP1 + 1;
        else
            start_with_UP2 = start_with_UP2 + 1;
        end
    end
    
    % % % %     if isempty(UP1_start_idx)
    % % % %         UP2_only = UP2_only +1;
    % % % %     elseif isempty(UP2_start_idx)
    % % % %         UP1_only = UP1_only +1;
    % % % %     end
    
    for jj = 1 : length(Down_start_idx)
        Down_temp(jj,1) = start_ts(i,1) + Down_start_idx(jj) ;
        Down_temp(jj,2) = start_ts(i,1) + Down_end_idx(jj) ;
    end
    
    for ii = 1:length(UP1_start_idx)
        UP1_temp(ii,1) = start_ts(i,1) + UP1_start_idx(ii) ;
        UP1_temp(ii,2) = start_ts(i,1) + UP1_end_idx(ii) ;
    end
    
    for iii = 1:length(UP2_start_idx)
        UP2_temp(iii,1) = start_ts(i,1) + UP2_start_idx(iii) ;
        UP2_temp(iii,2) = start_ts(i,1) + UP2_end_idx(iii) ;
    end
    
    Down = [Down; Down_temp];
    UP1 = [UP1;UP1_temp];
    UP2 = [UP2;UP2_temp];
    
end

Down_dur = Down(:,2) - Down(:,1);

UP1_dur = UP1(:,2)-UP1(:,1);
idx3 = find(UP1_dur < DurationLimit);
UP1(idx3,:) = [];
UP1_dur(UP1_dur < DurationLimit) = [];

UP2_dur = UP2(:,2)-UP2(:,1);
idx4 = find(UP2_dur < DurationLimit);
UP2(idx4,:) = [];
UP2_dur(UP2_dur < DurationLimit) = [];

idx_Down = find(Down_dur < DurationLimit);
Down(idx_Down,:) = [];
Down_dur(idx_Down) = [];

Down_full = (Down*10) + sleep_start;
UP1_full = (UP1*10) + sleep_start;
UP2_full = (UP2*10) + sleep_start;

% Added by Soroush later

Down_epochs = Down;
UP1_epochs = UP1;
UP2_epochs = UP2;

Down_spikes = cell(length(Down_epochs),1);
for j = 1 : length(Down_epochs)
    
   Down_spikes{j} = times(times > Down_epochs(j,1) & times < Down_epochs (j,2));
   A_Down(j) = length(Down_spikes{j})*1000/Down_dur(j);
    
end
Down_spikes_mat = cell2mat(Down_spikes);
Down_spikes_tot = length(Down_spikes_mat);

UP1_spikes = cell(length(UP1_epochs),1);
for i = 1:length(UP1_epochs)
    
    UP1_spikes{i} = times(times > UP1_epochs(i,1) & times < UP1_epochs(i,2));
    A_UP1(i) = length(UP1_spikes{i})*1000/UP1_dur(i);
    
end
UP1_spikes_mat = cell2mat(UP1_spikes);
UP1_spikes_tot = length(UP1_spikes_mat);


UP2_spikes = cell(length(UP2_epochs),1);
for i = 1:length(UP2_epochs)
    
    UP2_spikes{i} = times(times > UP2_epochs(i,1) & times < UP2_epochs(i,2));
    A_UP2(i) = length(UP2_spikes{i})*1000/UP2_dur(i);
    
end
UP2_spikes_mat = cell2mat(UP2_spikes);
UP2_spikes_tot = length(UP2_spikes_mat);




figure;
min=0;

A_Down_UPs_Tot = [A_Down/numNeurons, A_UP1/numNeurons, A_UP2/numNeurons];
Down_UPs_Duration_Tot = [Down_dur; UP1_dur; UP2_dur];
edgesx1= [0:.4:max(A_Down_UPs_Tot)];
edgesx2= [0:.025:max(Down_UPs_Duration_Tot/1000)];

axax1=subplot(1,3,1);

H1 = histogram(Down_dur/1000, edgesx2);
MaxHistogram1 = max(H1.Values);
title('Distribution of Down States', 'fontsize', 12);
xlabel('Down State Duration (s)','fontsize', 12);
ylabel('Number of Down State','fontsize', 12);
Mean_Down_Dur = mean(Down_dur / 1000)

hold on;

axax2 = subplot(1,3,3);

H2 = histogram(UP1_dur/1000, edgesx2);
MaxHistogram2 = max(H2.Values);
title('Distribution of UP Type1', 'fontsize' ,12);
xlabel('UP Type1 Duration (s)','fontsize', 12);
ylabel('Number of UP Type1','fontsize', 12);
Mean_UP1_Dur = mean(UP1_dur/1000)

hold on;


axax3=subplot(1,3,2);

H3 = histogram(UP2_dur/1000, edgesx2);
MaxHistogram2= max(H3.Values);
title('Distribution of UP Type2', 'fontsize' ,12);
xlabel('UP Type2 Duration(s)','fontsize', 12);
ylabel('Number of UP Type2','fontsize', 12);
Mean_UP2_Dur = mean(UP2_dur/1000)
%ylim([min MaxHistogram1]);
linkaxes([axax1, axax2, axax3], 'xy');

saveas(gcf, ['NewDownandTwoUPsDuration' num2str(Epoch_Select) '_Th' num2str(ThPer) '.fig' ]);

%%%%


figure;

axcd1 = subplot(1,3,1);
H4 = histogram(A_Down/numNeurons, edgesx1);
MaxHistogram4 = max(H4.Values);
title('Firing rate of Down State', 'fontsize' ,12);
xlabel('Firing rate (Hz)', 'fontsize' ,12);
ylabel('Number of Down State', 'fontsize' ,12);
Mean_A_Down = mean(A_Down/numNeurons)
SEM_Down_FiringRate = std(A_Down/numNeurons)/sqrt(length(A_Down/numNeurons));

hold on;

axcd2 = subplot(1,3,3);
H5 = histogram(A_UP1/numNeurons, edgesx1);
MaxHistogram3 = max(H5.Values);
title('Firing rate of UP Type1', 'fontsize' ,12);
xlabel('Firing rate (Hz)','fontsize', 12);
ylabel('Number of UP Type1','fontsize', 12);
Mean_A_UP1 = mean(A_UP1/numNeurons)
%%%ylim([min max(MaxHistogram3,MaxHistogram4)]);
SEM_UP1_FiringRate = std(A_UP1/numNeurons)/ sqrt(length(A_UP1/numNeurons));




hold on;

axcd3 = subplot(1,3,2);

H6 = histogram(A_UP2/numNeurons, edgesx1);
MaxHistogram4= max(H6.Values);
title('Firing rate of UP Type2', 'fontsize' ,12);
xlabel('Firing rate (Hz)','fontsize', 12);
ylabel('Number of UP Type2','fontsize', 12);
Mean_A_UP2 = mean(A_UP2/numNeurons)
%%%%%%ylim([min MaxHistogram4]);
SEM_UP2_FiringRate = std(A_UP2/numNeurons)/ sqrt(length(A_UP2/numNeurons));

linkaxes([axcd1, axcd2, axcd3], 'xy');

saveas(gcf, ['NewDownandTwoUPsFiringRates' num2str(Epoch_Select) '_Th' num2str(ThPer) '.fig' ]);



%%
%%%%% Finding UP states  - UP type1 and UP type 2 combined! %%%%%


UP = {};

for i = 1:numEpochs

    start_ts(i,1) = HMmodels(i).StartTime;
    
    S_UP = Vit{i};
    S_UP(S_UP == 3)=2; 
    idx1 = find(S_UP == 2);
    
    a = zeros(1,length(S_UP));
    a(idx1) = 1;
    b = diff([0 a 0]);
    
    UP_start = find(b == 1);
    UP_end = find(b == -1);
    UP_end = UP_end-1;
    
    
    UP{i} = [UP_start' UP_end'];
    UP{i} = UP{i}+start_ts(i,1);
        
end




UP_epochs = cell2mat(UP');

UP_epochs_full_ts = ((UP_epochs*10) + sleep_start);

UP_dur = (UP_epochs(:,2)-UP_epochs(:,1));



%%
%%%%%% This is the older script for 3 states or more in which I used eval  
%%%%%% and sprintf for different state %%%%

% % % % % 
% % % % % for j=1:states
% % % % %     
% % % % %     eval(sprintf('UP%d= [];',j))
% % % % %     eval(sprintf('start_with_UP%d = 0;',j))
% % % % %     eval(sprintf('UP%d_only = 0;', j))
% % % % %     
% % % % % end
% % % % % 
% % % % % for i = 1:numEpochs
% % % % %     
% % % % %     start_ts = HMmodels(i).StartTime;
% % % % %     
% % % % %     state_seq = Vit{i};
% % % % %     
% % % % %     if isempty(state_seq)
% % % % %         continue
% % % % %     end
% % % % %     
% % % % %     for j=1:states
% % % % %         
% % % % %         eval(sprintf('state%d = state_seq == %d;', j,j))
% % % % %         eval(sprintf('UP%d_diff = diff([0 state%d 0]);', j,j))
% % % % %         eval(sprintf('UP%d_start_idx = find(UP%d_diff == 1);', j, j))
% % % % %         eval(sprintf('UP%d_end_idx = find(UP%d_diff == -1) - 1;', j, j))
% % % % %         eval(sprintf('UP%d_temp = [];', j))
% % % % %         
% % % % %     end
% % % % %     
% % % % %     
% % % % %     
% % % % %     for j=1:states
% % % % %         
% % % % %         for ii = 1: eval(sprintf('length(UP%d_start_idx)', j))
% % % % %             
% % % % %             eval(sprintf('UP%d_temp(ii,1) = start_ts + UP%d_start_idx(ii) - 1;', j , j))
% % % % %             eval(sprintf('UP%d_temp(ii,2) = start_ts + UP%d_end_idx(ii) - 1;' ,j ,j))
% % % % %             
% % % % %         end
% % % % %     end
% % % % %     
% % % % %     
% % % % %     for j=1:states
% % % % %         
% % % % %         eval(sprintf('UP%d = [UP%d ; UP%d_temp];', j, j, j))
% % % % %         
% % % % %     end
% % % % % end
% % % % % 
% % % % % for j=1:states
% % % % %     
% % % % %     eval(sprintf('UP%d_dur = UP%d(:,2) - UP%d(:,1);', j, j, j))
% % % % %     eval(sprintf('idx%d = find( UP%d_dur < 10);', j,j))
% % % % %     eval(sprintf('UP%d(idx%d,:) = [];', j, j))
% % % % %     eval(sprintf('UP%d_dur(UP%d_dur < 10) = [];',j ,j ))
% % % % %     eval(sprintf('UP%d_full = ((UP%d*10) + sleep_start);' ,j,j ))
% % % % %     
% % % % % end
% % % % % 
% % % % % %%%%%% Added for UP figures%%%%%%%%
% % % % % 
% % % % % 
% % % % % for j=1:states
% % % % %     
% % % % %     eval(sprintf('UP%d_epochs = UP%d;', j ,j))
% % % % %     eval(sprintf('UP%d_spikes = cell(length(UP%d_epochs),1);' ,j,j))
% % % % %     
% % % % %     if eval(sprintf('~length(UP%d_epochs) == 0;', j))
% % % % %         
% % % % %         for iii = 1:length(eval(sprintf('UP%d_epochs(:,1)' ,j)))
% % % % %             
% % % % %             eval(sprintf('UP%d_spikes{iii} = times(times > UP%d_epochs(iii,1) & times < UP%d_epochs(iii,2));' ,j ,j ,j))
% % % % %             eval(sprintf('A_UP%d(iii) = length(UP%d_spikes{iii})*1000/UP%d_dur(iii);', j ,j, j))
% % % % %             eval(sprintf('UP%d_spikes_mat = cell2mat(UP%d_spikes);',j,j))
% % % % %             eval(sprintf('UP%d_spikes_tot = length(UP%d_spikes_mat);',j,j))
% % % % %             
% % % % %             
% % % % %         end
% % % % %         
% % % % %     else
% % % % %         eval(sprintf('A_UP%d =[];', j))
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%
% % % % % cd(curr_output_dir);
% % % % % 
% % % % % A_UP_Tot = 0;
% % % % % UP_Duration_Tot = [];
% % % % % 
% % % % % for ii=1:states
% % % % %     
% % % % %     eval(sprintf('Max_UP%d = A_UP%d/numNeurons;' ,ii , ii))
% % % % %     eval(sprintf('A_UP_Tot = [A_UP_Tot, A_UP%d/numNeurons];', ii))
% % % % %     eval(sprintf('UP_Duration_Tot = [UP_Duration_Tot , UP%d_dur''];', ii))
% % % % %     
% % % % %     eval(sprintf('Mean_FiringRate_UP%d = mean(A_UP%d/numNeurons);', ii , ii))
% % % % %     eval(sprintf('Firing_Mean_All(ii) = Mean_FiringRate_UP%d;', ii))
% % % % %     
% % % % % end
% % % % % 
% % % % % Firing_Mean_All(isnan(Firing_Mean_All)) = 0 ;
% % % % % [Max_Firing_state, Max_State_ID] = max(Firing_Mean_All);
% % % % % 
% % % % % [Sorted_state, ID_sort] = sort(Firing_Mean_All, 'descend');
% % % % % 
% % % % % 
% % % % % for ll = 1 : length(ID_sort)
% % % % %     
% % % % %     ID_ID = ID_sort(ll);
% % % % %     
% % % % %     eval(sprintf('UP%d_sorted = UP%d;' , ll, ID_ID))
% % % % %     eval(sprintf('UP%d_sorted_full = UP%d_full;' , ll , ID_ID))
% % % % %     eval(sprintf('A_UP%d_sort = A_UP%d;' , ll , ID_ID))
% % % % %     eval(sprintf('UP%d_sorted_dur = UP%d_dur;', ll, ID_ID))
% % % % %     
% % % % %     
% % % % % end
% % % % % 
% % % % % 
% % % % % for ll = 1 : length(ID_sort)
% % % % %     
% % % % %     eval(sprintf('UP%d = UP%d_sorted;' , ll, ll))
% % % % %     eval(sprintf('UP%d_full = UP%d_sorted_full;' , ll , ll))
% % % % %     
% % % % % end
% % % % % 
% % % % % edgesx1= [0:.1:max(A_UP_Tot)];
% % % % % edgesx2= [0:.025:max(UP_Duration_Tot/1000)];
% % % % % edgesx3 = [0:.0125:max(UP_Duration_Tot/1000)];
% % % % % min = 0 ;
% % % % % 
% % % % % 
% % % % % 
% % % % % figure;
% % % % % 
% % % % % for j = 1 : states
% % % % %     
% % % % %     eval(sprintf('ax%d = subplot(2,2,%d);', j, j))
% % % % %     eval(sprintf('H%d= histogram(UP%d_sorted_dur/1000, edgesx2);',j,j))
% % % % %     eval(sprintf('MaxHistogram%d = max(H%d.Values);', j ,j ))
% % % % %     title(sprintf('UP Type%d',j), 'fontsize' ,10);
% % % % %     xlabel('Duration(s)','fontsize', 10);
% % % % %     ylabel(sprintf('Number of UP Type%d',j), 'fontsize', 10);
% % % % %     eval(sprintf('Mean_UP%d_sorted_dur = mean(UP%d_sorted_dur/1000)', j ,j));
% % % % %     P(j) = eval(sprintf('Mean_UP%d_sorted_dur',j));
% % % % %     LG = legend(['Mean Duration = ' eraseBetween(num2str(P(j)),5,length(num2str(P(j)))) ' ms']);
% % % % %     LG.FontSize = 7;
% % % % %     hold on;
% % % % % end
% % % % % 
% % % % % figtitle('Duration of UP different types')
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % linkaxes([ax1,ax2,ax3],'y');
% % % % % 
% % % % % save(savefile,'UP1','UP1_full','UP2','UP2_full','UP3','UP3_full','filename','start_with_UP1','start_with_UP2','UP1_only','UP2_only')
% % % % % 
% % % % % 
% % % % % 
% % % % % saveas(gcf, sprintf('%dUPsandDownDistribution_Th60Motionless.fig', states-1));
% % % % % 
% % % % % 
% % % % % %%
% % % % % 
% % % % % 
% % % % % figure;
% % % % % 
% % % % % for j = 1 : states
% % % % %     eval(sprintf('axd%d = subplot(2,2,%d);', j, j))
% % % % %     eval(sprintf('H_UP%d= histogram(A_UP%d_sort/numNeurons, edgesx1);',j,j))
% % % % %     eval(sprintf('MaxHistogram%d = max(H_UP%d.Values);', j ,j ))
% % % % %     title(sprintf('UP Type%d',j), 'fontsize' ,10);
% % % % %     xlabel('Firing rate (Hz)','fontsize', 10);
% % % % %     ylabel(sprintf('Number of UP Type%d',j), 'fontsize', 10);
% % % % %     eval(sprintf('Mean_A_UP%d_sort = mean(A_UP%d_sort/numNeurons)', j ,j));
% % % % %     eval(sprintf('ylim([min MaxHistogram%d]);',j));
% % % % %     eval(sprintf('SEM_UP%d_FiringRate = std(A_UP%d_sort/numNeurons)/ sqrt(length(A_UP%d_sort/numNeurons))', j,j,j))
% % % % %     P(j) = eval(sprintf('Mean_A_UP%d_sort',j));
% % % % %     LG = legend(['Mean Firing Rate = ' eraseBetween(num2str(P(j)),5,length(num2str(P(j)))) ' ms']);
% % % % %     LG.FontSize = 7;
% % % % %     hold on;
% % % % % end
% % % % % 
% % % % % figtitle('Firing Rate of UP different types')
% % % % % 
% % % % % saveas(gcf, sprintf('%dUPsandDownssFiringRate_Th60Motionless.fig', states-1));
% % % % % 

cd(sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\New_Result\\%s', dataset));
save(savefile,'UP_epochs','UP_epochs_full_ts','UP_dur','Down','Down_full','Down_dur','A_Down','UP1','UP1_full','UP1_dur','A_UP1','A_UP2','UP2','UP2_dur','UP2_full','filename','start_with_UP1','start_with_UP2','UP1_only','UP2_only')