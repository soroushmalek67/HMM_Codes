% Select which sleep epoch to analyze (1 = sleep1; 2 = sleep2; 3 = sleep3)
clear all;

ThPer = 60;
Epoch_Select = 3;
Dataset = '7165_11p';
DurationLimit = 10;
Run_Num = 5; %for the save filename
% Select which HMM file to load
cd(sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s', Dataset));

filename = sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s/HMM_all_train_2states_10rep_run70_sleep3Th%d.mat', Dataset,ThPer);
%filename = sprintf('Z:\\leanna_UPstate\\Template Matching\\Data\\%s\\HMM_all_train_2states_10rep_run2_sleep3.mat', Dataset);

load(filename,'numEpochs','Vit','HMmodels');
load('ts.mat');

switch Epoch_Select
    case 1
        load('sleep1_HMM_data.mat')
        sleep_start = e.epochs.sleep1(1);
        savefile = ['sleep1_UP_epochs_run',num2str(Run_Num),'.mat'];
    case 2
        load('sleep2_HMM_data.mat')
        sleep_start = e.epochs.sleep2(1);
        savefile = ['sleep2_UP_epochs_run',num2str(Run_Num),'.mat'];
    case 3
        load('sleep3_HMM_data.mat')
        
        if exist('e','var') 
           sleep_start = e.epochs.sleep3(1);
         else
           sleep_start = epochs.sleep3(1);
        end;
        
        savefile = ['sleep3_UP_epochs_run50',num2str(Run_Num),'_Th' num2str(ThPer) '.mat'];
    otherwise
        error('Wrong Value for Epoch_select')
end

numNeurons=max(ids);
UP = {};

for i = 1:numEpochs

    start_ts = HMmodels(i).StartTime;
    
    S = Vit{i};
    idx1 = find(S == 2);
    
    a = zeros(1,length(S));
    a(idx1) = 1;
    b = diff([0 a 0]);
    
    UP_start = find(b == 1);
    UP_end = find(b == -1);
    UP_end = UP_end-1;
    
% % % % % % %     if S(end) == 2;
% % % % % % %         UP_start(end) = [];     
% % % % % % %         UP_end(end) = [];
% % % % % % %     end
    
    UP{i} = [UP_start' UP_end'];
    UP{i} = UP{i}+start_ts;
        
end




UP_epochs = cell2mat(UP');

UP_epochs_full_ts = ((UP_epochs*10) + sleep_start);

UP_dur = (UP_epochs(:,2)-UP_epochs(:,1));

%%%%%% Added by Soroush later %%%

idx0 = find(UP_dur < DurationLimit);

UP_dur(idx0,:) = [];
UP_epochs(idx0,:) = [] ;
%%%%%% Added by Soroush later %%%

idx2 = find(UP_dur >= 700 & UP_dur <= 1300);
samp_size = round(length(idx2)*0.7);
idx2_rand = randsample(idx2,samp_size);
idx2_rand = sort(idx2_rand);
UP_train_rand = UP_epochs(idx2_rand,:);



UP_train = {};

for i = 1:length(UP_train_rand);
        idx3 = find(times > UP_train_rand(i,1) & times < UP_train_rand(i,2));
        tt = times(idx3)-UP_train_rand(i,1);
        T = UP_train_rand(i,2)-UP_train_rand(i,1)+1;
        n = ids(idx3);
        UP_train{i} = sleepOBS(tt,n,T);
end

idx4 = find(UP_dur >= 500 & UP_dur <= 1500);
UP_test = UP_epochs(idx4,:);

% Down States 

Down = {};

for i = 1:numEpochs

    start_ts = HMmodels(i).StartTime;
    
    S = Vit{i};
    Down_idx1 = find(S == 1);
    
    da = zeros(1,length(S));
    da(Down_idx1) = 1;
    db = diff([0 da 0]);
    
    Down_start = find(db == 1);
    Down_end = find(db == -1);
    Down_end = Down_end-1;
    
% % % % % %     if S(end) == 1;
% % % % % %         Down_start(end) = [];     
% % % % % %         Down_end(end) = [];
% % % % % %     end
% % % % % %     
    Down{i} = [Down_start' Down_end'];
    Down{i} = Down{i}+start_ts;
        
end

Down_epochs = cell2mat(Down');

Down_epochs_full_ts = ((Down_epochs*10) + sleep_start);

Down_dur = (Down_epochs(:,2)-Down_epochs(:,1));
% determine the number of spikes in DOWN and UP states
% total spikes in all epochs
epochs = zeros(length(HMmodels),2);
for i = 1:length(HMmodels)
    epochs(i,1) = HMmodels(i).StartTime;
    epochs(i,2) = HMmodels(i).EndTime;
end

all_epochs = cell(size(epochs,1),1);
for i = 1:size(epochs,1)
    all_epochs{i} = times(times > epochs(i,1) & times < epochs(i,2));
end
all_epochs_mat = cell2mat(all_epochs);
all_spikes_tot = length(all_epochs_mat);

% total spikes in all UP states + in each UP state
UP_spikes = cell(length(UP_epochs),1);

for i = 1:length(UP_epochs)
    UP_spikes{i} = times(times > UP_epochs(i,1) & times < UP_epochs(i,2));
    A_UP(i) = length(UP_spikes{i})*1000/(UP_dur(i)*numNeurons);
end
UP_spikes_mat = cell2mat(UP_spikes);
UP_spikes_tot = length(UP_spikes_mat);
 
% total spikes in each Down state

Down_spikes = cell(length(Down_epochs),1);
for i=1:length(Down_epochs)
    Down_spikes{i}= times(times > Down_epochs(i,1) & times < Down_epochs(i,2));
    A_Down(i) = length(Down_spikes{i})*1000/(numNeurons*Down_dur(i));
end
% total spikes in all DOWN states
DOWN_spikes_tot = all_spikes_tot - UP_spikes_tot;

DOWN_spikes_per = (DOWN_spikes_tot / all_spikes_tot) * 100
UP_spikes_per = (UP_spikes_tot / all_spikes_tot) * 100


% spikes during UP and DOWN states in each UP/DOWN epoch
UP_spikes_each = cell(size(epochs,1),1);
for i = 1:(size(epochs,1))
    UP_spikes_each{i} = UP_spikes_mat(UP_spikes_mat > epochs(i,1) & UP_spikes_mat < epochs(i,2));
end

UP_spikes_each_per = [];
for i = 1:(size(all_epochs,1))
    UP_spikes_each_per(i,1) = (length(UP_spikes_each{i}) / length(all_epochs{i}) * 100);
end
DOWN_spikes_each_per = 100 - UP_spikes_each_per;

UP_epochs_idx = [];
for i = 1:length(UP)
    UP_epochs_idx = [UP_epochs_idx; i*ones(size(UP{i}, 1), 1)];
end

figure;

min=0;
A_UPandDown_Tot = [A_UP, A_Down];
UPandDown_Duration_Tot = [UP_dur; Down_dur];
edgesx1= [0:.05:max(A_UPandDown_Tot)];
edgesx2= [0:.025:max(UPandDown_Duration_Tot/1000)];

subplot(1,2,1);

H1 = histogram(Down_dur/1000, edgesx2);
MaxHistogram1 = max(H1.Values);
title('Distribution of Down states', 'fontsize' ,14);
xlabel('Down Duration(s)','fontsize', 14);
ylabel('Number of Down states','fontsize', 14);
Mean_Down_Dur = mean(Down_dur/1000)
SEM_Down_Dur = std(Down_dur/1000)/sqrt(length(Down_dur/1000))
ylim([min MaxHistogram1]);

subplot(1,2,2);

H2= histogram(UP_dur/1000, edgesx2);
MaxHistogram2= max(H2.Values);
title('Distribution of UP states', 'fontsize' ,14);
xlabel('UP Duration(s)','fontsize', 14);
ylabel('Number of UP states','fontsize', 14);
Mean_UP_Dur = mean(UP_dur/1000)
SEM_UP_Dur = std(UP_dur/1000)/sqrt(length(UP_dur/1000))
ylim([min max(MaxHistogram1,MaxHistogram2)]);

% % 
% % saveas(gcf, sprintf('UPandDownDistributionSleep%d_%s_Th%d.fig', Epoch_Select,Dataset,ThPer));

figure;
subplot(1,2,1);
H3 = histogram(A_Down, edgesx1);
MaxHistogram3= max(H3.Values);
title('Firing rate of Down state', 'fontsize' ,14);
xlabel('Firing rate (Hz)','fontsize', 14);
ylabel('Number of Down states','fontsize', 14);
A_Down(isnan(A_Down)) = [] ;
Mean_Down_FiringRate = mean(A_Down)
SEM_Down_FiringRate = std(A_Down)/ sqrt(length(A_Down))
ylim([min MaxHistogram3]);

hold on;

subplot(1,2,2);
H4 = histogram(A_UP, edgesx1);
MaxHistogram4= max(H4.Values);
title('Firing rate of UP states', 'fontsize' ,14);
xlabel('Firing rate (Hz)','fontsize', 14);
ylabel('Number of UP states','fontsize', 14);
A_UP(A_UP==0) = [] ;
Mean_UP_FiringRate = nanmean(A_UP)
SEM_UP_FiringRate = std(A_UP)/ sqrt(length(A_UP))
ylim([min max(MaxHistogram3,MaxHistogram4)]);
% % 
% % saveas(gcf, sprintf('UPandDownSpikesNumbersleep%d_%s_Th%d.fig', Epoch_Select,Dataset,ThPer));
% % 
% % save(savefile,'UP_epochs','UP_epochs_full_ts','UP_train','UP_test','UP_spikes',...
% % 'filename','epochs','all_spikes_tot','UP_spikes_tot','DOWN_spikes_tot','A_UP','A_Down',...
% % 'UP_spikes_per','Down_epochs','DOWN_spikes_per','UP_spikes_each','all_epochs', 'UP_dur', 'Down_dur',...
% % 'UP_spikes_each_per','DOWN_spikes_each_per','UP_epochs_idx')
