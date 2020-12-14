% HiddenMarkovModels.m


clear all;
ThPer = 60; 
States = 2 ; %2-20 
dataset = '8482_15p';
Repetitions = 10; %10
Run_Num = 10; %for the save filename
Epoch_Select = 3; % select which sleep epoch to analyze (1 = sleep1; 2 =
% sleep2; 3 = sleep3)
%UPdir = 'Z:\soroush_UPstate\Template Matching\Data\7165_16p'
output_dir = sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s', dataset);
% % % % data_file = sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s/sleep3_UP_epochs_run505.mat', dataset);
data_file = sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s/sleep3_UP_epochs_run505_Th%dMotionless.mat', dataset, ThPer);

cd(output_dir);
% Load data and UP/DOWN periods:
%***************************************
switch Epoch_Select
    case 1
        load('sleep1_HMM_data.mat') 
        load(data_file, 'UP_epochs')
        sleep = 1;
    case 2
        load('sleep2_HMM_data.mat')
        load(data_file, 'UP_epochs')
        sleep = 2;
    case 3
        load('sleep3_HMM_data.mat')
        load(data_file, 'UP_epochs')
        sleep = 3;
    otherwise
        error('Wrong Value for Epoch_select')
end

filename = ['E:\HMM - UP&Down\Soroush\Data\New_Result\',num2str(dataset),'\HMM_UPonly_all_train_',num2str(States),'states_',num2str(Repetitions),'_rep_run',num2str(Run_Num),'_sleep',num2str(sleep),'_Th',num2str(ThPer),'Motionless.mat'];

% load('sleep3_HMM_data.mat')
% filename = 'data0round.txt';
% [ids times]=textread(filename,'%f %f');

UP_train = {};
for i = 1:size(UP_epochs,1)
    idx = find(times > UP_epochs(i,1) & times < UP_epochs(i,2));
    tt = times(idx)-UP_epochs(i,1);
    T = UP_epochs(i,2)-UP_epochs(i,1)+1;
    n = ids(idx);
    UP_train{i} = sleepOBS(tt,n,T);
end

numNeurons=max(ids);
numEpochs=length(UP_epochs(:,1));

Vit  = cell(1,numEpochs);
Vit_all = cell(1,numEpochs);
BB = cell(1,numEpochs);

tic

[HMM,EMIS_all,TRANS_all,Likelihood]=sleepHMMTraining_new(UP_train,States,numNeurons,Repetitions,200);  

toc
tic

for ep=1:numEpochs
    
    TrainStart=UP_epochs(ep,1);
    TrainEnd=UP_epochs(ep,2);

    T=TrainEnd-TrainStart+1;

    % data in Train set:
    i=find(times<TrainEnd & times>TrainStart);
    tt=times(i);
    n=ids(i);
    tt=tt-TrainStart;
    
    % Build observations
    Obs=sleepOBS(tt,n,T);    
    

    %Training:
    %**********************************************************************
        % sleepHMMTraining(Obs,States,numNeurons,Repetitions,Iterations)
   
%     [HMM,EMIS_all,TRANS_all,Likelihood]=sleepHMMTraining_new(Obs,States,numNeurons,Repetitions,200);  
%     
% 
    TRANS=HMM.TRANSITION;
    EMIS=HMM.EMISSION;
    
    HMmodels(ep).TRANS= TRANS;
    HMmodels(ep).EMIS=EMIS;
    HMmodels(ep).TRANS_all=TRANS_all;
    HMmodels(ep).EMIS_all=EMIS_all;
    HMmodels(ep).Likelihood = Likelihood;
    HMmodels(ep).tt = tt;
    HMmodels(ep).n = n;
    HMmodels(ep).T = T;
    HMmodels(ep).Obs = Obs;
    HMmodels(ep).StartTime = TrainStart;
    HMmodels(ep).EndTime = TrainEnd;
    
    %Decoding stage (Viterbi)
    %**********************************************************************
    try
        S = hmmviterbi(Obs, TRANS, EMIS); % S: most likely sequence of states

        S_all = {};
        for ii = 1:Repetitions
            S_all{ii} = hmmviterbi(Obs, TRANS_all{ii}, EMIS_all{ii});
        end

        Vit{ep}=S;

        Vit_all{ep} = S_all;
        
    catch e 
        continue;
    end
    % Figures:
    %**********************************************************************
%     set(figure,'Position',[250 100 1000 400],'Color','w');
%     
%     subplot(3,1,2:3)
% 
%     plot(tt,n/numNeurons,'k.')
%     hold on
%     plot(1.2*S-1.3)
%     
%     tx1=2*floor(T/7);
%     tx2=tx1+3*floor(T/7);
% %     tx1 = 0;
% %     tx2 = TrainEnd;
% %     set(gca,'XColor','w','YColor','w','YLim',[-.3 1.5],'XLim',[tx1 tx2])
%     set(gca,'YLim',[-.3 1.5],'XLim',[tx1 tx2])
%     title(['epoch ' num2str(ep)],'FontSize',16)
% %     xlim([0 TrainEnd])
%    
%     pause(0);
    
end

toc

UP_test = UP_epochs;

% save HMMSTATES Vit Vit_all
% save MarkovModels HMmodels

save(filename, 'Vit', 'Vit_all', 'BB', 'HMmodels', 'States', 'numNeurons','Repetitions', 'numNeurons', 'numEpochs', 'data_file', 'UP_train','UP_test')
