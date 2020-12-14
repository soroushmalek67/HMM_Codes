% HiddenMarkovModels.m



clear all;

ThPer = 60;
States = 3;
Repetitions = 1 ;
Run_Num = 10 ; %for the save filename
Epoch_Select = 3 ; % select which sleep epoch to analyze (1 = sleep1; 2 =
% sleep2; 3 = sleep3)

UpDownFile = ['UpDownPeriodsSleep3_allprct_Th' num2str(ThPer)]; % select which file to load

% Load data and UP/DOWN periods:
%***************************************
load(UpDownFile)
switch Epoch_Select
    case 1
        load('sleep1_HMM_data.mat')
        sleep = 1;
    case 2
        load('sleep2_HMM_data.mat')
        sleep = 2;
    case 3
        load('sleep3_HMM_data.mat')
        sleep = 3;
    otherwise
        error('Wrong Value for Epoch_select')
end

% load('sleep3_HMM_data.mat')
% filename = 'data0round.txt';
% [ids times]=textread(filename,'%f %f');

% load('sleep3_mtnless.mat')
% load UpDownPeriods

if epochs(1,2) == 0
    epochs(1,:) = [];
end

if epochs(end,1) == epochs (end,2)
    epochs(end,:) = [];
end

% % % % % epochs = epochs(2:end, :);

numNeurons=max(ids);
numEpochs=length(epochs(:,1));

Vit = cell(1,numEpochs);
Vit_all = cell(1,numEpochs);
AA = cell(1,numEpochs);
train_all = cell(1,numEpochs);

for e = 1:numEpochs
    TrainStart = epochs(e,1);
    TrainEnd = epochs(e,2);
    
    T = TrainEnd-TrainStart+1;
    
    % data in Train set:
    i=find(times<TrainEnd & times>TrainStart);
    tt=times(i);
    n = ones(length(i),1);
    tt=tt-TrainStart;
    
    % Build observations
    train_all{e}=sleepOBS(tt,n,T);
end

tic
[HMM,EMIS_all,TRANS_all,Likelihood]=sleepHMMTraining_new(train_all,States,numNeurons,Repetitions,200);
toc

for ep=1:numEpochs
    tic
    TestStart=epochs(ep,1);
    TestEnd=epochs(ep,2);

    T=TestEnd-TestStart+1;

    % data in Train set:
    i=find(times<TestEnd & times>TestStart);
    tt=times(i);
    n = ones(length(i),1);
    tt=tt-TestStart;
    
    % Build observations
    Obs=sleepOBS(tt,n,T);    
    

    %Training:
    %**********************************************************************
        % sleepHMMTraining(Obs,States,numNeurons,Repetitions,Iterations)
   
%     [HMM,EMIS_all,TRANS_all,Likelihood]=sleepHMMTraining_new(Obs,States,numNeurons,Repetitions,200);  
    

    TRANS=HMM.TRANSITION;
    EMIS=HMM.EMISSION;
    
    HMmodels(ep).TRANS=TRANS;
    HMmodels(ep).EMIS=EMIS;
    HMmodels(ep).TRANS_all=TRANS_all;
    HMmodels(ep).EMIS_all=EMIS_all;
    HMmodels(ep).Likelihood = Likelihood;
    HMmodels(ep).tt = tt;
    HMmodels(ep).n = n;
    HMmodels(ep).T = T;
    HMmodels(ep).Obs = Obs;
    HMmodels(ep).StartTime = TestStart;
    HMmodels(ep).EndTime = TestEnd;
    
    %Decoding stage (Viterbi)
    %**********************************************************************
    nn=ids(i);
    S = hmmviterbi(Obs, TRANS, EMIS); % S: most likely sequence of states
    
    S_all = {};
    for ii = 1:Repetitions
        S_all{ii} = hmmviterbi(Obs, TRANS_all{ii}, EMIS_all{ii});
    end
    
    Vit{ep}=S;
    
    Vit_all{ep} = S_all;

    % Figures:
    %**********************************************************************
    set(0,'DefaultFigureWindowStyle','normal');
    H = set(figure,'Position',[250 100 1000 400],'Color','w');
    
    subplot(3,1,2:3)

    plot(tt,nn/numNeurons,'k.')
    hold on
    
    if S == 1 
    plot(1.2*S-1.3, 'g')
    elseif S == 2
    plot(1.2*S-1.3, 'r')
    elseif S == 3
    plot(S - 0.7, 'b')
    
    AA{ep} = S;
    tx1=2*floor(T/7);
    tx2=tx1+3*floor(T/7);
%     tx1 = 0;
%     tx2 = TrainEnd;
%     set(gca,'XColor','w','YColor','w','YLim',[-.3 1.5],'XLim',[tx1 tx2])
    set(gca,'YLim',[-.3 1.5],'XLim',[tx1 tx2])
    title(['epoch ' num2str(ep) '\Th' num2str(ThPer)],'FontSize',16)
     xlabel('Time (ms)','fontsize', 14);
%     xlim([0 TrainEnd])
   saveas(gcf, ['EpochsDetectedWithAllOnes' num2str(ep) 'Sleep' num2str(sleep) 'Th' num2str(ThPer) '_' num2str(States) 'States.fig']);
   set(H,'WindowStyle','docked');
    pause(0);
    toc
end

% save HMMSTATES Vit Vit_all
% save MarkovModels HMmodels

filename = ['HMM_all_train_With_AllOnes',num2str(States),'states_',num2str(Repetitions),'rep_run',num2str(Run_Num),'_sleep',num2str(sleep),'Th',num2str(ThPer),'.mat'];
save(filename, 'Vit', 'Vit_all', 'S','AA','HMmodels', 'States', 'Repetitions', 'numNeurons', 'numEpochs', 'UpDownFile')

