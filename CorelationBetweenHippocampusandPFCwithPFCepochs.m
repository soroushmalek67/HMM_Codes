%%%% Correlation between PFC and Hippocampus during PFC epochs %%%%
%%
clear all;

sigma = 0.01;
sz = 0.05;

Dataset = '8482_15p';
 
 tfiles_dir_Hippo = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus.PFCepochs\\tfiles', Dataset);
 epochs_dir = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/Hippocampus.PFCepochs', Dataset);
 epochs_dir_epochs = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC', Dataset);
 tfiles_dir_Cortex = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\PFC\\tfiles', Dataset);


cd(tfiles_dir_Hippo)
tfile = FindFiles('*.t');
S_Hippo = LoadSpikes(tfile);
Q_Hippo = MakeQfromS(S_Hippo,10);

 cd(epochs_dir)
 load('ts.mat')
 sleep = e.epochs.sleep3;
 savefile = 'sleep3_HMM_data.mat';
 
Q_cut_Hippo = Restrict(Q_Hippo,sleep(1),sleep(2));
QD_cut_Hippo = Data(Q_cut_Hippo);
QD_range_Hippo = Range(Q_cut_Hippo);

% extract neuron IDs and spike times from Q matrix in Hippocampus
[row_Hippo col_Hippo] = find(QD_cut_Hippo==1);
us_ids_Hippo = col_Hippo;
us_times_Hippo = QD_range_Hippo(row_Hippo);
us_times_Hippo = (us_times_Hippo-sleep(1))/10;
us_times_Hippo = floor(us_times_Hippo);

% sort times and neuron IDs
[times_Hippo,ix_Hippo] = sort(us_times_Hippo);
ids_Hippo = us_ids_Hippo(ix_Hippo);
numNeurons_Hippo=max(ids_Hippo);
%%%% Finding each Neurons activity individually %%%%
cd(epochs_dir_epochs);
UpDownFile = 'UpDownPeriodsSleep3_allprct_Th84';
load(UpDownFile);

Neurons_Hippo = cell(numNeurons_Hippo,1);

Neuronstime = cell(numNeurons_Hippo,1);
for iiii=1: numNeurons_Hippo
    Neurons_ids = find(ids_Hippo==iiii);

    Neurons_Hippo_Activity{iiii} = ids_Hippo(Neurons_ids);
    Neurons_Hippo_Activity_Time{iiii} =  times_Hippo(Neurons_ids);
    
end
%%
%%%%%% Building the function of activity for each neurons in Hippocampus %%%%

Activity = cell(numNeurons_Hippo, length(epochs));
Activity_Function = cell(numNeurons_Hippo, length(epochs));
Activity_ones = cell(numNeurons_Hippo, length(epochs));

for jj = 1: length(epochs)
    for jjj = 1: numNeurons_Hippo
       
       Activity{jjj, jj} = Neurons_Hippo_Activity_Time{jjj} (Neurons_Hippo_Activity_Time{jjj} > epochs(jj,1)& Neurons_Hippo_Activity_Time{jjj} < epochs(jj,2));
       Activity_Function{jjj,jj} = zeros (epochs(jj,2)- epochs(jj,1),1); 
       Activity_oness = Activity{jjj,jj};
       Activity_Function{jjj,jj}(Activity_oness-epochs(jj,1))=1;
       %Activity_Function =  Activity_Function(cellfun(Activity_oness == 1, jjj,jj));
       %Activity_Function{jjj,jj} = Activity_Function{Activity{jjj,jj}==1};
    end
end

%%
 %%%%% Finding the Convulotion for the neurons activity in Hippocampus %%%%% 

sigma_one_Neuron = .01;
sz_one_Neuron = .05;
x = linspace (-sz_one_Neuron/2, sz_one_Neuron/2); 
gaussFilter  = exp (-x.^2/(2*sigma_one_Neuron^2)); 
gaussFilter = gaussFilter / sum(gaussFilter); 


for ii = 1: length (epochs)
    for jj = 1: numNeurons_Hippo
 
     yfilt_one_neuron_Hippo{jj,ii} = filter(gaussFilter,1,(Activity_Function{jj,ii}));
     yfilt_notshifted_one_Hippo{jj,ii} = conv ((Activity_Function{jj,ii}), gaussFilter, 'same');

    end
end


% % % % % % % % figure;
% % % % % % % % plot(yfilt_one_neuron_Hippo{10,1})
% % % % % % % % figure;
% % % % % % % % plot(yfilt_notshifted_one_Hippo{10,1})
 
%%
%%%%%% Building the function of activity for all of neurons in each epochs in Hippocampus %%%%%
 
Activity_Function_Total = cell(length(epochs),1);

for ep=1:length(epochs)
    
    TestStart=epochs(ep,1);
    TestEnd=epochs(ep,2);

    T=TestEnd-TestStart+1;

    % data in Train set:
    i=find(times_Hippo<TestEnd & times_Hippo>TestStart);
    ttt=times_Hippo(i);
    n=ids_Hippo(i);
    ttt=ttt-TestStart;
    
    Activity_Function_Total{ep,1} = zeros(epochs(ep,2)- epochs(ep,1),1); 
    Activity_Function_Total{ep,1}(ttt)=1;
    
    for iirr=1: length(ttt)-1
        
        if ttt(iirr,1) == ttt(iirr+1,1)
        tttrrr = ttt(iirr,1);
        Activity_Function_Total{ep,1}(tttrrr)=2;
        
        end
    end
    
    
        for iirr=1: length(ttt)-2
        
        if ttt(iirr,1) == ttt(iirr+1,1) == ttt(iirr+2,1)
        tttrrrr = ttt(iirr,1);
        Activity_Function_Total{ep,1}(tttrrrr)=3;
        
        end
    end
    
end


%%
%%%%% Buidling Populatin Firing rate for Hippocampus with the same method we used in HMM %%%%%%


center_Hippo = cell (length(epochs), 1);
numSp_Hippo = cell (length(epochs),1);


for epp = 1 : length(epochs)

 
    
Start = epochs (epp,1);
VeryEnd = epochs (epp,2);

End = VeryEnd - Start;


i_HMM_Hippo=find(times_Hippo<VeryEnd & times_Hippo>Start);
tt_HMM_Hippo=times_Hippo(i_HMM_Hippo);
n_HMM_Hippo =ids_Hippo(i_HMM_Hippo);
tt_HMM_Hippo=tt_HMM_Hippo -Start;


bin=50; % chose bin size for population firing rate

numWind=floor(End/bin);
center_Hippo{epp,1}=zeros(1,numWind);
numSp_Hippo{epp,1}=zeros(1,numWind);


for w=1:numWind
   Wstart=(w-1)*bin;
   Wend=Wstart+bin;
   center_Hippo{epp,1}(w)=(Wend-Wstart)/2+Wstart;
   Sp_Hippo=find(tt_HMM_Hippo<Wend & tt_HMM_Hippo>Wstart); %%% Finding how many neurons fire during this bin
   numSp_Hippo{epp,1}(w) =1000*length(Sp_Hippo)/(numNeurons_Hippo*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
   
   if mod(w,50)==0
   
  
   
   end
end



% % % % % % % % 
% % % % % % % % set(figure,'Position',[250 100 760 760],'Color','w');
% % % % % % % % orient tall;
% % % % % % % % 
% % % % % % % % %subplot(3,1,1)
% % % % % % % % 
% % % % % % % % plot(center_Hippo{epp,1}/1000,numSp_Hippo{epp,1},'k')
% % % % % % % % title(['Population firing rate for Hippocampus for epochs ' num2str(epp)],'FontSize',14)
% % % % % % % % ylabel('spikes/s','FontSize',14)
% % % % % % % % set(gca,'FontSize',12)
% % % % % % % % %box off
% % % % % % % % text(.8,.9,['Window=' num2str(bin) 'ms'],'Units','Normalized','FontSize',11)
% % % % % % % % 
% % % % % % % % 
 end


%%  
 %%%%% Finding the Convulotion for the all the activities of neurons %%%%% 


x = linspace (-sz/2, sz/2); 
gaussFilter  = exp (-x.^2/(2*sigma^2)); 
gaussFilter = gaussFilter / sum(gaussFilter); 

for  ii = 1: length(epochs)
    
yfilt_Hippo{ii,1} = filter(gaussFilter,1,(Activity_Function_Total{ii,1}));
yfilt_notshifted_Hippo{ii,1} = conv ((Activity_Function_Total{ii,1}), gaussFilter, 'same');

end

% % % % % % % % % figure;
% % % % % % % % % plot(yfilt_Hippo{1,1})
% % % % % % % % % figure;
% % % % % % % % % plot(yfilt_notshifted_Hippo{1,1})

%%
%%%%%%%%%%%%%%%%%%%%%% Do all the steps for PFC %%%%%%%%%%%%%%%%%%%%

cd(tfiles_dir_Cortex)
tfile = FindFiles('*.t');
S= LoadSpikes(tfile);
Q = MakeQfromS(S,10);

 cd(epochs_dir_epochs)
 load('ts.mat')
 sleep = e.epochs.sleep3;
 savefile = 'sleep3_HMM_data.mat';
 
Q_cut = Restrict(Q,sleep(1),sleep(2));
QD_cut = Data(Q_cut);
QD_range = Range(Q_cut);

% extract neuron IDs and spike times from Q matrix in PFC
[row col] = find(QD_cut==1);
us_ids = col;
us_times = QD_range(row);
us_times = (us_times-sleep(1))/10;
us_times = floor(us_times);

% sort times and neuron IDs
[times,ix] = sort(us_times);
ids = us_ids(ix);
numNeurons=max(ids);
%%%% Finding each Neurons activity individually %%%%
cd(epochs_dir_epochs);
UpDownFile = 'UpDownPeriodsSleep3_allprct_Th84';
load(UpDownFile);

Neurons_PFC = cell(numNeurons,1);

Neuronstime_PFC = cell(numNeurons,1);
for iii=1: numNeurons
    Neurons_ids_PFC = find(ids==iii);

    Neurons_PFC_Activity{iii} = ids(Neurons_ids_PFC);
    Neurons_PFC_Activity_Time{iii} =  times(Neurons_ids_PFC);
    
end

%%
%%%%%% Building the function of activity for each neurons in PFC %%%%

Activity_PFC = cell(numNeurons, length(epochs));
Activity_Function_PFC = cell(numNeurons, length(epochs));
Activity_ones_PFC = cell(numNeurons, length(epochs));
for jj = 1: length(epochs)
    for jjj = 1: numNeurons
       
       Activity_PFC{jjj, jj} = Neurons_PFC_Activity_Time{jjj} (Neurons_PFC_Activity_Time{jjj} > epochs(jj,1)& Neurons_PFC_Activity_Time{jjj} < epochs(jj,2));
       Activity_Function_PFC{jjj,jj} = zeros (epochs(jj,2)- epochs(jj,1),1); 
       Activity_oness_PFC = Activity_PFC{jjj,jj};
       Activity_Function_PFC{jjj,jj}(Activity_oness_PFC-epochs(jj,1))=1;
       
     
       %Activity_Function =  Activity_Function(cellfun(Activity_oness == 1, jjj,jj));
       %Activity_Function{jjj,jj} = Activity_Function{Activity{jjj,jj}==1};
    end
end

%%
 %%%%% Finding the Convulotion for the neuron activity in PFC %%%%% 

sigma_one_Neuron = .01;
sz_one_Neuron = .05;
x = linspace (-sz_one_Neuron/2, sz_one_Neuron/2); 
gaussFilter  = exp (-x.^2/(2*sigma_one_Neuron^2)); 
gaussFilter = gaussFilter / sum(gaussFilter); 

for ii=1 : length(epochs)
    for jj = 1 : numNeurons
        
     yfilt_one_neuron_PFC{jj,ii} = filter(gaussFilter,1,(Activity_Function_PFC{jj,ii}));
     yfilt_notshifted_one_PFC{jj,ii} =conv ((Activity_Function_PFC{jj,ii}), gaussFilter, 'same');

    end
end

% % % % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % % plot(yfilt_one_neuron_PFC{1,1})
% % % % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % % plot(yfilt_notshifted_one_PFC{1,1})
 
%%
%%%%%% Building the function of activity for all of neurons in each epochs in PFC %%%%%
 
Activity_Function_PFC_Total = cell(length(epochs),1);

for ep=1:length(epochs)
    
    TestStart=epochs(ep,1);
    TestEnd=epochs(ep,2);

    T=TestEnd-TestStart+1;

    % data in Train set:
    i=find(times<TestEnd & times>TestStart);
    tt=times(i);
    n=ids(i);
    tt=tt-TestStart;
    
    Activity_Function_PFC_Total{ep,1} = zeros(epochs(ep,2)- epochs(ep,1),1); 
    Activity_Function_PFC_Total{ep,1}(tt)=1;
    
    
    for ir=1: length(tt)-1
        
        if tt(ir,1) == tt(ir+1,1)
        ttrr = tt(ir,1);
        Activity_Function_PFC_Total{ep,1}(ttrr)=2;
        
        end
    end
    
    
    
    
    for ir=1: length(tt)-2
        
        if tt(ir,1) == tt(ir+1,1) == tt(ir+2,1)
        tttrrr = tt(ir,1);
        Activity_Function_PFC_Total{ep,1}(tttrrr)=3;
        
        end
    end
    
     %%% Activity_Function_All_Epochs_PFC = Activity_Function_PFC_Total{ep,1}
    
end
    

%%
%%%%% Buidling Populatin Firing rate with the same method we used in HMM %%%%%%



center = cell (length(epochs), 1);
numSp = cell (length(epochs),1);


for epp = 1 : length(epochs)

    
    
    
    
Start = epochs (epp,1);
VeryEnd = epochs (epp,2);

End = VeryEnd - Start;


i_HMM=find(times<VeryEnd & times>Start);
tt_HMM=times(i_HMM);
n_HMM =ids(i_HMM);
tt_HMM=tt_HMM-Start;


bin=50; % chose bin size for population firing rate

numWind=floor(End/bin);
center{epp,1}=zeros(1,numWind);
numSp{epp,1}=zeros(1,numWind);


for w=1:numWind
   Wstart=(w-1)*bin;
   Wend=Wstart+bin;
   center{epp,1}(w)=(Wend-Wstart)/2+Wstart;
   Sp=find(tt_HMM<Wend & tt_HMM>Wstart); %%% Finding how many neurons fire during this bin
   numSp{epp,1}(w) =1000*length(Sp)/(numNeurons*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
   
   if mod(w,50)==0
  
  
   
   end
end



% % % % % % % 
% % % % % % % set(figure,'Position',[250 100 760 760],'Color','w');
% % % % % % % orient tall;
% % % % % % % 
% % % % % % % %subplot(3,1,1)
% % % % % % % 
% % % % % % % plot(center{epp,1}/1000,numSp{epp,1},'k')
% % % % % % % title('Population firing rate','FontSize',14)
% % % % % % % ylabel('spikes/s','FontSize',14)
% % % % % % % set(gca,'FontSize',12)
% % % % % % % %box off
% % % % % % % text(.8,.9,['Window=' num2str(bin) 'ms'],'Units','Normalized','FontSize',11)


end

%%
 %%%%% Finding the Convulotion for the all the activities of neurons in PFC %%%%% 


x = linspace (-sz/2, sz/2); 
gaussFilter  = exp (-x.^2/(2*sigma^2)); 
gaussFilter = gaussFilter / sum(gaussFilter); 

for ii = 1 : length(epochs)
    
yfilt_PFC{ii,1} = filter(gaussFilter,1,(Activity_Function_PFC_Total{ii,1}));
yfilt_notshifted_PFC{ii,1} = conv ((Activity_Function_PFC_Total{ii,1}), gaussFilter, 'same');

end

% % % % % % % % figure;
% % % % % % % % plot(yfilt_PFC{1,1})
% % % % % % % % figure;
% % % % % % % % plot(yfilt_notshifted_PFC{1,1})

corrPFCHippo_Sample = corr((yfilt_notshifted_Hippo{1,1}),(yfilt_notshifted_PFC{1,1}));
corrcoefPFCHippo_sample = corrcoef((yfilt_notshifted_Hippo{1,1}), (yfilt_notshifted_PFC{1,1}));

% corrPFCHippo_Sample2 =  corr((yfilt_notshifted_Hippo{16,1}), (yfilt_notshifted_PFC{16,1}));

% % % % % % 
% % % % % % for ii = 1 : length(epochs)
% % % % % %     
% % % % % %  H(ii)= figure;
% % % % % %  
% % % % % %  ax1 = subplot(2,1,1);
% % % % % %  plot(yfilt_notshifted_PFC{ii,1});
% % % % % %  title(['Population Activity for epoch' num2str(ii) ' in PFC '],'FontSize',12);
% % % % % %  xlabel('Time (ms)','fontsize', 14);
% % % % % %  ylabel('Amplitude Ratio','fontsize', 12);
% % % % % % 
% % % % % %  ax2 = subplot(2,1,2);
% % % % % %  plot(yfilt_notshifted_Hippo{ii,1});
% % % % % %  title(['Population Activity for epoch' num2str(ii) ' in Hippocampus '],'FontSize',12);
% % % % % %  xlabel('Time (ms)','fontsize', 14);
% % % % % %  ylabel('Amplitude Ratio','fontsize', 12);
% % % % % %  
% % % % % %  linkaxes([ax1, ax2],'y');
% % % % % %  
% % % % % %  saveas(gcf, ['PopulationActivityforEpochs ' num2str(ii) '.fig']);
% % % % % %  
% % % % % %  set(H(ii),'WindowStyle','docked');
% % % % % %  
% % % % % % end



%%%%%%%%%%%%%%%%% Cross Correlation for PFC and Hippocampus activities %%%%
% % % % % % % % % % for i=1:length(epochs)
% % % % % % % % % %   
% % % % % % % % % % [acor, lag] = xcorr((Activity_Function_PFC_Total{i,1}),(Activity_Function_Total{i,1}), epochs(i,2)-epochs(i,1),'coeff');
% % % % % % % % % % [~,I] = max(abs(acor));
% % % % % % % % % % lagDiff = lag(I);
% % % % % % % % % % %%stem (lag, acor)
% % % % % % % % % % figure; plot(lag,acor)
% % % % % % % % % %    
% % % % % % % % % % end
%%%%%%%%%%%%%%%%%% Cross Correlation %%%%%

% % % % [acor, lag] = xcorr((Activity_Function_PFC_Total{2,1}),(Activity_Function_PFC_Total{2,1}), 'coeff');
% % % % [~,I] = max(abs(acor));
% % % % lagDiff = lag(I)
% % % % stem (lag, acor)
% % % % figure; plot(lag,acor)

% % % % a3 = gca;
% % % % a3.XTick = sort([-3000:1000:3000 lagDiff]);
%%
% % % % % % % % % % % % % % for ii =1 : length(epochs)
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % [acor, lag] = xcorr(yfilt_notshifted_PFC{ii,1},yfilt_notshifted_Hippo{ii,1}, 1000,'biased');
% % % % % % % % % % % % % % [~,I] = max(abs(acor));
% % % % % % % % % % % % % % lagDiff = lag(I);
% % % % % % % % % % % % % % %stem (lag, acor)
% % % % % % % % % % % % % % set(0,'DefaultFigureWindowStyle','normal');
% % % % % % % % % % % % % % H(ii)= figure; 
% % % % % % % % % % % % % % plot(lag,acor);
% % % % % % % % % % % % % % DelayBetween (ii) = finddelay(yfilt_notshifted_PFC{ii,1}, yfilt_notshifted_Hippo{ii,1}, 400); 
% % % % % % % % % % % % % % D2 (ii) = -(lag(acor == max(acor)));
% % % % % % % % % % % % % % title(['Cross Correlation epoch of PFC and Hippocampus in PFC epoch ' num2str(ii)],'FontSize',12);
% % % % % % % % % % % % % % xlabel('Delay(ms)','fontsize', 14);
% % % % % % % % % % % % % % ylabel('Amplitude Ratio','fontsize', 12);
% % % % % % % % % % % % % % %legend('Max Correlation Coefficient', 'Delay Between')
% % % % % % % % % % % % % % %MaxText = ['Max= ', num2str(max(abs(acor)) ) newline  'Delay= ', num2str(DelayBetween) ' ms' ];
% % % % % % % % % % % % % % %text((epochs(ii,2)-epochs(ii,1))/4 - 1000, 0.1, MaxText , 'HorizontalAlignment','center');
% % % % % % % % % % % % % % saveas(gcf, ['CrossCorrelationBiasedPFCandHippocampusforEpoch' num2str(ii) '.fig']);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % set(H(ii),'WindowStyle','docked');
% % % % % % % % % % % % % % end

%%
for ii =1 : length(epochs)
    
[acorr, lagg] = xcorr(yfilt_notshifted_Hippo{ii,1},yfilt_notshifted_PFC{ii,1}, 202,'unbiased');
[Delayy, Locs] = findpeaks(acorr, lagg, 'SortStr','descend');
Location(ii,1) = Locs(1);
    
[acor, lag] = xcorr(yfilt_notshifted_Hippo{ii,1},yfilt_notshifted_PFC{ii,1}, 1000,'unbiased');
[~,I] = max(abs(acor));
lagDiff = lag(I);
%stem (lag, acor)
set(0,'DefaultFigureWindowStyle','normal');
H(ii)= figure; 
plot(lag,acor);
DelayBetween2 (ii) = finddelay(yfilt_notshifted_Hippo{ii,1}, yfilt_notshifted_PFC{ii,1},400); 
D22 (ii) = -(lag(acor == max(acor)));
 
findpeaks(acor,lag)
% % [pks,locs] = findpeaks(acor,lag);
text(Locs+.02,Delayy,num2str((1:numel(Delayy))'));

title(['Cross Correlation epoch of Hippocampus and PFC in PFC epoch ' num2str(ii)],'FontSize',12);
xlabel('Delay(ms)','fontsize', 12);
ylabel('Amplitude Ratio','fontsize', 12);

%legend('Max Correlation Coefficient', 'Delay Between')
%MaxText = ['Max= ', num2str(max(abs(acor)) ) newline  'Delay= ', num2str(DelayBetween) ' ms' ];
%text((epochs(ii,2)-epochs(ii,1))/4 - 1000, 0.1, MaxText , 'HorizontalAlignment','center');
saveas(gcf, ['CrossCrollationUnbiasedHippocampusandPFCforEpoch' num2str(ii) '.fig']);


set(H(ii),'WindowStyle','docked');
end
%%

figure;
histbins = [-200 : 20 :200]; 
histogram (Location, histbins);
title('Distribution of Lags between Hippocampus and PFC for all epochs of sleep3', 'fontsize' ,10);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, 'DistributionofLagsBetweenHippocampusandPFC.fig');

%%
%%% Finding all of UP States in PFC and Hippocampus and their cross-correlation %%%%

cd (epochs_dir);

load('sleep3_UP_epochs_run65.mat');

 yfilt_notshifted_UP_PFC = cell (length(UP_epochs), 1);
 yfilt_notshifted_UP_Hippo = cell (length(UP_epochs), 1);

for ii = 1 : length(epochs)  
    
for iii = 1: length(UP_epochs)
    

    
    if UP_epochs (iii,1) > epochs (ii,1) & UP_epochs (iii,2) < epochs (ii,2)
        
        UP_epochs_starts (iii) = UP_epochs (iii,1) - epochs(ii,1)  ;
        UP_epochs_ends (iii) = UP_epochs (iii,2) - epochs(ii,1) -1  ;
        
        yfilt_notshifted_UP_PFC{iii,1} = yfilt_notshifted_PFC{ii,1}( UP_epochs_starts (iii): UP_epochs_ends (iii),1);
        yfilt_notshifted_UP_Hippo{iii,1} = yfilt_notshifted_Hippo{ii,1}( UP_epochs_starts (iii): UP_epochs_ends (iii),1);

    end
    
    
    
end 
end


for iii = 1 :length(UP_epochs)

     [acor_Up, lag_UP] = xcorr(yfilt_notshifted_UP_Hippo{iii,1},yfilt_notshifted_UP_PFC{iii,1}, 200, 'unbiased');
     [Delay_UP, Locs_UP] = findpeaks(acor_Up, lag_UP, 'SortStr','descend');
     if ~ length(Locs_UP) == 0
        Location_UP(iii,1) = Locs_UP(1);
     else 
        Location_UP(iii,1) = 0;
     end

end

figure;
histbins_UP = [-200 : 20 :200]; 
histogram (Location_UP, histbins_UP);

title('Distribution of Lags between Hippocampus and PFC for all UPs of sleep3',  'fontsize' ,10);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUPs.fig']);

%%
%%% Finding all of UP States in PFC with a window 100 ms after and before %%%%

UP_Window = 20; %ms


Activity_Function_PFC_UP_Total = cell(length(UP_epochs),1);

for ep_UP=1:length(UP_epochs)
    
    TestStart=UP_epochs(ep_UP,1)- UP_Window;
    TestEnd=UP_epochs(ep_UP,2) + UP_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    i_UP=find(times<TestEnd & times>TestStart);
    tt_UP=times(i_UP);
    n_UP =ids(i_UP);
    tt_UP=tt_UP-TestStart;
    
    Activity_Function_PFC_UP_Total{ep_UP,1} = zeros(UP_epochs(ep_UP,2)- UP_epochs(ep_UP,1)+2*UP_Window,1); 
    Activity_Function_PFC_UP_Total{ep_UP,1}(tt_UP)=1;
    
    
    for ir=1: length(tt_UP)-1
        
        if tt_UP(ir,1) == tt_UP(ir+1,1)
        ttrr = tt_UP(ir,1);
        Activity_Function_PFC_UP_Total{ep_UP,1}(ttrr)=2;
        
        end
    end
    
    
    
    
    for ir=1: length(tt_UP)-2
        
        if tt_UP(ir,1) == tt_UP(ir+1,1) == tt_UP(ir+2,1)
        tttrrr = tt_UP(ir,1);
        Activity_Function_PFC_UP_Total{ep_UP,1}(tttrrr)=3;
        
        end
    end
    
     %%% Activity_Function_All_Epochs_PFC = Activity_Function_PFC_Total{ep,1}
    
end

%%
%%% Finding all of UP States in Hippocampus with a window 100 ms after and before %%%%



Activity_Function_Hippo_UP_Total = cell(length(UP_epochs),1);

for epp_UP=1:length(UP_epochs)
    
    TestStart=UP_epochs(epp_UP,1)- UP_Window;
    TestEnd=UP_epochs(epp_UP,2) + UP_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    ii_UP=find(times_Hippo<TestEnd & times_Hippo>TestStart);
    ttt_UP =times_Hippo(ii_UP);
    nn_UP=ids_Hippo(ii_UP);
    ttt_UP=ttt_UP-TestStart;
    
    Activity_Function_Hippo_UP_Total{epp_UP,1} = zeros(UP_epochs(epp_UP,2)- UP_epochs(epp_UP,1) +2*UP_Window,1); 
    Activity_Function_Hippo_UP_Total{epp_UP,1}(ttt_UP)=1;
    
    for iirr=1: length(ttt_UP)-1
        
        if ttt_UP(iirr,1) == ttt_UP(iirr+1,1)
        tttrrr = ttt_UP(iirr,1);
        Activity_Function_Hippo_UP_Total{epp_UP,1}(tttrrr)=2;
        
        end
    end
    
    
        for iirr=1: length(ttt_UP)-2
        
        if ttt_UP(iirr,1) == ttt_UP(iirr+1,1) == ttt_UP(iirr+2,1)
        tttrrrr = ttt_UP(iirr,1);
        Activity_Function_Hippo_UP_Total{epp_UP,1}(tttrrrr)=3;
        
        end
    end
    
end



%%
%%%  Cross-correlation of all of UP States in PFC and Hippocampus looking over 100 ms after and before them %%%%


for ii = 1 : length(UP_epochs)
    
yfilt_PFC_UP_win{ii,1} = filter(gaussFilter,1,(Activity_Function_PFC_UP_Total{ii,1}));
yfilt_notshifted_PFC_UP_win{ii,1} = conv ((Activity_Function_PFC_UP_Total{ii,1}), gaussFilter, 'same');

yfilt_Hippo_UP_win{ii,1} = filter(gaussFilter,1,(Activity_Function_Hippo_UP_Total{ii,1}));
yfilt_notshifted_Hippo_UP_win{ii,1} = conv ((Activity_Function_Hippo_UP_Total{ii,1}), gaussFilter, 'same');

end



for iii = 1 :length(UP_epochs)

     [acor_UP_win, lag_UP_win] = xcorr(yfilt_notshifted_Hippo_UP_win{iii,1},yfilt_notshifted_PFC_UP_win{iii,1}, 200, 'unbiased');
     [Delay_UP_win, Locs_UP_win] = findpeaks(acor_UP_win, lag_UP_win, 'SortStr','descend');
     if ~ length(Locs_UP_win) == 0
        Location_UP_win(iii,1) = Locs_UP_win(1);
     else 
        Location_UP_win(iii,1) = 0;
     end

end

figure;
histbins_UP = [-200 : 20 :200]; 
histogram (Location_UP_win, histbins_UP);

title(['Distribution of Lags between Hippocampus and PFC for all UPs of sleep3 with window of ' num2str(UP_Window)], 'fontsize' ,8);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUPswithwindowof ' num2str(UP_Window) '.fig']);


%%% Finding all of UP1 States in Hippocampus and PFC and their cross-correlation %%%%

load('sleep3_UP1_UP2_epochs_run_extra3.mat');


yfilt_notshifted_UP1_PFC = cell (length(UP1), 1);
 yfilt_notshifted_UP1_Hippo = cell (length(UP1), 1);

for ii = 1 : length(epochs)  
    
for iii = 1: length(UP1)
    

    
    if UP1(iii,1) > epochs (ii,1) & UP1(iii,2) < epochs (ii,2)
        
        UP1_starts (iii) = UP1 (iii,1) - epochs(ii,1)  ;
        UP1_ends (iii) = UP1 (iii,2) - epochs(ii,1) -1  ;
        
        yfilt_notshifted_UP1_PFC{iii,1} = yfilt_notshifted_PFC{ii,1}( UP1_starts (iii): UP1_ends (iii),1);
        yfilt_notshifted_UP1_Hippo{iii,1} = yfilt_notshifted_Hippo{ii,1}( UP1_starts (iii): UP1_ends (iii),1);

    end
    
    
    
end 
end


for iii = 1 :length(UP1)

     [acor_UP1, lag_UP1] = xcorr(yfilt_notshifted_UP1_Hippo{iii,1},yfilt_notshifted_UP1_PFC{iii,1}, 200, 'unbiased');
     [Delay_UP1, Locs_UP1] = findpeaks(acor_UP1, lag_UP1, 'SortStr','descend');
     if ~ length(Locs_UP1) == 0
        Location_UP1(iii,1) = Locs_UP1(1);
     else 
        Location_UP1(iii,1) = 0;
     end

end

figure;
histbins_UP1 = [-200 : 20 :200]; 
histogram (Location_UP1, histbins_UP1);

title('Distribution of Lags between Hippocampus and PFC for all UP1s of sleep3',  'fontsize' ,10);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUP1s.fig']);


%%
%%% Finding all of UP1 States in PFC with a window 20 ms after and before %%%%

UP1_Window = 20; %ms


Activity_Function_PFC_UP1_Total = cell(length(UP1),1);

for ep_UP1 =1:length(UP1)
    
    TestStart=UP1(ep_UP1,1)- UP1_Window;
    TestEnd=UP1(ep_UP1,2) + UP1_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    i_UP1=find(times<TestEnd & times>TestStart);
    tt_UP1=times(i_UP1);
    n_UP1 =ids(i_UP1);
    tt_UP1=tt_UP1-TestStart;
    
    Activity_Function_PFC_UP1_Total{ep_UP1,1} = zeros(UP1(ep_UP1,2)- UP1(ep_UP1,1)+2*UP1_Window,1); 
    Activity_Function_PFC_UP1_Total{ep_UP1,1}(tt_UP1)=1;
    
    
    for ir=1: length(tt_UP1)-1
        
        if tt_UP1(ir,1) == tt_UP1(ir+1,1)
        ttrr = tt_UP1(ir,1);
        Activity_Function_PFC_UP1_Total{ep_UP1,1}(ttrr)=2;
        
        end
    end
    
    
    
    
    for ir=1: length(tt_UP1)-2
        
        if tt_UP1(ir,1) == tt_UP1(ir+1,1) == tt_UP1(ir+2,1)
        tttrrr = tt_UP1(ir,1);
        Activity_Function_PFC_UP1_Total{ep_UP1,1}(tttrrr)=3;
        
        end
    end
    
     %%% Activity_Function_All_Epochs_PFC = Activity_Function_PFC_Total{ep,1}
    
end

%%
%%% Finding all of UP1 States in Hippocampus with a window 20 ms after and before %%%%



Activity_Function_Hippo_UP1_Total = cell(length(UP1),1);

for epp_UP1=1:length(UP1)
    
    TestStart=UP1(epp_UP1,1)- UP1_Window;
    TestEnd=UP1(epp_UP1,2) + UP1_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    ii_UP1=find(times_Hippo<TestEnd & times_Hippo>TestStart);
    ttt_UP1 =times_Hippo(ii_UP1);
    nn_UP1=ids_Hippo(ii_UP1);
    ttt_UP1=ttt_UP1-TestStart;
    
    Activity_Function_Hippo_UP1_Total{epp_UP1,1} = zeros(UP1(epp_UP1,2)- UP1(epp_UP1,1) +2*UP1_Window,1); 
    Activity_Function_Hippo_UP1_Total{epp_UP1,1}(ttt_UP1)=1;
    
    for iirr=1: length(ttt_UP1)-1
        
        if ttt_UP1(iirr,1) == ttt_UP1(iirr+1,1)
        tttrrr = ttt_UP1(iirr,1);
        Activity_Function_Hippo_UP1_Total{epp_UP1,1}(tttrrr)=2;
        
        end
    end
    
    
        for iirr=1: length(ttt_UP1)-2
        
        if ttt_UP1(iirr,1) == ttt_UP1(iirr+1,1) == ttt_UP1(iirr+2,1)
        tttrrrr = ttt_UP1(iirr,1);
        Activity_Function_Hippo_UP1_Total{epp_UP1,1}(tttrrrr)=3;
        
        end
    end
    
end

%%
%%%  Cross-correlation of all of UP1 States in Hippocampus and PFC looking over 20 ms after and before them %%%%


for ii = 1 : length(UP1)
    
yfilt_PFC_UP1_win{ii,1} = filter(gaussFilter,1,(Activity_Function_PFC_UP1_Total{ii,1}));
yfilt_notshifted_PFC_UP1_win{ii,1} = conv ((Activity_Function_PFC_UP1_Total{ii,1}), gaussFilter, 'same');

yfilt_Hippo_UP1_win{ii,1} = filter(gaussFilter,1,(Activity_Function_Hippo_UP1_Total{ii,1}));
yfilt_notshifted_Hippo_UP1_win{ii,1} = conv ((Activity_Function_Hippo_UP1_Total{ii,1}), gaussFilter, 'same');

end



for iii = 1 :length(UP1)

     [acor_UP1_win, lag_UP1_win] = xcorr(yfilt_notshifted_Hippo_UP1_win{iii,1},yfilt_notshifted_PFC_UP1_win{iii,1}, 200, 'unbiased');
     [Delay_UP1_win, Locs_UP1_win] = findpeaks(acor_UP1_win, lag_UP1_win, 'SortStr','descend');
     if ~ length(Locs_UP1_win) == 0
        Location_UP1_win(iii,1) = Locs_UP1_win(1);
     else 
        Location_UP1_win(iii,1) = 0;
     end

end

figure;
histbins_UP1 = [-200 : 20 :200]; 
histogram (Location_UP1_win, histbins_UP1);

title(['Distribution of Lags between Hippocampus and PFC for all UP1s of sleep3 with window of ' num2str(UP_Window)], 'fontsize' ,8);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUP1swithwindowof ' num2str(UP_Window) '.fig']);


%%

%%% Finding all of UP2 States in Hippocampus and PFC and their cross-correlation %%%%

load('sleep3_UP1_UP2_epochs_run_extra3.mat');


yfilt_notshifted_UP2_PFC = cell (length(UP2), 1);
 yfilt_notshifted_UP2_Hippo = cell (length(UP2), 1);

for ii = 1 : length(epochs)  
    
for iii = 1: length(UP2)
    

    
    if UP2(iii,1) > epochs (ii,1) & UP2(iii,2) < epochs (ii,2)
        
        UP2_starts (iii) = UP2 (iii,1) - epochs(ii,1)  ;
        UP2_ends (iii) = UP2 (iii,2) - epochs(ii,1) -1  ;
        
        yfilt_notshifted_UP2_PFC{iii,1} = yfilt_notshifted_PFC{ii,1}( UP2_starts (iii): UP2_ends (iii),1);
        yfilt_notshifted_UP2_Hippo{iii,1} = yfilt_notshifted_Hippo{ii,1}( UP2_starts (iii): UP2_ends (iii),1);

    end
    
    
    
end 
end


for iii = 1 :length(UP2)

     [acor_UP2, lag_UP2] = xcorr(yfilt_notshifted_UP2_Hippo{iii,1},yfilt_notshifted_UP2_PFC{iii,1}, 200, 'unbiased');
     [Delay_UP2, Locs_UP2] = findpeaks(acor_UP2, lag_UP2, 'SortStr','descend');
     if ~ length(Locs_UP2) == 0
        Location_UP2(iii,1) = Locs_UP2(1);
     else 
        Location_UP2(iii,1) = 0;
     end

end

figure;
histbins_UP2 = [-200 : 20 :200]; 
histogram (Location_UP2, histbins_UP2);

title('Distribution of Lags between Hippocampus and PFC for all UP2s of sleep3',  'fontsize' ,10);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUP2s.fig']);


%%
%%% Finding all of UP2 States in PFC with a window 20 ms after and before %%%%

UP2_Window = 20; %ms


Activity_Function_PFC_UP2_Total = cell(length(UP2),1);

for ep_UP2 =1:length(UP2)
    
    TestStart=UP2(ep_UP2,1)- UP2_Window;
    TestEnd=UP2(ep_UP2,2) + UP2_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    i_UP2=find(times<TestEnd & times>TestStart);
    tt_UP2=times(i_UP2);
    n_UP2 =ids(i_UP2);
    tt_UP2=tt_UP2-TestStart;
    
    Activity_Function_PFC_UP2_Total{ep_UP2,1} = zeros(UP2(ep_UP2,2)- UP2(ep_UP2,1)+2*UP2_Window,1); 
    Activity_Function_PFC_UP2_Total{ep_UP2,1}(tt_UP2)=1;
    
    
    for ir=1: length(tt_UP2)-1
        
        if tt_UP2(ir,1) == tt_UP2(ir+1,1)
        ttrr = tt_UP2(ir,1);
        Activity_Function_PFC_UP2_Total{ep_UP2,1}(ttrr)=2;
        
        end
    end
    
    
    
    
    for ir=1: length(tt_UP2)-2
        
        if tt_UP2(ir,1) == tt_UP2(ir+1,1) == tt_UP2(ir+2,1)
        tttrrr = tt_UP2(ir,1);
        Activity_Function_PFC_UP2_Total{ep_UP2,1}(tttrrr)=3;
        
        end
    end
    
     %%% Activity_Function_All_Epochs_PFC = Activity_Function_PFC_Total{ep,1}
    
end

%%
%%% Finding all of UP2 States in Hippocampus with a window 20 ms after and before %%%%



Activity_Function_Hippo_UP2_Total = cell(length(UP2),1);

for epp_UP2=1:length(UP2)
    
    TestStart=UP2(epp_UP2,1)- UP2_Window;
    TestEnd=UP2(epp_UP2,2) + UP2_Window;

    T=TestEnd-TestStart+1;

    % data in Train set:
    ii_UP2=find(times_Hippo<TestEnd & times_Hippo>TestStart);
    ttt_UP2 =times_Hippo(ii_UP2);
    nn_UP2=ids_Hippo(ii_UP2);
    ttt_UP2=ttt_UP2-TestStart;
    
    Activity_Function_Hippo_UP2_Total{epp_UP2,1} = zeros(UP2(epp_UP2,2)- UP2(epp_UP2,1) +2*UP2_Window,1); 
    Activity_Function_Hippo_UP2_Total{epp_UP2,1}(ttt_UP2)=1;
    
    for iirr=1: length(ttt_UP2)-1
        
        if ttt_UP2(iirr,1) == ttt_UP2(iirr+1,1)
        tttrrr = ttt_UP2(iirr,1);
        Activity_Function_Hippo_UP2_Total{epp_UP2,1}(tttrrr)=2;
        
        end
    end
    
    
        for iirr=1: length(ttt_UP2)-2
        
        if ttt_UP2(iirr,1) == ttt_UP2(iirr+1,1) == ttt_UP2(iirr+2,1)
        tttrrrr = ttt_UP2(iirr,1);
        Activity_Function_Hippo_UP2_Total{epp_UP2,1}(tttrrrr)=3;
        
        end
    end
    
end

%%
%%%  Cross-correlation of all of UP2 States in PFC and Hippocampus looking over 20 ms after and before them %%%%


for ii = 1 : length(UP2)
    
yfilt_PFC_UP2_win{ii,1} = filter(gaussFilter,1,(Activity_Function_PFC_UP2_Total{ii,1}));
yfilt_notshifted_PFC_UP2_win{ii,1} = conv ((Activity_Function_PFC_UP2_Total{ii,1}), gaussFilter, 'same');

yfilt_Hippo_UP2_win{ii,1} = filter(gaussFilter,1,(Activity_Function_Hippo_UP2_Total{ii,1}));
yfilt_notshifted_Hippo_UP2_win{ii,1} = conv ((Activity_Function_Hippo_UP2_Total{ii,1}), gaussFilter, 'same');

end



for iii = 1 :length(UP2)

     [acor_UP2_win, lag_UP2_win] = xcorr(yfilt_notshifted_Hippo_UP2_win{iii,1}, yfilt_notshifted_PFC_UP2_win{iii,1},200, 'unbiased');
     [Delay_UP2_win, Locs_UP2_win] = findpeaks(acor_UP2_win, lag_UP2_win, 'SortStr','descend');
     if ~ length(Locs_UP2_win) == 0
        Location_UP2_win(iii,1) = Locs_UP2_win(1);
     else 
        Location_UP2_win(iii,1) = 0;
     end

end

figure;
histbins_UP2 = [-200 : 20 :200]; 
histogram (Location_UP2_win, histbins_UP2);

title(['Distribution of Lags between Hippocampus and PFC for all UP2s of sleep3 with window of ' num2str(UP_Window)], 'fontsize' ,8);
xlabel('Lags Duration(ms)','fontsize', 14);
ylabel('Number of Lags','fontsize', 14);

saveas(gcf, ['DistributionofLagsBetweenHippocampusandPFCAllUP2swithwindowof ' num2str(UP_Window) '.fig']);

