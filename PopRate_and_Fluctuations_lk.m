 
% PopRate_and_Fluctuations.m


clear all;

% Select which sleep epoch to analyze (1 = sleep1; 2 = sleep2; 3 = sleep3)
Epoch_Select = 3;

% UPDATE 03/03/2019 by Soroush:
% Select which method to use for defining control period:
% 1 = original method using the first 300 seconds of the entire sleep epoch
% 2 = updated method that uses the mean + Th*std of the awake periods
% (determined by the video tracker data)
% 3 = uses the 97.5th percentile of the awake periods
% 4 = uses 2*mean of the awake periods
% 5 = uses mean + Th*std of the entire signal (the entire recording period)
% 6 = uses the Th percentile of the entire signal (entire recording period)
Control_Select = 6;

% Load data:
%**************************************
load('ts.mat')

switch Epoch_Select
    case 1
        load('sleep1_HMM_data.mat')
        sleep = 1;
        % % % % %
        ts_sleep = e.epochs.sleep1;
                a = (e.epochs.rest < ts_sleep(2) & e.epochs.rest > ts_sleep(1));
                ts_rest(:,1) = e.epochs.rest(a(:,1));
                ts_rest(:,2) = e.epochs.rest(a(:,1),2);
                ts_wake(1,1) = ts_sleep(1);
                ts_wake(1,2) = ts_rest(1,1);
                for i = 2:length(ts_rest)
                    ts_wake(i,1) = ts_rest(i-1,2);
                    ts_wake(i,2) = ts_rest(i,1);
                end
                ts_wake(length(ts_rest)+1,1) = ts_rest(end,2);
                ts_wake(end,2) = ts_sleep(2);
                ts_wake = (ts_wake - ts_sleep(1))/10;
                dur = ts_wake(:,2) - ts_wake(:,1);
                ts_wake(dur < 3000,:) = [];
        
        ts_start = 0;
        ts_end = (ts_sleep(2) - ts_sleep(1))/10;
        % % % % %
    case 2
        load('sleep2_HMM_data.mat')
        sleep = 2;
        
                ts_sleep = e.epochs.sleep2;
                a = (e.epochs.rest < ts_sleep(2) & e.epochs.rest>ts_sleep(1));
                ts_rest(:,1) = e.epochs.rest(a(:,1));
                ts_rest(:,2) = e.epochs.rest(a(:,1),2);
                ts_wake(1,1) = ts_sleep(1);
                ts_wake(1,2) = ts_rest(1,1);
                for i = 2:length(ts_rest)
                    ts_wake(i,1) = ts_rest(i-1,2);
                    ts_wake(i,2) = ts_rest(i,1);
                end
                ts_wake(length(ts_rest)+1,1) = ts_rest(end,2);
                ts_wake(end,2) = ts_sleep(2);
                ts_wake = (ts_wake - ts_sleep(1))/10;
                dur = ts_wake(:,2) - ts_wake(:,1);
                ts_wake(dur < 3000,:) = [];
        
        ts_start = 0;
        ts_end = (ts_sleep(2) - ts_sleep(1))/10;
        % % % % %
    case 3
        load('sleep3_HMM_data.mat')
        sleep = 3;

        if exist('e')
            ts_sleep = e.epochs.sleep3;
        else
            ts_sleep = epochs.sleep3;
        end
        a = (e.epochs.rest < ts_sleep(2) & e.epochs.rest>ts_sleep(1));
        ts_rest(:,1) = e.epochs.rest(a(:,1));
        ts_rest(:,2) = e.epochs.rest(a(:,1),2);
        ts_wake(1,1) = ts_sleep(1);
        ts_wake(1,2) = ts_rest(1,1);
        for i = 2:length(ts_rest)
            ts_wake(i,1) = ts_rest(i-1,2);
            ts_wake(i,2) = ts_rest(i,1);
        end
        ts_wake(length(ts_rest)+1,1) = ts_rest(end,2);
        ts_wake(end,2) = ts_sleep(2);
        ts_wake = (ts_wake - ts_sleep(1))/10;
        dur = ts_wake(:,2) - ts_wake(:,1);
        ts_wake(dur < 3000,:) = [];
        
        ts_start = 0;
        ts_end = (ts_sleep(2) - ts_sleep(1))/10;
        
    otherwise
        error('Wrong Value for Epoch_select')
end

%  load('sleep3_HMM_data.mat')
% filename = 'data0round.txt';
% [ids times]=textread(filename,'%f %f');


VeryEnd=max(times);
Start=0;
End=VeryEnd;


% spike times and Ids:
i=find(times<End & times>Start);
tt=times(i);
n=ids(i);
tt=tt-Start;

numNeurons=max(ids);

% Population Firing Rate:
%**************************************

bin=200; % chose bin size for population firing rate

numWind=floor(VeryEnd/bin);
center=zeros(1,numWind);
numSp=zeros(1,numWind);

display('calculating rate function...')
h = waitbar(0,'calculating rate function...');

for w=1:numWind
    Wstart=(w-1)*bin;
    Wend=Wstart+bin;
    center(w)=(Wend-Wstart)/2+Wstart;
    Sp=find(tt<Wend & tt>Wstart); %%% Finding how many neurons fire during this bin
    numSp(w)=1000*length(Sp)/(numNeurons*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
    
    if mod(w,50)==0
        waitbar(w / numWind)
    end
end

close(h)

% Fluctuations:
%******************************************

wind=3000; % chose bin size for population firing rate
Nwind=floor(VeryEnd/wind);
fluct=zeros(1,Nwind);
cent=zeros(1,Nwind);

display('calculating fluctuation...')



for w=1:Nwind
    Wstart=(w-1)*wind;
    win_start(w) = Wstart;
    Wend=Wstart+wind;
    win_end(w) = Wend;
    cent(w)=(Wend-Wstart)/2+Wstart;
    indx=find(center<Wend & center>Wstart);
    temp=numSp(indx);
    if temp>0
        fluct(w)=std(temp)/mean(temp);
    else
        fluct(w)=fluct(w-1);
    end
end

%filter:

filterlength=5;
myfilter=ones(1,filterlength)/filterlength;
filtfluct=filtfilt(myfilter,1,fluct);

% Compare fluctuations with control period:
%******************************************

switch Control_Select
    case 1
        con = 'org';
        % define control period:
        T=300000; % first T ms of the recording session
        FirstWindows=floor(T/wind);
        Th=2; % Define threshold for fluctuations
        Level=Th*mean(filtfluct(1:FirstWindows));
        % Level = 0.28;
        
    case 2
        con = 'std';
        conwin = cell(length(ts_wake),1);
        confiltfluct = cell(length(ts_wake),1);
        for i = 1:length(ts_wake);
            conidx = find(win_start <= ts_wake(i,2) & win_start >= ts_wake (i,1));
            conidx(end) = [];
            conwin{i} = conidx;
            confiltfluct{i} = filtfluct(conidx);
        end
        convec = cell2mat(confiltfluct');
         Th = 3; %%%% It was 3, I am just trying to find the result for
                    % other number as well!
% % % % % %         Th = 1;
        
        Level = mean(convec) + Th*std(convec);
        
    case 3
        con = 'prct';
        conwin = cell(length(ts_wake),1);
        confiltfluct = cell(length(ts_wake),1);
        for i = 1:length(ts_wake);
            conidx = find(win_start <= ts_wake(i,2) & win_start >= ts_wake (i,1));
            conidx(end) = [];
            conwin{i} = conidx;
            confiltfluct{i} = filtfluct(conidx);
        end
        convec = cell2mat(confiltfluct');
        Th = 97.5;
        Level = prctile(convec, Th);
        
    case 4
        con = 'mean';
        conwin = cell(length(ts_wake),1);
        confiltfluct = cell(length(ts_wake),1);
        for i = 1:length(ts_wake);
            conidx = find(win_start <= ts_wake(i,2) & win_start >= ts_wake (i,1));
            conidx(end) = [];
            conwin{i} = conidx;
            confiltfluct{i} = filtfluct(conidx);
        end
        convec = cell2mat(confiltfluct');
        Th = 3;
        Level = Th*mean(convec);
        
    case 5
        con = 'all';
        Th = 1;
        Level = mean(filtfluct) + Th*std(filtfluct);
        
    case 6
        con = 'allprct';
        Th = 60;
        Level = prctile(filtfluct, Th);
        
    otherwise
        error('Wrong Value for Control_Select')
end


% Finding UP and Down epochs:
%******************************************

ep=find(filtfluct>Level);
t=zeros(1,length(filtfluct));
t(ep)=1;

% edit : use diff function instead of circshift to find times
% of threshold crossing
b = diff([0 t 0]);
idx = find(b == 1);
IN = cent(idx);

idx2 = find(b == -1);
OUT = cent(idx2-1);

% edit : remove epochs that are less than 9000 seconds in
% length
idx3 = find(OUT-IN < 9000);
OUT(idx3) = [];
IN(idx3) = [];


% a=t;
% b=circshift(a,[0 -1]);
% c=find(a~=b);
% trans=c(1:end-1);
% transitions=cent(c); % times of threshold crossing

% t=[];
% IN=[];
% OUT=[];

% for u=transitions(1:end-1)
%    t=find(cent==u);
%    if filtfluct(t-2)<filtfluct(t+2)
%        IN=[IN u];
%    elseif filtfluct(t-2)>filtfluct(t+2)
%        OUT=[OUT u];
%    end
% end

%save transition times:
if length(IN)==length(OUT)
    epochs=[IN' OUT'];
elseif length(IN)>length(OUT)
    OUT=[OUT VeryEnd];
    epochs=[IN' OUT'];
end

% Save the detected UP/DOWN state periods:
%*********************************************
filename = ['UpDownPeriodsSleep',num2str(sleep),'_',con,'_Th',num2str(Th),'.mat'];
% % % % % save(filename,'epochs','Level','Th')

%       save UpDownPeriods epochs Level
%*********************************************


% Plot figure
%*****************************************************

set(figure,'Position',[250 100 760 760],'Color','w');
orient tall;

subplot(3,1,1)

plot(center/1000,numSp,'k')
title('Population firing rate','FontSize',14)
ylabel('spikes/s','FontSize',14)
set(gca,'FontSize',12)
box off
text(.8,.9,['Window=' num2str(bin) 'ms'],'Units','Normalized','FontSize',11)

subplot(3,1,2)

plot(cent/1000,fluct,'k')
title('Population firing rate fluctuations','FontSize',14)
ylabel('SD/mean','FontSize',14)
set(gca,'FontSize',12)
box off
text(.8,.9,['Window=' num2str(wind) 'ms'],'Units','Normalized','FontSize',11)


subplot(3,1,3)

rectangle('Position',[cent(1)/1000, 0, (cent(end)-cent(1))/1000, Level],'FaceColor',[.8 .8 .8],'EdgeColor','none')
hold on
plot([cent(1)/1000 cent(end)/1000],[Level Level],'Color',[.6 .6 .6])
text(cent(end)/1000,Level*1.1,'Threshold','FontSize',11)
plot(cent/1000,filtfluct,'k')

MAX=max(filtfluct);

xlabel('Time (s)','FontSize',14)
ylabel('SD/mean','FontSize',14)
set(gca,'FontSize',12,'YLim',[0 MAX*1.15])
box off
title('fluctuations (filtered)','FontSize',14)

for e=1:size(epochs,1)
    plot(epochs(e,:)/1000,[MAX*1.05 MAX*1.05],'b','LineWidth',4)
end

if Control_Select == 1
    plot(cent(1:FirstWindows)/1000,filtfluct(1:FirstWindows),'Color',[.97 .97 .97])
    text(cent(floor(FirstWindows/3))/1000,filtfluct(floor(FirstWindows/3))*.6,'Control','Color',[.97 .97 .97])
elseif Control_Select == 5
    plot([ts_start ts_end/1000],[MAX*1.1 MAX*1.1],'r','LineWidth',4)
elseif Control_Select == 6
    plot([ts_start ts_end/1000],[MAX*1.1 MAX*1.1],'r','LineWidth',4)
else
    for e=1:size(ts_wake,1)
        plot(ts_wake(e,:)/1000,[MAX*1.1 MAX*1.1],'r','LineWidth',4)
    end
end


% % % %  saveas(gcf, ['UpDownPeriodsSleep',num2str(sleep),'_',con,'_Th',num2str(Th),'.fig']);
figure;
motionless(1,:) = [ts_start,ts_wake(1,1)];
for m = 1 : length(ts_wake) - 1
motionless(m+1,:) = [ts_wake(m,2),ts_wake(m+1,1)];
end
motionless(length(motionless)+1,:) = [ts_wake(end,2), ts_end];

for e=1:size(motionless,1)
a = rectangle('Position',[motionless(e,1)/1000, 0, (motionless(e,2) - motionless(e,1))/1000 ,5 ],'FaceColor',[.8 .8 .8],'EdgeColor','none', 'Visible', 'on')
end
hold on
plot([cent(1)/1000 cent(end)/1000],[Level Level],'Color',[.6 .6 .6])
text(cent(end)/1000,Level,'Threshold','FontSize',8)
ac = plot(cent/1000,filtfluct,'k')
text(.8,1.05,['Window=' num2str(wind*5) 'ms'],'Units','Normalized','FontSize',11)
MAX=max(filtfluct);
xlabel('Time (s)','FontSize',14)
ylabel('SD/mean','FontSize',14)
set(gca,'FontSize',12,'YLim',[0 MAX*1.15])
box off
title('fluctuations (filtered)','FontSize',14)
for e=1:size(epochs,1)
ad = plot(epochs(e,:)/1000,[MAX*1.05 MAX*1.05],'b','LineWidth',4);
end
if Control_Select == 1
plot(cent(1:FirstWindows)/1000,filtfluct(1:FirstWindows),'Color',[.97 .97 .97])
text(cent(floor(FirstWindows/3))/1000,filtfluct(floor(FirstWindows/3))*.6,'Control','Color',[.97 .97 .97])
elseif Control_Select == 5
plot([ts_start ts_end/1000],[MAX*1.1 MAX*1.1],'r','LineWidth',4)
end
legend ([ad], 'Epochs')
