% function [] = UP_replay_hist_sleep_norm(UPdir, TMdir)
clear all;

bin=10; % chose bin size for population firing rate
Limited_Lag = 200 / bin ;

UP_Window = 20; UP1_Window = 20 ; UP2_Window = 20 ;
Dataset = '8482_15p';
Template = 8;
Compression = 7;
filename = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC/sleep3_UP1_UP2_epochs_run_extra3.mat', Dataset);
filenamee = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC/sleep3_UP_epochs_run505.mat', Dataset);
% load detected UP states

epochs_dir_epochs = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC', Dataset);
tfiles_dir_Cortex = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\PFC\\tfiles', Dataset);

cd(epochs_dir_epochs);
load('ts.mat')
sleep = e.epochs.sleep3;
savefile = 'sleep3_HMM_data.mat';
UpDownFile = 'UpDownPeriodsSleep3_allprct_Th84';
load(UpDownFile);

cd(tfiles_dir_Cortex)
tfile = FindFiles('*.t');
S= LoadSpikes(tfile);
Q = MakeQfromS(S,10);

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

load(filename, 'UP1_full', 'UP2_full' , 'UP1' , 'UP2')
UP1_full = UP1_full/10000;
UP2_full = UP2_full/10000;

load(filenamee, 'UP_epochs_full_ts' , 'UP_epochs')
UP_epochs_full = UP_epochs_full_ts/10000;

% load template matching results
TMdir = sprintf ('\\\\arcturus.uleth.ca\\workspace\\leanna\\TMOut\\TM-lk_v1_cluster\\TM-%s_Run001-%dX-tmpl%d\\all-variables.mat', Dataset, Compression, Template);


%%%outputDir2 = [TMdir,'\HMM_UP_replay_analysis\sleep3_norm_restrict_run1\all-variables.mat'];
load(TMdir, 'tmpl','QS_binsize')

% make folders for saving figures
outputDir = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC', Dataset);
outputDir3 = sprintf('E:/HMM - UP&Down/Soroush/Data/%s/PFC/Template%dCompression%d', Dataset, Template, Compression);

if ~exist(outputDir3,'dir')
    mkdir(outputDir3)
end

% 1/2 the template length
tmpl_length = size(Data(tmpl.Q),1)/2*QS_binsize/1000;
% restrict UP states to those that are longer than half the length of the
% template
ind1 = UP1_full(:,2)-UP1_full(:,1) >= tmpl_length;
UP1_rest = UP1_full(ind1,:);
ind2 = UP2_full(:,2)-UP2_full(:,1) >= tmpl_length;
UP2_rest = UP2_full(ind2,:);

cd(epochs_dir_epochs);

for i = 4:5
    
  % find peaks in template matching results
    Ct = tmpl.ZCCS(:,2);
    [p,loc] = findpeaks(Ct,'Minpeakheight',i,'Minpeakdistance',(round((size(Data(tmpl.Q),1))/2)));
    % to prevent multiple peaks from being detected in a cluster of Ct
    % points that pass threshold i, use a minpeakdistance that is half the
    % bin size of the template matrix (originally, I tried using the full
    % bin size of the template matrix, but this got rid of too many peaks)
    if isempty(p)
        break
    end
    t = tmpl.ZCCS(:,1);
    % shift Ct results by 1/2 the template bin width
    t = t + (size(Data(tmpl.Q),1)/2 * QS_binsize/1000);
    % create new variable with timestamps of peak locations
    pks = t(loc);
    
    %     figure; hold on
    %     plot(t,Ct,'k'); widefig(gca)
    %     plot(t(loc),Ct(loc),'r*')
    %     lim = get(gca,'xlim');
    %     line([lim(1) lim(2)],[i i],'Color','g')
    
    % find index of peaks in the template matching results that occur during individual
    % UP states
    idx1 = cell(length(UP1_rest),1);
    for n = 1:length(UP1_rest);
        idx1{n,1} = find(pks > UP1_rest(n,1) & pks < UP1_rest(n,2));
    end
    
    idx2 = cell(length(UP2_rest),1);
    for n = 1:length(UP2_rest);
        idx2{n,1} = find(pks > UP2_rest(n,1) & pks < UP2_rest(n,2));
    end
    
    % get timestamp of peaks within UP states
    UP1_pks = cell(length(idx1),1);
    for n = 1:length(idx1);
        UP1_pks{n,1} = pks(idx1{n,1});
    end
    
    UP2_pks = cell(length(idx2),1);
    for n = 1:length(idx2);
        UP2_pks{n,1} = pks(idx2{n,1});
    end
    
    % calculate where peaks occur within UP states (normalize the UP state
    % length to 1 ms)
    UP1_pks_pos = cell(length(UP1_pks),1);
    for n = 1:length(UP1_pks);
        UP1_pks_pos{n,1} = (UP1_pks{n,1}-UP1_rest(n,1))/(UP1_rest(n,2)-UP1_rest(n,1));
    end
    
    UP2_pks_pos = cell(length(UP2_pks),1);
    for n = 1:length(UP2_pks);
        UP2_pks_pos{n,1} = (UP2_pks{n,1}-UP2_rest(n,1))/(UP2_rest(n,2)-UP2_rest(n,1));
    end
    
     %%
    %%% The Part for onset and offset %%%%
    
    %%% For UP1 %%%
    
    UP1_pks_pos_nonnorm = cell(length(UP1_pks),1);
    UP1_pks_pos_zeros = cell(length(UP1_pks),1) ;
    
    for n = 1:length(UP1_pks);
        UP1_rest_dur(n) = UP1_rest(n,2)-UP1_rest(n,1);
        UP1_pks_pos_nonnorm{n,1} = (UP1_pks{n,1}-UP1_rest(n,1));
        UP1_pks_pos_zeros{n,1} = UP1_pks_pos_nonnorm{n,1} ;
    end
    
    [UP1_rest_dur_sorted index_sorted] = sort(UP1_rest_dur, 'descend');
    UP1_rest_dur_sorted = UP1_rest_dur_sorted' ;
    
    for iii = 1 : length(UP1_pks_pos_zeros)
        if isempty(UP1_pks_pos_zeros{iii,1})
            UP1_pks_pos_zeros{iii,1} = 0;
        end
    end
    
    %%%%UP1_pks_pos(cellfun('isnan',UP1_pks_pos))={0} ;
    UP1_pks_pos_withzeros = cell2mat(UP1_pks_pos_zeros);
    UP1_pks_pos_withzeros = UP1_pks_pos_withzeros(index_sorted);
    UP1_pks_pos_withzeros (UP1_pks_pos_withzeros == 0 ) = nan;
    UP1_pks_pos_nonnorm_sorted = UP1_pks_pos_nonnorm(index_sorted);
    %%%UP1_pks_pos(isnan(UP1_pks_pos)) = 0;
    
    x_UP1_dur = zeros(length(UP1_pks), round(UP1_rest_dur_sorted(n)*100) + 1 );
    % %     for n = 1 :  length(UP1_pks)
    % %     x_UP1_dur(n,1:floor(UP1_rest_dur_sorted*100) + 1) = [0:0.01:UP1_rest_dur_sorted];
    % %     y_UP1_dur(1,n) = n ;
    % %     end
    
    %f(x_UP1_dur, y_UP1_dur)
    figure ;
    plot(UP1_rest_dur_sorted(1:end,1), 1:length(UP1_pks));
    hold on;
    for n = length(UP1_pks) : -1 : 1
        Y_UP1_dur = [] ;
        X_UP1_dur_only = [0:0.00001: UP1_rest_dur_sorted(n,1)];
        Y_UP1_dur(round(X_UP1_dur_only*100000) +1 ) = n;
        plot(X_UP1_dur_only, Y_UP1_dur);
        hold on;
    end
% % % % % %     hold on;
% % % % % %     plot(UP1_pks_pos_withzeros(1:end,1), 1:length(UP1_pks), '*');
% % % % % %     hold on;
    for  OO= 1 : length(UP1_pks);
        
        if ~isempty(UP1_pks_pos_nonnorm_sorted{OO,1});
            
            plot([UP1_pks_pos_nonnorm_sorted{OO,1}(1:end,1)], OO,'b*' );
            
            hold on;
            
        else
            
        end;
        
    end;
    
    title (['UP1 Replay Peak Location - Z score' num2str(i)],'fontsize', 12);
    xlabel('Duration of UP1 (s)','fontsize', 12);
    ylabel('Number of UP1','fontsize' ,12);
    ylim([0 length(UP1_rest_dur)]);
    
    saveas(gcf, ['UP1 Durations Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);
    
    %%%% For UP2 %%%%
    
    UP2_pks_pos_nonnorm = cell(length(UP2_pks),1);
    UP2_pks_pos_zeros = cell(length(UP2_pks),1) ;
    
    for n = 1:length(UP2_pks);
        UP2_rest_dur(n) = UP2_rest(n,2)-UP2_rest(n,1);
        UP2_pks_pos_nonnorm{n,1} = (UP2_pks{n,1}-UP2_rest(n,1));
        UP2_pks_pos_zeros{n,1} = UP2_pks_pos_nonnorm{n,1} ;
    end
    
    [UP2_rest_dur_sorted index_sorted_UP2] = sort(UP2_rest_dur, 'descend');
    UP2_rest_dur_sorted = UP2_rest_dur_sorted' ;
    
    for iii = 1 : length(UP2_pks_pos_zeros)
        if isempty(UP2_pks_pos_zeros{iii,1})
            UP2_pks_pos_zeros{iii,1} = 0;
        end
    end
    
    %%%%UP2_pks_pos(cellfun('isnan',UP2_pks_pos))={0} ;
    UP2_pks_pos_withzeros = cell2mat(UP2_pks_pos_zeros);
    UP2_pks_pos_withzeros = UP2_pks_pos_withzeros(index_sorted_UP2);
    UP2_pks_pos_withzeros (UP2_pks_pos_withzeros == 0 ) = nan;
    UP2_pks_pos_nonnorm_sorted = UP2_pks_pos_nonnorm(index_sorted_UP2);

    %%%UP2_pks_pos(isnan(UP2_pks_pos)) = 0;
    
    x_UP2_dur = zeros(length(UP2_pks), round(UP2_rest_dur_sorted(n)*100) + 1 );
    % %     for n = 1 :  length(UP2_pks)
    % %     x_UP2_dur(n,1:floor(UP2_rest_dur_sorted*100) + 1) = [0:0.01:UP2_rest_dur_sorted];
    % %     y_UP2_dur(1,n) = n ;
    % %     end
    
    %f(x_UP2_dur, y_UP2_dur)
    
    figure ;
    plot(UP2_rest_dur_sorted(1:end,1), 1:length(UP2_pks));
    hold on;
    for n = length(UP2_pks) : -1 : 1
        Y_UP2_dur = [] ;
        X_UP2_dur_only = [0:0.00001: UP2_rest_dur_sorted(n,1)];
        Y_UP2_dur(round(X_UP2_dur_only*100000) +1 ) = n;
        plot(X_UP2_dur_only, Y_UP2_dur);
        hold on;
    end
% % % % % % % % % % %     hold on;
% % % % % % % % % % %     plot(UP2_pks_pos_withzeros(1:end,1), 1:length(UP2_pks), '*');
% % % % % % % % % % %     hold on;
    
    for  OO= 1 : length(UP2_pks);
        
        if ~isempty(UP2_pks_pos_nonnorm_sorted{OO,1});
            
            plot([UP2_pks_pos_nonnorm_sorted{OO,1}(1:end,1)], OO, 'b*');
            
            hold on;
            
        else
            
        end;
        
    end;
    
    title (['UP2 Replay Peak Location - Z score ' num2str(i)],'fontsize', 12);
    xlabel('Duration of UP2 (s)','fontsize', 12);
    ylabel('Number of UP2','fontsize' ,12);
    ylim([0 length(UP2_rest_dur)]);

    saveas(gcf, ['UP2 Durations Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);

    
  %%%%%%%%%%%%%%%% for UPs %%%%%%%

    
    idx_UP = cell(length(UP_epochs_full),1);
    for ii = 1:length(UP_epochs_full);
        idx_UP{ii,1} = find(pks > UP_epochs_full(ii,1) & pks < UP_epochs_full(ii,2));
    end
    
    
    % get timestamp of peaks within UP states
    UP_pks = cell(length(idx_UP),1);
    for ii = 1:length(idx_UP);
        UP_pks{ii,1} = pks(idx_UP{ii,1});
    end
    
    
    for ii = 1 : length (UP_epochs_full)
        
        UP_dur(ii) = UP_epochs_full(ii,2) - UP_epochs_full(ii,1);
        UP_pks_pos_nonnorm{ii,1} = (UP_pks{ii,1}-UP_epochs_full(ii,1));
        UP_pks_pos_zeros{ii,1} = UP_pks_pos_nonnorm{ii,1} ;
        
    end
    
    
    [UP_rest_dur_sorted index_sorted_UP] = sort(UP_dur, 'descend');
    UP_rest_dur_sorted = UP_rest_dur_sorted' ;
    
    for iii = 1 : length(UP_pks_pos_zeros)
        if isempty(UP_pks_pos_zeros{iii,1})
            UP_pks_pos_zeros{iii,1} = 0;
        end
    end
    
    %%%%UP_pks_pos(cellfun('isnan',UP_pks_pos))={0} ;
    UP_pks_pos_withzeros = cell2mat(UP_pks_pos_zeros);
    UP_pks_pos_withzeros = UP_pks_pos_withzeros(index_sorted_UP);
    UP_pks_pos_withzeros (UP_pks_pos_withzeros == 0 ) = nan;
    UP_pks_pos_nonnorm_sorted = UP_pks_pos_nonnorm(index_sorted_UP);
    %%%UP1_pks_pos(isnan(UP1_pks_pos)) = 0;
    
    x_UP_dur = zeros(length(UP_pks), round(UP_rest_dur_sorted(n)*100) + 1 );

    figure ;
% % %     plot(UP_rest_dur_sorted(1:end,1), 1:length(UP_pks) , 'w');
    hold on;
    for n = length(UP_pks) : -1 : 1
        Y_UP_dur = [] ;
        X_UP_dur_only = [0:0.00001: UP_rest_dur_sorted(n)];
        Y_UP_dur(round(X_UP_dur_only*100000) +1 ) = n;
        plot(X_UP_dur_only, Y_UP_dur);
        hold on;
    end
    
        for  OO= 1 : length(UP_pks);
        
        if ~isempty(UP_pks_pos_nonnorm_sorted{OO,1});
            
            plot([UP_pks_pos_nonnorm_sorted{OO,1}(1:end,1)], OO, 'b*');
            
            hold on;
            
        else
            
        end;
        
    end;
    title (['UP Replay Peak Location - Z score' num2str(i)],'fontsize', 14);
    xlabel('Duration of UP (s)','fontsize', 12);
    ylabel('Number of UP','fontsize' ,12);
    ylim([0 length(UP_dur)]);
    
    saveas(gcf, ['UP Durations Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);

    
    %%%%%%%%%%%% UPs, UP1, and UP2 combined figures %%%%%%
        %%
    cd(epochs_dir_epochs);
    
    idx_UP = cell(length(UP_epochs_full),1);
    for ii = 1:length(UP_epochs_full);
        idx_UP{ii,1} = find(pks > UP_epochs_full(ii,1) & pks < UP_epochs_full(ii,2));
    end
    
    
    % get timestamp of peaks within UP states
    UP_pks = cell(length(idx_UP),1);
    for ii = 1:length(idx_UP);
        UP_pks{ii,1} = pks(idx_UP{ii,1});
    end
    
    
    for ii = 1 : length (UP_epochs_full)
        
        UP_dur(ii) = UP_epochs_full(ii,2) - UP_epochs_full(ii,1);
        UP_pks_pos_nonnorm{ii,1} = (UP_pks{ii,1}-UP_epochs_full(ii,1));
        UP_pks_pos_zeros{ii,1} = UP_pks_pos_nonnorm{ii,1} ;
        
    end
    
    
    [UP_rest_dur_sorted index_sorted_UP] = sort(UP_dur, 'descend');
    UP_rest_dur_sorted = UP_rest_dur_sorted' ;
    
    for iii = 1 : length(UP_pks_pos_zeros)
        if isempty(UP_pks_pos_zeros{iii,1})
            UP_pks_pos_zeros{iii,1} = 0;
        end
    end
    
    %%%%UP_pks_pos(cellfun('isnan',UP_pks_pos))={0} ;
    UP_pks_pos_withzeros = cell2mat(UP_pks_pos_zeros);
    UP_pks_pos_withzeros = UP_pks_pos_withzeros(index_sorted_UP);
    UP_pks_pos_withzeros (UP_pks_pos_withzeros == 0 ) = nan;
    UP_pks_pos_nonnorm_sorted = UP_pks_pos_nonnorm(index_sorted_UP);
    %%%UP1_pks_pos(isnan(UP1_pks_pos)) = 0;
    
    x_UP_dur = zeros(length(UP_pks), round(UP_rest_dur_sorted(n)*100) + 1 );
    % %     for n = 1 :  length(UP1_pks)
    % %     x_UP1_dur(n,1:floor(UP1_rest_dur_sorted*100) + 1) = [0:0.01:UP1_rest_dur_sorted];
    % %     y_UP1_dur(1,n) = n ;
    % %     end
    
    %f(x_UP1_dur, y_UP1_dur)
    figure ;
% % %     plot(UP_rest_dur_sorted(1:end,1), 1:length(UP_pks) , 'w');
    hold on;
    for n = length(UP_pks) : -1 : 1
        Y_UP_dur = [] ;
        X_UP_dur_only = [0:0.00001: UP_rest_dur_sorted(n)];
        Y_UP_dur(round(X_UP_dur_only*100000) +1 ) = n;
        plot(X_UP_dur_only, Y_UP_dur , 'w');
        hold on;
    end
    
    
    UP1_of_UP = cell(length(UP_epochs_full),1);
    UP2_of_UP = cell(length(UP_epochs_full),1);
    
    for iii = 1 : length (UP_epochs_full)
        aa = 1;
        bb = 1; 
        for ii = 1 :length(UP1_full)
            if  UP1_full(ii,1) >= UP_epochs_full(iii,1) & UP1_full(ii,2) <= UP_epochs_full(iii,2)
                UP1_of_UP{iii,1}(aa,1) = UP1_full(ii,1) - UP_epochs_full(iii,1);
                UP1_of_UP{iii,1}(aa,2) =  UP1_full(ii,2) - UP_epochs_full(iii,1);
                aa = aa + 1;
            else
                
            end
            
        end
        
        
        for ii = 1 :length(UP2_full)
            if  UP2_full(ii,1) >= UP_epochs_full(iii,1) & UP2_full(ii,2) <= UP_epochs_full(iii,2)
                UP2_of_UP{iii,1}(bb,1) = UP2_full(ii,1) - UP_epochs_full(iii,1);
                UP2_of_UP{iii,1}(bb,2) =  UP2_full(ii,2) - UP_epochs_full(iii,1);
            bb = bb+1;
            
            else
                
            end
            
        end
        
    end
    
    
    UP1_of_UP = UP1_of_UP(index_sorted_UP,1);
    
    
    for n = length(UP_pks) : -1 : 1
        
        if ~isempty(UP1_of_UP{n,1});
         if  length(UP1_of_UP{n,1}(:,1))>1
            for aa= 1 : 2
                Y_UP_UP1_dur = [] ;
        X_UP_UP1_dur_only = [UP1_of_UP{n,1}(aa,1):0.00001: UP1_of_UP{n,1}(aa,2)];
        if length(X_UP_UP1_dur_only) > 2
            OO = length(X_UP_UP1_dur_only);
            Y_UP_UP1_dur( 1 : OO ) = n;
            plot(X_UP_UP1_dur_only, Y_UP_UP1_dur, 'b');
          
        end
            end
         else
                Y_UP_UP1_dur = [] ;
        X_UP_UP1_dur_only = [UP1_of_UP{n,1}(:,1):0.00001: UP1_of_UP{n,1}(:,2)];
        if length(X_UP_UP1_dur_only) > 2
            OO = length(X_UP_UP1_dur_only);
            Y_UP_UP1_dur( 1 : OO ) = n;
            plot(X_UP_UP1_dur_only, Y_UP_UP1_dur, 'b');
        end
         
        hold on;
        end
        end
    end
    
    
    UP2_of_UP = UP2_of_UP(index_sorted_UP,1);
    
    
    for n = length(UP_pks) : -1 : 1

        if ~isempty(UP2_of_UP{n,1});
            if  length(UP2_of_UP{n,1}(:,1))>1
           for bb = 1 : 2
               Y_UP_UP2_dur = [] ;
        X_UP_UP2_dur_only = [UP2_of_UP{n,1}(bb,1):0.00001: UP2_of_UP{n,1}(bb,2)];
        if length(X_UP_UP2_dur_only) > 2
            OOO = length(X_UP_UP2_dur_only);
            Y_UP_UP2_dur( 1 : OOO ) = n;
            plot(X_UP_UP2_dur_only, Y_UP_UP2_dur, 'r');
        end
           end
            else
               Y_UP_UP2_dur = [] ;
        X_UP_UP2_dur_only = [UP2_of_UP{n,1}(:,1):0.00001: UP2_of_UP{n,1}(:,2)];
        if length(X_UP_UP2_dur_only) > 2
            OOO = length(X_UP_UP2_dur_only);
            Y_UP_UP2_dur( 1 : OOO ) = n;
            plot(X_UP_UP2_dur_only, Y_UP_UP2_dur, 'r');
        end
        hold on;
            end
            
        end
    end
    
    
% % % % % % %     plot(UP_pks_pos_withzeros(1:end,1), 1:length(UP_pks), '*');


    for  OO= 1 : length(UP_pks);
        
        if ~isempty(UP_pks_pos_nonnorm_sorted{OO,1});
            
            plot([UP_pks_pos_nonnorm_sorted{OO,1}(1:end,1)], OO, 'w*');
            
            hold on;
            
        else
            
        end;
        
    end;
    title (['UP Replay Peak Location - Z score' num2str(i)],'fontsize', 14);
    xlabel('Duration of UP (s)','fontsize', 12);
    ylabel('Number of UP','fontsize' ,12);
    ylim([0 length(UP_dur)]);
    
    saveas(gcf, ['UP - UP1 - UP2 Durations Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);


    
    % % % % % % %
    % % % % % % %     %% Barh Figure %%%%%
    % % % % % % %
    % % % % % % %     UP1_of_UP = zeros(length(UP_epochs_full),2);
    % % % % % % %     UP2_of_UP = zeros(length(UP_epochs_full),2);
    % % % % % % %
    % % % % % % %     figure ;
    % % % % % % %
    % % % % % % %     for iii = 1 : length (UP_epochs_full)
    % % % % % % %
    % % % % % % %         for ii = 1 :length(UP1_full)
    % % % % % % %             if  UP1_full(ii,1) >= UP_epochs_full(iii,1) & UP1_full(ii,2) <= UP_epochs_full(iii,2)
    % % % % % % %                 UP1_of_UP(iii,1) = UP1_full(ii,1) - UP_epochs_full(iii,1);
    % % % % % % %                 UP1_of_UP(iii,2) =  UP1_full(ii,2) - UP_epochs_full(iii,1);
    % % % % % % %             else
    % % % % % % %
    % % % % % % %             end
    % % % % % % %
    % % % % % % %         end
    % % % % % % %
    % % % % % % %
    % % % % % % %         for ii = 1 :length(UP2_full)
    % % % % % % %             if  UP2_full(ii,1) >= UP_epochs_full(iii,1) & UP2_full(ii,2) <= UP_epochs_full(iii,2)
    % % % % % % %                 UP2_of_UP(iii,1) = UP2_full(ii,1) - UP_epochs_full(iii,1);
    % % % % % % %                 UP2_of_UP(iii,2) =  UP2_full(ii,2) - UP_epochs_full(iii,1);
    % % % % % % %             else
    % % % % % % %
    % % % % % % %             end
    % % % % % % %
    % % % % % % %         end
    % % % % % % %
    % % % % % % %     end
    % % % % % % %
    % % % % % % %
    % % % % % % %     UP1_of_UP(:,1) = UP1_of_UP(index_sorted_UP,1);
    % % % % % % %     UP1_of_UP(:,2) = UP1_of_UP(index_sorted_UP,2);
    % % % % % % %
    % % % % % % %     for n = length(UP_pks) : -1 : 1
    % % % % % % %         X_UP_UP1_dur_only(iii) = [UP1_of_UP(n,1): UP1_of_UP(n,2)];
    % % % % % % %
    % % % % % % %     UP2_of_UP(:,1) = UP2_of_UP(index_sorted_UP,1);
    % % % % % % %     UP2_of_UP(:,2) = UP2_of_UP(index_sorted_UP,2);
    % % % % % % %
    % % % % % % %         X_UP_UP2_dur_only(iii) = [UP2_of_UP(n,1); UP2_of_UP(n,2)];
    % % % % % % %
    % % % % % % %     end
    % % % % % % %
    % % % % % % %
    % % % % % % %
    % % % % % % %     Y_UP_UP1_UP2 = [X_UP_UP1_dur_only, X_UP_UP2_dur_only];
    % % % % % % %     X_UP_UP1_UP2 = [1:10];
    % % % % % % %
    % % % % % % %     barh(X_UP_UP1_UP2, Y_UP_UP1_UP2(1:10),'stacked');
    % % % % % % %
    % % % % % % %     hold on
    % % % % % % %
    % % % % % % %     plot(UP_pks_pos_withzeros(1:end,1), 1:length(UP_pks), '*');
    % % % % % % %     title (['UP durations ' num2str(i)],'fontsize', 14);
    % % % % % % %     xlabel('Duration of UPs','fontsize', 12);
    % % % % % % %     ylabel('Number of UPs','fontsize' ,12);
    % % % % % % %     ylim([0 length(UP_dur)]);
    % % % % % % %
    %%
    % create histogram of when replay peaks occur within the UP states
    % transform cell array data into matrix
    %%
    %%%%% Buidling Firing rate for all of UPs of Hippocampus and PFC with bining the Q-matrix %%%%%%

    
    center_UP = cell (length(UP_pks), 1);
    numSp_UP_PFC = cell (length(UP_pks),1);
    
    for ep_UP = 1 : length(UP_pks)
        
        
        Start = UP_epochs (ep_UP,1);
        VeryEnd = UP_epochs (ep_UP,2);
        
        End = VeryEnd - Start;
        
        
        i_UP=find(times<VeryEnd & times>Start);
        tt_UP=times(i_UP);
        n_UP =ids(i_UP);
        tt_UP=tt_UP-Start;
        
        numWind_UP=floor(End/bin);
        center_UP{ep_UP,1}=zeros(1,numWind_UP);
        numSp_UP_PFC{ep_UP,1}=zeros(1,numWind_UP);
        
        
        for w=1:numWind_UP
            Wstart=(w-1)*bin;
            Wend=Wstart+bin;
            center_UP{ep_UP,1}(w)=(Wend-Wstart)/2+Wstart;
            Sp_UP=find(tt_UP<Wend & tt_UP>Wstart); %%% Finding how many neurons fire during this bin
            numSp_UP_PFC{ep_UP,1}(w) =1000*length(Sp_UP)/(numNeurons*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
            
            if mod(w,5)==0
                
                
                
            end
        end
        
        
        
    end
    
    
    numSp_UP_PFC_withNaN = numSp_UP_PFC;
    
    [a b] = cellfun(@size,numSp_UP_PFC_withNaN(:,1));
    aa = max(b);
    
    for jj = 1 : length(numSp_UP_PFC_withNaN)
        
        numSp_UP_PFC_withNaN{jj,1}(1,end+1:aa) = nan;
        
    end
    
    
    for oo = 1 : aa - 50
        
        Mean_UP(oo) = nanmean(cellfun(@(c) c(1,oo), numSp_UP_PFC_withNaN(:,1)));
        
    end
    
    figure ; plot(1*bin/1000:bin/1000:length(Mean_UP)*bin/1000,Mean_UP)
    
    hold on;
    
    plot(UP_pks_pos_withzeros(1:end,1), 1, 'b*');
    title (['Mean of UP firing rate - Replay Peak Location - Z score ' num2str(i)],'fontsize', 10);
    xlabel('Duration of UP (s)','fontsize', 12);
    ylabel('Firing rate (Hz)','fontsize' ,12);

    saveas(gcf, ['Mean of UP firing rate - Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);

    
    % % %     ylim([0 length(UP_dur)]);
    
% % % % % %     indx_UP_reactivation = find (~isnan(UP_pks_pos_withzeros)) ;
% % % % % %     UP_With_reactivation = cell(length(indx_UP_reactivation),1);
% % % % % %     
% % % % % %     for jj = 1 : 5 :length(indx_UP_reactivation)
% % % % % %         
% % % % % %         UP_With_reactivation{jj} = numSp_UP_PFC{indx_UP_reactivation(jj,1),1}(1,:);
% % % % % %         
% % % % % %         
% % % % % %         Mean_UP_With_reactivation(jj) = mean (UP_With_reactivation{jj,1}(1,:));
% % % % % %         
% % % % % %         figure; plot(1*bin/1000:bin/1000:length(UP_With_reactivation{jj})*bin/1000 , UP_With_reactivation{jj});
% % % % % %         
% % % % % %         hold on;
% % % % % %         
% % % % % %         plot(UP_pks_pos_withzeros(indx_UP_reactivation(jj,1))/QS_binsize, 1,'*');
% % % % % %         
% % % % % %     end
% % % % % %     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%
    %%%%%%%%%%%%% Buidling Firing rate for all of UP1s of PFC with bining the Q-matrix %%%%%%
    
    
    center_UP1 = cell (length(UP1), 1);
    numSp_UP1_PFC = cell (length(UP1),1);
    
    for ep_UP1 = 1 : length(UP1)
        
        
        Start = UP1 (ep_UP1,1);
        VeryEnd = UP1 (ep_UP1,2);
        
        End = VeryEnd - Start;
        
        
        i_UP1=find(times<VeryEnd & times>Start);
        tt_UP1=times(i_UP1);
        n_UP1 =ids(i_UP1);
        tt_UP1=tt_UP1-Start;
        
        
        numWind_UP1=floor(End/bin);
        center_UP1{ep_UP1,1}=zeros(1,numWind_UP1);
        numSp_UP1_PFC{ep_UP1,1}=zeros(1,numWind_UP1);
        
        
        for w=1:numWind_UP1
            Wstart=(w-1)*bin;
            Wend=Wstart+bin;
            center_UP1{ep_UP1,1}(w)=(Wend-Wstart)/2+Wstart;
            Sp_UP1=find(tt_UP1<Wend & tt_UP1>Wstart); %%% Finding how many neurons fire during this bin
            numSp_UP1_PFC{ep_UP1,1}(w) =1000*length(Sp_UP1)/(numNeurons*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
            
            if mod(w,5)==0
                
                
                
            end
        end
        
    end
    
    numSp_UP1_PFC_withNaN = numSp_UP1_PFC ;
    [a_UP1 b_UP1] = cellfun(@size,numSp_UP1_PFC_withNaN(:,1));
    aa_UP1 = max(b_UP1);
    
    for jj = 1 : length(numSp_UP1_PFC_withNaN)
        
        numSp_UP1_PFC_withNaN{jj,1}(1,end+1:aa) = nan;
        
    end
    
    
    for oo = 1 : aa_UP1 - 50
        
        Mean_UP1(oo) = nanmean(cellfun(@(c) c(1,oo), numSp_UP1_PFC_withNaN(:,1)));
        
    end
    
    figure ; plot(1*bin/1000:bin/1000:length(Mean_UP1)*bin/1000,Mean_UP1)
    
    hold on;
    
    plot(UP1_pks_pos_withzeros(1:end,1), 1, 'b*');
    title (['Mean of UP1 firing rate - Replay Peak Location - Z score' num2str(i)],'fontsize', 10);
    xlabel('Duration of UP1 (s)','fontsize', 12);
    ylabel('Firing rate (Hz)','fontsize' ,12);
    
    saveas(gcf, ['Mean of UP1 firing rate - Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    %%%%%%%%%%%%% Buidling Firing rate for all of UP2s of Hippocampus and PFC with bining the Q-matrix %%%%%%
    
    
    center_UP2 = cell (length(UP2), 1);
    numSp_UP2_PFC = cell (length(UP2),1);
    
    for ep_UP2 = 1 : length(UP2)
        
        
        Start = UP2 (ep_UP2,1);
        VeryEnd = UP2 (ep_UP2,2);
        
        End = VeryEnd - Start;
        
        
        i_UP2=find(times<VeryEnd & times>Start);
        tt_UP2=times(i_UP2);
        n_UP2 =ids(i_UP2);
        tt_UP2=tt_UP2-Start;
        
        
        numWind_UP2=floor(End/bin);
        center_UP2{ep_UP2,1}=zeros(1,numWind_UP2);
        numSp_UP2_PFC{ep_UP2,1}=zeros(1,numWind_UP2);
        
        
        for w=1:numWind_UP2
            Wstart=(w-1)*bin;
            Wend=Wstart+bin;
            center_UP2{ep_UP2,1}(w)=(Wend-Wstart)/2+Wstart;
            Sp_UP2=find(tt_UP2<Wend & tt_UP2>Wstart); %%% Finding how many neurons fire during this bin
            numSp_UP2_PFC{ep_UP2,1}(w) =1000*length(Sp_UP2)/(numNeurons*bin);  %%% using the number of  firing (length(sp) find the average firing rate)
            
            if mod(w,5)==0
                
                
            end
        end
        
        
    end
    
    [a_UP2 b_UP2] = cellfun(@size,numSp_UP2_PFC(:,1));
    aa_UP2 = max(b_UP2);
    
    for jj = 1 : length(numSp_UP2_PFC)
        
        numSp_UP2_PFC{jj,1}(1,end+1:aa) = nan;
        
    end
    
    
    for oo = 1 : aa_UP2 - 50
        
        Mean_UP2(oo) = nanmean(cellfun(@(c) c(1,oo), numSp_UP2_PFC(:,1)));
        
    end
    
    figure ; plot(1*bin/1000:bin/1000:length(Mean_UP2)*bin/1000, Mean_UP2)
    
    hold on;
    
    plot(UP2_pks_pos_withzeros(1:end,1), 1, 'b*');
    title (['Mean of UP2 firing rate - Replay Peak Location - Z score' num2str(i)],'fontsize', 10);
    xlabel('Duration of UP2s','fontsize', 12);
    ylabel('Firing rate','fontsize' ,12);
    
    saveas(gcf, ['Mean of UP2 firing rate - Replay Peak Location - Z score' num2str(i) 'Sleep3.fig']);
   
    
    %%
    UP1_pks_pos_vector = cell2mat(UP1_pks_pos);
    UP2_pks_pos_vector = cell2mat(UP2_pks_pos);
    
    if isempty(UP1_pks_pos_vector) && isempty(UP2_pks_pos_vector)
        continue
    end
    
    % plot and save figures
    for j = [0.1 0.05 0.02 0.01]
        
        edges = 0:j:1;
        [pks_pos_count1,bin1] = histc(UP1_pks_pos_vector,edges);
        % normalize by total number of UP states
        num_UP1 = length(UP1_rest);
        pks_pos_count_norm1 = pks_pos_count1/num_UP1;
        
        figure;
        ax1 = subplot(1,2,1);
        hold on;
        o=bar(edges,pks_pos_count_norm1,'histc');
        set(o,'facecolor',[0.5 0.5 0.5])
        xlim([-0.1 1.1])
        title(['UP1 Distribution of Replay Peak Location Above Z = ' num2str(i)],'fontsize',8)
        xlabel('Replay Peak Location within Normalized UP States','fontsize',8);
        ylabel('Number of Replay Peaks/Number of UP States','fontsize',8);
        
        k = length(edges)-1;
        temp = [outputDir3,'\UP1_Replay_Peaks_Location_Z',num2str(i),'_bins',num2str(k),'.fig'];
        saveas(o,temp);
        
        edges = 0:j:1;
        [pks_pos_count2,bin2] = histc(UP2_pks_pos_vector,edges);
        % normalize by total number of UP states
        num_UP2 = length(UP2_rest);
        pks_pos_count_norm2 = pks_pos_count2/num_UP2;
        
        
        ax2 = subplot(1,2,2);
        o=bar(edges,pks_pos_count_norm2,'histc');
        set(o,'facecolor',[0.5 0.5 0.5])
        xlim([-0.1 1.1])
        title(['UP2 Distribution of Replay Peak Location Above Z = ' num2str(i)],'fontsize',8)
        xlabel('Replay Peak Location within Normalized UP States','fontsize',8);
        ylabel('Number of Replay Peaks/Number of UP States','fontsize',8);
        linkaxes([ax1,ax2],'xy');
        ylim([0 .15]);
        
        k = length(edges)-1;
        temp = [outputDir3,'\UP2_Replay_Peaks_Location_Z',num2str(i),'_bins',num2str(k),'.fig'];
        saveas(o,temp);
        
        temp_mat = [outputDir3,'\Replay_Peaks_Location_Z',num2str(i),'_bins',num2str(k),'.mat'];
        save(temp_mat,'p','loc','pks','idx1','idx2','UP1_pks','UP2_pks',...
            'UP1_pks_pos','UP2_pks_pos','UP1_pks_pos_vector',...
            'UP2_pks_pos_vector','j','edges','pks_pos_count1','num_UP1',...
            'pks_pos_count_norm1','pks_pos_count2','num_UP2',...
            'pks_pos_count_norm2','bin1','bin2')
        
    end
    
end

temp_mat = [outputDir3,'\other_variables.mat'];
save(temp_mat,'TMdir','filename','outputDir3','UP1_full',...
    'UP2_full','tmpl_length','UP1_rest','UP2_rest','ind1','ind2','Ct','t')

% end