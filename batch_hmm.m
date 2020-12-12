%% Updated by Soroush On 2020/12/11 %%

datasets = {'7165_11p', '7165_16p', '7165_31p', '7165_36p', '8202_07p', '8202_08p', '8202_09p', '8482_14p', '8482_15p', '8482_16p'};
in_dir = '/home/Soroush/tatsunonas2_1/soroush/HMM Result - Paper/New_Result/';
out_dir = '/home/Soroush/tatsunonas2_1/soroush/upstate/analysis/';
sv_dir = '/home/Soroush/huxley/workspace/leanna/TMOut/TM-lk_v1_cluster/';
sv_dir2 = '/home/Soroush/tatsunonas2_1/karim/upstate/TMOut/TM-lk_v1_cluster/';
csc_files = {'CSC1.dat', 'CSC1.dat', 'CSC1.dat', 'CSC2.dat', 'CSC08.ncs', 'CSC02.ncs', 'CSC01.ncs', 'CSC08.ncs', 'CSC07.ncs', 'CSC07.ncs'};
templates = [8, 2, 2, 5, 1, 4, 4, 2, 8, 8];
compressions = [7, 8, 8, 6, 8, 10, 9, 7, 6, 5];
z_thresh = 5:6;
sleep = 3;
plot_figs = false;

if sleep == 3
	load('/home/Soroush/tatsunonas2_1/soroush/upstate/kernel_threshold/thresholds.mat');
else
	load('/home/Soroush/tatsunonas2_1/soroush/upstate/kernel_threshold/thresholds_sleep1.mat');
	% load('/home/Soroush/tatsunonas2_1/soroush/upstate/kernel_threshold/7165_16p/thresholds.mat');
	% threshs = zeros(2, 1);
	% threshs(2) = thresholds(2);
end

freq_range = [0, 5; 5, 10; 10, 20; 20, 60; 60, 100; 100, 300];
freq_names = {'1-5Hz', '5-10Hz', '10-20Hz', '20-60Hz', '60-100hz', '100-300Hz'};

tau_all = zeros(length(datasets), 2);
tau_comp_all = zeros(length(datasets), 2);
tau_task_all = zeros(length(datasets), 1);
pow_all = zeros(length(freq_names), 2, length(datasets));
d_all = zeros(length(datasets), 1);
fr_all = zeros(length(datasets), 3);
fr_sem = zeros(length(datasets), 3);
dur_all = zeros(length(datasets), 3);
dur_sem = zeros(length(datasets), 3);
replay_count_z = zeros(length(datasets), 2, length(z_thresh));
start_state = zeros(length(datasets), 2);
start_state_replay = zeros(length(datasets), 2);
finish_state = zeros(length(datasets), 2);
finish_state_replay = zeros(length(datasets), 2);
tr_lfp_all = cell(length(datasets), 6);
up_types_all = zeros(length(datasets), 5);
up_bin_num_all = zeros(length(datasets), 6, 2);
up_bin_dur_all = zeros(length(datasets), 6, 2);

for i = 1:length(datasets)
    tic;
	curr_in = fullfile(in_dir, datasets{i});
	curr_out = fullfile(out_dir, datasets{i}, sprintf('sleep%d', sleep));
    if ~exist(curr_out, 'dir')
        mkdir(curr_out);
    end
	[ts, fr, dur, ep, fr_n, b, hmm] = compute_hmm(curr_in, sleep, threshs(i), curr_out, plot_figs);
%     continue;

	[~, down_ix] = min(cell2mat(cellfun(@mean, fr, 'UniformOutput', false)));
	sv_fn = sprintf('%s/TM-%s_Run001-%dX-tmpl%d/all-variables.mat', sv_dir, datasets{i}, compressions(i), templates(i));
	if i == 3 || i == 4
		sv_fn = sprintf('%s/TM-%s_Run001-%dX-tmpl%d/all-variables.mat', sv_dir2, datasets{i}, compressions(i), templates(i));
	end
	
	load(fullfile(curr_in, 'ts.mat'));
	if ~exist('e', 'var')
		e = struct();
		e.epochs = epochs;
	end
	if ~isfield(e.epochs, 'rest')
		if exist(fullfile(curr_in, 'epochs.mat'), 'file')
			load(fullfile(curr_in, 'epochs.mat'), 'rest');
			e.epochs.rest = rest;
		else
			load(fullfile(curr_in, 'ts_full.mat'), 'e');
		end
	end
	sleep_ep = e.epochs.(sprintf('sleep%d', sleep));
	ep = ep .* 10 + sleep_ep(1);

% 	load(fullfile(curr_in, sprintf('sleep%d_HMM_data.mat', sleep)), 'times', 'ids');
% 	n = max(ids);
% 	tt = max(times);
% 	q = full(sparse(ids(times > 0), times(times > 0), 1, n, tt));
    S = LoadSpikes(FindFiles('*.t', 'StartingDirectory', fullfile(curr_in, 'tfiles')));
    n = length(S);
	QS_binsize = 1;
% 	q2 = reshape(q(:, 1:(QS_binsize*floor(size(q, 2)/QS_binsize))), size(q, 1), QS_binsize, []);
% 	q2 = squeeze(sum(q2, 2));
% 	tt2 = 1:QS_binsize:tt;
% 	tt2 = tt2(1:size(q2, 2));
% 	Q_cut = ctsd(sleep_ep(1), QS_binsize*10, q2');
    Q = MakeQfromS(S, 10*QS_binsize);
    Q_cut = Restrict(Q, sleep_ep(1), sleep_ep(2));
	
    
    if exist(fullfile(curr_out, 'state_vector_decorrelation.mat'), 'file')
        load(fullfile(curr_out, 'state_vector_decorrelation.mat'), 'UP1_ix', 'UP2_ix');
        up1_ix = UP1_ix;
        up2_ix = UP2_ix;
    else
        fprintf('Running SVD on dataset %s\n', datasets{i});
        [up1_ix, up2_ix] = statevecdecorr(Q_cut, QS_binsize, ts, curr_out, down_ix, hmm, fr_n, true);
    end
    Q_task1 = Restrict(Q, e.epochs.task1(1), e.epochs.task1(2));
    Q_task2 = Restrict(Q, e.epochs.task2(1), e.epochs.task2(2));
    if exist(fullfile(curr_out, 'task_fit.mat'), 'file')
        load(fullfile(curr_out, 'task_fit.mat'));
    else
        [task1_fit, ~, ~, task1_mean, ~, ~, task1_ix, ~, task1_offset] = calc_svd(Range(Q_task1), full(Data(Q_task1)), e.epochs.task1, QS_binsize, 150);
        [task2_fit, ~, ~, task2_mean, ~, ~, task2_ix, ~, task2_offset] = calc_svd(Range(Q_task2), full(Data(Q_task2)), e.epochs.task2, QS_binsize, 150);
        save(fullfile(curr_out, 'task_fit.mat'), 'task1_fit', 'task2_fit', 'task1_mean', 'task2_mean', 'task1_ix', 'task2_ix', 'task1_offset', 'task2_offset');
    end
    
    fr_all(i, 1) = mean(fr{down_ix});
    fr_all(i, 2) = mean(fr{up1_ix});
    fr_all(i, 3) = mean(fr{up2_ix});
    
    fr_sem(i, 1) = std(fr{down_ix}) / sqrt(length(fr{down_ix}));
    fr_sem(i, 2) = std(fr{up1_ix}) / sqrt(length(fr{up1_ix}));
    fr_sem(i, 3) = std(fr{up2_ix}) / sqrt(length(fr{up2_ix}));
    
    dur_all(i, 1) = mean(dur{down_ix});
    dur_all(i, 2) = mean(dur{up1_ix});
    dur_all(i, 3) = mean(dur{up2_ix});
    
    dur_sem(i, 1) = std(dur{down_ix}) / sqrt(length(dur{down_ix}));
    dur_sem(i, 2) = std(dur{up1_ix}) / sqrt(length(dur{up1_ix}));
    dur_sem(i, 3) = std(dur{up2_ix}) / sqrt(length(dur{up2_ix}));
    
    if plot_figs
        tmp = cellfun(@(x)(mean(x, 2)), fr_n, 'UniformOutput', false);
        tmp2 = [tmp{down_ix}(:), tmp{up1_ix}(:), tmp{up2_ix}(:)];
        figure;
        bar(tmp2);
        xlabel('Neuron');
        ylabel('Firing Rate');
        legend({'Down', 'UP-1', 'UP-2'});
        title('Firing Rate per neuron');
        saveas(gcf, fullfile(curr_out, 'firing_rate_neuron.fig'));
        close;

        fcn = @(b, x)(b(1) * exp(-x ./ (b(2))) + b(3));
        figure;
        ax1 = subplot(2, 1, 1);
        hold all;
        plot(task1_mean);
        plot(fcn(task1_fit, (1:length(task1_mean)) - task1_offset));
        xlabel('Time (ms)');
        ylabel('Correlation');
        title(sprintf('Task 1 State Vector Decorrelation (\\tau = %gms)', task1_fit(2)*QS_binsize));
        ax2 = subplot(2, 1, 2);
        hold all;
        plot(task2_mean);
        plot(fcn(task2_fit, (1:length(task2_mean)) - task2_offset));
        xlabel('Time (ms)');
        ylabel('Correlation');
        title(sprintf('Task 2 State Vector Decorrelation (\\tau = %gms)', task2_fit(2)*QS_binsize));
        saveas(gcf, fullfile(curr_out, 'state_vector_decorr_task.fig'));
        close;
    end
    
    load(fullfile(curr_out, 'state_vector_decorrelation.mat'), 'tau');
    tau_all(i, 1) = tau(up1_ix);
    tau_all(i, 2) = tau(up2_ix);
    tau_task = mean([task1_fit(2), task2_fit(2)]);
    tau_task_all(i) = tau_task;
    tau_comp = tau_task ./ tau;
    tau_comp_all(i, 1) = tau_comp(up1_ix);
    tau_comp_all(i, 2) = tau_comp(up2_ix);
	   
    [cr, fs] = ReadCR_tsd(fullfile(curr_in, csc_files{i}), e.epochs.sleep3(1), e.epochs.sleep3(2));
    cr_t = Range(cr);
    cr_d = Data(cr);

    up_ts = cell(0, 2);
    up_ts_comb = cell(0, 1);
    
    up_length = zeros(0, 1);
    up_range = zeros(0, 2);
    up_start = zeros(0, 1);
    up_end = zeros(0, 1);
    up_ep = cell(size(ep, 1), 2);
    up_ep_length = zeros(size(ep, 1), 2);
    up_ep_dur = zeros(size(ep, 1), 2);
    up_ep_ts = ep;
    up_bin = cell(6, 2);
    up_bin_length = zeros(6, 2);
    up_bin_dur = zeros(6, 2);
    up_all_bin = cell(6, 1);
    tr_spec = cell(6, 1); % 1 = DOWN->UP1, 2 = DOWN->UP2, 3: UP1->UP2, 4: UP2->UP1, 5: UP1->DOWN, 6: UP2->DOWN
    tr_lfp = cell(6, 1);
    up_types = zeros(5, 1); % UP1 only, UP2 only, UP1->UP2, UP2->UP1, Multiple Transitions
    up_pow = cell(2, 1);
    fp = 1:2:300;
    task2_end = e.epochs.task2(2);
    for j = 1:6
        ts1 = task2_end + (j-1)*10*60e4;
        ts2 = task2_end + j*10*60e4;
        up_bin{j, 1} = ts{up1_ix}(ts{up1_ix}(:, 1) >= ts1 & ts{up1_ix}(:, 2) <= ts2, :);
        up_bin{j, 2} = ts{up2_ix}(ts{up2_ix}(:, 1) >= ts1 & ts{up2_ix}(:, 2) <= ts2, :);
        up_bin_length(j, 1) = size(up_bin{j, 1}, 1);
        up_bin_length(j, 2) = size(up_bin{j, 2}, 1);
        up_bin_dur(j, 1) = sum(diff(up_bin{j, 1}, 1, 2)./10);
        up_bin_dur(j, 2) = sum(diff(up_bin{j, 1}, 1, 2)./10);
    end
    up_bin_num_all(i, :, :) = up_bin_length;
    up_bin_dur_all(i, :, :) = up_bin_dur;
    for j = 1:size(ep, 1)
        up_ep{j, 1} = ts{up1_ix}(ts{up1_ix}(:, 1) >= ep(j, 1) & ts{up1_ix}(:, 2) <= ep(j, 2), :);
        up_ep{j, 2} = ts{up2_ix}(ts{up2_ix}(:, 1) >= ep(j, 1) & ts{up2_ix}(:, 2) <= ep(j, 2), :);
        up_ep_length(j, 1) = size(up_ep{j, 1}, 1);
        up_ep_length(j, 2) = size(up_ep{j, 2}, 1);
        curr_down = ts{down_ix}(ts{down_ix}(:, 2) >= ep(j, 1) & ts{down_ix}(:, 1) <= ep(j, 2), :); 
        if isempty(curr_down)
            continue;
        end
        curr_up_ep = [[ep(j, 1), curr_down(1, 1)]; [curr_down(1:end-1, 2), curr_down(2:end, 1)]; [curr_down(end, 2), ep(j, 2)]];
        curr_up_ep = sortrows(curr_up_ep);
        for k = 1:6
            ts1 = task2_end + (k-1)*10*60e4;
            ts2 = task2_end + k*10*60e4;
            up_all_bin{k} = [up_all_bin{k}; curr_up_ep(curr_up_ep(:, 1) >= ts1 & curr_up_ep(:, 2) <= ts2, :)];
        end
        for ii = 1:size(curr_up_ep, 1)
            curr_up = curr_up_ep(ii, :);
            curr_up1_ix = ts{up1_ix}(:, 2) > curr_up(1) & ts{up1_ix}(:, 1) < curr_up(2);
            curr_up2_ix = ts{up2_ix}(:, 2) > curr_up(1) & ts{up2_ix}(:, 1) < curr_up(2);
            if nnz(curr_up1_ix) == 0 && nnz(curr_up2_ix) == 0
                continue;
            end
            tmp_up1 = ts{up1_ix}(curr_up1_ix, :);
            tmp_up2 = ts{up2_ix}(curr_up2_ix, :);
            ts_all = sortrows([tmp_up1; tmp_up2]);
            up_ts{end+1, 1} = tmp_up1;
            up_ts{end, 2} = tmp_up2;
            ts_all2 = sortrows([[tmp_up1, ones(size(tmp_up1, 1), 1)]; [tmp_up2, 2.*ones(size(tmp_up2, 1), 1)]], 1);
            trange = 2e4;
            win = 256;
            overlap = floor(win/2);
            f = 0:1:30;
            
            if size(ts_all2, 1) == 1
                up_types(ts_all2(1, 3)) = up_types(ts_all2(1, 3)) + 1;
            elseif size(ts_all2, 1) == 2
                up_types(ts_all2(1, 3)+2) = up_types(ts_all2(1, 3)+2) + 1;
            else
                up_types(5) = up_types(5) + 1;
            end
            
            for jj = 1:size(ts_all2, 1)
                curr_cr = Restrict(cr, ts_all2(jj, 1), ts_all2(jj, 2));
                if length(Data(curr_cr)) >= 300
                    if isempty(up_pow{ts_all2(jj, 3)})
                        up_pow{ts_all2(jj, 3)}(1, :) = pwelch(Data(curr_cr), [], [], fp, fs);
                    else
                        up_pow{ts_all2(jj, 3)}(end+1, :) = pwelch(Data(curr_cr), [], [], fp, fs);
                    end
                end
                if jj == 1 % DOWN to UP
%                     curr_cr = Restrict(cr, ts_all2(jj, 1) - trange, ts_all2(jj, 1) + trange);
%                     crd = full(Data(curr_cr));
                    crd = cr_d(cr_t >= ts_all2(jj, 1) - trange & cr_t <= ts_all(jj, 1) + trange);
                    [~, ~, t_spec, curr_spec] = spectrogram(crd, win, overlap, f, fs);
                    if ts_all2(jj, 3) == 1 % DOWN to UP1
                        ix = 1;
                    else % DOWN to UP2
                        ix = 2;
                    end
                    if isempty(tr_spec{ix})
                        tr_spec{ix} = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix} = di(:);
                    else
                        tr_spec{ix}(:, :, end+1) = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix}(:, end+1) = di;
                    end
                end
                crd = cr_d(cr_t >= ts_all2(jj, 2) - trange & cr_t <= ts_all(jj, 2) + trange);
                if jj == size(ts_all2, 1) % UP to DOWN
%                     curr_cr = Restrict(cr, ts_all2(jj, 2) - trange, ts_all2(jj, 2) + trange);
%                     crd = Data(curr_cr);
                    [~, ~, t_spec, curr_spec] = spectrogram(crd, win, overlap, f, fs);
                    if ts_all2(jj, 3) == 1 % UP1 to DOWN
                        ix = 5;
                    else % UP2 to DOWN
                        ix = 6;
                    end
                    if isempty(tr_spec{ix})
                        tr_spec{ix} = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix} = di(:);
                    else
                        tr_spec{ix}(:, :, end+1) = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix}(:, end+1) = di;
                    end
                else
%                     curr_cr = Restrict(cr, ts_all2(jj, 2) - trange, ts_all2(jj, 2) + trange);
%                     crd = Data(curr_cr);
                    [~, ~, t_spec, curr_spec] = spectrogram(crd, win, overlap, f, fs);
                    if ts_all2(jj, 3) == 1 % UP1 to UP2
                        ix = 3;
                    else % UP2 to UP1
                        ix = 4;
                    end
                    if isempty(tr_spec{ix})
                        tr_spec{ix} = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix} = di(:);
                    else
                        tr_spec{ix}(:, :, end+1) = curr_spec;
                        cr_len = length(crd);
                        di = interp1(linspace(-2, 2, cr_len), crd, linspace(-2, 2, 2000));
                        tr_lfp{ix}(:, end+1) = di;
                    end
                end
            end
            up_ts_comb{end+1} = ts_all2;
            up_range(end+1, :) = [ts_all(1, 1), ts_all(end, 2)];
            up_length(end+1) = diff(up_range(end, :));
            
            if ~isempty(tmp_up1) && tmp_up1(1, 1) == ts_all(1, 1)
                up_start(end+1) = 1;
            else
                up_start(end+1) = 2;
            end

            if ~isempty(tmp_up1) && tmp_up1(end, 2) == ts_all(end, 2)
                up_end(end+1) = 1;
            else
                up_end(end+1) = 2;
            end
        end
    end
    tr_lfp_all(i, :) = tr_lfp;
    up_types_all(i, :) = up_types;
    if plot_figs
        up_ep_num_ratio = up_ep_length ./ sum(up_ep_length, 2);
        figure;
        bar(up_ep_num_ratio, 'stacked');
        title('Ratio of UP-1 and UP-2 per epoch');
        legend('UP-1', 'UP-2');
        xlabel('Epoch');
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_ep_num_ratio_stacked.fig'));
        close;

        figure;
        plot(up_ep_num_ratio);
        title('Ratio of UP-1 and UP-2 per epoch');
        legend('UP-1', 'UP-2');
        xlabel('Epoch');
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_ep_num_ratio_line.fig'));
        close;

        figure;
        bar(up_ep_num_ratio, 'grouped');
        title('Ratio of UP-1 and UP-2 per epoch');
        legend('UP-1', 'UP-2');
        xlabel('Epoch');
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_ep_num_ratio_grouped.fig'));
        close;

        up_bin_num_ratio = up_bin_length ./ sum(up_bin_length, 2);
        figure;
        bar(up_bin_num_ratio, 'stacked');
        title('Ratio of UP-1 and UP-2 for 10 min bins');
        legend('UP-1', 'UP-2');
        xticks([1, 2, 3]);
        xticklabels({'0-10 min', '10-20 min', '20-30 min'});
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_bin_num_ratio_stacked.fig'));
        close;

        figure;
        plot(up_bin_num_ratio);
        title('Ratio of UP-1 and UP-2 per epoch');
        legend('UP-1', 'UP-2');
        xticks([1, 2, 3]);
        xticklabels({'0-10 min', '10-20 min', '20-30 min'});
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_bin_num_ratio_line.fig'));
        close;

        figure;
        bar(up_bin_num_ratio, 'grouped');
        title('Ratio of UP-1 and UP-2 per epoch');
        legend('UP-1', 'UP-2');
        xticks([1, 2, 3]);
        xticklabels({'0-10 min', '10-20 min', '20-30 min'});
        ylabel('Ratio');
        saveas(gcf, fullfile(curr_out, 'up_bin_num_ratio_grouped.fig'));
        close;

        transition_names = {'DOWN to UP-1', 'DOWN to UP-2', 'UP-1 to UP-2', 'UP-2 to UP-1', 'UP-1 to DOWN', 'UP-2 to DOWN'};
        figure;
        for j = 1:6
            subplot(3, 2, j);
            plot(linspace(-trange/1e4, trange/1e4, size(tr_lfp{j}, 1)), mean(tr_lfp{j}, 2));
            title(sprintf('%s (# transitions: %d)', transition_names{j}, size(tr_lfp{j}, 2)));
            if j >= 5
                xlabel('Time (s)');
            end
            if mod(j, 2) > 0
                ylabel('ADC values');
            end
        end
        saveas(gcf, fullfile(curr_out, 'transition_lfp_mean.fig')); 
        close;

        figure;
        for j = 1:6
            subplot(3, 2, j);
            imagesc(linspace(-trange/1e4, trange/1e4, size(tr_spec{j}, 2)), f, 10*log10(mean(tr_spec{j}, 3)));
            axis xy;
            title(sprintf('%s (# transitions: %d)', transition_names{j}, size(tr_lfp{j}, 2)));
            if j >= 5
                xlabel('Time (s)');
            end
            if mod(j, 2) > 0
                ylabel('Frequency (Hz)');
            end
        end
        saveas(gcf, fullfile(curr_out, 'transition_spec_mean.fig')); 
        close;

        figure;
        for j = 1:6
            subplot(3, 2, j);
            curr_spec = mean(tr_spec{j}, 3);
            hold all;
            curr_x = linspace(-trange/1e4, trange/1e4, size(curr_spec, 2));
            plot(curr_x, mean(curr_spec(f >= 1 & f < 5, :)));
            plot(curr_x, mean(curr_spec(f >= 5 & f < 10, :)));
            plot(curr_x, mean(curr_spec(f >= 10 & f < 20, :)));
            if j >= 5
                xlabel('Time (s)');
            end
            if mod(j, 2) > 0
                ylabel('Power');
            end
            if j == 2
                legend('1-5Hz', '5-10Hz', '10-20Hz', 'Location', 'NorthEastOutside');
            end
        end
        saveas(gcf, fullfile(curr_out, 'transition_spec_bands.fig'));
        close;

        figure;
        hold all;
        errorbar(fp, mean(log(up_pow{1})), 2.*std(log(up_pow{1}))/sqrt(size(up_pow{1}, 1)));
        errorbar(fp, mean(log(up_pow{2})), 2.*std(log(up_pow{2}))/sqrt(size(up_pow{2}, 1)));
        legend('UP-1', 'UP-2');
        title('UP-1 vs UP-2 power spectral density (2*SEM)');
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        saveas(gcf, fullfile(curr_out, 'up_power_spectrum.fig'));
        close;

        figure;
        hold all;
        bm = zeros(length(freq_names), 2);
        bs = zeros(length(freq_names), 2);
        for j = 1:length(freq_names)
            bm(j, 1) = mean(mean(log(up_pow{1}(:, fp >= freq_range(j, 1) & fp < freq_range(j, 2)))));
            bm(j, 2) = mean(mean(log(up_pow{2}(:, fp >= freq_range(j, 1) & fp < freq_range(j, 2)))));
            bs(j, 1) = std(mean(log(up_pow{1}(:, fp >= freq_range(j, 1) & fp < freq_range(j, 2)))));
            bs(j, 2) = std(mean(log(up_pow{2}(:, fp >= freq_range(j, 1) & fp < freq_range(j, 2)))));
        end
        bar(bm);
        bw = 2/3.5;
        for j = 1:2
            errorbar((1:length(freq_names)) - bw/2 + (2*j-1) * bw / 4, bm(:, j), bs(:, j), 'k.');
        end
        xticks(1:length(freq_names));
        xticklabels(freq_names);
        legend('UP-1', 'UP-2');
        ylabel('Log Average power');
        title('UP-1 vs UP-2 power in frequency bands');
        saveas(gcf, fullfile(curr_out, 'up_power_bands.fig'));
        close;

        pow_all(:, :, i) = bm;

    
        figure;
        up_start1 = nnz(up_start == 1)/length(up_start);
        up_start2 = nnz(up_start == 2)/length(up_start);
        bar([up_start1, up_start2], 'b');
        xlim([0.5, 2.5]);
        ylim([0, 1]);
        set(gca, 'XTickLabels', {'UP-1', 'UP-2'});
        title('Starting UP-state subtype');
        saveas(gcf, fullfile(curr_out, 'up_starting_state.fig'));
        close;
        start_state(i, 1) = up_start1;
        start_state(i, 2) = up_start2;

        figure;
        up_end1 = nnz(up_end == 1)/length(up_end);
        up_end2 = nnz(up_end == 2)/length(up_end);
        bar([up_end1, up_end2], 'b');
        xlim([0.5, 2.5]);
        ylim([0, 1]);
        set(gca, 'XTickLabels', {'UP-1', 'UP-2'});
        title('Ending UP-state subtype');
        saveas(gcf, fullfile(curr_out, 'up_ending_state.fig'));
        close;
        finish_state(i, 1) = up_end1;
        finish_state(i, 2) = up_end2;
    end
    
    save(sprintf('%s/up_stats.mat', curr_out), 'up_ts', 'up_ts_comb', 'up_range', 'up_length', 'up_start', 'up_end', 'up_ep', 'up_ep_length', 'tr_spec', 'tr_lfp', 'trange', 'f', 'up_pow', 'tau_task', 'tau_comp', 'up_ep_ts', 'up_bin', 'up_bin_length', 'up_types', 'up_all_bin');
	load(sv_fn, 'tmpl', 'QS_binsize');
	for k = 1:length(z_thresh)
		[replay_pos, replay_peak, replay_num, replay_abs, replay_max, replay_max_abs] = get_replay_stats(tmpl, QS_binsize, ts, z_thresh(k), curr_out);
		curr_r_num = [sum(replay_num{up1_ix}), sum(replay_num{up2_ix})];
		figure;
		bar(curr_r_num);
		title(sprintf('Z-score = %d', z_thresh(k)));
		set(gca, 'XTickLabels', {'UP-1', 'UP-2'});
		% ylim([0, 1]);
		xlim([0.5, 2.5]);
		saveas(gcf, sprintf('%s/replay_count_bar_z%d.fig', curr_out, z_thresh(k)));
		close;
        replay_count_z(i, :, k) = curr_r_num;

% 		tmp_pca = compute_mapping(cell2mat(b([up1_ix, up2_ix])')', 'PCA');
        [coeff, tmp_pca, latent, ~, explained] = pca(cell2mat(b([up1_ix, up2_ix])')');
		pca_bin = tmp_pca;
		tmp_len = size(ts{up1_ix}, 1);
		up1_replay_ix = double(replay_num{up1_ix} > 0);
		up2_replay_ix = double(replay_num{up2_ix} > 0);
		figure;
		hold all;
		scatter(tmp_pca(1:tmp_len, 1), tmp_pca(1:tmp_len, 2), 40 .* (5 .* up1_replay_ix + 1), 'b', '.');
		scatter(tmp_pca(tmp_len+1:end, 1), tmp_pca(tmp_len+1:end, 2), 40 .* (5 .* up2_replay_ix + 1), 'r', '.');
		legend('UP-1', 'UP-1 replay', 'UP-2', 'UP-2 replay');
		xlabel('PC1');
		ylabel('PC2');
		title(sprintf('PCA Binary Z=%d', z_thresh(k)));
		saveas(gcf, sprintf('%s/pca_binary_z%d.fig', curr_out, z_thresh(k)));
		close;

% 		tmp_pca = compute_mapping(cell2mat(fr_n([up1_ix, up2_ix])')', 'PCA');
        [coeff, tmp_pca, latent, ~, explained] = pca(cell2mat(fr_n([up1_ix, up2_ix])')');
		pca_fr = tmp_pca;
		tmp_len = size(ts{up1_ix}, 1);
		up1_replay_ix = double(replay_num{up1_ix} > 0);
		up2_replay_ix = double(replay_num{up2_ix} > 0);
		figure;
		hold all;
		scatter(tmp_pca(1:tmp_len, 1), tmp_pca(1:tmp_len, 2), 40 .* (5 .* up1_replay_ix + 1), 'b', '.');
		scatter(tmp_pca(tmp_len+1:end, 1), tmp_pca(tmp_len+1:end, 2), 40 .* (5 .* up2_replay_ix + 1), 'r', '.');
		legend('UP-1', 'UP-2');
		xlabel('PC1');
		ylabel('PC2');
		title(sprintf('PCA Firing Rate (Unnormalized) Z=%d', z_thresh(k)));
		saveas(gcf, sprintf('%s/pca_fr_z%d.fig', curr_out, z_thresh(k)));
		close;

% 		tmp_pca = compute_mapping(zscore(cell2mat(fr_n([up1_ix, up2_ix])')'), 'PCA');
        tmp_z = zscore(cell2mat(fr_n([up1_ix, up2_ix])')', 0, 1);
        grp = [ones(size(fr_n{up1_ix}, 2), 1); 2.*ones(size(fr_n{up2_ix}, 2), 1)];
        [coeff, tmp_pca, latent, ~, explained] = pca(tmp_z);
		pca_fr_norm = tmp_pca;
        explained_fr_norm = explained;
        latent_fr_norm = latent;
        coeff_fr_norm = coeff;
		tmp_len = size(ts{up1_ix}, 1);
		up1_replay_ix = double(replay_num{up1_ix} > 0);
		up2_replay_ix = double(replay_num{up2_ix} > 0);
		figure;
		hold all;
% 		scatter(tmp_pca(1:tmp_len, 1), tmp_pca(1:tmp_len, 2), 40 .* (5 .* up1_replay_ix + 1), 'b', '.');
% 		scatter(tmp_pca(tmp_len+1:end, 1), tmp_pca(tmp_len+1:end, 2), 40 .* (5 .* up2_replay_ix + 1), 'r', '.');
        scatter(tmp_pca(find(~up1_replay_ix), 1), tmp_pca(find(~up1_replay_ix), 2), 20, [0.4, 0.4, 1], '.');
        s1 = scatter(tmp_pca(find(up1_replay_ix), 1), tmp_pca(find(up1_replay_ix), 2), 120, 'b', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'b');
        scatter(tmp_pca(tmp_len+find(~up2_replay_ix), 1), tmp_pca(tmp_len+find(~up2_replay_ix), 2), 20, [1, 0.4, 0.4], '.');
        s2 = scatter(tmp_pca(tmp_len+find(up2_replay_ix), 1), tmp_pca(tmp_len+find(up2_replay_ix), 2), 120, 'r', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'r');
		legend([s1, s2], 'UP-1', 'UP-2');
		xlabel('PC1');
		ylabel('PC2');
		title(sprintf('PCA Firing Rate (normalized) Z=%d', z_thresh(k)));
		saveas(gcf, sprintf('%s/pca_fr_norm_z%d.fig', curr_out, z_thresh(k)));
		close;
        
        figure;
        hold all;
        plot(latent);
        title('PCA Eigenvalues Firing Rate (normalized)');
        saveas(gcf, sprintf('%s/pca_latent_fr_norm.fig', curr_out));
        close;
        
        figure;
        hold all;
        bar(explained);
        ylabel('% explained');
        xlabel('PCA #');
        ylim([0, 50]);
        title('PCA % Explained Firing Rate (normalized)');
        saveas(gcf, sprintf('%s/pca_explained_fr_norm.fig', curr_out));
        close;

        d = cluster_dist(tmp_z, grp);
        d_all(i) = d;
        
        lc = @(g)(sum(1-chi2cdf(mahal(tmp_z(grp ~= g, :), tmp_z(grp == g, :)), n))./nnz(grp == g));
        
        lc1 = lc(1);
        lc2 = lc(2);
        
        figure;
        bar([lc1, lc2]);
        xticklabels({'UP-1', 'UP-2'});
        xlabel('Upstate');
        ylabel('L-ratio');
        title('L_{ratio} for PCA Firing Rate (normalized) clusters');
        saveas(gcf, sprintf('%s/pca_l_ratio_fr_norm.fig', curr_out));
        close;
        
        figure;
        silhouette(tmp_z, grp);
        title('Silhouette for PCA Firing Rate (normalized)');
        saveas(gcf, sprintf('%s/pca_silhouette_fr_norm.fig', curr_out));
        close;
        
        figure;
        subplot(1, 3, 1);
        hold all;
        scatter(tmp_pca(find(~up1_replay_ix), 1), tmp_pca(find(~up1_replay_ix), 2), 20, [0.4, 0.4, 1], '.');
        s1 = scatter(tmp_pca(find(up1_replay_ix), 1), tmp_pca(find(up1_replay_ix), 2), 120, 'b', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'b');
        scatter(tmp_pca(tmp_len+find(~up2_replay_ix), 1), tmp_pca(tmp_len+find(~up2_replay_ix), 2), 20, [1, 0.4, 0.4], '.');
        s2 = scatter(tmp_pca(tmp_len+find(up2_replay_ix), 1), tmp_pca(tmp_len+find(up2_replay_ix), 2), 120, 'r', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'r');
        title('PC1 vs PC2');
        subplot(1, 3, 2);
        hold all;
        scatter(tmp_pca(find(~up1_replay_ix), 1), tmp_pca(find(~up1_replay_ix), 3), 20, [0.4, 0.4, 1], '.');
        s1 = scatter(tmp_pca(find(up1_replay_ix), 1), tmp_pca(find(up1_replay_ix), 3), 120, 'b', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'b');
        scatter(tmp_pca(tmp_len+find(~up2_replay_ix), 1), tmp_pca(tmp_len+find(~up2_replay_ix), 3), 20, [1, 0.4, 0.4], '.');
        s2 = scatter(tmp_pca(tmp_len+find(up2_replay_ix), 1), tmp_pca(tmp_len+find(up2_replay_ix), 3), 120, 'r', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'r');
        title('PC1 vs PC3');
        subplot(1, 3, 3);
        hold all;
        scatter(tmp_pca(find(~up1_replay_ix), 2), tmp_pca(find(~up1_replay_ix), 3), 20, [0.4, 0.4, 1], '.');
        s1 = scatter(tmp_pca(find(up1_replay_ix), 2), tmp_pca(find(up1_replay_ix), 3), 120, 'b', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'b');
        scatter(tmp_pca(tmp_len+find(~up2_replay_ix), 2), tmp_pca(tmp_len+find(~up2_replay_ix), 3), 20, [1, 0.4, 0.4], '.');
        s2 = scatter(tmp_pca(tmp_len+find(up2_replay_ix), 2), tmp_pca(tmp_len+find(up2_replay_ix), 3), 120, 'r', '.', 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'r');
        title('PC2 vs PC3');
        saveas(gcf, fullfile(curr_out, sprintf('pca_pc123_z%d.fig', z_thresh(k))));
        close;
        
        bins = -3:0.2:8;
        figure;
        hold all;
        subplot(3, 1, 1);
        histogram(replay_max{down_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Z-score');
        ylabel('Count');
        title({'Z-scored Replay Peak distribution', 'DOWN'});
        subplot(3, 1, 2);
        histogram(replay_max{up1_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Z-score');
        ylabel('Count');
        title('UP-1');
        subplot(3, 1, 3);
        histogram(replay_max{up2_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Z-score');
        ylabel('Count');
        title('UP-2');
        saveas(gcf, sprintf('%s/replay_zscore_dist.fig', curr_out));
        close;
        
        bins = -0.3:0.05:0.5;
        figure;
        hold all;
        subplot(3, 1, 1);
        histogram(replay_max_abs{down_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Z-score');
        ylabel('Count');
        title({'Z-scored Replay Peak distribution', 'DOWN'});
        subplot(3, 1, 2);
        histogram(replay_max_abs{up1_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Correlation');
        ylabel('Count');
        title('UP-1');
        subplot(3, 1, 3);
        histogram(replay_max_abs{up2_ix}, bins, 'Normalization', 'probability');
        xlabel('Replay Correlation');
        ylabel('Count');
        title('UP-2');
        saveas(gcf, sprintf('%s/replay_corr_dist.fig', curr_out));
        close;
        
        up_replay = cell(0, 1);
        for j = 1:size(ep, 1)
            curr_down = ts{down_ix}(ts{down_ix}(:, 2) >= ep(j, 1) & ts{down_ix}(:, 1) <= ep(j, 2), :); 
            if isempty(curr_down)
                continue;
            end
            curr_up_ep = [[ep(j, 1), curr_down(1, 1)]; [curr_down(1:end-1, 2), curr_down(2:end, 1)]; [curr_down(end, 2), ep(j, 2)]];
            curr_up_ep = sortrows(curr_up_ep);
            for ii = 1:size(curr_up_ep, 1)
                curr_up = curr_up_ep(ii, :);
                curr_up1_ix = ts{up1_ix}(:, 2) > curr_up(1) & ts{up1_ix}(:, 1) < curr_up(2);
                curr_up2_ix = ts{up2_ix}(:, 2) > curr_up(1) & ts{up2_ix}(:, 1) < curr_up(2);
                if nnz(curr_up1_ix) == 0 && nnz(curr_up2_ix) == 0
                    continue;
                end
                tmp_replay1 = replay_pos{up1_ix}(curr_up1_ix);
                tmp_replay2 = replay_pos{up2_ix}(curr_up2_ix);
                up_replay{end+1} = sort([cell2mat(tmp_replay1); cell2mat(tmp_replay2)]);
            end
        end

		try
            if plot_figs
                up_react_ix = find(~cellfun(@isempty, up_replay));
                [~, sort_react_ix] = sort(up_length(up_react_ix), 'descend');
                up_react_ix = up_react_ix(sort_react_ix);
                up_react = up_ts(up_react_ix, :);
                up_react_range = up_range(up_react_ix, :);
                up_replay_pos = up_replay(up_react_ix);
                figure;
                hold all;
                tmp_up1 = zeros(0, 3);
                tmp_up2 = zeros(0, 3);
                tmp_replay = zeros(0, 2);
                for j = 1:size(up_react, 1)
                    tmp_up1 = [tmp_up1; [up_react{j, 1} - up_react_range(j, 1), j*ones(size(up_react{j, 1}, 1), 1)]];
                    tmp_up2 = [tmp_up2; [up_react{j, 2} - up_react_range(j, 1), j*ones(size(up_react{j, 2}, 1), 1)]];
                    tmp_replay = [tmp_replay; [up_replay_pos{j} - up_react_range(j, 1), j*ones(size(up_replay_pos{j}, 1), 1)]];

                end
                tmp_x1 = [tmp_up1(:, 1:2)'./1e4; NaN(1, size(tmp_up1, 1))];
                tmp_y1 = [tmp_up1(:, 3)'.*ones(2, size(tmp_up1, 1)); NaN(1, size(tmp_up1, 1))];
                tmp_x2 = [tmp_up2(:, 1:2)'./1e4; NaN(1, size(tmp_up2, 1))];
                tmp_y2 = [tmp_up2(:, 3)'.*ones(2, size(tmp_up2, 1)); NaN(1, size(tmp_up2, 1))];
                h1 = plot(tmp_x1(:), tmp_y1(:), 'b');
                h2 = plot(tmp_x2(:), tmp_y2(:), 'r');
                plot(tmp_replay(:, 1)./1e4, tmp_replay(:, 2), 'k.', 'MarkerSize', 10);
                legend([h1(1), h2(1)], 'UP-1', 'UP-2');
                xlabel('Duration (s)');
                ylabel('UP #');
                title(sprintf('State Sequence Sorted - Replay Only (Z=%d)', z_thresh(k)));
                saveas(gcf, sprintf('%s/state_sequence_replay_only_sorted_z%d.fig', curr_out, z_thresh(k)));
                close;

                up_react_ix = find(~cellfun(@isempty, up_replay));
                up_react = up_ts(up_react_ix, :);
                up_react_range = up_range(up_react_ix, :);
                up_replay_pos = up_replay(up_react_ix);
                figure;
                hold all;
                tmp_up1 = zeros(0, 3);
                tmp_up2 = zeros(0, 3);
                tmp_replay = zeros(0, 2);
                for j = 1:size(up_react, 1)
                    tmp_up1 = [tmp_up1; [up_react{j, 1} - up_react_range(j, 1), j*ones(size(up_react{j, 1}, 1), 1)]];
                    tmp_up2 = [tmp_up2; [up_react{j, 2} - up_react_range(j, 1), j*ones(size(up_react{j, 2}, 1), 1)]];
                    tmp_replay = [tmp_replay; [up_replay_pos{j} - up_react_range(j, 1), j*ones(size(up_replay_pos{j}, 1), 1)]];

                end
                tmp_x1 = [tmp_up1(:, 1:2)'./1e4; NaN(1, size(tmp_up1, 1))];
                tmp_y1 = [tmp_up1(:, 3)'.*ones(2, size(tmp_up1, 1)); NaN(1, size(tmp_up1, 1))];
                tmp_x2 = [tmp_up2(:, 1:2)'./1e4; NaN(1, size(tmp_up2, 1))];
                tmp_y2 = [tmp_up2(:, 3)'.*ones(2, size(tmp_up2, 1)); NaN(1, size(tmp_up2, 1))];
                h1 = plot(tmp_x1(:), tmp_y1(:), 'b');
                h2 = plot(tmp_x2(:), tmp_y2(:), 'r');
                plot(tmp_replay(:, 1)./1e4, tmp_replay(:, 2), 'k.', 'MarkerSize', 10);
                legend([h1(1), h2(1)], 'UP-1', 'UP-2');
                xlabel('Duration (s)');
                ylabel('UP #');
                title(sprintf('State Sequence Chronological - Replay Only (Z=%d)', z_thresh(k)));
                saveas(gcf, sprintf('%s/state_sequence_replay_only_chrono_z%d.fig', curr_out, z_thresh(k)));
                close;

                up_react_ix = find(up_length > 0);
                [~, sort_react_ix] = sort(up_length(up_react_ix), 'descend');
                up_react_ix = up_react_ix(sort_react_ix);
                up_react = up_ts(up_react_ix, :);
                up_react_range = up_range(up_react_ix, :);
                up_replay_pos = up_replay(up_react_ix);
                figure;
                hold all;
                tmp_up1 = zeros(0, 3);
                tmp_up2 = zeros(0, 3);
                tmp_replay = zeros(0, 2);
                for j = 1:size(up_react, 1)
                    tmp_up1 = [tmp_up1; [up_react{j, 1} - up_react_range(j, 1), j*ones(size(up_react{j, 1}, 1), 1)]];
                    tmp_up2 = [tmp_up2; [up_react{j, 2} - up_react_range(j, 1), j*ones(size(up_react{j, 2}, 1), 1)]];
                    tmp_replay = [tmp_replay; [up_replay_pos{j} - up_react_range(j, 1), j*ones(size(up_replay_pos{j}, 1), 1)]];

                end
                tmp_x1 = [tmp_up1(:, 1:2)'./1e4; NaN(1, size(tmp_up1, 1))];
                tmp_y1 = [tmp_up1(:, 3)'.*ones(2, size(tmp_up1, 1)); NaN(1, size(tmp_up1, 1))];
                tmp_x2 = [tmp_up2(:, 1:2)'./1e4; NaN(1, size(tmp_up2, 1))];
                tmp_y2 = [tmp_up2(:, 3)'.*ones(2, size(tmp_up2, 1)); NaN(1, size(tmp_up2, 1))];
                h1 = plot(tmp_x1(:), tmp_y1(:), 'b');
                h2 = plot(tmp_x2(:), tmp_y2(:), 'r');
                plot(tmp_replay(:, 1)./1e4, tmp_replay(:, 2), 'k.', 'MarkerSize', 10);
                legend([h1(1), h2(1)], 'UP-1', 'UP-2');
                xlabel('Duration (s)');
                ylabel('UP #');
                title(sprintf('State Sequence Sorted - All (Z=%d)', z_thresh(k)));
                saveas(gcf, sprintf('%s/state_sequence_all_sorted_z%d.fig', curr_out, z_thresh(k)));
                close;

                up_react_ix = find(up_length > 0);
                up_react = up_ts(up_react_ix, :);
                up_react_range = up_range(up_react_ix, :);
                up_replay_pos = up_replay(up_react_ix);
                figure;
                hold all;
                tmp_up1 = zeros(0, 3);
                tmp_up2 = zeros(0, 3);
                tmp_replay = zeros(0, 2);
                for j = 1:size(up_react, 1)
                    tmp_up1 = [tmp_up1; [up_react{j, 1} - up_react_range(j, 1), j*ones(size(up_react{j, 1}, 1), 1)]];
                    tmp_up2 = [tmp_up2; [up_react{j, 2} - up_react_range(j, 1), j*ones(size(up_react{j, 2}, 1), 1)]];
                    tmp_replay = [tmp_replay; [up_replay_pos{j} - up_react_range(j, 1), j*ones(size(up_replay_pos{j}, 1), 1)]];
                end
                tmp_x1 = [tmp_up1(:, 1:2)'./1e4; NaN(1, size(tmp_up1, 1))];
                tmp_y1 = [tmp_up1(:, 3)'.*ones(2, size(tmp_up1, 1)); NaN(1, size(tmp_up1, 1))];
                tmp_x2 = [tmp_up2(:, 1:2)'./1e4; NaN(1, size(tmp_up2, 1))];
                tmp_y2 = [tmp_up2(:, 3)'.*ones(2, size(tmp_up2, 1)); NaN(1, size(tmp_up2, 1))];
                h1 = plot(tmp_x1(:), tmp_y1(:), 'b');
                h2 = plot(tmp_x2(:), tmp_y2(:), 'r');
                plot(tmp_replay(:, 1)./1e4, tmp_replay(:, 2), 'k.', 'MarkerSize', 10);
                legend([h1(1), h2(1)], 'UP-1', 'UP-2');
                xlabel('Duration (s)');
                ylabel('UP #');
                title(sprintf('State Sequence Chronological - All (Z=%d)', z_thresh(k)));
                saveas(gcf, sprintf('%s/state_sequence_all_chrono_z%d.fig', curr_out, z_thresh(k)));
                close;

                up_react_ix = find(~cellfun(@isempty, up_replay));
                figure;
                up_start1 = nnz(up_start(up_react_ix) == 1)/length(up_start(up_react_ix));
                up_start2 = nnz(up_start(up_react_ix) == 2)/length(up_start(up_react_ix));
                bar([up_start1, up_start2], 'b');
                xlim([0.5, 2.5]);
                ylim([0, 1]);
                set(gca, 'XTickLabels', {'UP-1', 'UP-2'});
                title('Starting UP-state subtype (reactivating states)');
                saveas(gcf, fullfile(curr_out, 'up_starting_state_replay_only.fig'));
                close;
                start_state_replay(i, 1) = up_start1;
                start_state_replay(i, 2) = up_start2;

                figure;
                up_end1 = nnz(up_end(up_react_ix) == 1)/length(up_end(up_react_ix));
                up_end2 = nnz(up_end(up_react_ix) == 2)/length(up_end(up_react_ix));
                bar([up_end1, up_end2], 'b');
                xlim([0.5, 2.5]);
                ylim([0, 1]);
                set(gca, 'XTickLabels', {'UP-1', 'UP-2'});
                title('Ending UP-state subtype (reactivating states)');
                saveas(gcf, fullfile(curr_out, 'up_ending_state_replay_only.fig'));
                close;
                finish_state_replay(i, 1) = up_end1;
                finish_state_replay(i, 2) = up_end2;
            end
		catch
			close all;
		end

		save(sprintf('%s/up_sequence_z%d.mat', curr_out, z_thresh(k)), 'up_replay', 'pca_bin', 'pca_fr', 'pca_fr_norm', 'explained_fr_norm', 'latent_fr_norm', 'coeff_fr_norm');
    end

    if plot_figs
        hist_bins = 0:0.25:5;
        figure;
        subplot(1, 3, 1);
        histogram(fr{down_ix}, hist_bins, 'Normalization', 'probability');
        title('DOWN');
        ylabel('Probability');
        ylim([0, 1]);
        subplot(1, 3, 2);
        histogram(fr{up1_ix}, hist_bins, 'Normalization', 'probability');
        title('UP-1');
        xlabel('Firing Rate (Hz)');
        ylim([0, 1]);
        subplot(1, 3, 3);
        histogram(fr{up2_ix}, hist_bins, 'Normalization', 'probability');
        title('UP-2');
        ylim([0, 1]);
        saveas(gcf, fullfile(curr_out, 'firing_rate_hist.fig'));
        close;

        figure;
        hold all;
        histogram(fr{down_ix}, hist_bins, 'Normalization', 'probability');
        histogram(fr{up1_ix}, hist_bins, 'Normalization', 'probability');
        histogram(fr{up2_ix}, hist_bins, 'Normalization', 'probability');
        legend('Down', 'UP-1', 'UP-2');
        ylim([0, 1]);
        title('UP subtype Firing Rate');
        xlabel('Firing Rate (Hz)')
        ylabel('Probability');
        saveas(gcf, fullfile(curr_out, 'firing_rate_combined.fig'));
        close;

        figure;
        subplot(1, 3, 1);
        histogram(dur{down_ix}, hist_bins, 'Normalization', 'probability');
        title('DOWN');
        ylabel('Probability');
        ylim([0, 1]);
        subplot(1, 3, 2);
        histogram(dur{up1_ix}, hist_bins, 'Normalization', 'probability');
        title('UP-1');
        xlabel('Duration (s)');
        ylim([0, 1]);
        subplot(1, 3, 3);
        histogram(dur{up2_ix}, hist_bins, 'Normalization', 'probability');
        title('UP-2');
        ylim([0, 1]);
        saveas(gcf, fullfile(curr_out, 'duration_hist.fig'));
        close;

        figure;
        hold all;
        histogram(dur{down_ix}, hist_bins, 'Normalization', 'probability');
        histogram(dur{up1_ix}, hist_bins, 'Normalization', 'probability');
        histogram(dur{up2_ix}, hist_bins, 'Normalization', 'probability');
        legend('Down', 'UP-1', 'UP-2');
        ylim([0, 1]);
        title('UP subtype Duration');
        xlabel('Duration (s)');
        ylabel('Probability');
        saveas(gcf, fullfile(curr_out, 'duration_combined.fig'));
        close;
    end

	
    clear e;
    fprintf('%s took %g seconds\n', datasets{i}, toc);
end

if plot_figs
    figure;
    hold all;
    bar(mean(tau_all));
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    errorbar(mean(tau_all), std(tau_all)/sqrt(length(datasets)), 'k.');
    title('Average \tau across datasets');
    ylabel('\tau (ms)');
    saveas(gcf, fullfile(out_dir, 'tau_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(tau_comp_all));
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    errorbar(mean(tau_comp_all), std(tau_comp_all)/sqrt(length(datasets)), 'k.');
    title('Average \tau compression');
    ylabel('Compression (x faster)');
    saveas(gcf, fullfile(out_dir, 'tau_compression.fig'));
    close;

    figure;
    hold all;
    bm = mean(pow_all, 3);
    bs = std(pow_all, 0, 3);
    bar(bm);
    bw = 2/3.5;
    for j = 1:2
        errorbar((1:length(freq_names)) - bw/2 + (2*j-1) * bw / 4, bm(:, j), bs(:, j), 'k.');
    end
    xticks(1:length(freq_names));
    xticklabels(freq_names);
    legend('UP-1', 'UP-2');
    ylabel('Log Average power');
    title('UP-1 vs UP-2 power in frequency bands');
    saveas(gcf, fullfile(out_dir, 'up_power_bands_average.fig'));
    close;

    datasets2 = strrep(datasets, '_', '-');
    figure;
    bar(categorical(datasets2), d_all);
    xlabel('Dataset');
    ylabel('Cluster distance (larger = better)');
    title('Cluster distances across datasets')
    saveas(gcf, fullfile(out_dir, 'cluster_distance_datasets.fig'));
    close;

    figure;
    hold all;
    bar(mean(d_all));
    errorbar(mean(d_all), std(d_all)/sqrt(length(datasets)), 'k.');
    title('Average cluster distance');
    ylabel('Cluster Distance (larger = better)');
    saveas(gcf, fullfile(out_dir, 'cluster_distance_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(fr_all));
    errorbar(mean(fr_all), std(fr_all)/sqrt(length(datasets)), 'k.');
    title('Average firing rate');
    ylabel('Firing Rate');
    xticks([1, 2, 3]);
    xticklabels({'DOWN', 'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'fr_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(dur_all));xlabel('seconds')
    errorbar(mean(dur_all), std(dur_all)/sqrt(length(datasets)), 'k.');
    title('Average duration');
    ylabel('Duration (s)');
    xticks([1, 2, 3]);
    xticklabels({'DOWN', 'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'dur_average.fig'));
    close;

    for i = 1:length(z_thresh)
        tmp = replay_count_z(:, :, i) ./ sum(replay_count_z(:, :, i), 2);
        figure;
        hold all;
        bar(nanmean(tmp));
        errorbar(nanmean(tmp), nanstd(tmp)/sqrt(length(datasets)), 'k.');
        title(sprintf('Average replay count z=%d', z_thresh(i)));
        ylabel('Replay count');
        xticks([1, 2]);
        xticklabels({'UP-1', 'UP-2'});
        saveas(gcf, fullfile(out_dir, sprintf('replay_count_average_z%d.fig', z_thresh(i))));
        close;
    end

    figure;
    hold all;
    bar(mean(start_state));
    errorbar(mean(start_state), std(start_state)./sqrt(length(datasets)), 'k.');
    title('Average starting state proportion');
    ylabel('Probability');
    yticks([0, 0.5, 1]);
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'start_state_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(finish_state));
    errorbar(mean(finish_state), std(finish_state)./sqrt(length(datasets)), 'k.');
    title('Average ending state proportion');
    ylabel('Probability');
    yticks([0, 0.5, 1]);
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'finish_state_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(start_state_replay));
    errorbar(mean(start_state_replay), std(start_state_replay)./sqrt(length(datasets)), 'k.');
    title('Average starting state replay only proportion');
    ylabel('Probability');
    yticks([0, 0.5, 1]);
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'start_state_replay_average.fig'));
    close;

    figure;
    hold all;
    bar(mean(finish_state_replay));
    errorbar(mean(finish_state_replay), std(finish_state_replay)./sqrt(length(datasets)), 'k.');
    title('Average ending state replay only proportion');
    ylabel('Probability');
    yticks([0, 0.5, 1]);
    xticks([1, 2]);
    xticklabels({'UP-1', 'UP-2'});
    saveas(gcf, fullfile(out_dir, 'finish_state_replay_average.fig'));
    close;

    tr_lfp_m = cellfun(@(x)(mean(x, 2)'), tr_lfp_all, 'UniformOutput', false);
    tr_lfp_m2 = cell(1, 6);
    for i = 1:6
        tr_lfp_m2{i} = mean(cell2mat(tr_lfp_m(:, i)));
    end
    figure;
    for i = 1:6
        subplot(3, 2, i);
        plot(linspace(-2, 2, 2000), tr_lfp_m2{i}); title(transition_names{i});
    end
    saveas(gcf, fullfile(out_dir, 'tr_lfp_average.fig'));
    
    tr_len = cell2mat(cellfun(@(x)(size(x, 2)), tr_lfp_all, 'UniformOutput', false))
    tr_prc = 100 .* (tr_len ./ sum(tr_len, 2))
    figure;
    hold all;
    bar(mean(tr_prc));
    errorbar(mean(tr_prc), std(tr_prc)./sqrt(length(datasets)), 'k.')
    xticks(1:6)
    xticklabels(transition_names)
    xlim([0.5, 6.5])
    ylabel('%')
    saveas(gcf, fullfile(out_dir, 'transition_percent_average.fig'));
    close;
    
    up_type_prc = 100 .* (up_types_all ./ sum(up_types_all, 2));
    figure;
    hold all;
    bar(mean(up_type_prc));
    errorbar(mean(up_type_prc), std(up_type_prc)./sqrt(length(datasets)), 'k.');
    xticks(1:5);
    xticklabels({'UP-1 only', 'UP-2 only', 'UP-1 -> UP-2', 'UP-2 -> UP-1', 'Multiple'});
    xlim([0.5, 5.5]);
    ylabel('%');
    saveas(gcf, fullfile(out_dir, 'up_type_percent_average.fig'));
    close;
end

save(fullfile(out_dir, 'combined_stats.mat'), 'tau_all', 'tau_comp_all', 'pow_all', 'd_all', 'fr_all', 'fr_sem', 'dur_all', 'dur_sem', 'tau_task_all', 'replay_count_z', 'start_state', 'start_state_replay', 'finish_state', 'finish_state_replay', 'tr_lfp_all', 'up_types_all', 'up_bin_dur_all', 'up_bin_num_all');