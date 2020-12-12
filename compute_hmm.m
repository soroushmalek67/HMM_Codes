function [ts_full, fr, dur, epochs, fr_n, bin, hmm] = compute_hmm(dataset, sleep, threshold, out_dir, plot_figs, num_states)
    if ~exist('num_states', 'var') || isempty(num_states)
        num_states = 3;
    end
	if ~exist('plot_figs', 'var') || isempty(plot_figs)
		plot_figs = true;
	end

	if sleep == 1 && ~isempty(strfind(dataset, '8482_15p'))
		load(fullfile(out_dir, sprintf('sleep%d_HMM_data.mat', sleep)), 'times', 'ids');
	else
		load(fullfile(dataset, sprintf('sleep%d_HMM_data.mat', sleep)), 'times', 'ids');
	end
	load(fullfile(dataset, 'ts.mat'));

	num_neurons = max(ids);

	if ~exist('e', 'var')
		e = struct();
		e.epochs = epochs;
	end
	if ~isfield(e.epochs, 'rest')
		if exist(fullfile(dataset, 'epochs.mat'), 'file')
			load(fullfile(dataset, 'epochs.mat'), 'rest');
			e.epochs.rest = rest;
		else
			load(fullfile(dataset, 'ts_full.mat'), 'e');
		end
	end
	epochs = poprate_epochs_kernel(times, ids, e, threshold, out_dir, plot_figs, sleep);
	num_epochs = size(epochs, 1);
	obs = sleepOBS(times(times > 0), ids(times > 0), max(times));

	train_all = cell(size(epochs, 1), 1);
	for i = 1:size(epochs, 1)
		train_all{i} = obs(epochs(i, 1):epochs(i, 2));
	end

	hmm_file = fullfile(out_dir, sprintf('hmm_%dstates.mat', num_states));
	if ~exist(hmm_file, 'file')
		[hmm, emis, trans, l] = sleepHMMTraining_new(train_all, num_states, num_neurons, 10, 200);
	else
		tmp = load(hmm_file, 'hmm', 'emis', 'trans', 'l');
		hmm = tmp.hmm;
		emis = tmp.emis;
		trans = tmp.trans;
		l = tmp.l;
    end

	S = cell(num_epochs, 1);
	ts = cell(num_states, 1);
	fr = cell(num_states, 1);
	fr_n = cell(num_states, 1);
    fr_ep = cell(num_states, size(epochs, 1));
	bin = cell(num_states, 1);
	dur = cell(num_states, 1);
	q = full(sparse(ids(times > 0), times(times > 0), 1, num_neurons, max(times)));
	for i = 1:size(epochs, 1)
		S{i} = hmmviterbi(train_all{i}, hmm.TRANSITION, hmm.EMISSION);
		ix = times > epochs(i, 1) & times < epochs(i, 2);
		tt = times(ix) - epochs(i, 1);
		n = ids(ix);
		if plot_figs
			figure;
			hold all;
			plot(tt, n/num_neurons, 'k.');
			plot(1.1*S{i} - 1.2);
			ylim([0, 2.5]);
			saveas(gcf, fullfile(out_dir, sprintf('State_Transition_Epoch%d.fig', i)));
			close;
		end
		for j = 1:num_states
			tmp_bin = S{i} == j;
			tmpd = diff([0, tmp_bin, 0]);
			st = find(tmpd == 1);
			en = find(tmpd == -1);
			tmp_ts = [st', en'] + epochs(i, 1);
			tmp_dur = diff(tmp_ts, 1, 2) ./ 1000;
			tmp_fr = zeros(size(tmp_ts, 1), 1);
			tmp_fr_n = zeros(num_neurons, size(tmp_ts, 1));
			tmp_bin = zeros(num_neurons, size(tmp_ts, 1));
			for k = 1:size(tmp_ts, 1)
				curr_ix = find(times >= tmp_ts(k, 1) & times <= tmp_ts(k, 2));
				tmp_fr(k) = length(curr_ix);
				% tmp_obs = obs(tmp_ts(k, 1):tmp_ts(k, 2));
				% tmp_bin(tmp_obs, k) = 1;
				% for l = 1:length(tmp_obs)
				% 	tmp_fr_n(tmp_obs(l), k) = tmp_fr_n(tmp_obs(l), k) + 1;
				% end
				curr_q = q(:, tmp_ts(k, 1):tmp_ts(k, 2));
				tmp_fr_n(:, k) = sum(curr_q, 2);
			end
			tmp_bin = double(tmp_fr_n > 0);
			tmp_fr_n = tmp_fr_n ./ tmp_dur';
			tmp_fr_n(:, tmp_dur == 0) = 0;
			tmp_fr = tmp_fr ./ (num_neurons .* tmp_dur);
			tmp_fr(tmp_dur == 0) = 0;
			ts{j} = [ts{j}; tmp_ts];
			fr{j} = [fr{j}; tmp_fr];
            fr_ep{j,i} = tmp_fr;
			fr_n{j} = [fr_n{j}, tmp_fr_n];
			bin{j} = [bin{j}, tmp_bin];
			dur{j} = [dur{j}; tmp_dur];
		end
    end
    
	sleep_ep = e.epochs.(sprintf('sleep%d', sleep));
	ts_full = cellfun(@(x)((x.*10)+sleep_ep(1)), ts, 'UniformOutput', false);
	save(hmm_file, 'hmm', 'emis', 'trans', 'l', 'S', 'ts', 'dur', 'fr', 'fr_ep', 'ts_full', 'fr_n', 'bin', 'train_all');
end