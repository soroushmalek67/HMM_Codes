function job = compute_hmm(dataset, sleep, threshold, out_dir, plot_figs, num_states)
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
        s = parcluster;
        job = createCommunicatingJob(s, 'Type', 'pool');
        job.NumWorkersRange = [10, 10];
        createTask(job, @sleepHMMTraining_new, 4, {train_all, num_states, num_neurons, 10, 200});
        submit(job);
%         [hmm, emis, trans, l] = sleepHMMTraining_new(train_all, num_states, num_neurons, 10, 200);
    else
        tmp = load(hmm_file, 'hmm', 'emis', 'trans', 'l');
        job = struct();
        job.State = 'finished';
        task = struct();
        task.OutputArguments = {tmp.hmm, tmp.emis, tmp.trans, tmp.l};
        job.Tasks = task;
    end
end