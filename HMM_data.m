% SET DIRECTORIES AND SELECT EPOCH
Dataset = '7165_11p'

tfiles_dir = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\New_Result\\%s\\tfiles', Dataset);
epochs_dir = sprintf('E:/HMM - UP&Down/Soroush/Data/New_Result/%s', Dataset);
% select which part of the session you want to analyze (1 = sleep1; 2 =
% sleep2; 3 = sleep3)
Epoch_select = 3;

% make Q matrix from tfiles
cd(tfiles_dir)
tfile = FindFiles('*.t');
S = LoadSpikes(tfile);
Q = MakeQfromS(S,10);

% load sleep times
cd(epochs_dir)
load('ts.mat')

switch Epoch_select
    case 1
        sleep = e.epochs.sleep1;
        savefile = 'sleep1_HMM_data.mat';
    case 2
        sleep = e.epochs.sleep2;
        savefile = 'sleep2_HMM_data.mat';
    case 3
        sleep = e.epochs.sleep3;
        savefile = 'sleep3_HMM_data.mat';
    otherwise
        error('Wrong Value for Epoch_select')
end

% restrict Q matrix to sleep
Q_cut = Restrict(Q,sleep(1),sleep(2));
QD_cut = Data(Q_cut);
QD_range = Range(Q_cut);

% extract neuron IDs and spike times from Q matrix
[row col] = find(QD_cut==1);
us_ids = col;
us_times = QD_range(row);
us_times = (us_times-sleep(1))/10;
us_times = floor(us_times);

% sort times and neuron IDs
[times,ix] = sort(us_times);
ids = us_ids(ix);

 save(savefile,'ids','times')