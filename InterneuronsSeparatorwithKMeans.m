clear all;

Dataset = '8482_16p' ;
iteration = 10 ;
tfiles_dir = sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus', Dataset);

%%%% addpath(genpath('C:\Users\saeedeh\Documents\MATLAB\Ramp_task'))


%% Main
[Interneurons,Pyramidals] = deal(cell(12,1));

% % % % % % % % % %     Root = [sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus\\ClusterSummaryData', Dataset)];
% % % % % % % % % %     cd(Root)
% % % % % % % % % %     files = dir(Root); names = {files.name};
% % % % % % % % % %     data = []; testSet = [];cte = 0;
% % % % % % % % % % %     data is a matrix of size 3 by 4 times nCells in which
% % % % % % % % % % % first row is PeakHalfWidthPts and second row is PeakToTroughPts
% % % % % % % % % % % while third row is showing cell number (same order as in spikes ans trial_events)
% % % % % % % % % % % Since there are four electrodes per tetrode, we have 4 waveform per cell
% % % % % % % % % %     for i = 1:length(names)
% % % % % % % % % %         if strcmp(names{i}(1),'T')
% % % % % % % % % %             cte = cte+1;load(names{i})
% % % % % % % % % %             data = [data [CI.PeakHalfWidthPts;CI.PeakToTroughPts;cte*ones(1,4)]];
% % % % % % % % % %         end
% % % % % % % % % % %         clear CI
% % % % % % % % % %     end
   

 %% My Selection Part %%%%
    [Interneurons,Pyramidals] = deal(cell(12,1));
    
    data = []; testSet = [];cte = 0;
    Root = [sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus\\ClusterSummaryData', Dataset)];
    cd(Root)
    files = dir(Root); names = {files.name};
    namess = FindFiles('*.mat');
    
    for ii = 1 : length(namess)
        cte = cte +1 ;
        load(namess{ii})
        data = [data [CI.PeakHalfWidthPts;CI.PeakToTroughPts;cte*ones(1,4)]];
    end




%% Clustering and Plot
    cd(Root); 
    Dataa = data';
    rng default 
    figure;
    Dataa(:,1:2) = Dataa(:,1:2)/32;   % converting the unit to ms
    h1 = plot(Dataa(:,1),Dataa(:,2),'.k','MarkerSize',6); hold on
    tic;
    opts = statset('Display', 'final');
    [idx,C] = kmeans(Dataa(:,1:2),2, 'Replicates',iteration, 'Options', opts);
    toc
% % % % % % %     [x,y] = ginput(2);
% % % % % % %     h2 = plot(x,y,'-r','LineWidth',2);
% % % % % % %     b = y(1)+(diff(y)/diff(x))*(Dataa(1,:)-x(1));

    h3 = plot(Dataa(idx==2,1),Dataa(idx==2,2),'or');hold on
    h4 = plot(Dataa(idx==1,1),Dataa(idx==1,2),'og');
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
    title([sprintf('Classification of cells for %s', Dataset )]);
% % %     legend('raw Dataa','classes border','Interneurons','Pyramidals');
    legend('raw Dataa', 'Cluster 1','Cluster 2','Centroids','Location','NW');
    xlabel('PeakHalfWidthPts');
    ylabel('PeakToTroughPts')
    cd(tfiles_dir)
    saveas(gcf, sprintf('NeuronTypeKmeaniteration %d', iteration));
    
    group1 = Dataa(idx==2,3); % Interneurons
    group2 = Dataa(idx==1,3); % Pyramydal neurons
    
    % NOTE: we have 4 data points for each neuron
    % Now we need to decide the cell type for those shared in two groups
     if length(group2) > length(group1)
    Interneurons = unique(group1);
    Pyramidals = unique(group2);counts2 = histc(group2(:), Pyramidals);
     else 
         Interneurons = unique(group2);
    Pyramidals = unique(group1);counts2 = histc(group1(:), Pyramidals);
     end
    shared12 = Interneurons(ismember(Interneurons,Pyramidals));   % shared cells
    % Remove a cell from corresponding group if occured less than two data points(out of four)
    % For those which have two data points in each group -> put it in pyramidals
    remove1 = shared12(counts2(ismember(Pyramidals,shared12))'>=2);
    remove2 = shared12(counts2(ismember(Pyramidals,shared12))'<2);

    Interneurons(ismember(Interneurons,remove1)) = [];
    Pyramidals(ismember(Pyramidals,remove2)) = [];
%     close all

%% Save
cd(sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus', Dataset))
save('NeuronsTypeKmeans','Interneurons','Pyramidals')


%% Separate Neurons to the related folders
%% Separate InterNeurons to the InterNeuron folders

cd(tfiles_dir);

load('NeuronsTypeKmeans.mat');

Rootoftfiles = [sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus\\tfiles', Dataset)];
cd(Rootoftfiles)
files = dir(Rootoftfiles);
names = {files.name};
namess = FindFiles('*.t');

cd(tfiles_dir);

for jj = 1:length(namess)
           for jjj = 1:length(Interneurons)
               InterNeuronNum = Interneurons(jjj);
                if isequal(jj,InterNeuronNum)
                   
                      if ~exist('FinalInterneuronsKmeans')
                          mkdir FinalInterneuronsKmeans;
                      end
                      
                      copyfile(namess{jj}, 'FinalInterneuronsKmeans');
                      
                end    
           end

end

%% Separate Pyramidals to the Pyramidal folders


% % % % % Rootoftfiles = [sprintf('E:\\HMM - UP&Down\\Soroush\\Data\\%s\\Hippocampus\\tfiles', Dataset)];
% % % % % cd(Rootoftfiles)
% % % % % files = dir(Rootoftfiles);
% % % % % names = {files.name};

% % % % % cd(tfiles_dir);

for jj = 1:length(namess)
           for jjj = 1:length(Pyramidals)
               PyramidalNumbers = Pyramidals(jjj);
                if isequal(jj,PyramidalNumbers)
                       
                      if ~exist ('FinalPyramidalsKmeans')
                          mkdir FinalPyramidalsKmeans;
                      end
                      
                      copyfile(namess{jj}, 'FinalPyramidalsKmeans');
                      
                end    
           end

end

