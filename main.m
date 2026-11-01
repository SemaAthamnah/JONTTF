% JONTTF: Best NMI across SSBM noise levels
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University

clear; clc;
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'funs'));

lambda1_values = [0.01,1000];
lambda2_values = 1;
p_values = 1;
sampleIndex = 1;
datasets = {'NB'}; %'NW'
opts = struct('maxIter', 1000, 'tol', 1e-5, ...
    'tolChange', 1e-6, 'verbose', false);

for d = 1:numel(datasets)
    dataset = datasets{d};
    dataFile = fullfile(projectDir, 'data', ['NTF_SSBM_' dataset '_1s.mat']);
    assert(isfile(dataFile), 'Missing dataset: %s', dataFile);

    fileInfo = whos('-file', dataFile);
    variableNames = {fileInfo.name};
    matches = regexp(variableNames, '^SSBM_(\d+)$', 'tokens', 'once');
    valid = ~cellfun('isempty', matches);
    settingIds = sort(cellfun(@(x) str2double(x{1}), matches(valid)));
    assert(~isempty(settingIds), 'No SSBM_i variables found in %s.', dataFile);

    allResults = table();
    fprintf('\nDataset %s: %d noise levels\n', dataset, numel(settingIds));
    for s = 1:numel(settingIds)
        settingId = settingIds(s);
        adjacencyName = sprintf('SSBM_%d', settingId);
        labelName = sprintf('GT_%d', settingId);
        assert(ismember(labelName, variableNames), 'Missing %s in %s.', labelName, dataFile);
        loaded = load(dataFile, adjacencyName, labelName);
        SSBM = loaded.(adjacencyName);
        GT = loaded.(labelName);

        assert(ndims(SSBM) == 4, 'Expected samples x layers x nodes x nodes.');
        assert(size(SSBM,1) >= sampleIndex, 'Sample index exceeds available samples.');
        numLayers = size(SSBM,2);
        numNodes = size(SSBM,3);
        assert(size(SSBM,4) == numNodes, 'Adjacency matrices must be square.');
        assert(size(GT,1) >= sampleIndex && size(GT,2) == numLayers && ...
            size(GT,3) == numNodes, 'Unexpected GT dimensions.');

        A_cell = cell(1,numLayers);
        for l = 1:numLayers
            A_cell{l} = double(reshape(SSBM(sampleIndex,l,:,:),numNodes,numNodes));
        end
        Y = double(reshape(GT(sampleIndex,1,:),[],1));
        if min(Y) == 0
            Y = Y + 1;
        end
        assert(all(isfinite(Y)) && all(Y >= 1) && all(Y == fix(Y)), ...
            'Ground-truth labels must be positive integers.');
        k = numel(unique(Y));

        fprintf('Noise level %d | sample %d | n=%d | L=%d | k=%d\n', ...
            settingId,sampleIndex,numNodes,numLayers,k);
        [~, best] = JONTTF_grid_search( ...
            A_cell,k,Y,lambda1_values,lambda2_values,p_values,opts);
        assert(isfinite(best.NMI_barU), ...
            'All parameter combinations failed at noise level %d.',settingId);
        allResults = [allResults; table(settingId, best.lambda1, ...
            best.lambda2, best.p, best.NMI_barU, ...
            'VariableNames', {'noiseLevel','lambda1','lambda2','p','NMI'})];
    end

    disp(allResults);
    figure('Name', ['JONTTF - ' dataset], 'Color', 'w');
    plot(allResults.noiseLevel, allResults.NMI, '-o', ...
        'LineWidth', 1.8, 'MarkerSize', 7);
    xlabel('Noise level');
    ylabel('NMI');
    title(['JONTTF best NMI on SSBM-' dataset]);
    ylim([0 1]);
    xticks(settingIds);
    grid on;
end
