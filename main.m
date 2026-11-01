% JONTTF: Single-sample SSBM experiment with grid search
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University

clear; clc;
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'funs'));

% Choose 'NB' for inter-community noise or 'NW' for intra-community noise.
dataset = 'NB';
dataFile = fullfile(projectDir, 'data', ['NTF_SSBM_' dataset '_1s.mat']);
assert(isfile(dataFile), 'Missing dataset: %s', dataFile);

% Select one available noise setting and the first sample.
fileInfo = whos('-file', dataFile);
variableNames = {fileInfo.name};
match = regexp(variableNames, '^SSBM_(\d+)$', 'tokens', 'once');
valid = ~cellfun('isempty', match);
settingIds = sort(cellfun(@(x) str2double(x{1}), match(valid)));
assert(~isempty(settingIds), 'No SSBM_i variables found in %s.', dataFile);
settingId = settingIds(1);
loaded = load(dataFile, sprintf('SSBM_%d', settingId), sprintf('GT_%d', settingId));
SSBM = loaded.(sprintf('SSBM_%d', settingId));
GT = loaded.(sprintf('GT_%d', settingId));

% Expected shapes: SSBM = samples x layers x nodes x nodes;
% GT = samples x layers x nodes.
assert(ndims(SSBM) == 4, 'Expected SSBM shape: samples x layers x nodes x nodes.');
numLayers = size(SSBM, 2);
numNodes = size(SSBM, 3);
assert(size(SSBM, 4) == numNodes, 'Adjacency matrices must be square.');
assert(size(GT, 1) >= 1 && size(GT, 2) == numLayers && ...
    size(GT, 3) == numNodes, 'Unexpected GT dimensions.');
A_cell = cell(1, numLayers);
for l = 1:numLayers
    A_cell{l} = reshape(SSBM(1,l,:,:), numNodes, numNodes);
end

% The SSBM benchmark stores zero-based labels; use positive labels for evaluation.
Y = double(reshape(GT(1,1,:), [], 1));
if min(Y) == 0
    Y = Y + 1;
end
assert(all(Y >= 1) && all(Y == fix(Y)), 'Ground truth must contain positive integer labels.');
k = numel(unique(Y));

% Adjust the search ranges as needed.
lambda1_values = [0.01 0.1 1 10 100 1000];
lambda2_values = 1;
p_values = 1;
opts = struct('maxIter', 1000, 'tol', 1e-5, ...
    'tolChange', 1e-6, 'verbose', false);

fprintf('Dataset: %s | Setting: %d | Sample: 1 | Layers: %d | Nodes: %d | k: %d\n', ...
    dataset, settingId, numLayers, numNodes, k);
[Results, best] = JONTTF_grid_search( ...
    A_cell, k, Y, lambda1_values, lambda2_values, p_values, opts);
disp(Results);
if isfinite(best.NMI_barU)
    fprintf('Best lambda1=%g, lambda2=%g, p=%g, NMI=%.4f\n', ...
        best.lambda1, best.lambda2, best.p, best.NMI_barU);
else
    warning('All grid-search combinations failed. Inspect Results.errorMessage.');
end
