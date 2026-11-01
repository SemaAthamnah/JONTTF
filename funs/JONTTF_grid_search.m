function [Results,best] = JONTTF_grid_search(A_cell,k,GT,lambda1_values,lambda2_values,p_values,opts)
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University
%
% Evaluate averaged U with ClusteringMeasure.
% Return Results and best in memory; create no files or directories.
% GT contains labels starting at 1. Labels are used without an offset.

if nargin < 7 || isempty(opts), opts = struct; end
assert(isstruct(opts) && isscalar(opts),'opts must be a scalar struct.');
validateattributes(lambda1_values,{'numeric'}, ...
    {'real','finite','vector','nonempty','positive'});
validateattributes(lambda2_values,{'numeric'}, ...
    {'real','finite','vector','nonempty','nonnegative'});
validateattributes(p_values,{'numeric'}, ...
    {'real','finite','vector','nonempty','>',0,'<=',1});
assert(exist('JONTTF','file') ~= 0, ...
    'JONTTF.m must be on the MATLAB path.');
assert(exist('ClusteringMeasure','file') ~= 0, ...
    'ClusteringMeasure.m must be on the MATLAB path.');

% One adjacency matrix per cell; pass this cell array to JONTTF.
assert(iscell(A_cell) && ~isempty(A_cell), ...
    'A_cell must be a nonempty cell array of square matrices.');
n = size(A_cell{1},1);
L = numel(A_cell);
for l = 1:L
    validateattributes(A_cell{l},{'numeric'}, ...
        {'real','finite','nonempty','2d','square'});
    assert(isequal(size(A_cell{l}),[n,n]), ...
        'All layers must have the same size.');
    A_cell{l} = full(double(A_cell{l}));
end
validateattributes(k,{'numeric'}, ...
    {'real','finite','scalar','integer','positive','<=',n});

% Accept either an n-element vector or the supplied 1 x 1 x n format.
if isvector(GT)
    GT = reshape(GT,1,1,[]);
end
assert(isnumeric(GT) && size(GT,1)==1 && ...
    size(GT,2)==1 && numel(GT)==n, ...
    'GT must be an n-element vector or a 1 x 1 x n array.');
Y = double(squeeze(GT(1,1,:)));
Y = Y(:);
validateattributes(Y,{'numeric'}, ...
    {'real','finite','integer','positive','vector','numel',n});

solverOpts = opts;
% Ground truth is used only for evaluation, never for fitting.
for field = {'outputDir','GT','GT_U'}
    if isfield(solverOpts,field{1})
        solverOpts = rmfield(solverOpts,field{1});
    end
end
[g1,g2,gp] = ndgrid(lambda1_values(:),lambda2_values(:),p_values(:));
N = numel(g1);
metricNames = {'ACC','NMI','Purity','F','P','R','RI'};
Results = table((1:N)',g1(:),g2(:),gp(:), ...
    'VariableNames',{'combination','lambda1','lambda2','p'});
for j = 1:numel(metricNames)
    Results.(metricNames{j}) = nan(N,1);
end
Results.NMI_barU = nan(N,1); % Compatibility: NMI of averaged U.
Results.seconds = nan(N,1);
Results.runStatus = repmat({'pending'},N,1);
Results.errorMessage = repmat({''},N,1);
best = struct('NMI_barU',-Inf,'combination',NaN);
fprintf('JONTTF grid: %d combinations.\n',N);
parameterNames = {'lambda1','lambda2','p'};
parameterValues = {lambda1_values,lambda2_values,p_values};
for j = 1:numel(parameterNames)
    values = unique(parameterValues{j});
    if numel(values)==1
        fprintf('Fixed %s = %g\n',parameterNames{j},values);
    end
end

for i = 1:N
    timer = tic;
    fprintf('[%d/%d] lambda1=%g lambda2=%g p=%g\n', ...
        i,N,g1(i),g2(i),gp(i));
    U = []; V = []; S = []; Sbar = []; info = [];
    result = []; Y_pre = [];
    parameters = struct('combination',i,'lambda1',g1(i), ...
        'lambda2',g2(i),'p',gp(i),'k',k,'solverOpts',solverOpts);
    try
        [U,V,S,Sbar,info] = JONTTF( ...
            A_cell,k,g1(i),g2(i),gp(i),solverOpts);
        factors = {U,V,S,Sbar};
        for j = 1:numel(factors)
            factor = factors{j};
            if iscell(factor), factor = cat(3,factor{:}); end
            assert(isnumeric(factor) && ~isempty(factor) && ...
                all(isfinite(factor(:))),'Core returned invalid factors.');
        end

        % Use the requested evaluation exactly for cell-array U.
        if iscell(U)
            U_tensor = cat(3,U{:,:});
            assert(isequal(size(U_tensor,1),n) && ...
                size(U_tensor,2)==k && size(U_tensor,3)==L, ...
                'U must represent an n x k x L tensor.');
            [~, Y_pre] = max(mean(cat(3, U{:,:}),3), [], 2);
        else
            assert(size(U,1)==n && size(U,2)==k && size(U,3)==L, ...
                'U must represent an n x k x L tensor.');
            [~, Y_pre] = max(mean(U,3), [], 2);
        end
        result = ClusteringMeasure(Y, Y_pre);
        for j = 1:numel(metricNames)
            name = metricNames{j};
            assert(isstruct(result) && isscalar(result) && isfield(result,name), ...
                'ClusteringMeasure must return a struct containing %s.',name);
            validateattributes(result.(name),{'numeric'}, ...
                {'real','finite','scalar'});
            Results.(name)(i) = result.(name);
        end
        Results.NMI_barU(i) = result.NMI;
        Results.runStatus{i} = 'ok';
        fprintf('  ACC=%.4f | NMI=%.4f | Purity=%.4f | F=%.4f | P=%.4f | R=%.4f | RI=%.4f\n', ...
            result.ACC,result.NMI,result.Purity,result.F,result.P,result.R,result.RI);
    catch ME
        Results.runStatus{i} = 'failed';
        Results.errorMessage{i} = getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'  FAILED: %s\n',ME.message);
    end
    Results.seconds(i) = toc(timer);
    if strcmp(Results.runStatus{i},'ok') && result.NMI > best.NMI_barU
        best = struct('combination',i,'lambda1',g1(i),'lambda2',g2(i), ...
            'p',gp(i),'NMI_barU',result.NMI,'result',result, ...
            'Y_pre',Y_pre,'U',{U},'V',{V},'S',{S},'Sbar',{Sbar}, ...
            'info',info,'parameters',parameters,'Y',Y, ...
            'seconds',Results.seconds(i));
    end
end

if isfinite(best.NMI_barU)
    fprintf('Best: lambda1=%g lambda2=%g p=%g | NMI=%.4f\n', ...
        best.lambda1,best.lambda2,best.p,best.NMI_barU);
    disp(best.result);
else
    warning('JONTTF_grid_search:AllFailed', ...
        'All combinations failed. Inspect Results.errorMessage.');
end
end
