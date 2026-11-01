function [U,V,S,Sbar,info] = JONTTF(A_cell,k,lambda1,lambda2,p,opts)
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University
%
if nargin < 6 || isempty(opts)
    opts = struct;
end

defaults = struct( ...
    'maxIter',1000,'tol',1e-5,'tolChange',1e-6, ...
    'mu1',1e-5,'mu2',1e-5,'mu3',1e-5,'mu4',1e-5,'mu5',1e-5, ...
    'eta',1.1,'maxPenalty',1e12,'verbose',true);

fields = fieldnames(defaults);
for i = 1:numel(fields)
    if ~isfield(opts,fields{i}) || isempty(opts.(fields{i}))
        opts.(fields{i}) = defaults.(fields{i});
    end
end

A = cat(3,A_cell{:});

n = size(A,1);
L = size(A,3);
sX = [n,k,L];
r = min(k,L);

Ap_cell=cell(L,1);
An_cell=cell(L,1);
for i=1:L
    Ap_cell{i}=max(A(:,:,i),0);
    An_cell{i}=max(-A(:,:,i),0);
end

Ap_hat = time2frequency(Ap_cell);
An_hat = time2frequency(An_cell);

U_hat = cell(1,L);
V_hat = cell(1,L);
S_hat = cell(1,L);
Sbar_hat = cell(1,L);

for l = 1:L
    U_hat{l} = eye(n,k);
    V_hat{l} = eye(n,k);
    S_hat{l} = eye(k);
    Sbar_hat{l} = zeros(k,k);
end
U = time2frequency(U_hat);
V = time2frequency(V_hat);
S = time2frequency(S_hat);
Sbar = time2frequency(Sbar_hat);

P = cell(1,L); 
Q = cell(1,L);
Y1 = cell(1,L); 
for l=1:L
    P{l}=zeros(n,k);
    Q{l}=zeros(n,k);
    Y1{l}=zeros(n,k);
end

M = P; J = Q; 
Y2 = Y1; Y3 = Y1; Y4 = Y1; Y5 = Y1;



mu1 = opts.mu1; mu2 = opts.mu2; mu3 = opts.mu3;
mu4 = opts.mu4; mu5 = opts.mu5;

info.residuals = zeros(opts.maxIter,5);
info.change = zeros(opts.maxIter,1);
info.penalties = zeros(opts.maxIter,5);
info.converged = false;
vec = @(C) reshape(cat(3,C{:}),[],1);

for iteration = 1:opts.maxIter
    previous = [vec(U);vec(V);vec(S);vec(Sbar); ...
            vec(Q);vec(P);vec(J);vec(M)];
    info.penalties(iteration,:) = [mu1,mu2,mu3,mu4,mu5];

    % update U
    Q_hat  = time2frequency(Q);
    J_hat  = time2frequency(J);
    Y1_hat = time2frequency(Y1);
    Y2_hat = time2frequency(Y2);
    Y5_hat = time2frequency(Y5);
    for l = 1:L
        B = 2*Ap_hat{l}*V_hat{l}*S_hat{l} ...
          + 2*lambda1*An_hat{l}*V_hat{l}*Sbar_hat{l} ...
          + (mu1*Q_hat{l}-Y1_hat{l} + mu2*J_hat{l}-Y2_hat{l} ...
          + mu5*V_hat{l}-Y5_hat{l});

        [Lambda,~,R] = svd(B,'econ');
        U_hat{l} = Lambda*R';
    end

    % update V
    P_hat  = time2frequency(P);
    M_hat  = time2frequency(M);
    Y3_hat = time2frequency(Y3);
    Y4_hat = time2frequency(Y4);
    for l = 1:L
        C = 2*Ap_hat{l}*U_hat{l}*S_hat{l} ...
          + 2*lambda1*An_hat{l}*U_hat{l}*Sbar_hat{l} ...
          + (mu3*P_hat{l}-Y3_hat{l} + mu4*M_hat{l}-Y4_hat{l} + mu5*U_hat{l}+Y5_hat{l});

        [Lambda,~,R] = svd(C,'econ');
        V_hat{l} = Lambda*R';
    end

    % update S and Sbar
    Zplus_hat = cell(1,L);
    Zminus_hat = cell(1,L);
    for l = 1:L
        Zplus_hat{l} = U_hat{l}'*Ap_hat{l}*V_hat{l};
        Zminus_hat{l} = U_hat{l}'*An_hat{l}*V_hat{l};
    end

    Zplus = frequency2time(Zplus_hat);
    Zminus = frequency2time(Zminus_hat);

    for l = 1:L
        S{l} = diag(max(diag(Zplus{l}),0));
        Z = Zminus{l};
        Z = max((Z+Z.')/2,0);
        Z(1:k+1:end) = 0;
        Sbar{l} = Z;
    end


    % update Q, P, J, and M
    U = frequency2time(U_hat);
    V = frequency2time(V_hat);

    Y1=frequency2time(Y1_hat);
    Y2=frequency2time(Y2_hat);
    Y3=frequency2time(Y3_hat);
    Y4=frequency2time(Y4_hat);

    for l=1:L
    Q{l} = max(U{l}+(Y1{l}/mu1),0);
    P{l} = max(V{l}+(Y3{l}/mu3),0);
    TJ{l} = U{l}+(Y2{l}/mu2);
    TM{l} = V{l}+(Y4{l}/mu4);
    end

    if lambda2 == 0
        J = TJ;
        M = TM;
    else
        rhoJ = (lambda2/(2*mu2))*ones(r,1);
        rhoM = (lambda2/(2*mu4))*ones(r,1);

        TJ_v = cat(3, TJ{:,:});
        TM_v = cat(3, TM{:,:});
        [jvec,~] = wshrinkObj_weight_lp(TJ_v(:),rhoJ,sX,0,3,p);
        [mvec,~] = wshrinkObj_weight_lp(TM_v(:),rhoM,sX,0,3,p);
        
        J_tensor = real(reshape(jvec,sX));
        M_tensor = real(reshape(mvec,sX));
        
        J = cell(1,L);
        M = cell(1,L);
        for l = 1:L
            J{l} = J_tensor(:,:,l);
            M{l} = M_tensor(:,:,l);
         end
    end

E1 = cellfun(@minus,U,Q,'UniformOutput',false);
E2 = cellfun(@minus,U,J,'UniformOutput',false);
E3 = cellfun(@minus,V,P,'UniformOutput',false);
E4 = cellfun(@minus,V,M,'UniformOutput',false);
E5 = cellfun(@minus,U,V,'UniformOutput',false);

Y1 = cellfun(@(y,e) y+mu1*e,Y1,E1,'UniformOutput',false);
Y2 = cellfun(@(y,e) y+mu2*e,Y2,E2,'UniformOutput',false);
Y3 = cellfun(@(y,e) y+mu3*e,Y3,E3,'UniformOutput',false);
Y4 = cellfun(@(y,e) y+mu4*e,Y4,E4,'UniformOutput',false);
Y5 = cellfun(@(y,e) y+mu5*e,Y5,E5,'UniformOutput',false);

    U_tensor = cat(3,U{:}); V_tensor = cat(3,V{:});
    scale = max([1,norm(U_tensor(:)),norm(V_tensor(:))]);
    
residuals = [norm(vec(E1)), ...
             norm(vec(E2)), ...
             norm(vec(E3)), ...
             norm(vec(E4)), ...
             norm(vec(E5))];

current = [vec(U);vec(V);vec(S);vec(Sbar); ...
           vec(Q);vec(P);vec(J);vec(M)];
change = norm(current-previous)/max(1,norm(previous));

    S_hat = time2frequency(S);
    Sbar_hat = time2frequency(Sbar);


  
    info.residuals(iteration,:) = residuals;
    info.change(iteration) = change;

multipliers = [vec(Y1);vec(Y2);vec(Y3);vec(Y4);vec(Y5)];


if opts.verbose && (iteration == 1 || mod(iteration,10) == 0)
    fprintf('Iter %4d | residual %.3e | U-V %.3e | change %.3e\n', ...
        iteration,max(residuals),residuals(5),change);
end

    if max(residuals) < opts.tol && change < opts.tolChange
        info.converged = true;
        break;
    end

    if iteration < opts.maxIter
        mu1 = min(opts.eta*mu1,opts.maxPenalty);
        mu2 = min(opts.eta*mu2,opts.maxPenalty);
        mu3 = min(opts.eta*mu3,opts.maxPenalty);
        mu4 = min(opts.eta*mu4,opts.maxPenalty);
        mu5 = min(opts.eta*mu5,opts.maxPenalty);
    end
end

info.iterations = iteration;
info.residuals = info.residuals(1:iteration,:);
info.residual_names = {'U-Q','U-J','V-P','V-M','U-V'};
info.change = info.change(1:iteration);
info.penalties = info.penalties(1:iteration,:);
info.penalty_names = {'mu1','mu2','mu3','mu4','mu5'};
info.final_penalties = [mu1,mu2,mu3,mu4,mu5];
info.initialization = 'SVD of summed layers; U and V initialized identically';
info.options = opts;

if info.converged
    info.status = 'converged';
else
    info.status = 'maximum iterations reached';
end

if opts.verbose
    fprintf('JONTTF: %s (%d iterations).\n',info.status,iteration);
end

end
