function [X_cell] = frequency2time(Y_cell)
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University
Y_tensor = cat(3, Y_cell{:,:});
X = shiftdim(Y_tensor, 0);
X_hat = ifft(X, [], 3);
% X_hat=X;
X_shift = shiftdim(X_hat, 0);
[~, ~, num_V] = size(X_shift);
for v = 1:num_V
        X_cell{v} = X_shift(:,:,v);
end
end

