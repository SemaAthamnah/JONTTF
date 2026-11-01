function [Y_cell] = time2frequency(X_cell)
% Author: Sema Athamnah
% Email: athamnah@msu.edu
% Affiliation: Department of Biomedical Engineering
%              Michigan State University
X_tensor = cat(3, X_cell{:,:});  
Y = shiftdim(X_tensor, 0);    
Y_hat = fft(Y, [], 3);
% Y_hat = Y;
Y_shift = shiftdim(Y_hat, 0);  
[~, ~, num_V] = size(Y_shift);
for v = 1:num_V
        Y_cell{v} = Y_shift(:,:,v);
end
end

