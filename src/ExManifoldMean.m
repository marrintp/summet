function [MU, eigenvals, tellapsed] = ExManifoldMean(data,check)
% SYNTAX:   [MU, eigenvals, tellapsed] = ExManifoldMean(data,check)
%
% INPUTS:  'data' is the collection of subspaces you wish to find a
%           representative for.  It should be entered as a num_pts by 1
%           cell array with each cell being an n by q_i matrix, where q_i 
%           CANNOT vary.  In other words, the datacubes need to be unrolled
%           already when using this function.
%
%           'check' determines whether or not the algorithm will verify
%           that your data points are in fact orthonormal matrices.  If
%           check == 1, it does the check, otherwise it does not.
%
% OUTPUTS:  'MU' is the extrinsic manifold mean, returned as a matrix that
%           corresponds to the dominant q_i-dimensional eigenvectors. 
%
%           'eigenvals' is a vector of the dominant q_i ordered eigenvalues 
%           of A = sum(Ui*Ui').
%
%           'tellapsed' is the time it took to create the mean.
%
% NOTES:    Algorithm implemented from description in:
%
%           A. Srivastava and E. Klassen. "Monte Carlo extrinsic estimators
%           of manifold valued-parameters." IEEE Transactions on Signal
%           Processing, 50(2):299-308, 2002.
%
% LAST EDITED: 06/21/14 by Tim Marrinan at Colorado State University.
%
% ---------------------------------------------------------------------
%
% Copyright (c) 2014 Timothy P. Marrinan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
%    1. Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%  
%    2. Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the 
%       distribution.
%  
%    3. Neither name of copyright holders nor the names of its contributors
%       may be used to endorse or promote products derived from this 
%       software without specific prior written permission.
%  
%  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% ---------------------------------------------------------------------

MU = 0;
tellapsed = 0;
etol = 10e-10;
num_pts = size(data,1);
numvids = num_pts;
[~, ncz] = size(data{1});
if check == 1
    for i = 1:num_pts
        Y = squeeze(data{i});
        M = Y'*Y;
        II = eye(ncz);
    
        if norm(M-II,'fro') >= etol % find closest o.n. matrix
            [U, ~, V] = svd(Y,0);
            data{i} = U*V';
        end
    end
end

%% Compute the extrinsic manifold mean
tic;
tstart = tic;

% Take the outerproduct of each subspace with itself.
outer_products = cell(numvids,1);
for j = 1 : numvids
    outer_products{j} = data{j}*data{j}';
end

% Sum outer product matrices to create A
A = outer_products{1};
if numvids > 1
    for k = 2 : numvids
        A = A + outer_products{k};
    end
end

[U, L] = eig(A);
L = abs(diag(L));

[eigenvals, I] = sort(L,'descend');
Usort = U;
for hh = 1 : length(I)
    Usort(:,hh) = U(:,I(hh));
end
MU = Usort(:,1:ncz);
tellapsed = toc(tstart);
