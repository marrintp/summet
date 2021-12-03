function [MU, eigenvals, tellapsed] = FlagMean(DIM, type, option, data)
% SYNTAX:   [MU, eigenvals, tellapsed] = FlagMean(DIM, type, option, data)
%
% INPUTS:   'DIM'  is a vector specifying exactly which eigenvectors from
%           created flag you would like the function to return. Can also
%           specify [ALL] to return the full flag.
%
%           'type' is a string that is either [EIG] or  [SVD], depending on 
%           whether you wish the mean to be created from ordered
%           eigenvectors of A = sum(Ui*Ui') or left singular vectors of 
%           A = [U1 U2 U3 ...].
%
%           'option' is a string that determines how many vectors from each
%           constituent subspace you would like to use in the creation of
%           the mean.  Options are:
%           [SVD] for all the left singular vectors of the thin SVD
%           [NRG] for the singular vectors that contain 95% of the variance
%           [QR] for the thin QR factorization
%           [CLOSEST] for the nearest point on the Grassmannian manifold
%           Gr(n,m) where the input data is an m by n matrix.
%           [AsIs] for the matrices pulled directly from 'data'. Use
%           this option if you have already set the subspaces up as
%           desired.
%
%           'data' is the collection of subspaces you wish to find a
%           representative for.  It should be entered as a num_pts by 1 
%           cell array with each cell being an n by q_i matrix, where q_i 
%           can vary. In other words, the datacubes need to be unrolled 
%           already when using this function.
%
% OUTPUTS:  'MU' is the flag mean, returned as a matrix that
%           corresponds to the eigenvectors that you specified in 'DIM'. 
%
%           'eigenvals' is a vector of ALL the ordered eigenvalues of A.
%
%           'tellapsed' is the time it took to create the flag.
%
% NOTES:    This function calls the function 'subspacer' which IS included.
%
%           Algorithm implemented from description in:
%
%           B.Draper, M. Kirby, J. Marks, T. Marrinan, and C. Peterson. "A
%           flag representation for finite collections of subspaces of
%           mixed dimensions." Linear Algebra and its Applications,
%           451:15-32,2014.
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

if ~exist('type','var')
    type = 'SVD';
end

if ~exist('option','var')
    type = 'AsIs';
end

% Set up data
video_matrices = data;
numvids = length(video_matrices);

% Find an orthonormal basis for each video matrix using 'option'.
if ~strcmp(option,'AsIs')
    for ii = 1 : numvids
        P = video_matrices{ii};
        P = subspacer(P,option);
        video_matrices{ii} = P;
    end
end

% Depending on the method specified in 'type' find the vectors that will be
% used to create the mean subpace and the eigenvalues of A.
if strcmp(type,'EIG')
    tic;
    tstart = tic;
    
    % Take the outerproduct of each subspace with itself.
    outer_products = cell(numvids,1);
    for j = 1 : numvids
        outer_products{j} = video_matrices{j}*video_matrices{j}';
    end

    % Sum outer product matrices to create A
    A = outer_products{1};
    if numvids > 1
        for k = 2 : numvids
            A = A + outer_products{k};
        end
    end
    
    [U L] = eig(A);
    L = abs(diag(L));
    
    [eigenvals I] = sort(L,'descend');
    Usort = U;
    for hh = 1 : length(I)
        Usort(:,hh) = U(:,I(hh));
    end
    
    tellapsed = toc(tstart);
end
if strcmp(type,'SVD')
    % Possibly faster alternative to taking an outer product.  Instead, take
    % the SVD of a matrix of the concatenation of the original matrices.
    tic;
    tstart = tic;
    
    clear j
    A_prime =[];
    for j = 1 : numvids
        A_prime = [A_prime video_matrices{j}];
    end
    % Take the SVD of A_prime to get the eigenvectors of A.  They will
    % already be ordered.
    flip = size(A_prime);
    if flip(1) >= flip(2)
        [Usort, S, ~] =svd(A_prime,0);
        eigenvals = diag(S).*diag(S);
    else
        [~, S, Usort] = svd(A_prime',0);
        eigenvals = diag(S).*diag(S);
    end
    
    tellapsed = toc(tstart);
end

% Give the option to return ALL vectors of Usort
if isa(DIM, 'char')
    MU = Usort;
else   
    num_dim = length(DIM);
    MU = zeros(size(Usort,1),num_dim);
    cutoff = size(Usort);
    num_dim = min(num_dim, cutoff(2));
    for qq = 1 : num_dim
        MU(:,qq) = Usort(:,DIM(qq));
    end
end

%% Called function
function [U] = subspacer(matrix,option)
% SYNTAX:   [U] = subspacer(matrix,option)
%
% INPUTS:   'matrix' is an m by n matrix that you would like to turn into
%           a subspace.
%
%           'option' is the way you would like to create your subspace.
%           Options are:
%           [SVD] to use all the left singular vectors of the thin SVD.
%           [NRG] to use the left singular vectors that contain 95 percent
%           of the variance.
%           [QR] to use the thing QR factorization.
%           [CLOSEST] to use the nearest point on the Grassmannian manifold
%           Gr(n,m).
%
% OUTPUTS:  'U' is the subspace that has been created.
%
% NOTES:    
%
% LAST EDITED: 05/27/14 by Tim Marrinan at Colorado State University.
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

if ~exist('option','var')
    option = 'SVD';
end

if strcmp(option,'SVD')
    [U , ~, ~] = svd(matrix,0);
    U = U(:,1:rank(matrix));
end

if strcmp(option,'NRG')
    [U , S, ~] = svd(matrix,0);
    % turn singular values into eigenvalues
    S = S*S;
    S = diag(S);
    tot = sum(S);
    % find the number of eigenvectors needed to capture 95 percent of
    % variance
    counter = 1;
    nrg = 0;
    while nrg < .95
        nrg = nrg + S(counter)/tot;
        counter = counter + 1;
    end
    U = U(:,1:counter-1);
end

if strcmp(option,'QR')
    [U R] = qr(matrix,0);
    U = U(:,1:rank(matrix));
end

if strcmp(option,'CLOSEST')
    [U, ~, V] = svd(matrix,0);
    U = U*V';
end
