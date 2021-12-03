function [MU,err,iters,tellapsed] = L2Median(data,start_pt,tol,maxiters,check)
%% SYNTAX:   [MU, err, iters, tellapsed] = L2median(data, start_pt, tol, maxiters, check)
%
% INPUTS:   'data' is the collection of subspaces you wish to find a
%           representative for.  It should be entered as a num_pts by 1
%           cell array with each cell being an n by q_i matrix, where q_i 
%           CANNOT vary.  In other words, the datacubes need to be unrolled
%           already when using this function
%
%           'start_pt' is a point from 'data' that you would like to use as
%           your initial base point.
%
%           'tol' how close MU needs to be at step n, relative to MU at
%           step n+1 to say that the algorithm has converged.
%
%           'maxiters' is the maximum allowable number of iterations.  The
%           algorithm will stop iterating if this threshold is reached and
%           return the current value of MU.
%
%           'check' determines whether or not the algorithm will verify
%           that your data points are in fact orthonormal matrices.  If
%           check == 1, it does the check, otherwise it does not.
%
% OUTPUTS:  'MU' is the Karcher mean, returned as a matrix. 
%
%           'err' is the distance between MU at the second to last 
%           iteration and MU at the last iteration.
%
%           'tellapsed' is the total time it took to compute MU.
%
% NOTES:    This function calls GrExp, GrLog, and GrDist which are NOT
%           included.
%
%           Algorithm implemented from description in:
%
%           P.T. Fletcher, S. Venkatasubramanian, and S. Joshi. "The
%           geometric median on Riemannian manifolds wiht application to
%           robust atlas estimation." NeuroImage, 45(1 Suppl):S143, 2009.
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
err = 999999;
iters = 1;
step_size = 1;

tic;
tstart = tic;

%test to make sure points are on Grassmannian
[nr, nc] = size(start_pt);
M = start_pt'*start_pt;
II = eye(nc);
if norm(M-II,'fro') >= etol
    display('Starting point is not orthonormal.');
    [U, ~, ~] = svd(start_pt,0);
    start_pt = U*V';
end

num_pts = size(data,1);
if num_pts == 1
    MU = data{1};
    err = 0;
    iters = 1;
    tellapsed = toc(tstart);
    return;
end
    
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

% Main loop of the unweighted Geometric Median for Riemannian manifolds
MU = start_pt;
while err > tol
    if mod(iters,100) == 0
        fprintf('On iteration %g, the error is %g. \n',iters,err)
    end
    SUMT = zeros(nr, nc);
    temp_dist = zeros(num_pts,2);
    for i = 1:num_pts
       temp_dist(i,1) = GrDist(MU, data{i}, size(MU,2), 'DGEO');
       if temp_dist(i,1) > 10^-5
            temp_dist(i,2) = 1/temp_dist(i,1);
            SUMT  = SUMT + GrLog(MU,data{i})/temp_dist(i,1);
       end
    end
    
    % If weights are equal, they cancel out in this sum
    if sum(temp_dist(:,2)) < 10^-5
        display('Trying to divide by zero,returning current MU without convergence.')
        return;
    end
    
    SUMT = SUMT/sum(temp_dist(:,2));
    NewMU = GrExp(MU, step_size*SUMT);
    [err] = GrDist(NewMU, MU, size(MU,2), 'DGEO');
    MU = NewMU;
    
    if err < tol
        tellapsed = toc(tstart);
        return;
    end
    iters = iters  + 1;
    if iters >= maxiters
        tellapsed = toc(tstart);
        return
    end
end