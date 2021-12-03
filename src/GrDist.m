function [distance]=GrDist(P, Q, N, dist_type)
% SYNTAX:   [distance]=GrDist(P, Q, N, dist_type)
%
% INPUTS:   'P' and 'Q' are matrices that span the two subspaces you wish 
%           to measure the distance between.  They must be of the same
%           dimension.
%
%           'N' is the number of angles to be used in our distance
%           calculation. Default value is N = size(P,2);
%
%           'dist_type' is a string that specifies which distance measure
%           you would like to use.  Options are:
%           [DGEO] for the Geodesic distance
%           [DPF] for the Projection F-norm
%           [DC] for the Chordal distance
%           [DCF] for the Chordal Frobenius
%           [DFS] for the Fubini-Studi
%           [MAX] for the maximum angle
%           [MIN] for the minimum angle
%           [ALL] for all the angles
%           [EUC] for the Frobenius norm of P-Q
%           Default value is 'DGEO'
%
% OUTPUTS:  'distance' is the measure between the two subspaces using the
%           distance specified.
%
% USE:      The point of this function is compute the distance between two
%           subspaces that are thought of as points on a Grassmannian
%           manifold.
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

% Set up the default values
distance = 500;

if ~exist('dist_type','var')
    dist_type = 'DGEO';
end

q = min(size(Q,2),size(P,2));
N=min(N,q);

%Euclidean
if strcmp(dist_type,'EUC')
    distance = norm(P-Q,'fro');
    return
end

if size(P,1)~=size(Q,1)
    disp('Halting due to dimension mismatch.')
    size(P,1)
    size(Q,1)
    return
end

% Compute the principle angles between my input matrices and calculate 
% the distance.
S = svd(P'*Q);
T=zeros(N,1);

for k=1:N
 if S(k)<10^(-8)
     T(k)=sqrt(2*(1-S(k)));
 else
     T(k)=real(acos(S(k))); 
 end
end

%Geodesic distance
if strcmp(dist_type,'DGEO')
    distance = norm(T);
end

%Projection f-norm
if strcmp(dist_type,'DPF')
    distance = norm(sin(T));
end

%chordal distance
if strcmp(dist_type,'DC')
    distance = norm(sin(T));
end

%chordal Frobenius
if strcmp(dist_type,'DCF')
    distance = norm(2*sin(T/2));
end

%Fubini-Studi
if strcmp(dist_type,'DFS')
    distance = acos(prod(S(1:N)));
end

%Min Angle
if strcmp(dist_type,'MIN')
    distance = T(1);
end

%Max Angle
if strcmp(dist_type,'MAX')
    distance = T(N);
end

%All the angles
if strcmp(dist_type,'ALL')
    distance = T;
end



