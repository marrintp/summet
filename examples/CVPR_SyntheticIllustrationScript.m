%% CVPR 2014 Synthetic Illustration Example
% This script file should recreate Figure 3 from the paper,
%
%   T. Marrinan, J. Ross Beveridge, B. Draper, M. Kirby, and C. Peterson.
%   "Finding the Subspace Mean or Median to Fit Your Need." CVPR, 2014.
% 
% The process does generate lines in 2D chosen from a two Gaussian
% distributions, so the results will not match the paper exactly, but they
% should be qualitatively the same.
%
% LAST EDITED: by T. Marrinan 06/22/2014
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

n = 30;
sigma = 0.2;
ratio = 0;
mean_type = 'ALL';

[data] = Plot2DMeans(n,sigma,ratio,mean_type);

ratio = 0.25;
[data] = Plot2DMeans(n,sigma,ratio,mean_type,data);

ratio = 0.5;
[data] = Plot2DMeans(n,sigma,ratio,mean_type,data);

ratio = 0.75;
[data] = Plot2DMeans(n,sigma,ratio,mean_type,data);

ratio = 1;
[data] = Plot2DMeans(n,sigma,ratio,mean_type,data);