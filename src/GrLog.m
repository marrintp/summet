function [TY] = GrLog(X,Y)
% SYNTAX:   [TY] = GrLog(X,Y)
%
% INPUTS:   'X' is the point about which the tangent space has been
%           computed.
%
%           'Y' is the point on the Grassmannian manifold that is to be
%           mapped to the tangent space of X.
%
% OUTPUTS:  'TY' is the representation of Y in the tangent space of X.
%
% NOTES:    Algorithm implemented from description in:
%
%           E. Begelfor and M. Werman.  "Affine invariance revisited."
%           CVPR, 2:2087 - 2094, 2006    
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

m = size(X,1);

% temp = (eye(m)-X*X')*Y*inv(X'*Y);
% The following line is a slightly faster way to compute temp.
temp = eye(m)*Y*inv(X'*Y) - X*(X'*Y)*inv(X'*Y);
[U,S,V] = svd(temp,0);
Theta = atan(S);

TY = U*Theta*V';