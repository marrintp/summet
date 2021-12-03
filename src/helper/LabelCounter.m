function [output] = LabelCounter(action_labels,clustresult,which)
% SYNTAX:   [output] = LabelCounter(action_labels,clustresult,which)
%
% INPUTS:   'action_labels' is the structure of action labels that
%           corresponds to the DARPA Mind's Eye video clips that are used
%           for the task at hand.
%
%           'clustresult' is the vector of indices that correspond to
%           cluster membership for each of the subspaces that have been
%           clustered.
%
%           'which' is a string that specifies the task that this function
%           is being called for.  Options are [kmeans] and [exemplar].
%
% OUTPUTS:  'output' is a cell array of structures that contain the name 
%           of the dominant label in a given cluster, the percentage of
%           subspaces in that cluster that share that dominant label, the
%           full list of labels in the cluster and the count for each
%           label.
%
% NOTES:
%
% LAST EDITED: 06/22/14 by Tim Marrinan at Colorado State University.
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

ncluster = length(unique(clustresult));
%% For k-means
if strcmp(which,'kmeans')
    labels = action_labels.labels;
    labelidxs = action_labels.idxs;
end

%% For 2345
if strcmp(which,'exemplar')
    labels = action_labels.labellist;
    labelidxs = action_labels.labelidxs;
end

output = cell(ncluster,1);
for clustid=1:ncluster
    subset = find(clustresult==clustid);
    %subsetlen = length(subset);
    temp = ismember(labelidxs,subset);
    namelist = labels(temp);
    shortlist = unique(namelist);
    counts = zeros(length(shortlist),1);
    for i = 1 : length(shortlist)
       counts(i,1) = sum(strcmp(namelist,shortlist{i}));
    end
    ids = ismember(counts,max(counts));
    output{clustid}.dom = shortlist(ids);
    output{clustid}.card = max(counts)/length(namelist);
    output{clustid}.names = shortlist;
    output{clustid}.counts = counts;
end