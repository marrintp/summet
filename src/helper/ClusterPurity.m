function [response] = ClusterPurity(output,action_labels,which,method)
% SYNTAX:   [response] = ClusterPurity(output,action_labels,which,method)
%
% INPUTS:   'output' is a cell array of structures that comes from the 
%           function 'LabelCounter'.  Each cell contains the name of the 
%           dominant label in a given cluster, the percentage of subspaces 
%           in that cluster that share that dominant label, the full list 
%           of labels in the cluster and the count for each label.
%
%           'action_labels' is the structure of action labels that
%           corresponds to the DARPA Mind's Eye video clips that are used
%           for the task at hand.
%
%           'which' is a string that specifies the task that this function
%           is being called for.  Options are [kmeans] and [exemplar].
%
%           'method' is a string that specifies how to determine the purity
%           of a cluster.  Options are [A] or [B].  See below for details.
%
% OUTPUTS:  'response' is a scalar that measures the average purity of the
%           clusters.
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

ncluster = length(output);
%% For kmeans dataset
if strcmp(which,'kmeans')
    idxs = action_labels.idxs;
    clear action_labels
end
%% For 2345 dataset
if strcmp(which,'exemplar')
    idxs = action_labels.labelidxs;
    clear action_labels
end

nlabels = length(idxs);
response = 0;
%% Options:
%% (A)   For each cluster find the probability that the label I picked as
%       dominant is the same as the label of a randomly selected video.
%       (1/weight) * output{clustid}.card = Prob(correct)
%
%       Then how many videos would provide the correct solution?
%       Prob(correct) * sum(output{clustid}.counts) = Correct options for
%       cluster clustid
%
%       Then sum the number of correct options and divide by total options
%       response = response/length(unique(idxs));
if strcmp(method,'A')
    for clustid = 1:ncluster
        weight = length(output{clustid}.dom); % How many dominant labels are there in a cluster
        if weight ~= 0
            response = response + (1/weight)*output{clustid}.card*sum(output{clustid}.counts); % How many videos were correctly labeled.
        end
    end
    response = response/length(unique(idxs)); % Correct options/total options
end

%% (B)  Average the probability of picking the correct label for each cluster
%       (1/weight) * output{clustid}.card = Prob(correct)
%
%       Divide by the number of clusters
%       response = response/ncluster
if strcmp(method, 'B')
    for clustid = 1:ncluster
        weight = length(output{clustid}.dom); % How many dominant labels are there in a cluster
        if weight ~= 0
            response = response + (1/weight)*output{clustid}.card; % How many videos were correctly labeled.
        end
    end
    response = response/ncluster; % This meaures the mean of the purities of each cluster, not how many videos we got right out of the total number.
end


