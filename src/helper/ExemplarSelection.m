function [exemplar_idxs,info] = ExemplarSelection(clust_idx,mean_type,TrackletPath,LabelPath,ExemplarInputPath,param)
% SYNTAX:   [exemplar_idxs,info] = ExemplarSelection(clust_idx,mean_type,param)
%
% INPUTS:   'clust_idx' is an integer that specifies where the cut was made
%           in the agglomerative clustering.  In other words, it tells you
%           how many clusters the set of videos has been broken into.  [1]
%           corresponds to '50' clusters, [2] to '100', and so on in steps
%           of 50.
%
%           'mean_type' is a string that specifies which type of subspace
%           average is to be used.  The options are [Flag], [Karcher],
%           [L2-median], and [Extrinsic].
%
%           'TrackletPath' is a string that specifies the full path and
%           file name of the DARPA Mind's Eye video clip data.
%
%           'LabelPath' is a string that specifies the full path and
%           file name of the DARPA Mind's Eye data action labels.
%
%           'ExemplarInputPath' is a string that specifies the full path 
%           and file name of the results of clustering the Mind's Eye
%           videos as subspaces using agglomerative clustering with Ward's
%           linkage.
%
%           'param' is a Matlab structure contains the parameters needed
%           for each of the various means, as well as the distance metric
%           to be used and the number of angles to use for the distance 
%           calculation
%
% OUTPUTS:  'exemplar_idxs' is a vector that contains the indices of the
%           prototypes (or exemplars) chosen for each cluster
%
%           'info' is a Matlab structure that contains information about
%           the procedure that was performed and the parameters used, as
%           well as the percentage of videos in a given cluster that match
%           the label of the exemplar chosen.  That information is stored
%           at info.response.
%
% NOTES:    This function calls 'GrDist', 'FlagMean', 'ExManifoldMean',
%           'L2Median', and 'KarcherMean', which ARE NOT included.  This
%           function also calls 'exemplar_matcher' and 'vid2matrix' which 
%           ARE included below.
%
%           This function uses the Matlab parallel processing toolbox. 
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

%% Open 12 workers
matlabpool('open',12);

%% Load tracklets
load(TrackletPath, '-mat');

%% Load labels
load(LabelPath, '-mat');

%% Load agglomerative clustering clustresult
load(ExemplarInputPath, '-mat');
clustresult = clustresult{clust_idx};

display('done loading data')

%% Unfold data
darpa_unfolded = cell(3,1);
b = size(tracklets);
darpa_unfolded{1} = cell(b(1),1);

for i = 1 : b(1)
    tracklets{i} = tracklets{i}.mean_subtracted_normalized_data;
    darpa_unfolded{1}{i} = vid2matrix(tracklets{i},1);
end
on_bases = cell(b(1),1);
gr_pts = cell(b(1),1);
for i = 1 : b(1)
    [U,~,V] = svd(darpa_unfolded{1}{i},0);
    on_bases{i,1} = U(:,1:rank(darpa_unfolded{1}{i}));
    gr_pts{i,1} = U*V';
end
clear darpa_unfolded tracklets
display('Done unfolding data')

%% Work on one cluster
ncluster = length(unique(clustresult));
exemplar_idxs = zeros(ncluster,1);
time = zeros(ncluster,1);
dist_time = zeros(ncluster,1);
iters = zeros(ncluster,1);
clust_size = zeros(ncluster,1);
for clustid=1:ncluster 
    subset = find(clustresult==clustid);
    subsetlen = length(subset);
    clust_size(clustid,1) = subsetlen;
    % Create mean
    if strcmp(mean_type,'Flag')
        to_mean = on_bases(subset);
        param.dim = 'ALL';
        param.type = 'SVD';
        param.option = 'AsIs';
        [mean_sub, ~, time(clustid,1)] = FlagMean(param.dim, param.type, param.option, to_mean);
        iters(clustid,1) = 1;
    elseif strcmp(mean_type,'Karcher')
        check = 0;
        to_mean = gr_pts(subset);
        [mean_sub, ~, iters(clustid,1), time(clustid,1)] = KarcherMean(to_mean, to_mean{randi(subsetlen,1)}, param.tol, param.maxiters,check);
    elseif strcmp(mean_type,'L2-median')
        check = 0;
        to_mean = gr_pts(subset);
        size(to_mean);
        [mean_sub,~,iters(clustid,1), time(clustid,1)] = L2Median(to_mean,to_mean{randi(subsetlen,1)},param.tol,param.maxiters,check);
    elseif strcmp(mean_type,'Extrinsic')
        to_mean = gr_pts(subset);
        d = size(to_mean{1});
        check = 0;
        [mean_sub, ~, time(clustid,1)] = ExManifoldMean(to_mean, check);
        mean_sub = mean_sub(:,1:d(2));
        iters(clustid,1) = 1;
    end
    
    % Measure distances
    tic;
    tstart = tic;
    temp_dist = zeros(subsetlen,1);
    if strcmp(mean_type,'Flag')
        d = size(to_mean{1},2);
        tmp_mean = mean_sub(:,d);
        d = param.angles;
        parfor i = 1 : subsetlen
            [temp_dist(i,1)] = GrDist(to_mean{i},tmp_mean,d,'DGEO');
        end
    elseif strcmp(mean_type,'Karcher')
        %d = size(to_mean{1},2);
        d = param.angles;
        parfor i = 1 : subsetlen
            [temp_dist(i,1)] = GrDist(to_mean{i},mean_sub,d,'DGEO');
        end    
    elseif strcmp(mean_type,'L2-median')
        %d = size(to_mean{1},2);
        d = param.angles;
        parfor i = 1 : subsetlen
            [temp_dist(i,1)] = GrDist(to_mean{i},mean_sub,d,'DGEO');
        end    
    elseif strcmp(mean_type, 'Extrinsic')
        %d = size(to_mean{1},2);
        d = param.angles;
        parfor i = 1 : subsetlen
            [temp_dist(i,1)] = GrDist(to_mean{i},mean_sub,d,'DGEO');
        end       
    end
    [~, tmp_idx] = min(temp_dist);
    exemplar_idxs(clustid,1) = subset(tmp_idx);
    dist_time(clustid,1) = toc(tstart);
end

[response] = exemplar_matcher(exemplar_idxs,output{clust_idx},smaller_action_labels);

info.clust_idx = clust_idx;
info.param = param;
info.mean_type = mean_type;
info.dist_type = 'DGEO, with just first angle';
info.time = time;
info.iters = iters;
info.dist_time = dist_time;
info.clust_size = clust_size;
info.response = response;

matlabpool('close');

function [mat] = vid2matrix(VID,unrolling)
% SYNTAX:   [mat] = vid2matrix(VID,unrolling)
%
% INPUTS:   'VID' is a three-way array
%
%           'unrolling' specifies how you would like your data cube turned 
%           into a matrix. Options are:
%           [1] for num_rows*num_columns by time (standard vers.) 
%           [2] for num_rows*time by num_columns
%           [3] for num_columns*time by num_rows 
%
% OUTPUTS:  'mat' is the unfolded matrix representation of the three-way
%           array
%
% NOTES:
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
% --------------------------------------------------------------------

d = size(VID);

if ~exist('unrolling','var')
    unrolling = 1;
end

if unrolling == 1   % rows*columns by time
    mat = reshape(VID,d(1)*d(2),d(3));
end
if unrolling == 2   % rows*time by columns
    mat = zeros(d(1)*d(3),d(2));
    for j = 1 : d(3)
        for k = 1 : d(2)
            mat((j-1)*d(1)+1:j*d(1),k) = VID(:,k,j);
        end
    end
end
if unrolling == 3   % columns*time by rows
    mat = zeros(d(2)*d(3),d(1));
    for j = 1 : d(3)
        for k = 1 : d(1)
            mat((j-1)*d(2)+1:j*d(2),k) = VID(k,:,j)';
        end
    end
end

function [response] = exemplar_matcher(exemplars,output,action_labels)
% SYNTAX:   [response] = exemplar_matcher(exemplars,output,action_labels)
%
% INPUTS:   'exemplars' is a vector that contains the indices of the
%           prototypes (or exemplars) chosen for each cluster.
%
%           'output' is a cell array of structures that comes from the 
%           function 'LabelCounter' (in this case it was saved to the 
%           ExemplarInput data file).  Each cell contains the name of the 
%           dominant label in a given cluster, the percentage of subspaces 
%           in that cluster that share that dominant label, the full list 
%           of labels in the cluster and the count for each label.
%
%           'action_labels' is the structure of action labels that
%           corresponds to the DARPA Mind's Eye video clips that are used
%           for the task at hand.
%
% OUTPUTS:  'response' is a scalar that represents the percentage of videos
%           in a given cluster whose label matches that of the video chosen
%           as an exemplar.
%
% NOTES:
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
% --------------------------------------------------------------------

ncluster = length(exemplars);
labels = action_labels.labellist;
idxs = action_labels.labelidxs;
response = 0;
for clustid = 1:ncluster
    weight = length(output{clustid}.dom); % How many dominant labels are there in a cluster
    ex_labels = labels(ismember(idxs,exemplars(clustid))); % What label is associated with the exemplar
    ex_weight = length(ex_labels); % How many labels are associated with the exemplar
    match = 0;
    for i = 1 : weight
        for j = 1 : ex_weight
            if strcmp(ex_labels{j,1},output{clustid}.dom{i,1})
                match = 1;
            end
        end
    end
    response = response + (1/weight)*match; % If there was more than one ...
    % ... dominant label per cluster, we weight the response by the
    % probability that we would have given the cluster that label (by
    % random chance).
    
    %response = response + match;
end
response = response/ncluster;