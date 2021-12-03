function [clustresult,param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename)
% SYNTAX: [clustresult,param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename)
%
% INPUTS:   'mean_type' is a string that contains the name of the mean you 
%           would like to use to compute the centers.  The options are 
%           [Flag], [Extrinsic], [Karcher], and [L2-median].
%
%           'param' is a structure whose fields correspond to the
%           parameters required for the mean type you have selected, as
%           well as the k-means clustering itself. Necessary parameters 
%           can be found in the documentation for the mean you are using, 
%           but should include, [param.angles], [param.tol], 
%           [param.k], [param.maxiters], [param.dist_type], 
%           [param.kmeans_maxiters].
%
%           'TrackletPath' is a string that contains the full path and
%           filename where the k-means tracklet data is located.
%
%           'LabelPath' is a string that contains the full path and
%           filename where the k-means label data is located.
%
%           'savename' is a string that conatins the full path and filename
%           to which you would like to periodically save the output.  This
%           k-means process can take a long time, so it is worthwhile to
%           check your output occasionally to see if it is performing as
%           expected.  If you do not wish to save your output, comment out
%           lines 235 and 251.
%
% OUTPUTS:  'clustresult' is a vector that contains the index of the
%           cluster to which each subspace was assigned.
%
%           'param' is a structure that contains fields corresponding to
%           the parameters that were passed into the function originally,
%           as well as timing and iteration information.
%
% NOTES:    This function calls 'GrDist', 'FlagMean', 'ExManifoldMean',
%           'L2Median', and 'KarcherMean', which ARE NOT included.
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

%% Add paths of called functions to MATLAB SearchPath
%addpath('/home/katrina/a/marrinan/PAL/Flag_Mean/Code')
%addpath('/home/katrina/a/marrinan/PAL/Flag_Mean/Code/NotCleanButNeeded')
%addpath(genpath('/home/katrina/a/marrinan/PAL/Flag_Mean/Relevant_Papers/WACV/WACV_Code'))

%% Load tracklets
load(TrackletPath,'-mat');
if strcmp(mean_type,'Flag')
    data = Data.on_bases;
else
    data = Data.gr_pts;
end

clear Data

%% Load labels
load(LabelPath, '-mat');

num_pts = size(data,1);
N = size(data{1},2);
param.mean_type = mean_type;
angles = param.angles;
dist_type = param.dist_type;
%% Initialize Random Centers

centers = cell(1,param.k);
R = randi(num_pts,param.k,1);
for i = 1:param.k
    centers{i} = data{R(i)};
end

%% Initialize each video to a random cluster to see how many move
clustresult = randi(param.k,num_pts,1);
group = cell(param.k,1);
old_group = cell(param.k,1);
movers = cell(param.k,1);
for ijk = 1 : param.k
    group{ijk} = find(clustresult == ijk);
end
param.conv = zeros(param.kmeans_maxiters,1);
redemption = 0;


%% Run k-means algorithm 'maxiters' times
param.mean_time = zeros(param.k,param.kmeans_maxiters);
for times = 1 : param.kmeans_maxiters
    if times > 1
        if param.conv(times-1) < param.k
            if redemption == 0
                redemption = 1;
            else
                break;
            end
        end
    end

%% Find the closest center to each data point
    distmat = zeros(param.k,num_pts);
    for m = 1 : param.k
        tmp_center = centers{1,m};
        parfor j = 1 : num_pts
            distmat(m,j) = GrDist(data{j},tmp_center, angles, dist_type);
        end
    end
    [~, I] = min(distmat);
    clustresult = I';
    
%% Create a mean for each cluster
    for clustid = 1 : param.k

% Check to see how many videos moved clusters
        old_group{clustid} = group{clustid};
        group{clustid} = find(clustresult==clustid);
        A1 = setdiff(old_group{clustid},group{clustid});
        movers{clustid} = length(A1);
% Flag Mean
        if strcmp(param.mean_type,'Flag')
            subset = ismember(clustresult,clustid);
            if sum(subset) ~= 0
                subsetlen = sum(subset);
                param.dim = 'ALL';
                param.type = 'SVD';
                param.option = 'AsIs';
                to_mean = data(subset);
                P = size(to_mean{1},2);
                N = min(N,P);
                [mean_sub, ~, param.mean_time(clustid,times)] = FlagMean(param.dim, param.type, param.option, to_mean);
                centers{clustid} = mean_sub(:,1:N);
            else
                % If a cluster is empty give it a random center
                R = randi(num_pts,1,1);
                centers{clustid} = data{R(1)};
                display(strcat('Cluster #',num2str(clustid),' has no members'))
            end
        end 
        
% Extrinsic Mean
        if strcmp(param.mean_type,'Extrinsic')
            subset = ismember(clustresult,clustid);
            if sum(subset) ~= 0
                subsetlen = sum(subset);
                check = 0;
                to_mean = data(subset);
                [mean_sub, ~, param.mean_time(clustid,times)] = ExManifoldMean(to_mean,check);
                centers{clustid} = mean_sub(:,1:N);
            else
                % If a cluster is empty give it a random center
                R = randi(num_pts,1,1);
                centers{clustid} = data{R(1)};
                display(strcat('Cluster #',num2str(clustid),' has no members'))
            end
        end         
        
% Karcher Mean
        if strcmp(param.mean_type, 'Karcher')
           subset = ismember(clustresult,clustid);
           if sum(subset) ~= 0
               subsetlen = sum(subset);
               to_mean = data(subset);
               check = 0;
               [mean_sub, ~, param.mean_iters(clustid,times), param.mean_time(clustid,times)] = KarcherMean(to_mean, to_mean{randi(subsetlen,1)}, param.tol, param.maxiters,check);
               centers{clustid} = mean_sub;
           else
                % If a cluster is empty, choose a random center
                R = randi(num_pts,1,1);
                centers{clustid} = data{R(1)};
                display(strcat('Cluster #',num2str(clustid),' has no members'))
           end
        end
    
% L2-Median
        if strcmp(param.mean_type, 'L2-median')
           subset = ismember(clustresult,clustid);
           if sum(subset) ~= 0
               subsetlen = sum(subset);
               to_mean = data(subset);
               check = 0;
               [mean_sub, ~, param.mean_iters(clustid,times), param.mean_time(clustid,times)] = L2Median(to_mean, to_mean{randi(subsetlen,1)}, param.tol, param.maxiters,check);
               centers{clustid} = mean_sub;
           else
                % If a cluster is empty, choose a random center
                R = randi(num_pts,1,1);
                centers{clustid} = data{R(1)};
                display(strcat('Cluster #',num2str(clustid),' has no members'))
           end
        end        
    end
    this = 0;
    for pppp = 1 : param.k
        this = this + movers{pppp};
    end
    param.conv(times) = this;
    save(savename,'clustresult','param','-v7.3');
    %display('Done saving data') 
end

%% Find the closest centers one last time
display('Last time')
    distmat = zeros(param.k,num_pts);
    for m = 1 : param.k
        tmp_center = centers{1,m};
        parfor j = 1 : num_pts
            distmat(m,j) = GrDist(data{j},tmp_center, angles, dist_type);
        end
    end
    [~, I] = min(distmat);
    clustresult = I';
matlabpool close
save(savename,'clustresult','param','-v7.3')
display('All Done!')



