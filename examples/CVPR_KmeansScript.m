%% CVPR K-means Example
% This script file should recreate Figure 5 from the paper,
%
%   T. Marrinan, J. Ross Beveridge, B. Draper, M. Kirby, and C. Peterson.
%   "Finding the Subspace Mean or Median to Fit Your Need." CVPR, 2014.
%
% A couple of notes:
%
% (1) The first two lines in this script need to point to the kmeans data
% points and the kmeans action labels that are provided as part of the
% Subspace Mean and Median Evaluation Toolkit.
%
% (2) There are some files that are saved during the execution of this 
% script. They default to being saved in your current directory, but if 
% you would like to change those locations, they are all set in the first 
% section. (The string saved to the variable 'DataPath' is one place, and 
% each time the variable 'savename' is set are the others. If you DON'T 
% want those files saved, you need to edit the function 'GrKmeans' as well.
%
% (3) The GrKmeans function uses the Matlab parallel processing toolbox.
% If you do not have access to that, you need to edit that function
% appropriately.
%
% (4) The legends are currently commented out from the plots.
%
% (5) This entire procedure takes a VERY long time (look at the timing plot
% from the paper).  Just a warning.
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

TrackletPath = '/home/katrina/a/marrinan/PAL/Flag_Mean/Code/CVPR/CodeForPackage/NecessaryData/kmeans_pts.mat';
LabelPath = '/home/katrina/a/marrinan/PAL/Flag_Mean/Code/CVPR/CodeForPackage/NecessaryData/kmeans_action_labels.mat';
DataPath = 'GrKmeans_data.mat';

%%
flag = cell(5,20);

karcher = cell(5,20);

extrinsic = cell(5,20);

L2 = cell(5,20);

for tt = 1 : 20
    tt
    tic;
    tbegin = tic;
    j = 1;
    for i = 5:5:25
        % Params for all...
        param.k = i;
        param.angles = 1;
        param.kmeans_maxiters = 20;
        param.tol = 0.1;
        param.maxiters = 100;
        param.dist_type = 'DGEO';
    
        % Flag
        tic;
        tstart = tic;
        mean_type = 'Flag';
        savename = strcat('GrKmeansFlag_',num2str(i),'.mat');
        [flag{j,tt}.clustresult,flag{j,tt}.param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename);
        display('Flag iteration done')
        flag_time = toc(tstart)
    
        % Karcher
        tic;
        tstart = tic;
        mean_type = 'Karcher';
        savename = strcat('GrKmeansKarcher_',num2str(i),'.mat');
        [karcher{j,tt}.clustresult,karcher{j,tt}.param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename);
        display('Karcher iteration done')
        karcher_time = toc(tstart)
    
        % Extrinsic
        tic;
        tstart = tic;
        mean_type = 'Extrinsic';
        savename = strcat('GrKmeansExtrinsic_',num2str(i),'.mat');
        [extrinsic{j,tt}.clustresult,extrinsic{j,tt}.param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename);
        display('Extrinsic iteration done')
        extrinsic_time = toc(tstart)
    
        % L2-median
        tic;
        tstart = tic;
        mean_type = 'L2-median';
        savename = strcat('GrKmeansL2_',num2str(i),'.mat');
        [L2{j,tt}.clustresult,L2{j,tt}.param] = GrKmeans(mean_type,param,TrackletPath,LabelPath,savename);
        display('L2-median iteration done')
        L2_time = toc(tstart)

        save(DataPath,'flag','karcher','L2','extrinsic','-v7.3')
        display('Finished one iteration of each mean.')
        j = j+1
        
    end
    LapTime = toc(tbegin)
end


%% Create CVPR K-means Box Plots

runs = 10;

% If the data has been saved someplace, rather than being open in your
% workspace, load it here.
% load(DataPath, '-mat')
% load(LabelPath, '-mat')

action_labels = kmeans_action_labels;
clear kmeans_action_labels;

Flag.response = zeros(5,runs);
Karcher.response = zeros(5,runs);
Extrinsic.response = zeros(5,runs);
Median.response = zeros(5,runs);


%% Iterate the label accuracy computations
for i = 1 : runs
    for j = 1 : 5
% Find the dominant class of each cluster and give the cluster that label
        flag{j,i}.output = LabelCounter(action_labels,flag{j,i}.clustresult,'kmeans');
        karcher{j,i}.output = LabelCounter(action_labels,karcher{j,i}.clustresult,'kmeans');
        L2{j,i}.output = LabelCounter(action_labels,L2{j,i}.clustresult,'kmeans');
        extrinsic{j,i}.output = LabelCounter(action_labels,extrinsic{j,i}.clustresult,'kmeans');
    end
end

%%
for i = 1 : runs
    for j = 1 : 5
% Compute the number of videos that are correctly labeled from cluster labels        
        Flag.response(j,i) = ClusterPurity(flag{j,i}.output,action_labels,'kmeans','B');
        Karcher.response(j,i) = ClusterPurity(karcher{j,i}.output,action_labels,'kmeans','B');
        Median.response(j,i) = ClusterPurity(L2{j,i}.output,action_labels,'kmeans','B');
        Extrinsic.response(j,i) = ClusterPurity(extrinsic{j,i}.output,action_labels,'kmeans','B');
        
        Flag.convergence(j,i) = length(find(flag{j,i}.param.conv));
        Karcher.convergence(j,i) = length(find(karcher{j,i}.param.conv));
        Median.convergence(j,i) = length(find(L2{j,i}.param.conv));
        Extrinsic.convergence(j,i) = length(find(extrinsic{j,i}.param.conv));
        
% Compute the total time to create all means for each run of k
        Flag.times(j,i) = sum(sum(flag{j,i}.param.mean_time));
        Karcher.times(j,i) = sum(sum(karcher{j,i}.param.mean_time));
        Median.times(j,i) = sum(sum(L2{j,i}.param.mean_time));
        Extrinsic.times(j,i) = sum(sum(extrinsic{j,i}.param.mean_time));
    end
end

%% Create the accuracy plot - using aboxplot
ntypes = 4;
nruns = runs;
nkvals = 5;

data = zeros(ntypes,nruns,nkvals);

data(1,:,:) = Flag.response';
data(2,:,:) = Karcher.response';
data(3,:,:) = Median.response';
data(4,:,:) = Extrinsic.response';

flag_color = [17,119,51]/255;
karcher_color = [68 119 170]/255;
L2_color = [204 102 119]/255;
extrinsic_color = [221 204 119]/255;

TempColor = [flag_color;  
             karcher_color;
             L2_color;  
             extrinsic_color];
        

h = fig('units','centimeters','width',17,'height',11.33,'font','Times New Roman');
aboxplot(data,'labels',[5,10,15,20,25],'colormap',TempColor)
xlabel('Number of clusters (k-value)')
ylabel('Average cluster purity')
%legend('Flag','Karcher (\epsilon = 0.1)','L2-median (\epsilon = 0.1)','Extrinsic','Random','Location','EastOutside')
% fn = 'plot';%filename
%     l = 1:length(fn);
%     fn = sprintf('%s_LEGEND',fn(l));%give it a separate filename
%     a = get(gcf,'children');% link to ledgend is now a(1)
%     b = get(gca,'children');% link to the data curve/s
%     set(b,'visible','off'); %hide data...
%     set(a(2),'visible','off'); %hide axes etc...
%     legfs = get(a(1),'Fontsize'); %get legend fontsize
%     set(a(1),'Fontsize',legfs+1); %make legend appear larger
%     pause

%% Create the total timing plot - using aboxplot
ntypes = 4;
nruns = runs;
nkvals = 5;

data = zeros(ntypes,nruns,nkvals);

data(1,:,:) = Flag.times';
data(2,:,:) = Karcher.times';
data(3,:,:) = Median.times';
data(4,:,:) = Extrinsic.times';

TempColor = [flag_color;  
             karcher_color;
             L2_color;  
             extrinsic_color];
        

h = fig('units','centimeters','width',17,'height',11.33,'font','Times New Roman');
aboxplot(data,'labels',[5,10,15,20,25],'colormap',TempColor)
xlabel('Number of clusters (k-value)')
ylabel('Total time to compute means (seconds)')
%legend('Flag','Karcher (\epsilon = 0.1)','L2-median (\epsilon = 0.1)','Extrinsic')
set(gca,'YScale','log')