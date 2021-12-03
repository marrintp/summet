%% CVPR Prototype Selection Example
% This script file should recreate Figure 4 from the paper,
%
%   T. Marrinan, J. Ross Beveridge, B. Draper, M. Kirby, and C. Peterson.
%   "Finding the Subspace Mean or Median to Fit Your Need." CVPR, 2014.
%
% A couple of notes:
%
% (1) The first three lines in this script need to point to the DARPA 
% Mind's Eye video clip data, the associated labels,  and data set that 
% contains the results of clustering those videos using agglomerative 
% clustering with Ward's linkage.  These files are all provided as 
% part of the Subspace Mean and Median Evaluation Toolkit.
%
% (2) There are some files that are saved during the execution of this 
% script. They default to being saved in your current directory, but if 
% you would like to change those locations, they are all set in the section 
% titled 'Save Data'. (The string saved to the variable 'DataPath' 
% specifies the location and file name. If you DON'T want those files 
% saved, you need to edit the function 'GrKmeans' as well.
%
% (3) The ExemplarSelection function uses the Matlab parallel processing 
% toolbox. If you do not have access to that, you need to edit that 
% function appropriately.
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

TrackletPath = '/home/katrina/a/marrinan/PAL/Flag_Mean/Code/CVPR/CodeForPackage/NecessaryData/DARPA_tracklets_2345_09_05.mat';
LabelPath = '/home/katrina/a/marrinan/PAL/Flag_Mean/Code/CVPR/CodeForPackage/NecessaryData/smaller_action_labels_2345.mat';
ExemplarInputPath = '/home/katrina/a/marrinan/PAL/Flag_Mean/Code/CVPR/CodeForPackage/NecessaryData/ExemplarInputs.mat';
DataPath = 'GrPrototypeSelection_data.mat';
len = 7;

%% Perform prototype (exemplar) selection experiment

% Flag Mean
flag.ex_ids = cell(len,1);
flag.info = cell(len,1);
for i = 1 : len
    clust_idx = i
    mean_type = 'Flag';
    param.angles = 1;
    [flag.ex_ids{i},flag.info{i}] = ExemplarSelection(clust_idx,mean_type,TrackletPath,LabelPath,ExemplarInputPath,param);
end
disp('Flag mean done')

% Extrinsic mean
extrinsic.ex_ids = cell(len,1);
extrinsic.info = cell(len,1);
for i = 1 : len
    clust_idx = i;
    mean_type = 'Extrinsic';
    param.angles = 1;
    [extrinsic.ex_ids{i},extrinsic.info{i}] = ExemplarSelection(clust_idx,mean_type,TrackletPath,LabelPath,ExemplarInputPath,param);
end
disp('Extrinsic mean done')

% Karcher Mean
karcher.ex_ids = cell(len,1);
karcher.info = cell(len,1);
for i = 1 : len
    tic;
    tstart = tic;
    clust_idx = i
    mean_type = 'Karcher';
    param.angles = 1;
    param.tol = 0.01;
    param.maxiters = 100;
    [karcher.ex_ids{i},karcher.info{i}] = ExemplarSelection(clust_idx,mean_type,TrackletPath,LabelPath,ExemplarInputPath,param);
    laptime = toc(tstart)
end
disp('Karcher mean done')

% L2-median
L2.ex_ids = cell(len,1);
L2.info = cell(len,1);
for i = 1 : len
    tic;
    tstart = tic;
    clust_idx = i
    mean_type = 'L2-median';
    param.angles = 1;
    param.tol = 0.01;
    param.maxiters = 100;
    [L2.ex_ids{i},L2.info{i}] = ExemplarSelection(clust_idx,mean_type,TrackletPath,LabelPath,ExemplarInputPath,param);
    laptime = toc(tstart)
end
disp('L2-median done')

%% Save data
save(DataPath,'flag','karcher','L2','extrinsic','-v7.3')


%% CVPR Plot Exemplar Selection Data

% Load data if need be.
% If the data has been saved someplace, rather than being open in your
% workspace, load it here.
% load(DataPath, '-mat')
% load(LabelPath, '-mat')

%% Plot accuracy data points
h = fig('units','centimeters','width',17,'height',11.33,'font','Times New Roman');
hold on
%len = 7;
axis([50 50*len  30 55]);
Xtick = 50:50:1200;
FlagVec = zeros(len,1);
KarcherVec = zeros(len,1);
L2Vec = zeros(len,1);
ExtrinsicVec = zeros(len,1);

for i = 1 : len
    FlagVec(i) = flag.info{i}.response;
    KarcherVec(i) = karcher.info{i}.response;
    L2Vec(i) = L2.info{i}.response;
    ExtrinsicVec(i) = extrinsic.info{i}.response;
end
flag_color = [17,119,51]/255;
karcher_color = [68 119 170]/255;
L2_color = [204 102 119]/255;
extrinsic_color = [221 204 119]/255;

h1 = plot(Xtick(1:len),FlagVec(1:len)*100,'color',flag_color,'LineWidth',3);
h2 = plot(Xtick(1:len),KarcherVec(1:len)*100,'color',karcher_color,'LineWidth',3);
h3 = plot(Xtick(1:len),L2Vec(1:len)*100,'color',L2_color,'LineWidth',3);
h4 = plot(Xtick(1:len),ExtrinsicVec(1:len)*100,'color',extrinsic_color,'LineWidth',3);

%legend([h1,h2,h3,h4],{'Flag','Karcher (\epsilon = 0.01)','L2-median (\epsilon = 0.01)','Extrinsic'},'Location','SouthEast');
xlabel('Number of clusters')
ylabel('Percent of exemplars that match their cluster label')
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

%% Plot timing data points - Total time per cluster, semilog y-axis
h = fig('units','centimeters','width',17,'height',11.33,'font','Times New Roman');
%len = 7;
Xtick = 50:50:1200;
FlagVec = zeros(len,1);
KarcherVec = zeros(len,1);
L2Vec = zeros(len,1);
ExtrinsicVec = zeros(len,1);
EuclideanVec = zeros(len,1);

for i = 1 : len
    FlagVec(i) = sum(flag.info{i}.time);
    KarcherVec(i) = sum(karcher.info{i}.time);
    L2Vec(i) = sum(L2.info{i}.time);
    ExtrinsicVec(i) = sum(extrinsic.info{i}.time);
    EuclideanVec(i) = sum(euclidean.info{i}.time);
end

flag_color = [17,119,51]/255;
karcher_color = [68 119 170]/255;
L2_color = [204 102 119]/255;
extrinsic_color = [221 204 119]/255;

h1 = semilogy(Xtick(1:len),FlagVec(1:len),'color',flag_color,'LineWidth',3);
hold on
h2 = semilogy(Xtick(1:len),KarcherVec(1:len),'color',karcher_color,'LineWidth',3);
h3 = semilogy(Xtick(1:len),L2Vec(1:len),'color',L2_color,'LineWidth',3);
h4 = semilogy(Xtick(1:len),ExtrinsicVec(1:len),'color',extrinsic_color,'LineWidth',3);
 
%legend([h1,h2,h3,h4],{'Flag','Karcher (\epsilon = 0.01)','L2-median (\epsilon = 0.01)','Extrinsic'},'Location','SouthEast');
xlabel('Number of clusters')
ylabel('Time to compute the means of all clusters (seconds)')

% fn = 'plot';%filename
%     l = 1:length(fn);
%     fn = sprintf('%s_LEGEND',fn(l));%give it a separate filename
%     a = get(gcf,'children');% link to ledgend is now a(1)
%     b = get(gca,'children');% link to the data curve/s
%     set(b,'visible','off'); %hide data...
%     set(a(2),'visible','off'); %hide axes etc...
%     legfs = get(a(1),'Fontsize'); %get legend fontsize
%     set(a(1),'Fontsize',legfs+1); %make legend appear larger
%     pause;