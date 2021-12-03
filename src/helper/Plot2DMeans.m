function [data] = Plot2DMeans(n,sigma,ratio,mean_type,prev_data)
% SYNTAX:   [data] = Plot2DMeans(n,sigma,ratio,mean_type,prev_data)
%
% INPUTS:   'n' is the number of points that will be sampled from the first
%           normal distributon.  Typically about the point [0,1]';
%
%           'sigma' is the standard deviation in the sampling.
%
%           'ratio' is the ratio of points that should be sampled from the
%           second distribution relative to 'n'.
%
%           'mean_type' is a string that specifies which means you would
%           like included in the plots.  Currently only accepts [ALL].
%
%           'prev_data' is an optional argument that when used should just
%           be the output 'data' from previous call of this function.
%           Allows the user to use the same random points for multiple
%           tries.  If using this argument, 'n' and 'sigma' must match the
%           previous 'n' and 'sigma,' but 'ratio' can be different.  Will
%           compute the means again if this option is used. 
%
% OUTPUTS:  'data' is a structure that contains the fields,
%
%           'data.n{1}'/'data.n{2}' - number of points in dist. 1/2
%
%           'data.sigma' - sigma
%
%           'data.ratio' - ratio
%
%           'data.start_id' - is the index of the pt to start the karcher
%           mean and L2-median so that you can start with the same point if
%           desired.
%
%           'data.mu{1}'/'data.mu{2}' - coordinates of the means
%
%           'data.pts{1}'/'data.pts{2}' - the coordinates of points in
%           dist. 1/2
%
%           'data.Z' - the coordinates of the points in a cell array for
%           the mean functions
%
%           'data.flagmean' - flag mean
%
%           'data.karchermean' - karcher mean
%
%           'data.gr_median' - L2-median
%
% NOTES:    This funtion calls 'FlagMean', 'L2Median', and 'KarcherMean' 
%           which ARE NOT included, and it also calls 'plot_vec' which IS 
%           included below.
%
%           The use of the extrinsic manifold mean is currently commented
%           out, because in 2D the solution is the same as the first vector
%           of the flag mean.  If you would like to see both, uncomment
%           lines 192 and 197.
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

if exist('prev_data','var')
    fresh = 0;
    data = prev_data;
    data.n{1} = n;
    data.n{2} = ceil(n*ratio);
    data.sigma = sigma;
    data.ratio = ratio;
else
    fresh = 1;
    data.n{1} = n;
    data.n{2} = ceil(n*ratio);
    data.sigma = sigma;
    data.ratio = ratio;
    data.start_id = randi(data.n{1}+data.n{2},1);
end

%% Create data points if none currently exist.
if fresh
    data.mu{1} = [0 1]';
    data.mu{2} = [1 0]';
    data.pts{1} = zeros(2,data.n{1});
    data.pts{2} = zeros(2,data.n{2});
    data.Z = cell(data.n{1}+data.n{2},1);
    for i = 1 : data.n{1}
        tmp = normrnd(data.mu{1},data.sigma,[2,1]);
        data.pts{1}(:,i) = tmp/norm(tmp);
        data.Z{i} = data.pts{1}(:,i);
    end
    for j = 1 : data.n{2}
        tmp = normrnd(data.mu{2},data.sigma,[2,1]);
        data.pts{2}(:,j) = tmp/norm(tmp);
        data.Z{data.n{1}+j} = data.pts{2}(:,j);
    end
end

%% Create data points if reusing old set.
if ~fresh
    diff = data.n{2} - size(data.pts{2},2);
    if diff > 0
        % If diff > 0, need to add more points to the 2nd distribution
        data.pts{2} = [data.pts{2}, zeros(2,diff)];
        for j = size(data.pts{2},2)-diff+1 : data.n{2}
            tmp = normrnd(data.mu{2},data.sigma,[2,1]);
            data.pts{2}(:,j) = tmp/norm(tmp);
            data.Z{data.n{1}+j} = data.pts{2}(:,j);
        end
    else
        % If diff <= 0, there are too many points in the 2nd distribution
        % already, so just use the first ones.
        clear data.Z;
        data.Z = cell(data.n{1}+data.n{2},1);
        for i = 1 : data.n{1}
            data.Z{i} = data.pts{1}(:,i);
        end
        for j = 1 : data.n{2}
            data.Z{data.n{1}+j} = data.pts{2}(:,j);
        end
    end
end

%% Compute means of all data points, excluding true means.
if strcmp(mean_type,'ALL') || strcmp(mean_type,'legend') 
    % Flag mean
    [flag_mean, ~, ~] = FlagMean('ALL', 'SVD', 'AsIs', data.Z);
    
    % Karcher mean
    [karcher_mean, ~, ~,~] = KarcherMean(data.Z, data.Z{data.start_id}, 10^-3, 100,1);
    
    % L2-median
    [gr_median,~,~,~] = L2Median(data.Z,data.Z{data.start_id},10^-3,100,1);
    
    % Extrinsic mean
    [ex_mean, ~, ~] = ExManifoldMean(data.Z,1);
    
    % GMEB center
    [~, optimal] = GMEB_dual_subgrad(data.Z, 1, 10^-4, ...
    'line_search', 'true');
    gmeb_mean = optimal.u;
end

%% Plot data points
%h=fig('units','centimeters','width',11,'height',11,'font','Helvetica','fontsize',16);
h=fig('units','centimeters','width',11,'height',11);
hold on
axis([-1 1 -1 1]);
set(gca,'YTick',[-1 0 1])
set(gca,'XTick',[-1 0 1])
for i = 1 : data.n{1}
    plot_vec(data.pts{1}(:,i),'k:');
end
for j = 1 : data.n{2}
    plot_vec(data.pts{2}(:,j),'k:');
end

%% Plot means
%h1 = plot_vec(mu,'k--');

flag_color = [17,119,51]/255;
karcher_color = [68 119 170]/255;
L2_color = [204 102 119]/255;
extrinsic_color = [221 204 119]/255;

%m2 = plot_vec(ex_mean(:,1),'c->','LineWidth',2);
m2 = plot_vec(gmeb_mean,'color',extrinsic_color,'LineWidth',3);
h2 = plot_vec(flag_mean(:,1),'color',flag_color,'LineWidth',3);
h3 = plot_vec(karcher_mean,'color',karcher_color,'LineWidth',3);
m1 = plot_vec(gr_median,'color',L2_color,'LineWidth',3);

%legend([h2,h3,m1],{'Flag/Extrinsic mean','Karcher mean','L2-median'},'Location','SouthEast');
%legend([m2,h2,h3,m1],{'Extrinsic mean','Flag mean','Karcher mean','L2-median'},'Location','SouthEast');
legend([m2,h2,h3,m1],{'GMEB center','Flag/Extrinsic mean','Karcher mean','L2-median'},'Location','SouthEast');

if strcmp(mean_type,'legend')
    fn = 'plot';%filename
    l = 1:length(fn);
    fn = sprintf('%s_LEGEND',fn(l));%give it a separate filename
    a = get(gcf,'children');% link to ledgend is now a(1)
    b = get(gca,'children');% link to the data curve/s
    set(b,'visible','off'); %hide data...
    set(a(2),'visible','off'); %hide axes etc...
    legfs = get(a(1),'Fontsize'); %get legend fontsize
    set(a(1),'Fontsize',legfs+1); %make legend appear larger
end


data.flagmean = flag_mean;
data.karchermean = karcher_mean;
data.gr_median = gr_median;
data.exmean = ex_mean;
data.gmeb = gmeb_mean;

%% Called function
function [handle] = plot_vec(vec,varargin)
% SYNTAX:   [handle] = plot_vec(vec,varargin)
%
% NOTES:    This function just plots unit vectors in 2D or 3D so I am not 
%           going to document it.
%
%           Vectors need to be in R^2 or R^3.
%
%           'varagin' should contain normal Matlab plotting options.
%
% LAST EDITED: by T. Marrinan 06/21/2013
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
b = size(vec);
line = -1:0.1:1;
d = length(line);
newvec = zeros(b(1),d);
for i = 1 : d
    newvec(:,i) = line(i)*vec;
end
newvec = newvec';
if b(1) == 2
    handle = plot(newvec(:,1),newvec(:,2),varargin{:});
end


if b(1) == 3
    handle = plot3(newvec(:,1),newvec(:,2),newvec(:,3),varargin{:});
end






