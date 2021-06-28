%% bulit FDC
%%% built connectivity distance and FDC based on functional connectivity and geodesic diatance
%%% Brainspace toolbox is necessary: https://github.com/MICA-MNI/BrainSpace 
%%% parcellation: Schaefer-200; number of nodes: 200; number of subjects: 50
%%% NOTE: you need to modify lines 11 to match your data directories
%%% self-defined parameter for number of parcells and subjects: line 16 and 17


clear;clc; close all
% load FC and geodesic distance, parameter seetings
Path     = strcat('yourpath\code available');
load(strcat(Path,'\data\FC_z_schaefer200_all.mat'));
load(strcat(Path,'\data\gd_schaefer200_all.mat'));
FC_all   = FC_z_all;
FC_group = mean(FC_all,3);
num_parc = 200;
num_sub  = 50;

% load surface template
load(strcat(Path,'\template\fsa5_midsurface_LR.mat'));
parcFSNew        = dlmread(strcat(Path,'\template\schaefer-200_mics.csv'));
[surf_l, surf_r] = split_surfaces(G);

%% bulid CD map
CD_all = zeros(num_parc,num_sub);

% bulit CD for each individual
for ii         = 1:num_sub
    % threshold FC by row
    FC_sub     = FC_all(:,:,ii);
    FC_thre    = zeros(num_parc,num_parc);
    for jj     = 1:num_parc
        FC_row = FC_sub(jj,:);
        p      = prctile(FC_row,90);
        FCbin  = ones(size(FC_row));
        FCbin(FC_row < p) = 0;
        FC_thre(jj,:)     = FCbin;
    end
    % fill inter-hemisphere connections in GD Matrix
    GD_sub                = gd_all(:,:,ii);
    GD_fill               = (GD_sub(1:100,1:100)+GD_sub(101:200,101:200))/2;
    GD_sub(1:100,101:200) = GD_fill;
    GD_sub(101:200,1:100) = GD_fill;
    % build CD by combining FC and GD
    GD_FCthr     = FC_thre.*GD_sub;
    CD_map       = mean(GD_FCthr,2);
    CD_all(:,ii) = CD_map;
end

% vasualization: add node 1 and node 102 as the mid-wall mask 
group_CD           = mean( CD_all,2);
CD_view            = zeros(num_parc+2,1);
CD_view(2:101,:)   = group_CD(1:100,:);
CD_view(103:202,:) = group_CD(101:200,:);
plot_hemispheres(CD_view,{surf_l, surf_r}, ...
    'parcellation', parcFSNew, ...
    'labeltext',{'CD '});

%% bulid group-level FDC using group-level FC and group-level CD
FC_thre    = zeros(num_parc,num_parc);
    for jj = 1:num_parc
        FC_row             = FC_group(jj,:);
        p                  = prctile(FC_row,90);
        FC_row(FC_row < p) = 0;
        FC_thre(jj,:)      = FC_row;
    end
    
% correlate FC with CD
FC       = FC_thre';
[rh,ph]  = corr(FC,group_CD);
FDCthre  = rh(:,1);

% vasualization: add node 1 and node 102 as the mid-wall mask 
FDC_view = zeros(num_parc+2,1);
FDC_view(2:101,:)   = FDCthre(1:100,:);
FDC_view(103:202,:) = FDCthre(101:200,:);
plot_hemispheres(FDC_view,{surf_l, surf_r}, ...
    'parcellation', parcFSNew, ...
    'labeltext',{'FDC'});
%% bulid group-level FDC using individual-level FC and individual-level CD
% calculate individual-level FDC
FDCthre_all = zeros(num_parc,num_sub);
for ii      = 1:num_sub
    FC_thre = zeros(num_parc,num_parc);
    FC_sub  = FC_z_all(:,:,ii);
    for jj  = 1:num_parc
        FC_row             = FC_sub(jj,:);
        p                  = prctile(FC_row,90);
        FC_row(FC_row < p) = 0;
        FC_thre(jj,:)      = FC_row;
    end
    % correlate FC with CD for each individual
    FC                     = FC_thre';
    [rh,ph]                = corr(FC,CD_all(:,ii));
    FDCthre_all(:,ii)      = rh;         
end

% calculate group-level FDC by averaging individual-level FDC
FDCthre_mean        = mean(FDCthre_all,2);
% vasualization: add node 1 and node 102 as the mid-wall mask 
FDC_view            = zeros(num_parc+2,1);
FDC_view(2:101,:)   = FDCthre_mean(1:100,:);
FDC_view(103:202,:) = FDCthre_mean(101:200,:);
plot_hemispheres(FDC_view,{surf_l, surf_r}, ...
    'parcellation', parcFSNew, ...
    'labeltext',{'FDC'});
