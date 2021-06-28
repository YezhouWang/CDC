%% gadients estimation, using qT1 mpc data as example
%%% BrainSpace toolbox is necessary: https://github.com/MICA-MNI/BrainSpace 
%%% For more information about BrainSpace, please see our tutorial:
%%% https://brainspace.readthedocs.io/en/latest/index.html
%%% parcellation: Schaefer-200; number of nodes: 200
%%% NOTE: you need to modify lines 12 to match your data directories
%%% self-defined parameter for number of parcells: line 18
%%% line 19,20,21 need to be modified to match your parcellation

clear;clc; close all
% load surface atlas, template, and qT1 mpc matrix (please see bulid_mpc.m)
Path             = strcat('yourpath\code available');
load(strcat(Path,'\template\fsa5_midsurface_LR.mat'));
parcFS           = dlmread(strcat(Path,'\template\schaefer-200_mics.csv'));
[surf_l, surf_r] = split_surfaces(G);
load(strcat(Path,'\data\qT1_MPC_schaefer200_all.mat'));
% define parcellations and mid masks
parcell_Num      = 200;
mask_label_1     = 1;
mask_label_2     = 102;
gradient_Num     = 13;
% get the group-level mpc matrix and remove the mask
group_mpc=mean(mpc_matrix_all,3);
group_mpc([mask_label_1 mask_label_2], : ) = [];
group_mpc(:, [mask_label_1 mask_label_2])  = [];
% calculate the template gradient from group-level mpc matrix
mpc_gradient    = GradientMaps('n_components', gradient_Num);
mpc_gradient    = mpc_gradient.fit(group_mpc);
Gref            = mpc_gradient.gradients{1};

%calculate the gradients of individual-level mpc matrices
Gradient_all_aligned  = zeros(parcell_Num+2,gradient_Num,size(mpc_matrix_all,3));
for i                 = 1:size(mpc_matrix_all,3)
    % remove mask for each individual
    sub_mpc_matrix                                  = mpc_matrix_all(:,:,i);
    sub_mpc_matrix([mask_label_1 mask_label_2], : ) = [];
    sub_mpc_matrix(:, [mask_label_1 mask_label_2])  = [];
    
    % calculate gradients for each individual and aligh them to template
    mpc_gradient2    = GradientMaps('n_components', gradient_Num,'alignment','pa');
    mpc_gradient2    = mpc_gradient2.fit(sub_mpc_matrix,'reference',Gref);
    gradient         = mpc_gradient2.aligned{1}(:,:);  
    
    % vasualization: add node 1 and node 102 as the mid-wall mask 
    Gradient_sub                = zeros(parcell_Num+2,gradient_Num);
    Gradient_sub(2:101,:)       = gradient(1:100,:);
    Gradient_sub(103:202,:)     = gradient(101:200,:);
    Gradient_all_aligned(:,:,i) = Gradient_sub;
end

% calculate the group-level gradient of all individual-level gradients
Gradient_average_align=mean(Gradient_all_aligned,3);

% visualization of gradient 1 and gradient 2
% gradients could be flipped for consistency with gradients of other scales
Gs      = Gradient_average_align(:,1:2);
Gs_flip = Gs.*-1;
plot_hemispheres(Gs_flip,{surf_l, surf_r}, ...
    'parcellation', parcFS, ...
    'labeltext',{'Gradient 1','Gradient 2'});

