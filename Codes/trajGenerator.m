%% Trajectory generator
% This function generates synthetic trajectories for a squared 2D space of
% size l * l (sqm)

%% inputs:
% l - the length of the side of the squared space
% n_select - the number of trajectories to be taken as the input
% measurements for the later regression; here by default we keep n_select
% as 4, one trajectory horizontal, vertical, diagonal and sinoidally
% curved.


%% outputs:
% selected_pts_nonRepeated - the points along the selected trajectories
% which will be taken as the location for the input points to GP for
% regression. Any repeated points should have already been removed. 
% n_select - as an output rather here, as it will be needed later
% traj_select_Index - the index of the Trajectories that have been selected
% Trajectory - the cell array that stores all the trajectories that has been generated 

function [selected_pts_nonRepeated, n_select, Trajectory, traj_select_Index] = trajGenerator(l)
l = floor(l);  % make sure l is an integer

%% Define the trajectories
% generate trajectories - imitate the sensor path
% line trajectories
% horiztal lines
n_hori = floor(l/1) + 1;
for p = 1:n_hori
    Trajectory{p} = 0:1:l;
    n_pathL = size(Trajectory{p},2);
    Trajectory{p} = [Trajectory{p};zeros(3,n_pathL)];
    Trajectory{p}(2,:) = (p-1)*ones(1,n_pathL); 
    for i = 1:n_pathL
        x_temp = Trajectory{p}(1:2,i);
    end 
end

% vertical lines
n_vert = floor(l/1) + 1;
for p = 1:n_vert
    n_pathL = length(0:1:l);
    Trajectory{n_hori+p} = zeros(4,n_pathL);
    Trajectory{n_hori+p}(2,:) = 0:1:20;
    Trajectory{n_hori+p}(1,:) = (p-1)*ones(1,n_pathL); 
    for i = 1:n_pathL
        x_temp = Trajectory{n_hori+p}(1:2,i);
    end 
end

% diagonals
Index_start = size(Trajectory,2);
n_pathL = length(0:1:l);
for p = 1:2
    Trajectory{Index_start+p} = zeros(4,n_pathL);
    if p == 1
        Trajectory{Index_start+p}(1:2,:) = [0:1:l;0:1:l];
    else
        Trajectory{Index_start+p}(1,:) = 0:1:l;
        Trajectory{Index_start+p}(2,:) = l:-1:0;
    end
    for i = 1:n_pathL
        x_temp = Trajectory{Index_start+p}(1:2,i);
    end 
end

% curved trajectory 
Index_start = size(Trajectory,2);
n_pathL = length(0:0.5:l);
for p = 1:floor(l/2)
    Trajectory{Index_start+p} = zeros(4,n_pathL);
    Trajectory{Index_start+p}(2,:) = [0:0.5:l];
    Trajectory{Index_start+p}(1,:) = (l/2)*(1+sin(0.1*p*Trajectory{Index_start+p}(2,:)));
    for i = 1:n_pathL
        x_temp = Trajectory{Index_start+p}(1:2,i);
    end 
end

%% select n_select trajectories to use - randomly use which one
n_traj = size(Trajectory,2);
rng('shuffle');
n_select = 4;
traj_select_Index = zeros(4,1);
traj_select_Index(1) = randi([1 21],1,1); % generate one horizontal traj
traj_select_Index(2) = randi([22 42],1,1); % generate one vertical traj
traj_select_Index(3) = randi([43 44],1,1); % generate one diagonal traj
traj_select_Index(4) = randi([44 54],1,1); % generate one curved traj

% randperm(n_traj,n_select); % generate non-repeated integer values

% keep track of the input space locations of these measurements along the selected trajectories
% keep a record of all data points along the selected trajectories
n_pts = 21*3+41;
selected_pts_all = zeros(n_pts,3);
index_counter = 1;
for i = 1:n_select
    traj_pts = size(Trajectory{traj_select_Index(i)},2);
    selected_pts_all(index_counter:index_counter+traj_pts-1,1) = traj_select_Index(i) * ones(traj_pts,1);
    selected_pts_all(index_counter:index_counter+traj_pts-1,2:3) = Trajectory{traj_select_Index(i)}(1:2,:).';
    index_counter =  index_counter + traj_pts;
end

% remove the repeated points and get a list of unique locations 
selected_pts_nonRepeated = selected_pts_all;
index_counter = 1;
for i = 1:n_select
    if i == 1
        index_counter = index_counter + size(Trajectory{traj_select_Index(i)},2);
        repeat_counter = index_counter;
    else
        for p = index_counter:n_pts
            pt_to_check = selected_pts_all(p,2:3).';
            x_index = find(selected_pts_all(1:p-1,2) == pt_to_check(1));
            y_index = find(selected_pts_all(1:p-1,3) == pt_to_check(2));
            if (isempty(x_index) && isempty(y_index))
                for m = 1:length(x_index)
                    repeat_index = find(y_index == x_index(m));
                end
            else 
                repeat_index = [];
            end
            if ~isempty(repeat_index)
                selected_pts_nonRepeated(repeat_counter,:) = [];
                repeat_counter = repeat_counter - 1;
            end
            repeat_counter = repeat_counter + 1;
        end
    end
end
end