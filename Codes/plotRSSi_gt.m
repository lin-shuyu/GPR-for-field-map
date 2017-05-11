%% plot the variation of RSSi over 20*20 (m) space
clear all
close all

% define the input space
step = 0.05;
x1 = 0:step:20;
x2 = 0:step:20;
x_AP = [6;6];
e = 5;

% get the function values using the RSSi defined in RSSi.m
n_pts = length(x1);
y = zeros(n_pts);
f = zeros(n_pts);
for i = 1:n_pts
    for j = 1:n_pts
        x_temp = [x1(i);x2(j)];
        [f(j,i),y(j,i)] = RSSi(x_temp,x_AP,e);
    end
end

% plot the 3D plot
figure(1)
mesh(x1,x2,f)
title('Synthetic model of the variation of RSSi over space - Noise free')
xlabel('Position in x (m)')
ylabel('Position in y (m)')

figure(2)
mesh(x1,x2,y)
title('Synthetic model of the variation of RSSi over space - with Gaussian white noise sd=5dBm')
xlabel('Position in x (m)')
ylabel('Position in y (m)')

% to show a slice of the 3D plot to show the comparison between noise-free
% and the function value with noise
figure(3)
Index_slice = x_AP(1)/step - 5;
plot(x1,f(Index_slice,:))
hold on
plot(x1,y(Index_slice,:))
title('A slice of the 3D model - comparison between noise-free and noisy cases')

%% generate trajectories - imitate the sensor path
% line trajectories
% horiztal lines
n_hori = floor(20/1) + 1;
for p = 1:n_hori
    Trajectory{p} = 0:1:20;
    n_pathL = size(Trajectory{p},2);
    Trajectory{p} = [Trajectory{p};zeros(3,n_pathL)];
    Trajectory{p}(2,:) = (p-1)*ones(1,n_pathL); 
    for i = 1:n_pathL
        x_temp = Trajectory{p}(1:2,i);
        [Trajectory{p}(3,i),Trajectory{p}(4,i)] = RSSi(x_temp,x_AP,e);
    end 
end

% vertical lines
n_vert = floor(20/1) + 1;
for p = 1:n_vert
    n_pathL = length(0:1:20);
    Trajectory{n_hori+p} = zeros(4,n_pathL);
    Trajectory{n_hori+p}(2,:) = 0:1:20;
    Trajectory{n_hori+p}(1,:) = (p-1)*ones(1,n_pathL); 
    for i = 1:n_pathL
        x_temp = Trajectory{n_hori+p}(1:2,i);
        [Trajectory{n_hori+p}(3,i),Trajectory{n_hori+p}(4,i)] = RSSi(x_temp,x_AP,e);
    end 
end

% diagonals
Index_start = size(Trajectory,2);
n_pathL = length(0:1:20);
for p = 1:2
    Trajectory{Index_start+p} = zeros(4,n_pathL);
    if p == 1
        Trajectory{Index_start+p}(1:2,:) = [0:1:20;0:1:20];
    else
        Trajectory{Index_start+p}(1,:) = 0:1:20;
        Trajectory{Index_start+p}(2,:) = 20:-1:0;
    end
    for i = 1:n_pathL
        x_temp = Trajectory{Index_start+p}(1:2,i);
        [Trajectory{Index_start+p}(3,i),Trajectory{Index_start+p}(4,i)] = RSSi(x_temp,x_AP,e);
    end 
end

% curved trajectory 
Index_start = size(Trajectory,2);
n_pathL = length(0:0.5:20);
for p = 1:10
    Trajectory{Index_start+p} = zeros(4,n_pathL);
    Trajectory{Index_start+p}(2,:) = [0:0.5:20];
    Trajectory{Index_start+p}(1,:) = 10*(1+sin(0.1*p*Trajectory{Index_start+p}(2,:)));
    for i = 1:n_pathL
        x_temp = Trajectory{Index_start+p}(1:2,i);
        [Trajectory{Index_start+p}(3,i),Trajectory{Index_start+p}(4,i)] = RSSi(x_temp,x_AP,e);
    end 
end

%% GPR
% take measurements along selected trajectories to be the training points
n_traj = size(Trajectory,2);
% select n_select trajectories to use - randomly use which one
rng('shuffle');
n_select = 5;
traj_select_Index = randi([1,n_traj],n_select,1);

for i = 1:n_select
    if i == 1
        x_train = Trajectory{traj_select_Index(i)}(1:2,:).';
        y_train = Trajectory{traj_select_Index(i)}(4,:).';
        y_train_gt = Trajectory{traj_select_Index(i)}(3,:).';
        n_sample = size(Trajectory{traj_select_Index(i)},2);
    else
        n_add = size(Trajectory{traj_select_Index(i)},2);
        for j = 1:n_add
            data_to_add = Trajectory{traj_select_Index(i)}(1:2,j).';
            index_Repeat = find(x_train == data_to_add);
            if isempty(index_Repeat) 
                n_sample = n_sample+1;
                x_train(n_sample,:) = data_to_add;
                y_train(n_sample,:) = Trajectory{traj_select_Index(i)}(4,j);
                y_train_gt(n_sample,:) = Trajectory{traj_select_Index(i)}(3,j);
            elseif (length(index_Repeat) ~=2) || (index_Repeat(1) ~= (index_Repeat(2)-n_sample))
                n_sample = n_sample+1;
                x_train(n_sample,:) = data_to_add;
                y_train(n_sample,:) = Trajectory{traj_select_Index(i)}(4,j);
                y_train_gt(n_sample,:) = Trajectory{traj_select_Index(i)}(3,j);
                
    %         else
    %             y_train(index_Repeat(1),:) = mean([y_train(index_Repeat(1),:),Trajectory{traj_select_Index(i)}(4,j)]);
            end
        end
    end
end

% remove the singularity point, which will cause trouble
index_Singular = find(y_train == Inf);
x_train(index_Singular,:) = [];
y_train(index_Singular) = [];
y_train_gt(index_Singular) = [];

% cov matrix
theta = [1,2,5];   % [sigma_f, l, sigma_n]
K = covMatrixSE(x_train,x_train,theta);

% prediction input space locations
pred_grid = 0:0.5:20;
n_grid = length(pred_grid);
x_pred = zeros(n_grid^2,2);
for i = 1:n_grid
    x_pred((i-1)*n_grid+1:i*n_grid,1) = pred_grid.';
    x_pred((i-1)*n_grid+1:i*n_grid,2) = pred_grid(i)*ones(n_grid,1);
end
n_pred = n_grid^2;

% prediction cov matrices
K_star = covMatrixSE(x_train,x_pred,theta);
K_star2 = covMatrixSE(x_pred,x_pred,theta);

% make prediciton 
K_y = K + (theta(3)^2)*eye(length(y_train));
K_inv = pinv(K_y);
GPR_mean = K_star.' * K_inv * y_train;
GPR_cov = K_star2 - K_star.' * K_inv * K_star;

%% plot the 3d map
Z = reshape(GPR_mean,[n_grid,n_grid]);
figure(4)
mesh(pred_grid,pred_grid,Z);