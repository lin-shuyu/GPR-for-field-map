%% GPML - Chapter 2 Regression Ex1
% To generate the random candidate functions from a prior GP model
clear all
close all

% the number of sample functions to be generated in this demo
n_prior = 3;
n_post = 3;

%% prior GP model
% define the step sizes to generate the equi-distant input values
step = zeros(n_prior+2,1);
step(1:2) = [0.2;0.1];
step(3:end) = 0.01*ones(n_prior,1);

for n = 1:n_prior+2
    % input values - x
    x = -5:step(n):5;
    x = x(:);
    n_pts = length(x);
    
    % mean function 
    m = zeros(n_pts,1);
    cov = computeCov(x,x);
    eps = 10^(-5); % to ensure cov is positive definite

    cov =  cov +  eps * eye(n_pts);
    s2 = ones(n_pts,1);

    % generate samples according to the distribution defined by GP - N(m,cov)
    L = chol(cov,'lower'); % Cholesky decomposition such that L * L' = cov
    u = diag(normrnd(zeros(n_pts),eye(n_pts)));
    f = m + L * u;
    
    % plot the result 
    figure(2)
    hold on
    % plot the 2*sd
    yupper = m + 2*sqrt(s2);
    ylower = m - 2*sqrt(s2);

    if n==1
        set(fill([x; flipud(x)],[yupper; flipud(ylower)],[0.9 0.9 0.9]),'EdgeColor',[0.9 0.9 0.9])
        alpha(0.7)
        title('Prior Distribution of functions (Zero mean, Squared Exp cov function)')
    end

    ylim([-3,3])
    if n<=2
        plot(x,f,'LineStyle','none','Marker','.')
    else
        plot(x,f,'LineStyle','-')
    end
end



%% posterior distribution with obervation added 
% generate some artificial measurments
n_obs = 5;
obs = zeros(n_obs,2);
obs(:,1) = [-4;-3;-1;0;2]; % column vector
obs(:,2) = [-2;0;1;2;-1];  % column

%  compute the cov matrix with obs and prior
k = computeCov(obs(:,1),obs(:,1));
k_star = computeCov(obs(:,1),x);
k_star2 = computeCov(x,x);

% infer the posterior function values by conditioning the joint Gaussian
% prior on the observations f_star
cov_post = k_star2 - k_star.' * pinv(k) * k_star;
cov_post = cov_post + eps * eye(n_pts);
m_post = k_star.' * pinv(k) * obs(:,2);
L_post = chol(cov_post,'lower'); % Cholesky decomposition such that L * L' = cov

for n = 1:n_post
    u_post = diag(normrnd(zeros(n_pts),eye(n_pts)));
    f_star = m_post + L_post * u_post;

    % plot the posterior distribution
    figure(3)
    hold on
    % plot the 2*sd
    s2_post = diag(cov_post);
    yupper = m_post + 2*sqrt(s2_post);
    ylower = m_post - 2*sqrt(s2_post);
    if n==1
        set(fill([x; flipud(x)],[yupper; flipud(ylower)],[0.9 0.9 0.9]),'EdgeColor',[0.9 0.9 0.9])
        alpha(0.7)
        figure(3)
        plot(obs(:,1),obs(:,2),'LineStyle','none','Marker','+','MarkerEdgeColor','k')
        title('Posterior Distribution of functions with 5 observations added to the prior GP model')
    end

    ylim([-3,3])
    plot(x,f_star,'LineStyle','-') 
end

%% define a function which outputs the covariance matrix given 
% the two input space vectors

% X and X_star should be devised in a way that 
% each row vector is an input point x


function cov = computeCov(X,X_star)
% get the dimension of the cov matrix
n_train = size(X,1);
n_test = size(X_star,1);

% compute K(X,X_star)
k = @squaredExp;
cov = zeros(n_train,n_test);
for i = 1:n_test
    for j = 1:n_train
        cov(j,i) = k(X(j,:),X_star(i,:));
    end 
end
end


%% define the Squared Exp covariance function
function y = squaredExp(x1,x2)
sigma_f = 1;
l = 1;
y = sigma_f^2 * exp(-(norm(x1 - x2))^2/(2*l^2));
end
