%% this function computes the Squared Exponential Cov martix given the
% two input space vectors

%% inputs: 
% x - the observed input space vector (column vector with n elements)
% x_star - the predictive input space vector (column vector with n_star elements)
% theta - the hyperparameters of the SE GP. It contains (sigma_f, l, sigma_n)

%% outputs:
% K - a cov matrix (dimension: n*n_star)

function K = covMatrixSE(x,x_star,theta)
n = size(x,1);
n_star = size(x_star,1);

% covariance function
k = @squaredExp;
K = zeros(n,n_star);
for i = 1:n_star
    for j = 1:n
        K(j,i) = k(x(j,:),x_star(i,:),theta(1),theta(2));
%         % where to add on the noise term (?)
%         if (n == n_star && i == j)
%             K(j,i) = K(j,i) + theta(3)^2;
%         end
    end
end


%% define the Squared Exp covariance function
function y = squaredExp(x1,x2,sigma_f,l)
y = sigma_f^2 * exp(-(norm(x1 - x2))^2/(2*l^2));
end
end