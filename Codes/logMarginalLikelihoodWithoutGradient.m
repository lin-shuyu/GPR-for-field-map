%% compute the log marginal likelihood & the gradient of it
%% inputs:
% theta - hyparameteres for GP to be optimised over
% can only have one input to be optimised over
% K and y have to be internalised
% K - cov matrix from covMatrixSE using the obervation data points (x,y)
% y - obervation values


%% outputs:
% f - log Marginal Likelihood value
% g - gradient of f with respect to theta


function f = logMarginalLikelihoodWithoutGradient(theta)
myVars = {'x_train','y'};
S = load('trainingData.mat',myVars{:});

x_train = S.x_train;
y = S.y;
n = size(y,1);

K = covMatrixSE(x_train,x_train,theta);
K_y = K + theta(3) * eye(n);
f =  0.5*y.'*K_y*y + 0.5*log(det(K_y)) + n/0.5*log(2*pi);

% compute the gradient to speed up (to be complete)
end