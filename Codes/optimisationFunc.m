%% to learn about constrained optimisation function in matlab fmincon()
% let's use a convex function y = exp(-((x - mu)/(2*sigma))^2)
% so mu and sigma are the parameters

%% define the function - a 2D Gaussian - which we'll run the optimisation 
% function over with various plots to show it 
% define the parameters
parameter_Optimal = struct('mean', [1 1], 'var', [1 0; 0 2]);
x1 = -2:0.25:4;
x2 = x1;
[X1 X2] = meshgrid(x1,x2);
y = exp(-1/2*((1/parameter_Optimal.var(1,1))*(X1-parameter_Optimal.mean(1)).^2+(1/parameter_Optimal.var(2,2))*(X2-parameter_Optimal.mean(2)).^2));

% a surface plot
figure(1)
surf(X1,X2,y)
colorbar

% a perspective plot using mesh() with colormap
figure(2)
mesh(X1,X2,y)
xlabel('x (ph)')
ylabel('y (ph)')
set(gca,'ydir','reverse');
colorbar

% a contour plot (black and white)
figure(3)
v = 0:0.1:1;
[C,h] = contour(X1,X2,y,v);
axis square
clabel(C,h)
xlabel('x (ph)')
ylabel('y (ph)')
set(gca,'ydir','reverse');
set(gca,'XAxisLocation','top');
colormap([0 0 0]);



%% Optimisation
fun = @(x) -exp(-1/2*((1/parameter_Optimal.var(1,1))*(x(1)-parameter_Optimal.mean(1)).^2+(1/parameter_Optimal.var(2,2))*(x(2)-parameter_Optimal.mean(2)).^2));
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0;0];
ub = [];
x0 = [0.5;0.5];
nonlcon = [];
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm','interior-point');
[x, fval, exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);


