##clear all;
close all;
clc;

addpath ".."
minimization_algorithms;
res_dir = "results";

% SETUP

global p = 100;

tol = 0.00000001;
alphamethod = 'parabolic';
betamethod = 'none';
newtonmethod = 'bfgs';
iterlimit = 400;

x0=[0.5,0]';
s = 0.0005;

function f=banane(x)
  global p;
  x1 = x(1,:);
  x2 = x(2,:);
  f = (x1 - 0.5).^2 + p*(x1.^2 - x2).^2;
end

function gr=gr(x)
  global p;
  x1 = x(1,:);
  x2 = x(2,:);
  gr1 = 2*(x1-0.5) + 4*p*(x1.^2 - x2).*x1;
  gr2 = -2*p*(x1.^2 - x2);
  gr = [gr1; gr2];
end

% VISUALIZATION

##x1 = linspace(-2,2,100);
##x2 = linspace(-2,4,100);
x1 = linspace(0.3,0.6,1001);
x2 = linspace(-0.05,0.5,1001);
[X1, X2] = meshgrid(x1, x2);
X_shape = size(X1);

F = banane([X1(:)';X2(:)']);
F = reshape(F, X_shape);

##H = figure;
##hold on;
##colorbar;
##
##d0 = -gr(x0);
##quiver(x0(1), x0(2), s*d0(1), s*d0(2), 'linewidth', 3);

% SOLUTION

[xmin, fmin, nbiter, iters, CONVCRIT] = minimize(x0, @banane, @gr, 'tol', tol, 'alphamethod', alphamethod, 'betamethod', betamethod, 'iterlimit', iterlimit, 'newtonmethod', newtonmethod);

H = figure;
hold on;
plot([iters.x](1,:), [iters.x](2,:), 'bo-');
labels = {};
for i = 1:length(iters)
  labels(end+1) = ['  x_{', num2str(i-1), '}'];
end
labels(2:end-1) = " ";
text([iters.x](1,:), [iters.x](2,:), labels, 'horizontalalignment',"left", 'verticalalignment',"bottom", 'fontsize',35);
contour(x1, x2, F, logspace(log10(fmin),log10(max(max(F))+1),15));
##axis equal;

##title(["x_f: ", num2str(xmin'), "\ncost_f: ", num2str(fmin), "\nNb-iter: ", num2str(nbiter)]);
grid on;
saveas(H, [res_dir, filesep, "tp1-p1-bfgs-parabolic-path"], "png");

disp(["xmin: ", num2str(xmin(1)),", ", num2str(xmin(2))]);
disp(["fmin: ", num2str(fmin)]);
disp(["nbiter: ", num2str(nbiter)]);
disp(["Stop criteria: ", CONVCRIT]);

result.name = "Parabolic, \\textit{BFGS}";
result.solution = xmin;
result.cost = fmin;
result.iters = nbiter;
results(end+1) = result;
