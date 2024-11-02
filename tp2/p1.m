clear all;
close all;
##clc;

addpath ".."
conjugate_gradient_search;

tol = 0.001;
alphamethod = 'aramijo';
##alphamethod = 'parabolic';
betamethod = 'none';
##betamethod = 'fletcher';
##betamethod = 'biere';
iterlimit = 7;

x0=[0.5; 1];

global eta = 0.01;

function f1=f1(x)
  f1 = 2*x(1,:).^2 + 3*x(1,:).*x(2,:) + 2*x(2,:).^2;
end

function p=p(x)
  x(1,:) += 0.5;
  x(2,:) += 0.5;
  x(1,x(1,:)<0) = 0;
  x(2,x(2,:)<0) = 0;
  p = sum(x.*x, 1);
end

function f=f(x)
  global eta;
  f = f1(x) + 1/eta*p(x);
end

function g=grf1(x)
  g1 = 4*x(1,:) + 3*x(2,:);
  g2 = 3*x(1,:) + 4*x(2,:);
  g = [g1; g2];
end

function g=grp(x)
  x(1,:) += 0.5;
  x(2,:) += 0.5;
  x(1,x(1,:)<0) = 0;
  x(2,x(2,:)<0) = 0;
  g1 = 2*x(1,:);
  g2 = 2*x(2,:);
  g = [g1; g2];
end

function g=gr(x)
  global eta;
  g = grf1(x) + 1/eta*grp(x);
end

[xmin, fmin, nbiter, iters, CONVCRIT] = steepest(x0, @f, @gr, 'tol', tol, 'alphamethod', alphamethod, 'betamethod', betamethod, 'iterlimit', iterlimit);

disp(["xmin: ", num2str(xmin')]);
disp(["fmin: ", num2str(fmin)]);
disp(["nbiter: ", num2str(nbiter)]);
disp(["Stop criteria: ", CONVCRIT]);

x1 = linspace(-1,.6,101);
x2 = linspace(-1,1.2,101);

[XX1, XX2] = meshgrid(x1, x2);
xx1 = XX1(:)';
xx2 = XX2(:)';
F = f([xx1;xx2]);
F = reshape(F, size(XX1));

minF = min(min(F));
maxF = max(max(F));

% 2 dimensional plot

figure;
hold on;
plot([iters.x](1,:), [iters.x](2,:), 'bo-');
labels = {};
for i = 0:nbiter
  labels(end+1) = ['  x_{', num2str(i), '}'];
end
text([iters.x](1,:), [iters.x](2,:), labels, 'horizontalalignment',"left", 'verticalalignment',"bottom", 'fontsize',15);

contour(x1, x2, F, logspace(-4,log10(maxF-minF),25)+minF);
colorbar;

x0 = [-0.75,-0.75]';
gr0 = gr(x0);
s = 0.0005;
quiver(x0(1), x0(2), s*gr0(1), s*gr0(2), 'r')
##axis equal;

% 3 dimensional plot

figure;
hold on;
surface(x1, x2, F, 'linestyle', 'none');
plot3([iters.x](1,:), [iters.x](2,:), [iters.f](1,:), 'bo-');
labels = {};
for i = 0:nbiter
  labels(end+1) = ['  x_{', num2str(i), '}'];
end
text([iters.x](1,:), [iters.x](2,:), [iters.f], labels, 'horizontalalignment',"left", 'verticalalignment',"bottom", 'fontsize',15);
colorbar;
