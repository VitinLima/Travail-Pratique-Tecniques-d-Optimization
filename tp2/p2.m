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

x0=[0.5];

global eta = 0.01;

function f2=f2(x)
  f2 = x.*x + 10*sin(x);
end

function p=p(x)
  xa = x - 2;
  xb = x - 10;
  xa(xa>0) = 0;
  xb(xb<0) = 0;
  p = xa.^2+xb.^2;
end

function f=f(x)
  global eta;
  f = f2(x) + 1/eta*p(x);
end

function g=grf1(x)
  g = 2*x + 10*cos(x);
end

function g=grp(x)
  xa = x - 2;
  xb = x - 10;
  xa(xa>0) = 0;
  xb(xb<0) = 0;
  g = 2*xa+2*xb;
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

x = linspace(0,12,101);

% 2 dimensional plot

figure;
hold on;
plot(x, f(x), 'displayname', "fp");
plot(x, f2(x), 'displayname', "f");
plot([iters.x], [iters.f], 'bo-');
labels = {};
for i = 0:nbiter
  labels(end+1) = ['  x_{', num2str(i), '}'];
end
text([iters.x], [iters.f], labels, 'horizontalalignment',"left", 'verticalalignment',"bottom", 'fontsize',15);
