clear all;
close all;
##clc;

addpath ".."
conjugate_gradient_search;

tol = 0.00001;
alphamethod = 'aramijo';
##alphamethod = 'parabolic';
betamethod = 'fletcher';
##betamethod = 'biere';
iterlimit = 400;

x0=[0.5];

global eta = 0.01;

function f0=f0(x)
  f0 = 1 + x + 1/3*x.*x.*x;
end

function p=p(x)
  p = min([x; zeros(1,length(x))], [], 1);
  p .*= p;
end

function f=f(x)
  global eta;
  f = f0(x) + 1/eta*p(x);
end

function g=grf0(x)
  g = 1 + x.*x;
end

function g=grp(x)
  g = 2*x;
  g(x>0) = 0;
end

function g=gr(x)
  global eta;
  g = grf0(x) + 1/eta*grp(x);
end

x = linspace(-0.3,1,1001);

##eta = 1;
##figure;
##hold on;
##F = f1(x);
##plot(x, F, 'displayname', "f0");
##for i = 1:4
##  F = f(x);
##  plot(x, F, 'displayname', ["fp, e=", num2str(eta)]);
##  eta /= 10;
##end
##legend;

##figure;
##plot(x, gr(x));

[xmin, fmin, nbiter, iters, CONVCRIT] = steepest(x0, @f, @gr, 'tol', tol, 'alphamethod', alphamethod, 'betamethod', betamethod, 'iterlimit', iterlimit);

disp(["xmin: ", num2str(xmin)]);
disp(["fmin: ", num2str(fmin)]);
disp(["nbiter: ", num2str(nbiter)]);

labels = {};
for i = 1:nbiter
  labels(end+1) = num2str(i);
end

figure;
hold on;
plot([iters.x], [iters.f], '*', 'markersize', 10);
plot(x, grp(x));

text([iters.x], [iters.f], labels, 'horizontalalignment',"left", 'verticalalignment',"bottom", 'fontsize',15);