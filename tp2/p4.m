clear all;
close all;
##clc;

addpath ".."
conjugate_gradient_search;

tol = 0.0001;
alphamethod = 'aramijo';
##alphamethod = 'parabolic';
betamethod = 'none';
##betamethod = 'fletcher';
##betamethod = 'biere';
iterlimit = 400;

x0=[10 -1 0.1 4*1]';

global eta = 0.01;
global xp = [1 1 1 1]';
global A = [3 1 0 2];
global b = 1;

function f2=f2(x)
  global xp;
  f2 = 1/2*(x-xp)'*(x-xp);
end

function p=p(x)
  global A;
  global b;
  p = (A*x - b).^2;
end

function f=f(x)
  global eta;
  f = f2(x) + 1/eta*p(x);
end

function g=grf1(x)
  global xp;
  g = x - xp;
end

function g=grp(x)
  global A;
  global b;
  g = 2*(A*x - b)*A';
end

function g=gr(x)
  global eta;
  g = grf1(x) + 1/eta*grp(x);
end

[xmin, fmin, nbiter, iters, SC] = steepest(x0, @f, @gr, 'tol', tol, 'alphamethod', alphamethod, 'betamethod', betamethod, 'iterlimit', iterlimit);

disp(["xmin: ", num2str(xmin')]);
disp(["fmin: ", num2str(fmin)]);
disp(["nbiter: ", num2str(nbiter)]);
disp(["Stop criteria: ", SC]);

xan = xp+A'*inv(A*A')*(b-A*xp);
fan = 1/2*(xan-xp)'*(xan-xp);
dann = sqrt(sum((xan-xmin).^2));

disp(["Analytical result: ", num2str(xan')]);
disp(["Analytical distance: ", num2str(fan)]);
disp(["Distance between analytical and numerical results: ", num2str(dann)]);
