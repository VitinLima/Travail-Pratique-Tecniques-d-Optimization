clear all;
close all;
##clc;

addpath "..";
dfp_algorithm;

tol = 1e-6;
##alphamethod = 'aramijo';
alphamethod = 'parabolic';
betamethod = 'none';
##betamethod = 'fletcher';
##betamethod = 'biere';
iterlimit = 1000;

x0=[2, 1, 2]';

global eta = 0.01;
global A = [5 3 1; 3 2 1; 1 1 4];

function f3=f3(x)
  global A;
  f3 = x'*A*x;
end

function p=p(x)
  x1 = x(1);
  x2 = x(2);
  x3 = x(3);

  p1 = (x1 + x2 - x3).^2;
  p1(x1 + x2 < x3) = 0;

  p2 = (3*x1 - 2*x2 - x3).^2;

  p = p1+p2;
end

function f=f(x)
  global eta;
  f = f3(x) + 1/eta*p(x);
end

function g=grf3(x)
  global A;
  gx1 = A(1,:)*x + x'*A(:,1);
  gx2 = A(2,:)*x + x'*A(:,2);
  gx3 = A(3,:)*x + x'*A(:,3);
  g = [gx1, gx2, gx3]';
end

function g=grp(x)
  x1 = x(1);
  x2 = x(2);
  x3 = x(3);

  g1 = 2*(x1 + x2 - x3)*[1 1 -1]';
  g1(x1 + x2 < x3) = 0;

  g2 = 2*(3*x1 - 2*x2 - x3)*[3 -2 -1]';

  g = g1 + g2;
end

function g=gr(x)
  global eta;
  g = grf3(x) + 1/eta*grp(x);
end

[xmin, fmin, nbiter, iters, SC] = dfp_algorithm(x0, @f, @gr, 'tol', tol, 'alphamethod', alphamethod, 'betamethod', betamethod, 'iterlimit', iterlimit);

disp(["xmin: ", num2str(xmin')]);
disp(["fmin: ", num2str(fmin)]);
disp(["nbiter: ", num2str(nbiter)]);
disp(["Stop criteria: ", SC]);

##xan = xp+A'*inv(A*A')*(b-A*xp);
##fan = 1/2*(xan-xp)'*(xan-xp);
##dann = sqrt(sum((xan-xmin).^2));
##
##disp(["Analytical result: ", num2str(xan')]);
##disp(["Analytical distance: ", num2str(fan)]);
##disp(["Distance between analytical and numerical results: ", num2str(dann)]);
