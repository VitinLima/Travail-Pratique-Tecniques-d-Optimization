clear all;
close all;
clc;

addpath ".."
conjugate_gradient_search;
res_dir = "results";

% SETUP

global p = 100;

tol = 0.001;
alphamethod = 'aramijo';
betamethod = 'none'; % Gradient à pas optimal
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

x1 = linspace(-2,2,100);
x2 = linspace(-2,4,100);
[X1, X2] = meshgrid(x1, x2);
X_shape = size(X1);

F = banane([X1(:)';X2(:)']);
F = reshape(F, X_shape);

H = figure;
hold on;
surface(X1, X2, F, 'linestyle', 'none');
grid on;

view(25, 45);
saveas(H, [res_dir, filesep, 'Banane'], 'png');
