grk = gr(x0);
ngrk = norm(grk);
dk = -grk/ngrk;

a = linspace(-1,1,1001);
fa = zeros(1,length(a));
for i = 1:length(a)
  fa(i) = f(x0+a(i)*dk);
end

close all;
hold on;

plot(a, fa);
plot([-0.2, 0.2], f(x0)+0.2*[ngrk, -ngrk]);
