clear all;
close all;
clc;

result.name = [];
result.solution = [];
result.cost = [];
results = {}

aramijo;
parabolic;
aramijo_fletcher;
aramijo_ribiere;
parabolic_fletcher;
parabolic_ribiere;
dpf_parabolic;
bfgs_parabolic;

for i = 1:length(results)
  result = cell2mat(results(i));
  if i==length(results)
    disp(["\t", result.name, " & [", strjoin(strsplit(num2str(result.solution')),', '), "] & ", num2str(result.cost), " & ", num2str(result.iters)]);
  else
    disp(["\t", result.name, " & [", strjoin(strsplit(num2str(result.solution')),', '), "] & ", num2str(result.cost), " & ", num2str(result.iters), " \\\\"]);
  endif
end
