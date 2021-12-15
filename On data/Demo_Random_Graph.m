close all
clear
clc

addpath('functions/')

%% Generate random graph and data
n = 50; 
V = ones(n,1);
r = 0.72;
p = 1/n^r;
number_of_samples = 50*n;
[G,D] = Generate_Data('Erdos-Renyi',number_of_samples,n,p);
m = clique_number(G);
if isDiamondfree(G)
    disp("G is diamond-free")
else
    disp("G is not diamond-free")
end
fprintf("Clique number: %d \n\n",m);

%% Markov boundary discovery
alpha_Mb = 2/n^2;
Mb = ComputeMb_TC(D, alpha_Mb); % Learning  Markov boundaries using TC algorithm

%% Calling RSL functions
alpha = 0.01;
[G_RSLD, tests_RSLD, SC_RSLD] = RSL_D(D, Mb, alpha, alpha_Mb);
[G_RSLW, tests_RSLW, SC_RSLW] = RSL_W(D, V, Mb, alpha, alpha_Mb, m);

%% Accuracy of the learned graph
report_accuracy("RSL_D",G,G_RSLD,tests_RSLD);
report_accuracy("RSL_W",G,G_RSLW,tests_RSLW);