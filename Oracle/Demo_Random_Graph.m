close all
clear
clc

addpath('functions/')

%% Generate random graph and data
n = 100;
r = 0.72;
p = 1/n^r;
number_of_samples = 50*n;
G = Generate_Graph('Erdos-Renyi',n,p);
m = clique_number(G);
V = ones(n,1);
if isDiamondfree(G)
    disp("G is diamond-free")
else
    disp("G is not diamond-free")
end
fprintf("Clique number: %d \n\n",m);

%% Markov boundary discovery
Mb = ComputeMb_oracle(G);

%% Calling RSL functions
[G_RSLD, tests_RSLD, SC_RSLD] = RSL_D_oracle(G, Mb);
[G_RSLW, tests_RSLW, SC_RSLW] = RSL_W_oracle(G,V,Mb,m);

%% Accuracy of the learned graph
report_accuracy("RSL_D",G,G_RSLD,tests_RSLD);
report_accuracy("RSL_W",G,G_RSLW,tests_RSLW);