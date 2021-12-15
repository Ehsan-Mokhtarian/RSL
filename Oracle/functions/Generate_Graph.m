function A = Generate_Graph(graph_type,n,p)


if strcmp(graph_type, 'Erdos-Renyi')
    A = rand(n)<p;
    A = triu(A,+1);
elseif strcmp(graph_type, 'alarm')
    f = matfile('alarm.mat');
    A = f.A;
elseif strcmp(graph_type, 'andes')
    f = matfile('andes.mat');
    A = f.A;
elseif strcmp(graph_type, 'asia')
    f = matfile('asia.mat');
    A = f.A;
elseif strcmp(graph_type, 'barley')
    f = matfile('barley.mat');
    A = f.A;
elseif strcmp(graph_type, 'cancer')
    f = matfile('cancer.mat');
    A = f.A;
elseif strcmp(graph_type, 'carpo')
    f = matfile('carpo.mat');
    A = f.A;
elseif strcmp(graph_type, 'diabetes')
    f = matfile('diabetes.mat');
    A = f.A;
elseif strcmp(graph_type, 'hailfinder')
    f = matfile('hailfinder.mat');
    A = f.A;
elseif strcmp(graph_type, 'hepar2')
    f = matfile('hepar2.mat');
    A = f.A;
elseif strcmp(graph_type, 'insurance')
    f = matfile('insurance.mat');
    A = f.A;
elseif strcmp(graph_type, 'Link')
    f = matfile('Link.mat');
    A = f.A;
elseif strcmp(graph_type, 'mehra')
    f = matfile('mehra.mat');
    A = f.A;
elseif strcmp(graph_type, 'mildew')
    f = matfile('mildew.mat');
    A = f.A;
elseif strcmp(graph_type, 'munin1')
    f = matfile('munin1.mat');
    A = f.A;
elseif strcmp(graph_type, 'munin2')
    f = matfile('munin2.mat');
    A = f.A;
elseif strcmp(graph_type, 'munin3')
    f = matfile('munin3.mat');
    A = f.A;
elseif strcmp(graph_type, 'munin4')
    f = matfile('munin4.mat');
    A = f.A;
elseif strcmp(graph_type, 'pathfinder')
    f = matfile('pathfinder.mat');
    A = f.A;
elseif strcmp(graph_type, 'pigs')
    f = matfile('pigs.mat');
    A = f.A;
elseif strcmp(graph_type, 'sachs')
    f = matfile('sachs.mat');
    A = f.A;
elseif strcmp(graph_type, 'water')
    f = matfile('water.mat');
    A = f.A;
elseif strcmp(graph_type, 'win95pts')
    f = matfile('win95pts.mat');
    A = f.A;
else
    disp('Error: Structure not found!')
    A = zeros(n);
end
end

