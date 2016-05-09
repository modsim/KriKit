clear
clc
close all
%%
mex determineParetoSet_Mex.cpp determineParetoPoints.cpp


vec = 1:10;
dataset = [vec',sort(vec,'descend')'];
dataset = [dataset;10,10;1,9;6,7];
dataset2 = [10,10];
% figure()
% plot(dataset(:,1),dataset(:,2),'*')

% paretoOptimal_nDimensional([])
paretoOptimal_nDimensional(dataset)
paretoOptimal_nDimensional(dataset2)

determineParetoSet(dataset)
determineParetoSet(dataset2)
determineParetoSet([])
determineParetoSet
determineParetoSet(ones(3,1))

dataSetRandom = rand(1e4,3);
tic
paretoSetCheck = paretoOptimal_nDimensional(dataSetRandom);
time1 = toc
tic
paretoSetTest = determineParetoSet(dataSetRandom);
time2 = toc;
time1/time2

max(max(abs(paretoSetCheck-paretoSetTest)))
