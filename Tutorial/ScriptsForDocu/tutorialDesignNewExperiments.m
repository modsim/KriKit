%% Tutorial: How to use KriKit for Optimization
%
% Ingore the next three line:
clc
clear
close all

%% Test Function
y1 = @(x)4*x(:,1).^2 + 4*x(:,2).^2;
y2 = @(x)(x(:,1)-5).^2 + (x(:,2)-5).^2;

nDataFull = 100;
XFull = createNDGRID([0,0],[5,3],nDataFull);
Y1Full = y1(XFull);
Y2Full = y2(XFull);

figure
plot(Y1Full,Y2Full,'*')

figure
surf(unique(XFull(:,1)),unique(XFull(:,2)),reshape(Y1Full,nDataFull,nDataFull)')
shading(gca,'interp')
xlabel('x1')
ylabel('x2')
set(gca,'FontSize',20)

figure
surf(unique(XFull(:,1)),unique(XFull(:,2)),reshape(Y2Full,nDataFull,nDataFull)')
shading(gca,'interp')
campos([26.2542   28.9225  259.8676])
xlabel('x1')
ylabel('x2')
set(gca,'FontSize',20)

% XTest = createNDGRID([0,0],[5,3],3);
XTest = [0,0;...
         2.5,0;...
         5,0;...
         0,1.5;...
         2.5,1.5;...
         0,3];
Y1Test = y1(XTest);
Y2Test = y2(XTest);

figure
plot(Y1Test,Y2Test,'*')

%% Contruct Kriging Models
krigingObj = BayesianOptimizationClass;
krigingObj.addKrigingObject(1,'Y1')
krigingObj.addKrigingObject(1,'Y2')
nObj = 2;

for iObj = 1:nObj
    krigingObj.KrigingObjects{iObj}.setInputData(XTest)
end
krigingObj.KrigingObjects{1}.setOutputData(Y1Test)
krigingObj.KrigingObjects{2}.setOutputData(Y2Test)

for iObj = 1:nObj
    krigingObj.KrigingObjects{iObj}.setCovariogramModelChoice(6)
    krigingObj.KrigingObjects{iObj}.setUseMatlabRegressionGP(true)
    krigingObj.KrigingObjects{iObj}.generateRegressionGPModel
end

%% Single Objective Optimization
figure
hold
surf(unique(XFull(:,1)),unique(XFull(:,2)),reshape(Y1Full,nDataFull,nDataFull)')
shading(gca,'interp')
campos([26.2542   28.9225  259.8676])
xlabel('x1')
ylabel('x2')
set(gca,'FontSize',20)
plot3(XTest(:,1),XTest(:,2),Y1Test,'ko','MarkerFaceColor','r');
campos([-33.0915   -8.0370  640.7167])

nIteration = 10;
for iIteration = 1:nIteration
    krigingObj = BayesianOptimizationClass;
    krigingObj.addKrigingObject(1,'Y1')
    krigingObj.KrigingObjects{1}.setInputData(XTest)
    krigingObj.KrigingObjects{1}.setOutputData(Y1Test)

    krigingObj.KrigingObjects{1}.setCovariogramModelChoice(6)
    krigingObj.KrigingObjects{1}.setUseMatlabRegressionGP(true)
    krigingObj.KrigingObjects{1}.generateRegressionGPModel
    
    % Design 10 new experiments based on Markov Chain and Expected
    % Improvement
    krigingObj.setnNewSamples(10)
    % Calculate a Markov-Chain with 1e4 links
    krigingObj.setnMCMCLinks(1e4)
    % Ignore first 1e3 links( algorithm might need some time until it converges)
    krigingObj.setnCutLinks(1e3) 
    newSamplePoint = krigingObj.calcNewSamplesViaMCMC(1,'DRAM');
end
