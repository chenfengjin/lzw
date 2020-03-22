%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA101
% Project Title: Implementation of Real-Coded Genetic Algorithm in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
tic
clc;
clear;
close all;
%% Problem Definition

CostFunction=@(x) -1*case39_YC(x);     % Cost Function

nVar=6 * 3;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size 这个变量没用了

VarMin=0;         % Lower Bound of Variables
VarMax= 120;         % Upper Bound of Variables

VarMin2=1;         % Lower Bound of Variables
VarMax2=22;         % Upper Bound of Variables
%% GA Parameters

MaxIt=100;     % Maximum Number of Iterations

nPop=1;       % Population Size

pc=0.7;                 % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parnets)
gamma=0.4;              % Extra Range Factor for Crossover %交叉的范围

pm=0.3;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants
mu=0.1;         % Mutation Rate  %多少个报价发生变异

ANSWER=questdlg('Choose selection method:','Genetic Algorith',...
    'Roulette Wheel','Tournament','Random','Roulette Wheel');

UseRouletteWheelSelection=strcmp(ANSWER,'Roulette Wheel');
UseTournamentSelection=strcmp(ANSWER,'Tournament');
UseRandomSelection=strcmp(ANSWER,'Random');

if UseRouletteWheelSelection
    beta=8; % Selection Pressure
end

if UseTournamentSelection
    TournamentSize=3;   % Tournamnet Size
end

pause(0.01); % Due to a bug in older versions of MATLAB

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    % Initialize Position
    pop(i).Position(1:5)=sort(unifrnd(VarMin,VarMax,1,15)); %均匀分布
    pop(i).Position(6)=unifrnd(VarMin2,VarMax2,3);
    % Evaluation
%     pop(i).Cost=1;
%     pop(1).Position(1)=16.1114;
%     pop(1).Position(2)=54.8892;
    pop(i).Cost=CostFunction(pop(i).Position);
    
end

% Sort Population  从大到小排序
Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);%sort是个排序函数
pop=pop(SortOrder);

% Store Best Solution
BestSol=pop(1);
disp(BestSol);
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);
BestPrice=zeros(MaxIt,nVar);
% Store Cost
WorstCost=pop(end).Cost;

%% Main Loop

it=1;

%% Main Loop
if (exist("mat.mat","file")) ==2
   load('mat.mat')
   fprintf('find existing result,begin at %d iteration',it)
   pause(5)
end

while(it<=MaxIt)
    
    % Calculate Selection Probabilities
    if UseRouletteWheelSelection
        P=exp(-beta*Costs/WorstCost);
        P=P/sum(P);
    end
    
    % Crossover
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
        if UseRouletteWheelSelection
            i1=RouletteWheelSelection(P);
            i2=RouletteWheelSelection(P);
        end
        if UseTournamentSelection
            i1=TournamentSelection(pop,TournamentSize);
            i2=TournamentSelection(pop,TournamentSize);
        end
        if UseRandomSelection
            i1=randi([1 nPop]);   %从50个里面随机选一个数
            i2=randi([1 nPop]);
        end

        % Select Parents
        p1=pop(i1); %从pop里取出来
        p2=pop(i2);   
        
        % Apply Crossover
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position,gamma,VarMin,VarMax,VarMin2,VarMax2);
        
        % Evaluate Offsprings
%         popc(k,1).Cost=1;
%         popc(k,2).Cost=1;
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);%变成一维
    
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        popm(k).Position=Mutate(p.Position,mu,VarMin,VarMax,VarMin2,VarMax2);
        
        % Evaluate Mutant
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
    % Create Merged Population
    pop=[pop
         popc
         popm]; %#ok
     
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Truncation
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BestPrice(it,:)=BestSol.Position;
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    plot(-BestCost,'LineWidth',2);
    pause(1);
    save('mat.mat')
    it=it+1;
end

%% Results
figure;
semilogy(-BestCost,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
xlabel('遗传代数/代');
ylabel('收益/元');
grid on;

figure
plot(BestPrice(1,1:5),'LineWidth',2);
hold on
plot(BestPrice(1,6),'LineWidth',2);
xlabel('Iteration');
ylabel('Price');
grid on;
toc