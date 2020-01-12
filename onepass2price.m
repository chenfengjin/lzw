%% 说明：本程序为电能量市场和调频市场依次出清，计算分析抽水蓄能效益??
%case39 33??7为蓄能，31??8为汽机，30??5为火??32??6为风??
function Income_PSH=onepass2price(quoted_prices)
try
        disp(quoted_prices)
        quoted_price_energy=quoted_prices(1:5);
        quoted_price_frequency=quoted_prices(6);
    %% 输入原始数据
    % 选用39节点标准算例
    mpc=case39lzw_tiao;
    % 导入、确定参考节??
    [~, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);%赋??
    %% 采用传??进来的报价数??
    % gencost = [
    %        40000  40000  0   100   250  200    312.5  300     375   400     437.5   500   500;
    %        20000  20000  0   170   100  300    200    400     300   550     400     650   500;
    %        0      0      0    quoted_price_energy   100    quoted_price_energy     125      quoted_price_energy     150     quoted_price_energy    175       quoted_price_energy  200;
    %        40000  40000  0   150   200  250    240    350     280   500     320     600   400;
    %        0      0      0    quoted_price_energy   100    quoted_price_energy     125      quoted_price_energy    150     quoted_price_energy    175       quoted_price_energy  200;
    %        20000  20000  0   190   80   320    160    450     240   600     320     700   400;
    %     ];
    gencost = [
        40000  40000  0   100   250  110    312.5  300     375   500     437.5   600   500;
        20000  20000  0   100   100  200    200    550     300   750     400     800   500;
        0      0      0    quoted_price_energy(1)   100    quoted_price_energy(2)     125      quoted_price_energy(3)     150     quoted_price_energy(4)    175       quoted_price_energy(5)  200;
        40000  40000  0   100   200  120    240    300     280   320     320     400   400;
        0      0      0    quoted_price_energy(1)   100    quoted_price_energy(2)     125      quoted_price_energy(3)    150     quoted_price_energy(4)    175       quoted_price_energy(5)  200;
        20000  20000  0   100   80   200    160    400     240   550     320     600   400;
        ];
    % [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);%matpower通用函数，创建参考节点，PV 节点??PQ 节点的节点向??
    % 导入负荷数据
    Pforecast0=0.2*xlsread('qifuwo',1,'A1:A24'); %导入24个时刻的负荷数据

    %% 系统常规设置
    Horizon = 96;%电能量市场时间尺??
    T= 24;%调频市场时间尺度
    t=1:Horizon;
    t0=1:T;
    Pforecast=interp1(t0,Pforecast0,1+(t-1)*0.25,'PCHIP');  %差??得到96个时刻的负荷
    % Pforecast=Pforecast0';
    % Max_Pforecast=max(Pforecast);%负荷????值，用于调参
    % Min_Pforecast=min(Pforecast);%负荷????值，用于调参

    %系统备用
    Rup   =0.01*Pforecast;%负荷??/20
    Rdown =0.01*Pforecast;
    % SRup  =0.005*Pforecast;%时段t上调旋转备用要求  
    % SRdown=0.005*Pforecast;%时段t下调旋转备用要求

    %节点??
    N_bus=size(bus,1);
    %线路??
    N_branch=size(branch,1);

    %% 发电??
    %1.机组节点编号
    %????机组节点编号
    BusN_G=gen(:,1);%列向量，发电机组????节点编号
    %火电节点编号数组
    NN=cell(3,1);
    NN(1)={[30,35]}; %现在繁琐，后续有大量????输入时体现出????
    %汽机节点编号数组
    NN(2)={[31,38]};
    %抽蓄节点编号数组
    NN(3)={[33,37]};
    %火电机组节点编号
    genis_THE=ismember(BusN_G,NN{1}); 
    BusN_THE=BusN_G(genis_THE);
    %汽机节点编号
    genis_LNG=ismember(BusN_G,NN{2});
    BusN_LNG=BusN_G(genis_LNG);
    %抽蓄机组节点编号
    genis_PSH=ismember(BusN_G,NN{3});
    BusN_PSH=BusN_G(genis_PSH);
    % %风电机组节点编号
    % BusN_wind=[32,36];



    %机组数量
    N_ALL=size(BusN_G,1);
    N_THE=size(BusN_THE,1);
    N_LNG=size(BusN_LNG,1);
    N_PSH=size(BusN_PSH,1);

    %机组出力上下??
    Pmax_ALL=gen(:,9);%????机组的出力上限：额定有功功率
    Pmin_ALL=gen(:,10);%下限
    Pmaxt_ALL=repmat(Pmax_ALL,1,Horizon);
    Pmint_ALL=repmat(Pmin_ALL,1,Horizon);
    % Sum_Pmax_ALL=sum(Pmax_ALL);%????出力之和，用于调??
    % Sum_Pmin_ALL=sum(Pmin_ALL);%????出力之和，用于调??
    % Pmax_THE=max(Pmax_ALL(genis_THE));
    % Pmax_LNG=max(Pmax_ALL(genis_LNG));
    % Pmax_PSH=max(Pmax_ALL(genis_PSH));

    %机组爬坡??
    Ramp_U_ALL=gen(:,17);
    Ramp_D_ALL=gen(:,18);


    %????????/停机时间
    Ton_min_ALL=gen(:,20);
    % Toff_min_ALL=gen(:,21);
    Ton_min_THE=Ton_min_ALL(genis_THE);
    Toff_min_THE=Ton_min_ALL(genis_THE);
    %抽蓄、气电无此约??

    %抽水蓄能
    %离散抽水功率
    % P_pump1_PSH=200;
    P_pump1_PSH=100;
    %定义抽水蓄能机组库容??
    V_min=3.782*10^5;
    V_max=1.949*10^6;
    V_initial=9.782*10^5;
    V_mint=repmat(V_min,N_PSH,Horizon);
    V_maxt=repmat(V_max,N_PSH,Horizon);
    %定义抽水状??下水电转化系?? 
    Factor_pump=150;
    %定义发电状??下水电转化系??
    Factor_gen=200;
    % %风电
    % P_wind=xlsread('gen_wind.xlsx');
    % A=Pforecast-sum(P_wind);
    % A_MAX=max(A);
    % A_MIN=min(A);
    % figure
    % P2=plot(value(sum(P_wind))');
    % figure
    % P=plot(value(A)');

    %报价(五段报价，与成本有关)
    %后面看看能不能编到前面来
    Len.Cost=size(gencost,2);

    %状??变量 状??和出??
    U_ALL=binvar(N_ALL,Horizon,'full');
    P_ALL=sdpvar(N_ALL,Horizon,'full');

    U_THE=U_ALL(genis_THE,:);
    %抽水蓄能机组
    U_pump_PSH=binvar(N_PSH,Horizon,'full');%抽水状??
    P_pump_PSH=sdpvar(N_PSH,Horizon,'full');%抽水功率
    ADD=binvar(N_PSH,Horizon,'full');%抽水功率附加变量
    V=sdpvar(N_PSH,Horizon,'full');%库容

    %机组启停
    % yg(:,1)=U_ALL(:,1);%若t=1时刻机组????，则该时刻启动次数为1
    % yg(:,2:Horizon)=diff(U_ALL,1,2);
    % L_on=gencost(:,1);%每个机组每一时刻启动费用
    % % L_off=gencost(:,2);%每个机组每一时刻停机费用
    % yg_obj=abs(yg);

    %% 两个市场的关??
    Nfu=N_ALL;%调频服务提供者数??
    Pfu=sdpvar(Nfu,T,'full');
    Pfu2=sdpvar(Nfu,Horizon,'full');

    I1=1:4:93;
    I2=2:4:94;
    I3=3:4:95;
    I4=4:4:96;
    Pfu2(:,I1)=Pfu;
    Pfu2(:,I2)=Pfu;
    Pfu2(:,I3)=Pfu;
    Pfu2(:,I4)=Pfu;

    %% 联合出清·电能量市??

    %约束
    Constraints=[];

    %1.系统
    % %(1)系统负荷平衡约束
    Constraints=[Constraints,sum(P_ALL)-sum(P_pump_PSH)==Pforecast]; 
    % %(1)系统负荷平衡约束
    % Constraints=[Constraints,sum(P_ALL)-sum(P_pump_PSH)+sum(P_sc1)-sum(P_sc2)==Pforecast]; 

    % Constraints=[Constraints,P_sc1>=0];
    % Constraints=[Constraints,P_sc2>=0];
    %(2)系统正备用容量约??
    Constraints=[Constraints,sum(bsxfun(@times,U_ALL,Pmax_ALL))>=Pforecast+Rup];
    %(3)系统负备用容量约??
    Constraints=[Constraints,sum(bsxfun(@times,U_ALL,Pmin_ALL))-sum(bsxfun(@times,U_pump_PSH,100))<=Pforecast-Rdown];
    % %(4)系统旋转备用要求 
    % Constraints=[Constraints,sum(bsxfun(@min,bsxfun(@plus,-P_ALL,Pmax_ALL),Ramp_U_ALL))>=SRup];
    % Constraints=[Constraints,sum(bsxfun(@min,max(bsxfun(@minus,P_ALL,Pmin_ALL),0),Ramp_D_ALL))>=SRdown];
    % 
    %2.机组
    %(1)机组出力上下限约??
    Constraints=[Constraints,(bsxfun(@times,U_ALL,Pmin_ALL)<=P_ALL-Pfu2)&(P_ALL+Pfu2<=bsxfun(@times,U_ALL,Pmax_ALL))];
    % Constraints=[Constraints,P_sc3>=0];
    % Constraints=[Constraints,P_sc4>=0];

    Ramp_add=[0;0;Pmax_ALL(3);0;Pmax_ALL(5);0];
    Rampt_add=repmat(Ramp_add,1,Horizon-1);
    %(2)机组爬坡约束
    Constraints=[Constraints,diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,1:Horizon-1),Ramp_U_ALL)+bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,2:Horizon)),Pmax_ALL)];
    Constraints=[Constraints,-diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,2:Horizon),Ramp_D_ALL)-bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,1:Horizon-1)),Pmax_ALL)];

    %(3)火电机组特有约束
    %  1)机组????连续????时间约束
    for t=2:Horizon 
    for i=1:N_THE    
        k=max(t-Ton_min_THE(i),1):(t-1);
        LUP(i,t)=sum(U_THE(i,k));%时段t时已连续????时间（到t-1??
        k=max(t-Toff_min_THE(i),1):(t-1);
        LDOWN(i,t)=sum(1-U_THE(i,k));%时段t时已连续停机时间（到t-1??
        Constraints=[Constraints,LDOWN(i,t)>=(U_THE(i,t)-U_THE(i,t-1))*Toff_min_THE(i)];%????连续停机时间约束
        Constraints=[Constraints,LUP(i,t)>=(U_THE(i,t-1)-U_THE(i,t))*Ton_min_THE(i)];   %????连续????时间约束
    end
    end

    %??）抽水蓄能特有约??
    %   1）抽水功率约??
        Constraints=[Constraints,P_pump_PSH==P_pump1_PSH*U_pump_PSH-100*ADD+100];

    %   2）状态变量约??
        Constraints=[Constraints,max(U_ALL(genis_PSH,:))+max(U_pump_PSH)<=1];
        Constraints=[Constraints,U_pump_PSH+ADD>=1];      
    %   3）库容约??
        Constraints=[Constraints,V(:,1)==V_initial+Factor_pump*P_pump_PSH(:,1)-Factor_gen*P_ALL(genis_PSH,1)];
        for t=2:Horizon
        Constraints=[Constraints,V(:,t)==V(:,t-1)+Factor_pump*P_pump_PSH(:,t)-Factor_gen*P_ALL(genis_PSH,t)];
        end
        Constraints=[Constraints,(V_mint<=V)&(V<=V_maxt)];
        Constraints=[Constraints,V(:,Horizon)==V_initial];  %日调??
    % 
    %(5)线路潮流约束
    H=makePTDF(mpc);%节点对线路的功率转移分布因子,行数为支路数，列数为节点数??
    Hg=H(:,BusN_G);
    SLlf=sdpvar(N_branch,Horizon,'full');%行为区域内支路????
    SLlz=sdpvar(N_branch,Horizon,'full');
    Plmax=branch(:,6);%列向量，rataA长期允许功率
    Plmax(17)=90;%制??阻塞
    Plmax(28)=170;

    Plmaxt=repmat(Plmax,1,Horizon);
    Pd=bus(:,3);%列向量，母线负荷??
    Pnforecast=repmat(Pforecast,N_bus,1);
    Pdt=bsxfun(@times,Pnforecast./sum(Pd),Pd);%每一个节点每????时刻的负荷??

    Constraints=[Constraints,-Plmaxt<=Hg*P_ALL-H(:,BusN_PSH)*P_pump_PSH-H*Pdt-SLlz+SLlf];
    Constraints=[Constraints,Hg*P_ALL-H(:,BusN_PSH)*P_pump_PSH-H*Pdt-SLlz+SLlf<=Plmaxt];

    Constraints=[Constraints,SLlz>=0];
    Constraints=[Constraints,SLlf>=0];

    %% 调频辅助服务市场

    %调频????
    Cneed=0.04*max(Pforecast0-min(Pforecast0),max(Pforecast0)-Pforecast0);%次日24小时各时段调频容量需求（MW??max函数作用后，是按列排列的
    co_F=Cneed'./max(Cneed);

    %调频费用（报价）
    Cfucost=[24;36;quoted_price_frequency;22;quoted_price_frequency;34]*co_F;
    %，对次日24个时段进行调频里程报价（??MW??????单位0.1MW??

    %调节速率
    v=0.5*[3/100*Pmax_ALL(1);4.5/100*Pmax_ALL(2);Pmax_ALL(3);3/100*Pmax_ALL(4);Pmax_ALL(5);4.5/100*Pmax_ALL(6)];%发电单元标准调节速率
    ki=[0.6;1;1.8;0.6;1.8;1];%各台发电机的综合调频性能指标，不小于0.5
    Pi=ki/max(ki);%归一化??能指??
    Cfucost_Pi=bsxfun(@rdivide,Cfucost,Pi);%调频里程排序价格=调频里程报价/ki


    %调频容量约束
    Constraints=[Constraints,sum(Pfu)==Cneed'];%


    %状??变量关联
    for i=1:T %SCUC??6时段，调频辅助服务市场为24时段，将4时段启停状??相加形成fU_ALL
    fU_ALL(:,i)=sum((U_ALL(:,Horizon/T*i-3:Horizon/T*i)),2);
    end

    %调频容量申报上下??

    Constraints=[Constraints,0<=Pfu];
    Constraints=[Constraints,Pfu<=0.5* min(repmat(v,1,T)*3,0.15*repmat(Pmax_ALL,1,T)).*max(fU_ALL-3,0)];


    %% 目标函数
    Objective=0;

    % % 启停成本                                                     
    % C_start=sum(sum(bsxfun(@times,L_on,yg_obj)));
    % Objective=Objective+C_start;%总成本初始??为火电启动成本，抽蓄无启动成??

    % % 报价成本
    % for m=(Len.Cost-2):-2:3%每个申报区间的端点（对应MW??
    %     c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0)- max(bsxfun(@minus,P_ALL,gencost(:,m+2)),0);%取出各个机组每一时刻位于某一区间的出力大??
    %     Objective=Objective+sum(sum(bsxfun(@times,c,gencost(:,m+1)),2));
    % end

    for j=(Len.Cost-1):-2:6
        gencost(:,j)=gencost(:,j)-gencost(:,j-2);
    end

    for m=3:2:Len.Cost-2%每个申报区间的端点（对应MW??
        c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0);%取出各个机组每一时刻位于某一区间的出力大??
        Objective=Objective+sum(sum( bsxfun(@times,c,gencost(:,m+1))));
    end

    % Objective=Objective+sum(sum(bsxfun(@times,[400;600;0;500;0;700],P_ALL)));


    M=1e9;%网络约束松弛罚因??
    Objective=Objective+M*sum(sum(SLlz+SLlf));
    
    % C_f=sum(sum(Cfucost_Pi.*Pfu+0.1*bsxfun(@times,Pfu,1-Pi)+0.01*bsxfun(@times,Pfu,max(ki)-ki) ));%目标函数调频费用????,当发电单元排序价格相同时，优先出??P 值高的发电单元；当发电单??P 值相同时，优先出??k 值高的发电单??
    C_f=sum(sum(Cfucost_Pi.*Pfu));%目标函数调频费用????,当发电单元排序价格相同时，优先出??P 值高的发电单元；当发电单??P 值相同时，优先出??k 值高的发电单??

    Objective=Objective+C_f;

    % % Objective=Objective+M*(sum(sum(P_sc1+P_sc2+P_sc3+P_sc4)))+M*(sum(sum(P_sc5+P_sc6+P_sc7+P_sc8)));
    % Objective=Objective+M*(sum(sum(P_sc1+P_sc2+P_sc3+P_sc4)))+M*(sum(sum(P_sc7+P_sc8)));

    %求解
    ops=sdpsettings('solver','gurobi','verbos',1);
    % ops.showprogress=1;
    ops.gurobi.MIPGAP=0.0002;
    optimize(Constraints,Objective,ops);



    %固定0/1变量
    U_ALL=double(U_ALL);
    U_pump_PSH=double(U_pump_PSH);
    fU_ALL=double(fU_ALL);
    ADD=double(ADD);


    %% 安全约束经济调度
    %约束条件
    Constraints2=[];

    %% 调频市场
    %(1)功率平衡约束
    Constraints2=[Constraints2,sum(P_ALL)-sum(P_pump_PSH)==Pforecast]; %为求节点电价方便放在这里

    %调频容量约束
    Constraints2=[Constraints2,sum(Pfu)==Cneed'];%

    %调频容量申报上下??
    % Constraints2=[Constraints2,0.5* repmat(min(v*1,3/100*(Pmax_ALL-Pmin_ALL)),1,T).*max(fU_ALL-3,0)<=Pfu];
    Constraints2=[Constraints2,0<=Pfu];
    Constraints2=[Constraints2,Pfu<=0.5* min(repmat(v,1,T)*3,0.15*repmat(Pmax_ALL,1,T)).*max(fU_ALL-3,0)];


    %(2)系统旋转备用要求 
    % Constraints2=[Constraints2,sum(bsxfun(@min,bsxfun(@plus,-P_ALL,Pmax_ALL),Ramp_U_ALL))>=SRup];
    % Constraints2=[Constraints2,sum(bsxfun(@min,bsxfun(@max,(bsxfun(@minus,P_ALL,Pmin_ALL)),0),Ramp_D_ALL))>=SRdown];

    %(3)机组出力上下限约??
    % Constraints2=[Constraints2,(bsxfun(@times,U_ALL,Pmint_ALL)<=P_ALL+P_sc)&(P_ALL+P_sc<=bsxfun(@times,U_ALL,Pmaxt_ALL))];
    Constraints2=[Constraints2,(bsxfun(@times,U_ALL,Pmint_ALL)<=P_ALL-Pfu2)&(P_ALL+Pfu2<=bsxfun(@times,U_ALL,Pmaxt_ALL))];

    %(4)机组爬坡约束
    Constraints2=[Constraints2,diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,1:Horizon-1),Ramp_U_ALL)+bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,2:Horizon)),Pmax_ALL)];
    Constraints2=[Constraints2,-diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,2:Horizon),Ramp_D_ALL)-bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,1:Horizon-1)),Pmax_ALL)];    
        
    %(5)抽水蓄能特有约束
    %   1）抽水功率约??
    Constraints2=[Constraints2,P_pump_PSH==P_pump1_PSH*U_pump_PSH-100*ADD+100]; 
    %   2）库容约??
    Constraints2=[Constraints2,V(:,1)==V_initial+Factor_pump*P_pump_PSH(:,1)-Factor_gen*P_ALL(genis_PSH,1)];
        for t=2:Horizon
    Constraints2=[Constraints2,V(:,t)==V(:,t-1)+Factor_pump*P_pump_PSH(:,t)-Factor_gen*P_ALL(genis_PSH,t)];
        end
    Constraints2=[Constraints2,(V_mint<=V)&(V<=V_maxt)];
    Constraints2=[Constraints2,V(:,Horizon)==V_initial];  %日调??
    %     
    % AAA=zeros(N_branch,Horizon);
    % AA=zeros(N_branch,1);
    % (6)线路潮流约束
    for i=1:N_branch
    Constraints2=[Constraints2,-Plmaxt(i,:)<=Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:)];
    Constraints2=[Constraints2,Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:)<=Plmaxt(i,:)];
    % AAA(i,:)=value(Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:));
    % AA(i)=min(Plmax(i,:)-abs(AAA(i,:)));
    % A=min(AA);
    end

    Constraints2=[Constraints2,SLlz>=0];
    Constraints2=[Constraints2,SLlf>=0];




    %目标函数
    Objective2=0;%总成??

    for m=3:2:Len.Cost-2%每个申报区间的端点（对应MW??
        c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0);%取出各个机组每一时刻位于某一区间的出力大??
        Objective2=Objective2+sum(sum( bsxfun(@times,c,gencost(:,m+1))));
    end

    % Objective2=Objective2+sum(sum(bsxfun(@times,[400;600;0;500;0;700],P_ALL)));

    % Price_pump=10;
    % C_p_PSH_pump=sum(sum(Price_pump*P_pump_PSH));
    % 
    % Objective2=Objective2-C_p_PSH_pump;

    %网络约束
    Objective2=Objective2+M*sum(sum(SLlz+SLlf));
    

    C_f=sum(sum(Cfucost_Pi.*Pfu));%目标函数调频费用????,当发电单元排序价格相同时，优先出??P 值高的发电单元；当发电单??P 值相同时，优先出??k 值高的发电单??

    Objective2=Objective2+C_f;

    %求解
    ops2=sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1);
    % ops2.showprogress=1;
    ops2.gurobi.MIPGAP=0.0002;
    optimize(Constraints2,Objective2,ops2);


    %% 求解节点电价
    mup=zeros(N_branch,Horizon);%线路????正向潮流约束拉格朗日乘子
    mun=zeros(N_branch,Horizon);
    for j=1:N_branch
        findlz=find(value(SLlz(j,:))~=0);%找出线路正向潮流松弛变量中不为零的角标，即为正向潮流越限??
        findlf=find(value(SLlf(j,:))~=0);%找出线路反向潮流松弛变量中不为零的角标，即为反向潮流越限??
        mun(j,:)=dual(Constraints2(end-2*N_branch-3+2*j))';
        mup(j,:)=dual(Constraints2(end-2*N_branch-2+2*j))';
        mup(j,findlz)=M;%线路正向潮流越限时，该拉格朗日乘子为网络潮流约束松弛罚因??
        mun(j,findlf)=M;%线路反向潮流越限时，该拉格朗日乘子为网络潮流约束松弛罚因??
    end
    Price_F=-dual(Constraints2(2));
    lameda=dual(Constraints2(1));
    LMP0=zeros(N_bus,Horizon);%行数为区域内节点数，列数为时段数
    for k=1:N_bus
        LMP0(k,:)=-lameda-sum(bsxfun(@times,mup,H(:,k))-bsxfun(@times,mun,H(:,k)));
    end
    LMP=zeros(N_bus,T);
    W_ALL=zeros(N_ALL,T);%机组每小时发电机发电??
    W_THE=zeros(N_THE,T);%
    W_LNG=zeros(N_LNG,T);%
    W_PSH=zeros(N_PSH,T);%
    W_PSH_gen=zeros(N_PSH,T);
    W_PSH_pump=zeros(N_PSH,T);
    P_THE=P_ALL(genis_THE,:);
    P_LNG=P_ALL(genis_LNG,:);
    P_PSH=P_ALL(genis_PSH,:);

    for Period=1:T%每小时的平均节点电价为每小时????5分钟节点电价的算术平均??
        LMP(:,Period)=sum(LMP0(:,4*Period-3:4*Period),2)/4;
        W_ALL(:,Period)=sum(value(P_ALL(:,4*Period-3:4*Period)),2)/4;
        W_THE(:,Period)=sum(value(P_THE(:,4*Period-3:4*Period)),2)/4;
        W_LNG(:,Period)=sum(value(P_LNG(:,4*Period-3:4*Period)),2)/4;
        W_PSH_gen(:,Period)=sum(value(P_PSH(:,4*Period-3:4*Period)),2)/4;
        W_PSH_pump(:,Period)=sum(value(P_pump_PSH(:,4*Period-3:4*Period)),2)/4;
        W_PSH=W_PSH_gen-W_PSH_pump;
    end

    % Income_THE_e=sum(sum(LMP([30,35],:).*W_THE));
    % Income_LNG_e=sum(sum(LMP([31,38],:).*W_LNG));
    Income_PSH_e=sum(sum(LMP([33,37],:).*W_PSH));
    % Income_THE_r=sum(sum(repmat(Price_F,N_THE,1).*Pfu(genis_THE,:)));
    % Income_LNG_r=sum(sum(repmat(Price_F,N_LNG,1).*Pfu(genis_LNG,:)));
    Income_PSH_r=sum(sum(repmat(Price_F,N_PSH,1).*Pfu(genis_PSH,:)));
    % Income_THE=Income_THE_e+Income_THE_r;
    % Income_LNG=Income_LNG_e+Income_LNG_r;
    Income_PSH=Income_PSH_e+Income_PSH_r;




    %% 返回结果，计算抽水蓄能收??

    Income_PSH=double(Income_PSH);
catch
    Income_PSH=0
    disp('onepass error,return 0 as default')
end
yalmip('clear');
end
