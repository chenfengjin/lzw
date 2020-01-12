function Income_PSH=case39_YC(quoted_prices)
        disp(quoted_prices)
        quoted_price_energy=quoted_prices(1:5);
        quoted_price_frequency=quoted_prices(6);

%% 电力系统安全机组组合
%% 设置不同cases
nt=96;%调度时段数24个,爬坡设置为60。备用若96点 需要/4吗？
nf=24;
%% 输入电网相关数据
mpc=case39psh;%matpower通用函数，将数据从文件或者结构体中导入到数据矩阵中
[baseMVA, bus, gen, branch,gencost,load24] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch,mpc.gencost,mpc.load24);%赋值
[ref, pv, pq] = bustypes(bus,gen);%matpower通用函数，创建参考节点，PV 节点和 PQ 节点的节点向量

nb=size(bus,1);%电网节点个数
ng=size(gen,1);%发电机组个数
nbr=length(branch(:,1));%线路的个数

%提取负荷信息
% Pload0 = bus(:,3);%提取节点负荷,主要用各节点的负荷比

%分布系数矩阵
gbus=gen(:,1);%发电机组所在电网节点
gbus_THE=gen(1:2,1);
gbus_LNG=gen(3:4,1);
gbus_PSH=gen(5:6,1);
H = makePTDF(baseMVA, bus, branch, ref);%matpower通用函数，形成直流 PTDF 矩阵,即分布系数H 节点与线路的转移矩阵
Hg=H(:,gbus);%发电机组和线路的PTDF矩阵
Hg_PSH=H(:,gbus_PSH);
%提取发电机组参数：初始功率、初始状态、最大功率、最小功率、所在节点发电机台数、爬坡率、最小开机时间、最小关机时间、初始已开机时间、初始已关机时间
P0=gen(:,2);
P0_all=sum(P0,1);
I0=gen(:,8);
Pmax=gen(:,9);
Pmin=gen(:,10);
ramp=gen(:,11);%向上/下爬坡率设为一样的
Ton=gen(:,12);
Toff=gen(:,12);
Xon0=gen(:,13);
Xoff0=gen(:,14);

UT=max(0,min(nt,(Ton-Xon0).*I0));%开始调度后必须保持开机的时间长度
DT=max(0,min(nt,(Toff-Xoff0).*(1-I0)));%开始调度后必须保持关机的时间长度

%提取开机费用
su=gencost(:,1);

% 提取线路的传输功率
Plmax=branch(:,6);
Plmax(17)=90;%制造阻塞
Plmax(28)=170;
% Plmax=20000;%用于排查

%% 各种关系矩阵
Vbg=sparse(gbus,(1:ng)',1,nb,ng,ng);%发电机组-电网节点的关系矩阵

%% 输入时序负荷相关数据
Pdtt24=load24'; %导入24个时刻的负荷数据
t=1:nt;
t0=1:nf;
Pdtt96=interp1(t0,Pdtt24,1+(t-1)*0.25,'PCHIP');  %差值得到96个时刻的负荷

P_f=0.03*Pdtt24;

P_k=bus(:,3);
Pdtt=P_k./sum(P_k).*Pdtt96;
A=sum(Pdtt);
%% 优化模型
%设置变量
P=sdpvar(ng,nt);%机组有功出力
spinR_up=sdpvar(ng,nt);%上旋转备用
spinR_dn=sdpvar(ng,nt);%下旋转备用
SU=sdpvar(ng,nt);%开机费用

I=binvar(ng,nt);%ON/OFF status, 0-1
y=binvar(ng,nt);%startup indicator, 0-1
z=binvar(ng,nt);%shutdown indicator, 0-1

Pl=sdpvar(nbr,nt);%输电线路传输容量

%抽蓄变量
N_PSH=2;
I_gen_PSH=I(5:6,:);
P_gen_PSH=P(5:6,:);
I_pump_PSH=binvar(N_PSH,nt);%抽水状态
P_pump_PSH=sdpvar(N_PSH,nt);%抽水功率
% 定义功率-水量转换系数
%发电状态,
UpFactor=200;
DownFactor=160;
%定义抽水状态离散功率点
P_pump_Cos=200;

%抽水蓄能机组库容量约束,有三个水库，前两个水库有八个抽蓄机组，后一个有四个机组
N_ponds=1;
V=sdpvar(N_ponds,nt);
V_min=3.782*10^5;
V_max=1.949*10^6;
V_initial=9.782*10^5;

gencost(5,[4,6,8,10,12])=quoted_price_energy;
gencost(6,[4,6,8,10,12])=quoted_price_energy;
gencost(5,2)=quoted_price_frequency;
gencost(6,2)=quoted_price_frequency;
%报价系数差分处理
Len.Cost=13;
for k=(Len.Cost-1):-2:6
    gencost(:,k)=gencost(:,k)-gencost(:,k-2);
end

%调频市场相关数据

co_f=P_f./max(P_f); %系数
C_f=co_f.*gencost(:,2);%报价 待修改
v=gen(:,15);%调节速率
ki=gen(:,16);%综合调频性能指标
Pi=ki/max(ki);%归一化后综合调频性能指标
C_f_Pi=C_f./Pi;%调频里程排序价格

%两个市场的关联
Pfu=sdpvar(ng,nf);%机组调频容量，按nf时刻出清
Pfu2=sdpvar(ng,nt);%
I_f=binvar(ng,nf);%联合出清时需要如此设计！！bug


%能量关联系数，用于96时刻
I1=1:4:93;
I2=2:4:94;
I3=3:4:95;
I4=4:4:96;

Pfu2(:,I1)=Pfu;
Pfu2(:,I2)=Pfu;
Pfu2(:,I3)=Pfu;
Pfu2(:,I4)=Pfu;

%% 目标函数
% obj=beta1'*sum(P,2)+beta2'*sum(regR,2)+beta3'*sum(spinR,2)+sum(sum(SU))+sum(sum(SD));%目标函数,矩阵行向量求和
obj=0;


for t=1:1:nt
   obj=obj+sum(P(:,t).*gencost(:,4));
   for m=5:2:Len.Cost-2%每个申报区间的端点（对应MW）
   detaP=max(P(:,t)-gencost(:,m),0);%取出各个机组每一时刻位于某一区间的出力大小
   obj=obj+sum(detaP.*gencost(:,m+1));
   end
end

obj=obj/4;


obj=obj+sum(sum(SU));%目标函数,矩阵行向量求和
% 
% % c=gencost(:,4);
% % obj=obj+c'*sum(P,2);
% 
for t=1:1:nf
   obj=obj+sum(Pfu(:,t).*C_f_Pi(:,t));%调频里程是调频容量3倍
end
% 
% 
% M=1e18;%弛罚因子
% obj=obj+M*sum(sum(P_sc3+P_sc4+P_sc5+P_sc6));
%% 约束条件 能量市场
con=[];
for t=1:1:nt
    con=con+[(sum(P(:,t))-sum(P_pump_PSH(:,t))==sum(Pdtt(:,t)))];%功率平衡方程

    con=con+[sum(Pmax.*I(:,t))>=1.2*(sum(Pdtt(:,t)))];%南网系统正备用
    con=con+[sum(Pmin.*I(:,t))<=0.9*(sum(Pdtt(:,t)))];%南网系统负备用  
     
    con=con+[P(:,t)+spinR_up(:,t)+Pfu2(:,t)<=Pmax.*I(:,t)];%发电机组输出功率约束，正旋转备用
    con=con+[Pmin.*I(:,t)<=P(:,t)-spinR_dn(:,t)-Pfu2(:,t)];%发电机组输出功率约束，负旋转备用 
            
    con=con+[0<=spinR_up(:,t)<=10*ramp.*I(:,t)];%发电机组旋转备用爬坡率约束
    con=con+[0<=spinR_dn(:,t)<=10*ramp.*I(:,t)];%发电机组旋转备用爬坡率约束
    
    con=con+[y(:,t)+z(:,t)<=1];
    
    con=con+[P_pump_PSH(:,t)==P_pump_Cos*I_pump_PSH(:,t)];
    con=con+[max(I_gen_PSH(:,t))+max(I_pump_PSH(:,t))<=1];
    con=con+[V_min<=V(t)]+[V(t)<=V_max];%库容上下限约束
    con=con+[V(nt)==V_initial];  %日调节

    if t==1
        con=con+[I(:,t)-I0==y(:,t)-z(:,t)];
        con=con+[P(1:2,t)-P0(1:2)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%煤电机组输出功率上爬坡约束
        con=con+[P0(1:2)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%煤电机组输出功率下爬坡约束
        con=con+[P(3:6,t)-P0(3:6)<=15*ramp(3:6)];%抽蓄、汽电机组输出功率上爬坡约束（不受最小出力限制）
        con=con+[P0(3:6)-P(3:6,t)<=15*ramp(3:6)];%煤电机组输出功率下爬坡约束（不受最小出力限制）
        
       % 库容约束
        con=con+[V(t)==V_initial-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];%广蓄
      
    else
        con=con+[I(:,t)-I(:,t-1)==y(:,t)-z(:,t)];
        con=con+[P(1:2,t)-P(1:2,t-1)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%发电机组输出功率爬坡约束
        con=con+[P(1:2,t-1)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%发电机组输出功率爬坡约束
        con=con+[P(3:6,t)-P(3:6,t-1)<=15*ramp(3:6)];%抽蓄、汽电机组输出功率上爬坡约束（不受最小出力限制）
        con=con+[P(3:6,t-1)-P(3:6,t)<=15*ramp(3:6)];%抽蓄、汽电机组输出功率下爬坡约束（不受最小出力限制）
        con=con+[V(t)==V(t-1)-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];
    end

    con=con+[SU(:,t)==su.*y(:,t)];%开机费用约束
    
    con=con+[sum(spinR_up(:,t))>=0.05*sum(Pdtt(:,t))];%旋转备用总需求约束，5%总负荷
    con=con+[sum(spinR_dn(:,t))>=0.05*sum(Pdtt(:,t))];%旋转备用总需求约束，5%总负荷
    

    con=con+[Pl(:,t)==Hg*P(:,t)-H*Pdtt(:,t)-Hg_PSH*P_pump_PSH(:,t)];%线路潮流约束
    con=con+[(Pl(:,t)<=Plmax)];%线路潮流约束
    con=con+[(-Pl(:,t)<=Plmax)];%线路潮流约束

end
% 
% 调频市场时间尺度
for t=1:1:nf
    I_f(:,t)=sum(I(:,nt/nf*t-3:nt/nf*t),2);%两个市场的状态变量关联
    con=con+[sum(Pfu(:,t))==P_f(t)];%调频容量约束
    con=con+[Pfu(:,t)>=0];%调频容量下限
    con=con+[Pfu(:,t)<=3*v.*max(I_f(:,t)-3,0)];%调频容量上限1
    con=con+[Pfu(:,t)<=0.15*Pmax.*max(I_f(:,t)-3,0)];%调频容量上限2
end



% minimum on/off time constraints
% for i=1:1:2
%     if UT(i)>0   %UT为保持开机的最少时长
%         con=con+[UT(i)-sum(I(i,1:1:UT(i)))==0];
%      for j=UT(i)+1:1:nt-Ton(i)+1
%         con=con+[sum(I(i,j:1:j+Ton(i)-1))>=Ton(i)*y(i,j)];
%      end 
%     for j=nt-Ton(i)+2:1:nt
%         con=con+[sum(I(i,j:1:nt))-(nt-j+1)*y(i,j)>=0];
%     end
%     end
%     
%     if DT(i)>0   %DT为保持关机的最少时长      
%         con=con+[sum(I(i,1:1:DT(i)))==0];
%     for j=DT(i)+1:1:nt-Toff(i)+1
%         con=con+[Toff(i)-sum(I(i,j:1:j+Toff(i)-1))>=Toff(i)*z(i,j)];
%     end
%     for j=nt-Toff(i)+2:1:nt
%         con=con+[(nt-j+1)-sum(I(i,j:1:nt))-(nt-j+1)*z(i,j)>=0];
%     end
%     end
% end


%% 求解
%选择求解器
ses = sdpsettings('verbose','2','solver','gurobi','debug','1');
%设置求解值
ses.showprogress=1;
% ses.gurobi.MIPGAP=0.0002;
%计算
dd=solvesdp(con,obj,ses);

if dd.problem~=0
    error('模型不可行!');
end

%% 固定0/1变量
I=value(I);
I_f=value(I_f);
y=value(y);
z=value(z);
I_pump_PSH=value(I_pump_PSH);%抽水状态


%% 经济调度，计算节点电价
%% 目标函数
obj2=0;%目标函数

% c=gencost(:,4);用于快速测算
% obj2=obj2+c'*sum(P,2);

for t=1:1:nt
   obj2=obj2+sum(P(:,t).*gencost(:,4));
   for m=5:2:Len.Cost-2%每个申报区间的端点（对应MW）
   detaP2=max(P(:,t)-gencost(:,m),0);%取出各个机组每一时刻位于某一区间的出力大小
   obj2=obj2+sum(detaP2.*gencost(:,m+1));
   end
end

obj2=obj2/4;

for t=1:1:nf
   obj2=obj2+sum(Pfu(:,t).*C_f_Pi(:,t));%调频里程是调频容量3倍
end

%% 约束条件
con2=[];
for t=1:1:nt
    P_pump_PSH(:,t)=P_pump_Cos*I_pump_PSH(:,t);
    con2=con2+[(sum(P(:,t))-sum(P_pump_PSH(:,t))==sum(Pdtt(:,t))):['lambda' num2str(t)]];%功率平衡方程
    
    con2=con2+[P(:,t)+Pfu2(:,t)+spinR_up(:,t)<=Pmax.*I(:,t)];%发电机组输出功率约束，正旋转备用
    con2=con2+[Pmin.*I(:,t)<=P(:,t)-Pfu2(:,t)-spinR_dn(:,t)];%发电机组输出功率约束，负旋转备用
    
    con2=con2+[0<=spinR_up(:,t)<=10*ramp.*I(:,t)];%发电机组旋转备用爬坡率约束
    con2=con2+[0<=spinR_dn(:,t)<=10*ramp.*I(:,t)];%发电机组旋转备用爬坡率约束
    
    con2=con2+[V_min<=V(t)]+[V(t)<=V_max];%库容上下限约束
    con2=con2+[V(nt)==V_initial];  %日调节
    
    if t==1
        con2=con2+[P(1:2,t)-P0(1:2)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%煤电机组输出功率上爬坡约束
        con2=con2+[P0(1:2)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%煤电机组输出功率下爬坡约束
        con2=con2+[P(3:6,t)-P0(3:6)<=15*ramp(3:6)];%抽蓄、汽电机组输出功率上爬坡约束（不受最小出力限制）
        con2=con2+[P0(3:6)-P(3:6,t)<=15*ramp(3:6)];%煤电机组输出功率下爬坡约束（不受最小出力限制）
        
       % 库容约束
        con2=con2+[V(t)==V_initial-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];%广蓄
      
    else
        
        con2=con2+[P(1:2,t)-P(1:2,t-1)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%发电机组输出功率爬坡约束
        con2=con2+[P(1:2,t-1)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%发电机组输出功率爬坡约束
        con2=con2+[P(3:6,t)-P(3:6,t-1)<=15*ramp(3:6)];%抽蓄、汽电机组输出功率上爬坡约束（不受最小出力限制）
        con2=con2+[P(3:6,t-1)-P(3:6,t)<=15*ramp(3:6)];%煤电机组输出功率下爬坡约束（不受最小出力限制）
        con2=con2+[V(t)==V(t-1)-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];
    end

    con2=con2+[sum(spinR_up(:,t))>=0.05*sum(Pdtt(:,t))];%旋转备用总需求约束，8%总负荷
    con2=con2+[sum(spinR_dn(:,t))>=0.05*sum(Pdtt(:,t))];%旋转备用总需求约束，8%总负荷

    con2=con2+[Pl(:,t)==Hg*P(:,t)-H*Pdtt(:,t)-Hg_PSH*P_pump_PSH(:,t)];
    con2=con2+[(Pl(:,t)<=Plmax):['mup' num2str(t)]];%线路潮流约束
    con2=con2+[(-Pl(:,t)<=Plmax):['mun' num2str(t)]];%线路潮流约束
    

end

%调频市场时间尺度
for t=1:1:nf
    con2=con2+[(sum(Pfu(:,t))==P_f(t)):['lambdaf' num2str(t)]];%调频容量约束   
    con2=con2+[Pfu(:,t)>=0];%调频容量下限
    con2=con2+[Pfu(:,t)<=3*v.*max(I_f(:,t)-3,0)];%调频容量上限1
    con2=con2+[Pfu(:,t)<=0.15*Pmax.*max(I_f(:,t)-3,0)];%调频容量上限2
end


%选择求解器
ses2 = sdpsettings('verbose','2','solver','gurobi','debug','1');
%设置求解值
ses2.showprogress=1;
% ses2.gurobi.MIPGAP=0.0002;
%计算
dd2=solvesdp(con2,obj2,ses2);

if dd2.problem~=0
    error('模型不可行!');
end

%% 结果输出
lambda=zeros(1,nt);
lambdaf=zeros(1,nf);
mup=zeros(nbr,nt);
mun=zeros(nbr,nt);
LMP=zeros(nb,nt);
for t=1:1:nt
   lambda(t)=dual(con2(['lambda' num2str(t)])) ;
   mup(:,t)=dual(con2(['mup' num2str(t)])) ;
   mun(:,t)=dual(con2(['mun' num2str(t)])) ; 
   LMP(:,t)=-lambda(t)-H'*(mup(:,t)-mun(:,t));
end

for t=1:1:nf
   lambdaf(t)=-dual(con2(['lambdaf' num2str(t)])) ;
end

% P=double(P);%机组有功出力
Pfu=double(Pfu);
% Pfu2=double(Pfu2);
% detaP=double(detaP);
% Pl=double(Pl);
% spinR_up=double(spinR_up);%上旋转备用
% spinR_dn=double(spinR_dn);%下旋转备用
% SU=double(SU);%开机费用
% SD=double(SD);%关机费用
% I=double(I);%ON/OFF status, 0-1
% y=double(y);%startup indicator, 0-1
% z=double(z);%shutdown indicator, 0-1
P_pump_PSH=double(P_pump_PSH);
P_gen_PSH=double(P_gen_PSH);
% V=double(V);
% I_pump_PSH=double(I_pump_PSH);
% I_gen_PSH=double(I_gen_PSH);


% obj=double(obj);%初步总费用
% obj2=double(obj2);

% P_sc1=double(P_sc1);
% P_sc2=double(P_sc2);
% P_sc3=double(P_sc3);
% P_sc4=double(P_sc4);
%开机费用
% Fee_THE_SU=sum(sum(SU(1:2,:)));
% Fee_LNG_SU=sum(sum(SU(3:4,:)));


%能量市场收入初值
% Fee_THE=0;
% Fee_LNG=0;
Fee_PSH=0;
% 调频辅助服务市场收入初值
% Feef_THE=0;
% Feef_LNG=0;
Feef_PSH=0;

%机组收入
%辅助服务市场
for t=1:1:nf
% Feef_THE=Feef_THE+lambdaf(t)*sum(Pfu(1:2,t));
% Feef_LNG=Feef_LNG+lambdaf(t)*sum(Pfu(3:4,t));
Feef_PSH=Feef_PSH+lambdaf(t)*sum(Pfu(5:6,t));
end


%能量市场
for t=1:1:nt
% Fee_THE=Fee_THE+sum(LMP(gbus_THE,t).*P(1:2,t));
% Fee_LNG=Fee_LNG+sum(LMP(gbus_LNG,t).*P(3:4,t));
Fee_PSH=Fee_PSH+sum(LMP(gbus_PSH,t).*(P_gen_PSH(:,t)-P_pump_PSH(:,t)));
end

% LMP24=4*LMP;


% W_e=Fee_THE+Fee_LNG+Fee_PSH;
% W_f=Feef_THE+Feef_LNG+Feef_PSH;
% W_all=W_e+W_f;
% SU_all=sum(sum(SU));
% S=W_all+SU_all;


% %add order
% LMP_busPSH=sum(LMP24(gbus_PSH,:))/2;
% LMP_ave=sum(LMP24)/39;
% LMP_ave_24=zeros(1,24);
% for i=1:1:24
% LMP_ave_24(i)=sum(LMP_ave(4*i-3:4*i))/4;
% end


% P_THE=sum(P(1:2,:));
% P_LNG=sum(P(3:4,:));
% P_gen_PSH_all=sum(P_gen_PSH);
% P_pump_PSH_all=sum(P_pump_PSH);
% Pfu_THE=sum(Pfu(1:2,:));
% Pfu_LNG=sum(Pfu(3:4,:));
% Pfu_PSH=sum(Pfu(5:6,:));

Income_PSH=Fee_PSH+Feef_PSH;


%% 返回结果 计算抽水蓄能收益


Income_PSH=double(Income_PSH);

yalmip('clear');
end
