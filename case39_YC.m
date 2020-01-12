function Income_PSH=case39_YC(quoted_prices)
        disp(quoted_prices)
        quoted_price_energy=quoted_prices(1:5);
        quoted_price_frequency=quoted_prices(6);

%% ����ϵͳ��ȫ�������
%% ���ò�ͬcases
nt=96;%����ʱ����24��,��������Ϊ60��������96�� ��Ҫ/4��
nf=24;
%% ��������������
mpc=case39psh;%matpowerͨ�ú����������ݴ��ļ����߽ṹ���е��뵽���ݾ�����
[baseMVA, bus, gen, branch,gencost,load24] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch,mpc.gencost,mpc.load24);%��ֵ
[ref, pv, pq] = bustypes(bus,gen);%matpowerͨ�ú����������ο��ڵ㣬PV �ڵ�� PQ �ڵ�Ľڵ�����

nb=size(bus,1);%�����ڵ����
ng=size(gen,1);%����������
nbr=length(branch(:,1));%��·�ĸ���

%��ȡ������Ϣ
% Pload0 = bus(:,3);%��ȡ�ڵ㸺��,��Ҫ�ø��ڵ�ĸ��ɱ�

%�ֲ�ϵ������
gbus=gen(:,1);%����������ڵ����ڵ�
gbus_THE=gen(1:2,1);
gbus_LNG=gen(3:4,1);
gbus_PSH=gen(5:6,1);
H = makePTDF(baseMVA, bus, branch, ref);%matpowerͨ�ú������γ�ֱ�� PTDF ����,���ֲ�ϵ��H �ڵ�����·��ת�ƾ���
Hg=H(:,gbus);%����������·��PTDF����
Hg_PSH=H(:,gbus_PSH);
%��ȡ��������������ʼ���ʡ���ʼ״̬������ʡ���С���ʡ����ڽڵ㷢���̨���������ʡ���С����ʱ�䡢��С�ػ�ʱ�䡢��ʼ�ѿ���ʱ�䡢��ʼ�ѹػ�ʱ��
P0=gen(:,2);
P0_all=sum(P0,1);
I0=gen(:,8);
Pmax=gen(:,9);
Pmin=gen(:,10);
ramp=gen(:,11);%����/����������Ϊһ����
Ton=gen(:,12);
Toff=gen(:,12);
Xon0=gen(:,13);
Xoff0=gen(:,14);

UT=max(0,min(nt,(Ton-Xon0).*I0));%��ʼ���Ⱥ���뱣�ֿ�����ʱ�䳤��
DT=max(0,min(nt,(Toff-Xoff0).*(1-I0)));%��ʼ���Ⱥ���뱣�ֹػ���ʱ�䳤��

%��ȡ��������
su=gencost(:,1);

% ��ȡ��·�Ĵ��书��
Plmax=branch(:,6);
Plmax(17)=90;%��������
Plmax(28)=170;
% Plmax=20000;%�����Ų�

%% ���ֹ�ϵ����
Vbg=sparse(gbus,(1:ng)',1,nb,ng,ng);%�������-�����ڵ�Ĺ�ϵ����

%% ����ʱ�򸺺��������
Pdtt24=load24'; %����24��ʱ�̵ĸ�������
t=1:nt;
t0=1:nf;
Pdtt96=interp1(t0,Pdtt24,1+(t-1)*0.25,'PCHIP');  %��ֵ�õ�96��ʱ�̵ĸ���

P_f=0.03*Pdtt24;

P_k=bus(:,3);
Pdtt=P_k./sum(P_k).*Pdtt96;
A=sum(Pdtt);
%% �Ż�ģ��
%���ñ���
P=sdpvar(ng,nt);%�����й�����
spinR_up=sdpvar(ng,nt);%����ת����
spinR_dn=sdpvar(ng,nt);%����ת����
SU=sdpvar(ng,nt);%��������

I=binvar(ng,nt);%ON/OFF status, 0-1
y=binvar(ng,nt);%startup indicator, 0-1
z=binvar(ng,nt);%shutdown indicator, 0-1

Pl=sdpvar(nbr,nt);%�����·��������

%�������
N_PSH=2;
I_gen_PSH=I(5:6,:);
P_gen_PSH=P(5:6,:);
I_pump_PSH=binvar(N_PSH,nt);%��ˮ״̬
P_pump_PSH=sdpvar(N_PSH,nt);%��ˮ����
% ���幦��-ˮ��ת��ϵ��
%����״̬,
UpFactor=200;
DownFactor=160;
%�����ˮ״̬��ɢ���ʵ�
P_pump_Cos=200;

%��ˮ���ܻ��������Լ��,������ˮ�⣬ǰ����ˮ���а˸�������飬��һ�����ĸ�����
N_ponds=1;
V=sdpvar(N_ponds,nt);
V_min=3.782*10^5;
V_max=1.949*10^6;
V_initial=9.782*10^5;

gencost(5,[4,6,8,10,12])=quoted_price_energy;
gencost(6,[4,6,8,10,12])=quoted_price_energy;
gencost(5,2)=quoted_price_frequency;
gencost(6,2)=quoted_price_frequency;
%����ϵ����ִ���
Len.Cost=13;
for k=(Len.Cost-1):-2:6
    gencost(:,k)=gencost(:,k)-gencost(:,k-2);
end

%��Ƶ�г��������

co_f=P_f./max(P_f); %ϵ��
C_f=co_f.*gencost(:,2);%���� ���޸�
v=gen(:,15);%��������
ki=gen(:,16);%�ۺϵ�Ƶ����ָ��
Pi=ki/max(ki);%��һ�����ۺϵ�Ƶ����ָ��
C_f_Pi=C_f./Pi;%��Ƶ�������۸�

%�����г��Ĺ���
Pfu=sdpvar(ng,nf);%�����Ƶ��������nfʱ�̳���
Pfu2=sdpvar(ng,nt);%
I_f=binvar(ng,nf);%���ϳ���ʱ��Ҫ�����ƣ���bug


%��������ϵ��������96ʱ��
I1=1:4:93;
I2=2:4:94;
I3=3:4:95;
I4=4:4:96;

Pfu2(:,I1)=Pfu;
Pfu2(:,I2)=Pfu;
Pfu2(:,I3)=Pfu;
Pfu2(:,I4)=Pfu;

%% Ŀ�꺯��
% obj=beta1'*sum(P,2)+beta2'*sum(regR,2)+beta3'*sum(spinR,2)+sum(sum(SU))+sum(sum(SD));%Ŀ�꺯��,�������������
obj=0;


for t=1:1:nt
   obj=obj+sum(P(:,t).*gencost(:,4));
   for m=5:2:Len.Cost-2%ÿ���걨����Ķ˵㣨��ӦMW��
   detaP=max(P(:,t)-gencost(:,m),0);%ȡ����������ÿһʱ��λ��ĳһ����ĳ�����С
   obj=obj+sum(detaP.*gencost(:,m+1));
   end
end

obj=obj/4;


obj=obj+sum(sum(SU));%Ŀ�꺯��,�������������
% 
% % c=gencost(:,4);
% % obj=obj+c'*sum(P,2);
% 
for t=1:1:nf
   obj=obj+sum(Pfu(:,t).*C_f_Pi(:,t));%��Ƶ����ǵ�Ƶ����3��
end
% 
% 
% M=1e18;%�ڷ�����
% obj=obj+M*sum(sum(P_sc3+P_sc4+P_sc5+P_sc6));
%% Լ������ �����г�
con=[];
for t=1:1:nt
    con=con+[(sum(P(:,t))-sum(P_pump_PSH(:,t))==sum(Pdtt(:,t)))];%����ƽ�ⷽ��

    con=con+[sum(Pmax.*I(:,t))>=1.2*(sum(Pdtt(:,t)))];%����ϵͳ������
    con=con+[sum(Pmin.*I(:,t))<=0.9*(sum(Pdtt(:,t)))];%����ϵͳ������  
     
    con=con+[P(:,t)+spinR_up(:,t)+Pfu2(:,t)<=Pmax.*I(:,t)];%��������������Լ��������ת����
    con=con+[Pmin.*I(:,t)<=P(:,t)-spinR_dn(:,t)-Pfu2(:,t)];%��������������Լ��������ת���� 
            
    con=con+[0<=spinR_up(:,t)<=10*ramp.*I(:,t)];%���������ת����������Լ��
    con=con+[0<=spinR_dn(:,t)<=10*ramp.*I(:,t)];%���������ת����������Լ��
    
    con=con+[y(:,t)+z(:,t)<=1];
    
    con=con+[P_pump_PSH(:,t)==P_pump_Cos*I_pump_PSH(:,t)];
    con=con+[max(I_gen_PSH(:,t))+max(I_pump_PSH(:,t))<=1];
    con=con+[V_min<=V(t)]+[V(t)<=V_max];%����������Լ��
    con=con+[V(nt)==V_initial];  %�յ���

    if t==1
        con=con+[I(:,t)-I0==y(:,t)-z(:,t)];
        con=con+[P(1:2,t)-P0(1:2)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%ú������������������Լ��
        con=con+[P0(1:2)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%ú������������������Լ��
        con=con+[P(3:6,t)-P0(3:6)<=15*ramp(3:6)];%�����������������������Լ����������С�������ƣ�
        con=con+[P0(3:6)-P(3:6,t)<=15*ramp(3:6)];%ú������������������Լ����������С�������ƣ�
        
       % ����Լ��
        con=con+[V(t)==V_initial-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];%����
      
    else
        con=con+[I(:,t)-I(:,t-1)==y(:,t)-z(:,t)];
        con=con+[P(1:2,t)-P(1:2,t-1)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%������������������Լ��
        con=con+[P(1:2,t-1)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%������������������Լ��
        con=con+[P(3:6,t)-P(3:6,t-1)<=15*ramp(3:6)];%�����������������������Լ����������С�������ƣ�
        con=con+[P(3:6,t-1)-P(3:6,t)<=15*ramp(3:6)];%�����������������������Լ����������С�������ƣ�
        con=con+[V(t)==V(t-1)-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];
    end

    con=con+[SU(:,t)==su.*y(:,t)];%��������Լ��
    
    con=con+[sum(spinR_up(:,t))>=0.05*sum(Pdtt(:,t))];%��ת����������Լ����5%�ܸ���
    con=con+[sum(spinR_dn(:,t))>=0.05*sum(Pdtt(:,t))];%��ת����������Լ����5%�ܸ���
    

    con=con+[Pl(:,t)==Hg*P(:,t)-H*Pdtt(:,t)-Hg_PSH*P_pump_PSH(:,t)];%��·����Լ��
    con=con+[(Pl(:,t)<=Plmax)];%��·����Լ��
    con=con+[(-Pl(:,t)<=Plmax)];%��·����Լ��

end
% 
% ��Ƶ�г�ʱ��߶�
for t=1:1:nf
    I_f(:,t)=sum(I(:,nt/nf*t-3:nt/nf*t),2);%�����г���״̬��������
    con=con+[sum(Pfu(:,t))==P_f(t)];%��Ƶ����Լ��
    con=con+[Pfu(:,t)>=0];%��Ƶ��������
    con=con+[Pfu(:,t)<=3*v.*max(I_f(:,t)-3,0)];%��Ƶ��������1
    con=con+[Pfu(:,t)<=0.15*Pmax.*max(I_f(:,t)-3,0)];%��Ƶ��������2
end



% minimum on/off time constraints
% for i=1:1:2
%     if UT(i)>0   %UTΪ���ֿ���������ʱ��
%         con=con+[UT(i)-sum(I(i,1:1:UT(i)))==0];
%      for j=UT(i)+1:1:nt-Ton(i)+1
%         con=con+[sum(I(i,j:1:j+Ton(i)-1))>=Ton(i)*y(i,j)];
%      end 
%     for j=nt-Ton(i)+2:1:nt
%         con=con+[sum(I(i,j:1:nt))-(nt-j+1)*y(i,j)>=0];
%     end
%     end
%     
%     if DT(i)>0   %DTΪ���ֹػ�������ʱ��      
%         con=con+[sum(I(i,1:1:DT(i)))==0];
%     for j=DT(i)+1:1:nt-Toff(i)+1
%         con=con+[Toff(i)-sum(I(i,j:1:j+Toff(i)-1))>=Toff(i)*z(i,j)];
%     end
%     for j=nt-Toff(i)+2:1:nt
%         con=con+[(nt-j+1)-sum(I(i,j:1:nt))-(nt-j+1)*z(i,j)>=0];
%     end
%     end
% end


%% ���
%ѡ�������
ses = sdpsettings('verbose','2','solver','gurobi','debug','1');
%�������ֵ
ses.showprogress=1;
% ses.gurobi.MIPGAP=0.0002;
%����
dd=solvesdp(con,obj,ses);

if dd.problem~=0
    error('ģ�Ͳ�����!');
end

%% �̶�0/1����
I=value(I);
I_f=value(I_f);
y=value(y);
z=value(z);
I_pump_PSH=value(I_pump_PSH);%��ˮ״̬


%% ���õ��ȣ�����ڵ���
%% Ŀ�꺯��
obj2=0;%Ŀ�꺯��

% c=gencost(:,4);���ڿ��ٲ���
% obj2=obj2+c'*sum(P,2);

for t=1:1:nt
   obj2=obj2+sum(P(:,t).*gencost(:,4));
   for m=5:2:Len.Cost-2%ÿ���걨����Ķ˵㣨��ӦMW��
   detaP2=max(P(:,t)-gencost(:,m),0);%ȡ����������ÿһʱ��λ��ĳһ����ĳ�����С
   obj2=obj2+sum(detaP2.*gencost(:,m+1));
   end
end

obj2=obj2/4;

for t=1:1:nf
   obj2=obj2+sum(Pfu(:,t).*C_f_Pi(:,t));%��Ƶ����ǵ�Ƶ����3��
end

%% Լ������
con2=[];
for t=1:1:nt
    P_pump_PSH(:,t)=P_pump_Cos*I_pump_PSH(:,t);
    con2=con2+[(sum(P(:,t))-sum(P_pump_PSH(:,t))==sum(Pdtt(:,t))):['lambda' num2str(t)]];%����ƽ�ⷽ��
    
    con2=con2+[P(:,t)+Pfu2(:,t)+spinR_up(:,t)<=Pmax.*I(:,t)];%��������������Լ��������ת����
    con2=con2+[Pmin.*I(:,t)<=P(:,t)-Pfu2(:,t)-spinR_dn(:,t)];%��������������Լ��������ת����
    
    con2=con2+[0<=spinR_up(:,t)<=10*ramp.*I(:,t)];%���������ת����������Լ��
    con2=con2+[0<=spinR_dn(:,t)<=10*ramp.*I(:,t)];%���������ת����������Լ��
    
    con2=con2+[V_min<=V(t)]+[V(t)<=V_max];%����������Լ��
    con2=con2+[V(nt)==V_initial];  %�յ���
    
    if t==1
        con2=con2+[P(1:2,t)-P0(1:2)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%ú������������������Լ��
        con2=con2+[P0(1:2)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%ú������������������Լ��
        con2=con2+[P(3:6,t)-P0(3:6)<=15*ramp(3:6)];%�����������������������Լ����������С�������ƣ�
        con2=con2+[P0(3:6)-P(3:6,t)<=15*ramp(3:6)];%ú������������������Լ����������С�������ƣ�
        
       % ����Լ��
        con2=con2+[V(t)==V_initial-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];%����
      
    else
        
        con2=con2+[P(1:2,t)-P(1:2,t-1)<=15*ramp(1:2).*(1-y(1:2,t))+Pmin(1:2).*y(1:2,t)];%������������������Լ��
        con2=con2+[P(1:2,t-1)-P(1:2,t)<=15*ramp(1:2).*(1-z(1:2,t))+Pmin(1:2).*z(1:2,t)];%������������������Լ��
        con2=con2+[P(3:6,t)-P(3:6,t-1)<=15*ramp(3:6)];%�����������������������Լ����������С�������ƣ�
        con2=con2+[P(3:6,t-1)-P(3:6,t)<=15*ramp(3:6)];%ú������������������Լ����������С�������ƣ�
        con2=con2+[V(t)==V(t-1)-UpFactor*sum(P_gen_PSH(:,t))+DownFactor*sum(P_pump_PSH(:,t))];
    end

    con2=con2+[sum(spinR_up(:,t))>=0.05*sum(Pdtt(:,t))];%��ת����������Լ����8%�ܸ���
    con2=con2+[sum(spinR_dn(:,t))>=0.05*sum(Pdtt(:,t))];%��ת����������Լ����8%�ܸ���

    con2=con2+[Pl(:,t)==Hg*P(:,t)-H*Pdtt(:,t)-Hg_PSH*P_pump_PSH(:,t)];
    con2=con2+[(Pl(:,t)<=Plmax):['mup' num2str(t)]];%��·����Լ��
    con2=con2+[(-Pl(:,t)<=Plmax):['mun' num2str(t)]];%��·����Լ��
    

end

%��Ƶ�г�ʱ��߶�
for t=1:1:nf
    con2=con2+[(sum(Pfu(:,t))==P_f(t)):['lambdaf' num2str(t)]];%��Ƶ����Լ��   
    con2=con2+[Pfu(:,t)>=0];%��Ƶ��������
    con2=con2+[Pfu(:,t)<=3*v.*max(I_f(:,t)-3,0)];%��Ƶ��������1
    con2=con2+[Pfu(:,t)<=0.15*Pmax.*max(I_f(:,t)-3,0)];%��Ƶ��������2
end


%ѡ�������
ses2 = sdpsettings('verbose','2','solver','gurobi','debug','1');
%�������ֵ
ses2.showprogress=1;
% ses2.gurobi.MIPGAP=0.0002;
%����
dd2=solvesdp(con2,obj2,ses2);

if dd2.problem~=0
    error('ģ�Ͳ�����!');
end

%% ������
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

% P=double(P);%�����й�����
Pfu=double(Pfu);
% Pfu2=double(Pfu2);
% detaP=double(detaP);
% Pl=double(Pl);
% spinR_up=double(spinR_up);%����ת����
% spinR_dn=double(spinR_dn);%����ת����
% SU=double(SU);%��������
% SD=double(SD);%�ػ�����
% I=double(I);%ON/OFF status, 0-1
% y=double(y);%startup indicator, 0-1
% z=double(z);%shutdown indicator, 0-1
P_pump_PSH=double(P_pump_PSH);
P_gen_PSH=double(P_gen_PSH);
% V=double(V);
% I_pump_PSH=double(I_pump_PSH);
% I_gen_PSH=double(I_gen_PSH);


% obj=double(obj);%�����ܷ���
% obj2=double(obj2);

% P_sc1=double(P_sc1);
% P_sc2=double(P_sc2);
% P_sc3=double(P_sc3);
% P_sc4=double(P_sc4);
%��������
% Fee_THE_SU=sum(sum(SU(1:2,:)));
% Fee_LNG_SU=sum(sum(SU(3:4,:)));


%�����г������ֵ
% Fee_THE=0;
% Fee_LNG=0;
Fee_PSH=0;
% ��Ƶ���������г������ֵ
% Feef_THE=0;
% Feef_LNG=0;
Feef_PSH=0;

%��������
%���������г�
for t=1:1:nf
% Feef_THE=Feef_THE+lambdaf(t)*sum(Pfu(1:2,t));
% Feef_LNG=Feef_LNG+lambdaf(t)*sum(Pfu(3:4,t));
Feef_PSH=Feef_PSH+lambdaf(t)*sum(Pfu(5:6,t));
end


%�����г�
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


%% ���ؽ�� �����ˮ��������


Income_PSH=double(Income_PSH);

yalmip('clear');
end
