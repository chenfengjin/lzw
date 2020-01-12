%% ˵����������Ϊ�������г��͵�Ƶ�г����γ��壬���������ˮ����Ч��??
%case39 33??7Ϊ���ܣ�31??8Ϊ������30??5Ϊ��??32??6Ϊ��??
function Income_PSH=onepass2price(quoted_prices)
try
        disp(quoted_prices)
        quoted_price_energy=quoted_prices(1:5);
        quoted_price_frequency=quoted_prices(6);
    %% ����ԭʼ����
    % ѡ��39�ڵ��׼����
    mpc=case39lzw_tiao;
    % ���롢ȷ���ο���??
    [~, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);%��??
    %% ���ô�??�����ı�����??
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
    % [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);%matpowerͨ�ú����������ο��ڵ㣬PV �ڵ�??PQ �ڵ�Ľڵ���??
    % ���븺������
    Pforecast0=0.2*xlsread('qifuwo',1,'A1:A24'); %����24��ʱ�̵ĸ�������

    %% ϵͳ��������
    Horizon = 96;%�������г�ʱ���??
    T= 24;%��Ƶ�г�ʱ��߶�
    t=1:Horizon;
    t0=1:T;
    Pforecast=interp1(t0,Pforecast0,1+(t-1)*0.25,'PCHIP');  %��??�õ�96��ʱ�̵ĸ���
    % Pforecast=Pforecast0';
    % Max_Pforecast=max(Pforecast);%����????ֵ�����ڵ���
    % Min_Pforecast=min(Pforecast);%����????ֵ�����ڵ���

    %ϵͳ����
    Rup   =0.01*Pforecast;%����??/20
    Rdown =0.01*Pforecast;
    % SRup  =0.005*Pforecast;%ʱ��t�ϵ���ת����Ҫ��  
    % SRdown=0.005*Pforecast;%ʱ��t�µ���ת����Ҫ��

    %�ڵ�??
    N_bus=size(bus,1);
    %��·??
    N_branch=size(branch,1);

    %% ����??
    %1.����ڵ���
    %????����ڵ���
    BusN_G=gen(:,1);%���������������????�ڵ���
    %���ڵ�������
    NN=cell(3,1);
    NN(1)={[30,35]}; %���ڷ����������д���????����ʱ���ֳ�????
    %�����ڵ�������
    NN(2)={[31,38]};
    %����ڵ�������
    NN(3)={[33,37]};
    %������ڵ���
    genis_THE=ismember(BusN_G,NN{1}); 
    BusN_THE=BusN_G(genis_THE);
    %�����ڵ���
    genis_LNG=ismember(BusN_G,NN{2});
    BusN_LNG=BusN_G(genis_LNG);
    %�������ڵ���
    genis_PSH=ismember(BusN_G,NN{3});
    BusN_PSH=BusN_G(genis_PSH);
    % %������ڵ���
    % BusN_wind=[32,36];



    %��������
    N_ALL=size(BusN_G,1);
    N_THE=size(BusN_THE,1);
    N_LNG=size(BusN_LNG,1);
    N_PSH=size(BusN_PSH,1);

    %�����������??
    Pmax_ALL=gen(:,9);%????����ĳ������ޣ���й�����
    Pmin_ALL=gen(:,10);%����
    Pmaxt_ALL=repmat(Pmax_ALL,1,Horizon);
    Pmint_ALL=repmat(Pmin_ALL,1,Horizon);
    % Sum_Pmax_ALL=sum(Pmax_ALL);%????����֮�ͣ����ڵ�??
    % Sum_Pmin_ALL=sum(Pmin_ALL);%????����֮�ͣ����ڵ�??
    % Pmax_THE=max(Pmax_ALL(genis_THE));
    % Pmax_LNG=max(Pmax_ALL(genis_LNG));
    % Pmax_PSH=max(Pmax_ALL(genis_PSH));

    %��������??
    Ramp_U_ALL=gen(:,17);
    Ramp_D_ALL=gen(:,18);


    %????????/ͣ��ʱ��
    Ton_min_ALL=gen(:,20);
    % Toff_min_ALL=gen(:,21);
    Ton_min_THE=Ton_min_ALL(genis_THE);
    Toff_min_THE=Ton_min_ALL(genis_THE);
    %��������޴�Լ??

    %��ˮ����
    %��ɢ��ˮ����
    % P_pump1_PSH=200;
    P_pump1_PSH=100;
    %�����ˮ���ܻ������??
    V_min=3.782*10^5;
    V_max=1.949*10^6;
    V_initial=9.782*10^5;
    V_mint=repmat(V_min,N_PSH,Horizon);
    V_maxt=repmat(V_max,N_PSH,Horizon);
    %�����ˮ״??��ˮ��ת��ϵ?? 
    Factor_pump=150;
    %���巢��״??��ˮ��ת��ϵ??
    Factor_gen=200;
    % %���
    % P_wind=xlsread('gen_wind.xlsx');
    % A=Pforecast-sum(P_wind);
    % A_MAX=max(A);
    % A_MIN=min(A);
    % figure
    % P2=plot(value(sum(P_wind))');
    % figure
    % P=plot(value(A)');

    %����(��α��ۣ���ɱ��й�)
    %���濴���ܲ��ܱൽǰ����
    Len.Cost=size(gencost,2);

    %״??���� ״??�ͳ�??
    U_ALL=binvar(N_ALL,Horizon,'full');
    P_ALL=sdpvar(N_ALL,Horizon,'full');

    U_THE=U_ALL(genis_THE,:);
    %��ˮ���ܻ���
    U_pump_PSH=binvar(N_PSH,Horizon,'full');%��ˮ״??
    P_pump_PSH=sdpvar(N_PSH,Horizon,'full');%��ˮ����
    ADD=binvar(N_PSH,Horizon,'full');%��ˮ���ʸ��ӱ���
    V=sdpvar(N_PSH,Horizon,'full');%����

    %������ͣ
    % yg(:,1)=U_ALL(:,1);%��t=1ʱ�̻���????�����ʱ����������Ϊ1
    % yg(:,2:Horizon)=diff(U_ALL,1,2);
    % L_on=gencost(:,1);%ÿ������ÿһʱ����������
    % % L_off=gencost(:,2);%ÿ������ÿһʱ��ͣ������
    % yg_obj=abs(yg);

    %% �����г��Ĺ�??
    Nfu=N_ALL;%��Ƶ�����ṩ����??
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

    %% ���ϳ��塤��������??

    %Լ��
    Constraints=[];

    %1.ϵͳ
    % %(1)ϵͳ����ƽ��Լ��
    Constraints=[Constraints,sum(P_ALL)-sum(P_pump_PSH)==Pforecast]; 
    % %(1)ϵͳ����ƽ��Լ��
    % Constraints=[Constraints,sum(P_ALL)-sum(P_pump_PSH)+sum(P_sc1)-sum(P_sc2)==Pforecast]; 

    % Constraints=[Constraints,P_sc1>=0];
    % Constraints=[Constraints,P_sc2>=0];
    %(2)ϵͳ����������Լ??
    Constraints=[Constraints,sum(bsxfun(@times,U_ALL,Pmax_ALL))>=Pforecast+Rup];
    %(3)ϵͳ����������Լ??
    Constraints=[Constraints,sum(bsxfun(@times,U_ALL,Pmin_ALL))-sum(bsxfun(@times,U_pump_PSH,100))<=Pforecast-Rdown];
    % %(4)ϵͳ��ת����Ҫ�� 
    % Constraints=[Constraints,sum(bsxfun(@min,bsxfun(@plus,-P_ALL,Pmax_ALL),Ramp_U_ALL))>=SRup];
    % Constraints=[Constraints,sum(bsxfun(@min,max(bsxfun(@minus,P_ALL,Pmin_ALL),0),Ramp_D_ALL))>=SRdown];
    % 
    %2.����
    %(1)�������������Լ??
    Constraints=[Constraints,(bsxfun(@times,U_ALL,Pmin_ALL)<=P_ALL-Pfu2)&(P_ALL+Pfu2<=bsxfun(@times,U_ALL,Pmax_ALL))];
    % Constraints=[Constraints,P_sc3>=0];
    % Constraints=[Constraints,P_sc4>=0];

    Ramp_add=[0;0;Pmax_ALL(3);0;Pmax_ALL(5);0];
    Rampt_add=repmat(Ramp_add,1,Horizon-1);
    %(2)��������Լ��
    Constraints=[Constraints,diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,1:Horizon-1),Ramp_U_ALL)+bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,2:Horizon)),Pmax_ALL)];
    Constraints=[Constraints,-diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,2:Horizon),Ramp_D_ALL)-bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,1:Horizon-1)),Pmax_ALL)];

    %(3)����������Լ��
    %  1)����????����????ʱ��Լ��
    for t=2:Horizon 
    for i=1:N_THE    
        k=max(t-Ton_min_THE(i),1):(t-1);
        LUP(i,t)=sum(U_THE(i,k));%ʱ��tʱ������????ʱ�䣨��t-1??
        k=max(t-Toff_min_THE(i),1):(t-1);
        LDOWN(i,t)=sum(1-U_THE(i,k));%ʱ��tʱ������ͣ��ʱ�䣨��t-1??
        Constraints=[Constraints,LDOWN(i,t)>=(U_THE(i,t)-U_THE(i,t-1))*Toff_min_THE(i)];%????����ͣ��ʱ��Լ��
        Constraints=[Constraints,LUP(i,t)>=(U_THE(i,t-1)-U_THE(i,t))*Ton_min_THE(i)];   %????����????ʱ��Լ��
    end
    end

    %??����ˮ��������Լ??
    %   1����ˮ����Լ??
        Constraints=[Constraints,P_pump_PSH==P_pump1_PSH*U_pump_PSH-100*ADD+100];

    %   2��״̬����Լ??
        Constraints=[Constraints,max(U_ALL(genis_PSH,:))+max(U_pump_PSH)<=1];
        Constraints=[Constraints,U_pump_PSH+ADD>=1];      
    %   3������Լ??
        Constraints=[Constraints,V(:,1)==V_initial+Factor_pump*P_pump_PSH(:,1)-Factor_gen*P_ALL(genis_PSH,1)];
        for t=2:Horizon
        Constraints=[Constraints,V(:,t)==V(:,t-1)+Factor_pump*P_pump_PSH(:,t)-Factor_gen*P_ALL(genis_PSH,t)];
        end
        Constraints=[Constraints,(V_mint<=V)&(V<=V_maxt)];
        Constraints=[Constraints,V(:,Horizon)==V_initial];  %�յ�??
    % 
    %(5)��·����Լ��
    H=makePTDF(mpc);%�ڵ����·�Ĺ���ת�Ʒֲ�����,����Ϊ֧·��������Ϊ�ڵ���??
    Hg=H(:,BusN_G);
    SLlf=sdpvar(N_branch,Horizon,'full');%��Ϊ������֧·????
    SLlz=sdpvar(N_branch,Horizon,'full');
    Plmax=branch(:,6);%��������rataA����������
    Plmax(17)=90;%��??����
    Plmax(28)=170;

    Plmaxt=repmat(Plmax,1,Horizon);
    Pd=bus(:,3);%��������ĸ�߸���??
    Pnforecast=repmat(Pforecast,N_bus,1);
    Pdt=bsxfun(@times,Pnforecast./sum(Pd),Pd);%ÿһ���ڵ�ÿ????ʱ�̵ĸ���??

    Constraints=[Constraints,-Plmaxt<=Hg*P_ALL-H(:,BusN_PSH)*P_pump_PSH-H*Pdt-SLlz+SLlf];
    Constraints=[Constraints,Hg*P_ALL-H(:,BusN_PSH)*P_pump_PSH-H*Pdt-SLlz+SLlf<=Plmaxt];

    Constraints=[Constraints,SLlz>=0];
    Constraints=[Constraints,SLlf>=0];

    %% ��Ƶ���������г�

    %��Ƶ????
    Cneed=0.04*max(Pforecast0-min(Pforecast0),max(Pforecast0)-Pforecast0);%����24Сʱ��ʱ�ε�Ƶ��������MW??max�������ú��ǰ������е�
    co_F=Cneed'./max(Cneed);

    %��Ƶ���ã����ۣ�
    Cfucost=[24;36;quoted_price_frequency;22;quoted_price_frequency;34]*co_F;
    %���Դ���24��ʱ�ν��е�Ƶ��̱��ۣ�??MW??????��λ0.1MW??

    %��������
    v=0.5*[3/100*Pmax_ALL(1);4.5/100*Pmax_ALL(2);Pmax_ALL(3);3/100*Pmax_ALL(4);Pmax_ALL(5);4.5/100*Pmax_ALL(6)];%���絥Ԫ��׼��������
    ki=[0.6;1;1.8;0.6;1.8;1];%��̨��������ۺϵ�Ƶ����ָ�꣬��С��0.5
    Pi=ki/max(ki);%��һ��??��ָ??
    Cfucost_Pi=bsxfun(@rdivide,Cfucost,Pi);%��Ƶ�������۸�=��Ƶ��̱���/ki


    %��Ƶ����Լ��
    Constraints=[Constraints,sum(Pfu)==Cneed'];%


    %״??��������
    for i=1:T %SCUC??6ʱ�Σ���Ƶ���������г�Ϊ24ʱ�Σ���4ʱ����ͣ״??����γ�fU_ALL
    fU_ALL(:,i)=sum((U_ALL(:,Horizon/T*i-3:Horizon/T*i)),2);
    end

    %��Ƶ�����걨����??

    Constraints=[Constraints,0<=Pfu];
    Constraints=[Constraints,Pfu<=0.5* min(repmat(v,1,T)*3,0.15*repmat(Pmax_ALL,1,T)).*max(fU_ALL-3,0)];


    %% Ŀ�꺯��
    Objective=0;

    % % ��ͣ�ɱ�                                                     
    % C_start=sum(sum(bsxfun(@times,L_on,yg_obj)));
    % Objective=Objective+C_start;%�ܳɱ���ʼ??Ϊ��������ɱ���������������??

    % % ���۳ɱ�
    % for m=(Len.Cost-2):-2:3%ÿ���걨����Ķ˵㣨��ӦMW??
    %     c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0)- max(bsxfun(@minus,P_ALL,gencost(:,m+2)),0);%ȡ����������ÿһʱ��λ��ĳһ����ĳ�����??
    %     Objective=Objective+sum(sum(bsxfun(@times,c,gencost(:,m+1)),2));
    % end

    for j=(Len.Cost-1):-2:6
        gencost(:,j)=gencost(:,j)-gencost(:,j-2);
    end

    for m=3:2:Len.Cost-2%ÿ���걨����Ķ˵㣨��ӦMW??
        c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0);%ȡ����������ÿһʱ��λ��ĳһ����ĳ�����??
        Objective=Objective+sum(sum( bsxfun(@times,c,gencost(:,m+1))));
    end

    % Objective=Objective+sum(sum(bsxfun(@times,[400;600;0;500;0;700],P_ALL)));


    M=1e9;%����Լ���ɳڷ���??
    Objective=Objective+M*sum(sum(SLlz+SLlf));
    
    % C_f=sum(sum(Cfucost_Pi.*Pfu+0.1*bsxfun(@times,Pfu,1-Pi)+0.01*bsxfun(@times,Pfu,max(ki)-ki) ));%Ŀ�꺯����Ƶ����????,�����絥Ԫ����۸���ͬʱ�����ȳ�??P ֵ�ߵķ��絥Ԫ�������絥??P ֵ��ͬʱ�����ȳ�??k ֵ�ߵķ��絥??
    C_f=sum(sum(Cfucost_Pi.*Pfu));%Ŀ�꺯����Ƶ����????,�����絥Ԫ����۸���ͬʱ�����ȳ�??P ֵ�ߵķ��絥Ԫ�������絥??P ֵ��ͬʱ�����ȳ�??k ֵ�ߵķ��絥??

    Objective=Objective+C_f;

    % % Objective=Objective+M*(sum(sum(P_sc1+P_sc2+P_sc3+P_sc4)))+M*(sum(sum(P_sc5+P_sc6+P_sc7+P_sc8)));
    % Objective=Objective+M*(sum(sum(P_sc1+P_sc2+P_sc3+P_sc4)))+M*(sum(sum(P_sc7+P_sc8)));

    %���
    ops=sdpsettings('solver','gurobi','verbos',1);
    % ops.showprogress=1;
    ops.gurobi.MIPGAP=0.0002;
    optimize(Constraints,Objective,ops);



    %�̶�0/1����
    U_ALL=double(U_ALL);
    U_pump_PSH=double(U_pump_PSH);
    fU_ALL=double(fU_ALL);
    ADD=double(ADD);


    %% ��ȫԼ�����õ���
    %Լ������
    Constraints2=[];

    %% ��Ƶ�г�
    %(1)����ƽ��Լ��
    Constraints2=[Constraints2,sum(P_ALL)-sum(P_pump_PSH)==Pforecast]; %Ϊ��ڵ��۷����������

    %��Ƶ����Լ��
    Constraints2=[Constraints2,sum(Pfu)==Cneed'];%

    %��Ƶ�����걨����??
    % Constraints2=[Constraints2,0.5* repmat(min(v*1,3/100*(Pmax_ALL-Pmin_ALL)),1,T).*max(fU_ALL-3,0)<=Pfu];
    Constraints2=[Constraints2,0<=Pfu];
    Constraints2=[Constraints2,Pfu<=0.5* min(repmat(v,1,T)*3,0.15*repmat(Pmax_ALL,1,T)).*max(fU_ALL-3,0)];


    %(2)ϵͳ��ת����Ҫ�� 
    % Constraints2=[Constraints2,sum(bsxfun(@min,bsxfun(@plus,-P_ALL,Pmax_ALL),Ramp_U_ALL))>=SRup];
    % Constraints2=[Constraints2,sum(bsxfun(@min,bsxfun(@max,(bsxfun(@minus,P_ALL,Pmin_ALL)),0),Ramp_D_ALL))>=SRdown];

    %(3)�������������Լ??
    % Constraints2=[Constraints2,(bsxfun(@times,U_ALL,Pmint_ALL)<=P_ALL+P_sc)&(P_ALL+P_sc<=bsxfun(@times,U_ALL,Pmaxt_ALL))];
    Constraints2=[Constraints2,(bsxfun(@times,U_ALL,Pmint_ALL)<=P_ALL-Pfu2)&(P_ALL+Pfu2<=bsxfun(@times,U_ALL,Pmaxt_ALL))];

    %(4)��������Լ��
    Constraints2=[Constraints2,diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,1:Horizon-1),Ramp_U_ALL)+bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,2:Horizon)),Pmax_ALL)];
    Constraints2=[Constraints2,-diff(P_ALL,1,2)<=Rampt_add+bsxfun(@times,U_ALL(:,2:Horizon),Ramp_D_ALL)-bsxfun(@times,diff(U_ALL,1,2),Pmin_ALL)+bsxfun(@times,(1-U_ALL(:,1:Horizon-1)),Pmax_ALL)];    
        
    %(5)��ˮ��������Լ��
    %   1����ˮ����Լ??
    Constraints2=[Constraints2,P_pump_PSH==P_pump1_PSH*U_pump_PSH-100*ADD+100]; 
    %   2������Լ??
    Constraints2=[Constraints2,V(:,1)==V_initial+Factor_pump*P_pump_PSH(:,1)-Factor_gen*P_ALL(genis_PSH,1)];
        for t=2:Horizon
    Constraints2=[Constraints2,V(:,t)==V(:,t-1)+Factor_pump*P_pump_PSH(:,t)-Factor_gen*P_ALL(genis_PSH,t)];
        end
    Constraints2=[Constraints2,(V_mint<=V)&(V<=V_maxt)];
    Constraints2=[Constraints2,V(:,Horizon)==V_initial];  %�յ�??
    %     
    % AAA=zeros(N_branch,Horizon);
    % AA=zeros(N_branch,1);
    % (6)��·����Լ��
    for i=1:N_branch
    Constraints2=[Constraints2,-Plmaxt(i,:)<=Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:)];
    Constraints2=[Constraints2,Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:)<=Plmaxt(i,:)];
    % AAA(i,:)=value(Hg(i,:)*P_ALL-H(i,BusN_PSH)*P_pump_PSH-H(i,:)*Pdt-SLlz(i,:)+SLlf(i,:));
    % AA(i)=min(Plmax(i,:)-abs(AAA(i,:)));
    % A=min(AA);
    end

    Constraints2=[Constraints2,SLlz>=0];
    Constraints2=[Constraints2,SLlf>=0];




    %Ŀ�꺯��
    Objective2=0;%�ܳ�??

    for m=3:2:Len.Cost-2%ÿ���걨����Ķ˵㣨��ӦMW??
        c=max(bsxfun(@minus,P_ALL,gencost(:,m)),0);%ȡ����������ÿһʱ��λ��ĳһ����ĳ�����??
        Objective2=Objective2+sum(sum( bsxfun(@times,c,gencost(:,m+1))));
    end

    % Objective2=Objective2+sum(sum(bsxfun(@times,[400;600;0;500;0;700],P_ALL)));

    % Price_pump=10;
    % C_p_PSH_pump=sum(sum(Price_pump*P_pump_PSH));
    % 
    % Objective2=Objective2-C_p_PSH_pump;

    %����Լ��
    Objective2=Objective2+M*sum(sum(SLlz+SLlf));
    

    C_f=sum(sum(Cfucost_Pi.*Pfu));%Ŀ�꺯����Ƶ����????,�����絥Ԫ����۸���ͬʱ�����ȳ�??P ֵ�ߵķ��絥Ԫ�������絥??P ֵ��ͬʱ�����ȳ�??k ֵ�ߵķ��絥??

    Objective2=Objective2+C_f;

    %���
    ops2=sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1);
    % ops2.showprogress=1;
    ops2.gurobi.MIPGAP=0.0002;
    optimize(Constraints2,Objective2,ops2);


    %% ���ڵ���
    mup=zeros(N_branch,Horizon);%��·????������Լ���������ճ���
    mun=zeros(N_branch,Horizon);
    for j=1:N_branch
        findlz=find(value(SLlz(j,:))~=0);%�ҳ���·�������ɳڱ����в�Ϊ��ĽǱ꣬��Ϊ������Խ��??
        findlf=find(value(SLlf(j,:))~=0);%�ҳ���·�������ɳڱ����в�Ϊ��ĽǱ꣬��Ϊ������Խ��??
        mun(j,:)=dual(Constraints2(end-2*N_branch-3+2*j))';
        mup(j,:)=dual(Constraints2(end-2*N_branch-2+2*j))';
        mup(j,findlz)=M;%��·������Խ��ʱ�����������ճ���Ϊ���糱��Լ���ɳڷ���??
        mun(j,findlf)=M;%��·������Խ��ʱ�����������ճ���Ϊ���糱��Լ���ɳڷ���??
    end
    Price_F=-dual(Constraints2(2));
    lameda=dual(Constraints2(1));
    LMP0=zeros(N_bus,Horizon);%����Ϊ�����ڽڵ���������Ϊʱ����
    for k=1:N_bus
        LMP0(k,:)=-lameda-sum(bsxfun(@times,mup,H(:,k))-bsxfun(@times,mun,H(:,k)));
    end
    LMP=zeros(N_bus,T);
    W_ALL=zeros(N_ALL,T);%����ÿСʱ���������??
    W_THE=zeros(N_THE,T);%
    W_LNG=zeros(N_LNG,T);%
    W_PSH=zeros(N_PSH,T);%
    W_PSH_gen=zeros(N_PSH,T);
    W_PSH_pump=zeros(N_PSH,T);
    P_THE=P_ALL(genis_THE,:);
    P_LNG=P_ALL(genis_LNG,:);
    P_PSH=P_ALL(genis_PSH,:);

    for Period=1:T%ÿСʱ��ƽ���ڵ���ΪÿСʱ????5���ӽڵ��۵�����ƽ��??
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




    %% ���ؽ���������ˮ������??

    Income_PSH=double(Income_PSH);
catch
    Income_PSH=0
    disp('onepass error,return 0 as default')
end
yalmip('clear');
end
