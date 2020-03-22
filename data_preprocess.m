tic
clear all;
value=zeros(1,31);
for i=0:1:30
    value(i+1)=case39_YC_PSH_f([0,0,0,0,0,i]);
    disp(value(i+1))
    plot(value)
end
toc
