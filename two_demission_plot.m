tic
x_range=0:10:220;
y_range=0:1:20;
[X,Y]=meshgrid(x_range,y_range);
[rows,columns] = size(X);
Z=zeros(rows,columns);
for i= 1:rows
    for j=1:columns
        x=X(i,j);
        y=Y(i,j);
        Z(i,j)=GD_YC([x,x,x,x,x,y]);
        mesh(X,Y,Z);%Ïß¿òÍ¼
        pause(0.01)
    end
end
toc
