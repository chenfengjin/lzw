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

function [y1, y2]=Crossover(x1,x2,gamma,VarMin,VarMax,VarMin2,VarMax2)
    alpha=unifrnd(-gamma,1+gamma,size(x1));
    y1=alpha.*x1+(1-alpha).*x2;
    y2=alpha.*x2+(1-alpha).*x1;
    y1(1:5)=min(max(y1(1:5),VarMin),VarMax);
    y2(1:5)=min(max(y2(1:5),VarMin),VarMax);
    
    y1(16)=min(max(y1(16),VarMin2),VarMax2);
    y2(16)=min(max(y2(16),VarMin2),VarMax2);
    y1(17)=min(max(y1(17),VarMin2),VarMax2);
    y2(17)=min(max(y2(17),VarMin2),VarMax2);
    y1(18)=min(max(y1(18),VarMin2),VarMax2);
    y2(18)=min(max(y2(18),VarMin2),VarMax2);

    y1(1:5)=sort(y1(1:5))
    y2(1:5)=sort(y2(1:5))
    y1(6:10)=sort(y1(6:10))
    y2(6:10)=sort(y2(6:10))
    y1(11:15)=sort(y1(11:15))
    y2(11:15)=sort(y2(11:15))
end