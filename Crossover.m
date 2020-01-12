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
 
    
    y1(1)=max(y1(1),VarMin);
    y1(2)=max(y1(2),VarMin2);
    y1(1)=min(y1(1),VarMax);
    y2(1)=min(y1(2),VarMax2);
    
    y2(1)=max(y2(1),VarMin);
    y2(2)=max(y2(2),VarMin2);
    y2(1)=min(y2(1),VarMax);
    y2(2)=min(y2(2),VarMax2);

end