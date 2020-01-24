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

function y=Mutate(x,mu,VarMin,VarMax,VarMin2,VarMax2)

    nVar=numel(x);
    
    nmu=ceil(mu*nVar);
    
    j=randsample(nVar,nmu);
    sigma=zeros(nVar,1);
    sigma(1)=0.1*(VarMax-VarMin);
    sigma(2)=0.1*(VarMax2-VarMin2);
    
    while(1)
        y=x;
        y(j)=x(j)+sigma(j)*randn(size(j));
        
        y(1:5)=min(max(y(1:5),VarMin),VarMax)        
        y(6)=min(max(y(6),VarMin2),VarMax2)
        if all(sort(y(1:5)==y(1:5)))
            return 
        end
    end

end