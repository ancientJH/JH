function [xi, w] = GetQuadratureRule(ElementType, nQuad)
switch ElementType
    case 'P11D'
        switch nQuad
            case 2
                xi = [-1, 1]/sqrt(3);
                w = [1, 1];
        end
     case 'P21D'
        switch nQuad
            case 2
                xi = [-1, 1]/sqrt(3);
                w = [1, 1];
        end
        
    case 'P12D'
        switch nQuad
            case 3
                xi = [2/3, 1/6, 1/6; 1/6, 1/6, 2/3];
                w = [1, 1, 1]/6;%三个积分点，每个积分点的权重是1/6，面积坐标
        end 
        
    case 'P22D'
        switch nQuad
            case 4
                xi = [1/3, 1/5, 1/5, 3/5; 1/3, 1/5, 3/5, 1/5];
                w = [-27, 25, 25, 25]/96;
        end 
        
    case 'Q12D'
        switch nQuad
            case 4
                xi = [-1, 1, 1, -1; -1, -1, 1, 1]/sqrt(3);
                w = [1, 1, 1, 1];
        end                
end