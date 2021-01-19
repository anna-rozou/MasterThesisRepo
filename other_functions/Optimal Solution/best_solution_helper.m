function [matrix, gamma, matrix2] = best_solution_helper_18_11_20()

global M W G 
result = squeeze(num2cell(reshape(perms(1:2*M),[],2,M),2));

fprintf('-> We are inside best_solution_helper <-\n');

matrix = zeros(M,2,length(result));
gamma = zeros(M,2,length(result));
for k5=1:length(result)
    for k6=1:M
        r1 = result{k5,k6};
        matrix(k6,:,k5) = r1(:);
        
        gamma(k6,1,k5) = G(matrix(k6,1,k5),k6);
        gamma(k6,2,k5) = G(matrix(k6,2,k5),k6);
    end
end
matrix2 = matrix;
k5=1;
while k5<=length(matrix)
    %for k5=1:length(result)
    
    for k6=1:M
        if gamma(k6,1,k5) > gamma(k6,2,k5)* W(k6,2)/W(k6,1)
            cond = 0;
        else
            cond = 1;
            break;
        end
    end
    if cond == 1
        matrix(:,:,k5) = [];
        gamma(:,:,k5) = [];
    else
        k5=k5+1;
    end
end

check1 = 0;
for k5=1:length(matrix)
    for k6=1:M
        if gamma(k6,2,k5)>gamma(k6,1,k5)* W(k6,2)/W(k6,1)
            check1 = check1 +1;
        end
    end
end
if check1 == 0
    fprintf('\nCORRECT ALLOCATION !!!')
else
    fprintf('\nWRONG ALLOCATION !!!')
end

end





    
    