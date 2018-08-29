% WGG
% to change diag values to zero
function y = putDiagZeros(mat)
    y = mat;
    for i = 1:1:size(mat,1)
        for j = 1:1:size(mat,2)
            if i == j
               y(i,j) = 0; 
            end
        end
    end
end