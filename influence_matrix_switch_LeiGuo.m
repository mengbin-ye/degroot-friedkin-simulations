function [ C, c ] = influence_matrix_switch_LeiGuo(n)
%This function creates the influence matrix which is to be used. 
%Input:    n = the index of the graph
%Output:   C = the relative influence matrix
%          c = the dominant left eigenvector of C. Because C is
%          doubly-stochastic, c should be the vector of ones (normalised)

switch n
    case 1
        %Graph is a cyclic graph
        C = [0 0 0 1;
            1 0 0 0;
            0 1 0 0;
            0 0 1 0];

        
    case 2
                %Relative influence matrix (adjacency matrix)
        C = [0 0.5 0 0.5;
            0 0 1 0;
            1 0 0 0;
            0 0 1 0];

    case 3
         C = [0 1 0 0;
            0.9 0 0 0.1;
            0 0.9 0 0.1;
            0.9 0 0.1 0];
end


% Compute dominant left eigenvector of C
[a,b] = eig(C');
c = diag(b);    %Obtain eigenvalues of C' as a vector
[~,e] = max(c);   %Require index, e, indicating dominant eigenvalue of C
w = a(:,e); w = abs(w);   %Compute dominant left eigenvector, ensure positivity
w = w./norm(w,1);   %Normalise w
if (0.99 <= w'*ones(length(w),1)) && ( w'*ones(length(w),1) >= 1.01)
    disp('Wrong C Eigenvector Found')
end
c = w;
clear a b w e


end

