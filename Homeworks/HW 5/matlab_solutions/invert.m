function Ainv =  invert(A)
% function to invert a mtrix using the Gauss Jordan method
    [U, B1] = GaussianElimination(A,eye(size(A,1)));
    [I, Ainv]  = upwardGaussianElimination(U,B1);
end

function [U, B1] = GaussianElimination(A,B )
%Gaussian Elimination
% Input: A matrix to invert
%        B matrix
% Output:U upper triangular Matrix
%        B1 corresponding matrix

    % obtaining the sizes of the matrices
    [n,m] = size(A); 
    % checking for non square matrices
    % and incompatible dimensions
    if n~=m
        error('Non-square matrix \n');
    elseif m~=size(B,1)
        error('Dimension mismatch \n');
    end
    
    % building the augmented matrix
    M = [A, B];
    % loop over the columns
    for i = 1:n-1
        % finding the non zero pivots
        indNonZero = find(abs(M(i:end,i))>eps);
        % if no-pivots availables throw an error
        if isempty(indNonZero)
            error('Singular Matrix');
        end
        % compute the index of the first non zero pivot
        indj = i-1+indNonZero(1);
        % performing row swap if necessary
        if i ~= indj
           buffer = M(i,:);
           M(i,:) = M(indj,:);
           M(indj,:) = buffer;
        end
        % eliminating the entried below the diagonal
        % vectorized form to gain speed
        M(i+1:end,:) = M(i+1:end,:) - (M(i+1:end,i)*(M(i,:)/M(i,i)));
    end
    if abs(M(n,n)) < eps 
        error('Singular Matrix');
    end
    % splitting the matrix
    U = M(1:n,1:n);
    B1 = M(:,n+1:end);
end

function [I, X] = upwardGaussianElimination(U,B1)
    % function to compute a triangular solve
    % Input:   U,  upper triangular non-singular matrix
    %          B,  right-hand side
    % Ouput:   I, this should be almost an identity matrix
    %          X,  solution to, UX = B1
    
    [n,m] = size(U); 
    if n~=m
        error('Non-square matrix \n');
    elseif m~=size(B1,1)
        error('Dimension mismatch \n');
    end
    % allocating the answer
   % building the augmented matrix
    M = [U, B1];
    % loop over the columns
    for i = n:-1:1
        % normalizing the diagonal
        M(i,:)  =  M(i,:)/M(i,i);
        % eliminating the entries above the diagonal
        % vectorized form to gain speed
        M(1:i-1,:) = M(1:i-1,:) - M(1:i-1,i)*(M(i,:));
    end
    
    % splitting the matrix
    I = M(1:n,1:n);
    X = M(:,n+1:end);
end