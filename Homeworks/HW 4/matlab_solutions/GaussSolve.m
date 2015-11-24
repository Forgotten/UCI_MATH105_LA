function x =  SolveGauss(A,b)
    [U, b1] = GaussianElimination(A,b);
    x = BackwardSubstitution(U,b1);
end

function [U, b1] = GaussianElimination(A,b )
%Gaussian Elimination
% Input: A matric
%        b vector
% Output:U upper triangular Matrix
%        b1 corresponding vector

    % obtaining the sizes of the matrices
    [n,m] = size(A); 
    % checking for non square matrices
    % and incompatible dimensions
    if n~=m
        error('Non-square matrix \n');
    elseif m~=size(b,1)
        error('Dimension mismatch \n');
    end
    
    % building the augmented matrix
    M = [A, b];
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
    % splitting the matrix
    U = M(1:n,1:n);
    b1 = M(:,n+1);
end

function x = BackwardSubstitution(U,b1)
    % function to compute a triangular solve
    % Input:   U,  upper triangular non-singular matrix
    %          b1, right-hand side
    % Ouput:   x,  solution to, Ux = b1
    
    [n,m] = size(U); 
    if n~=m
        error('Non-square matrix \n');
    elseif m~=size(b1,1)
        error('Dimension mismatch \n');
    end
    % allocating the answer
    x = zeros(n,1);
    % solving the last element
    x(n) = b1(n)/U(n,n);
    % perform the backsubstitution
    for i=n-1:-1:1
       x(i) = (b1(i) - U(i,i+1:end)*x(i+1:end))/U(i,i);
    end
    
end