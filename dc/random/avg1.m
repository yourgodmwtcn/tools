function [A] = avg1(X)
    if size(X,1) == 1
        X = X';
        m = 1;
    else m = 0;
    end
    
    A = (X(1:end-1,:) + X(2:end,:))/2;
    if m, A = A'; end