function X = cheb2ultra( X )
% CHEB2ULTRA    Convert vector of Chebyshev coefficients to C^(3/2)
%
%    CHEB2ULTRA(X) applies the conversion to each column of X if X is
%    matrix.

% First convert the matrix of Chebyshev coefficients to a matrix of
% Legendre coefficients:
m = size(X, 1);
if ( m <= 10000 ) % Determined experimentally
    S = cheb2leg_mat( m );
    X = S * X;
else
    X = cheb2leg( X );
end

% Now, convert the matrix of Legendre coefficients to a matrix of
% ultraspherical coefficients:
S = leg2ultra_mat( m );
X = S * X;

end
