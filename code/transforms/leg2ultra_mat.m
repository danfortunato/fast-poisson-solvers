function S = leg2ultra_mat( n )
% Conversion matrix from Legendre coefficients to C^(3/2).
%
% Given coefficients in the Legendre basis the C^(3/2) coefficients
% can be computed via
%
%     c = rand(10, 1);    % Legendre coefficients
%     S = leg2ultra_mat( length(c) ); % conversion matrix
%     d = S * c;           % C^(3/2) coefficients
%
% Alex Townsend, 5th May 2016

lam = 1/2;
dg = lam./(lam + (2:n-1))';
v  = [1 ; lam./(lam+1) ; dg];
w  = [0 ; 0 ; -dg];
S  = spdiags( [v w], [0 2], n, n );

end
