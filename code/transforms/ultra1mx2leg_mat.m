function S = ultra1mx2leg_mat( n )
% Conversion matrix for (1-x^2)C^(3/2) to Legendre.
%
% Given coefficients in the (1-x^2)C^(3/2) basis the Legendre coefficients
% can be computed via
%
%     c = rand(10, 1);     % (1-x^2)C^(3/2) coefficients
%     S = ultra1mx2leg_mat( length(c) ); % conversion matrix
%     c_leg = S * c;       % Legendre coefficients
%
% Alex Townsend, 5th May 2016

d = ones(n, 1);
S = spdiags(((1:n).*(2:(n+1))./2./(3/2:n+1/2))', 0, n, n);
S = spdiags( [d,-d], [0,-2], n, n ) * S;

end
