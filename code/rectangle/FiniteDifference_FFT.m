function X = FiniteDifference_FFT( f, n )
% Fast Poisson solver using 2nd order central differences.
% 
%  Solves    -laplacian(u) = f   with zero Dirichlet on 
%               on [0 1]x[0 1]   with an nxn grid. 
%
% Alex Townsend, September 2014. Updated September 2015. 

% Evaluate rhs on grid:
x = linspace(0, 1, n); x([1,end])=[];
[xx, yy] = meshgrid(x);
F = f(xx, yy) ;

% Set up and solve matrix equation KX + XK = F:
h = 2/(n-2); s = sqrt(2/(n-1));
d = -(2-2*cos((1:n-2)'*pi/(n-1)))/h^2;     % Eigenvalues of K.
F = s*chebfun.dst(s*chebfun.dst(F,1).',1).';          % S \ F * S
Y = F./(d*ones(1,n-2) + ones(n-2,1)*d');% Divide by eigvals
X = s*chebfun.dst(s*chebfun.dst(Y,1).',1).';          % X = S \ Y * S

% Add back in zero dirichlet conditions
Z = zeros(1,n-2);
X = [0   Z   0 ;
     Z'  X   Z';
     0   Z   0 ];
X = real( X );

% Plot
% x = linspace(0, 1, n);
% [xx, yy] = meshgrid(x);
% surf(xx, yy, X,'edgealpha',0,'facecolor', 'interp'),
% view(0,90), axis square, set(gca, 'fontsize', 16)
end
