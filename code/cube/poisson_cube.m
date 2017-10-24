function X = poisson_cube( F, tol )
%POISSON_CUBE   Fast Poisson solver for the cube.
%   POISSON_CUBE( F ) solves laplacian(U) = F on [-1,1]x[-1,1]x[-1,1] with
%   zero Dirichlet boundary conditions. That is, U satisfies
%
%     U_{x,x} + U_{y,y} + U_{z,z} = F, on [-1,1]^3    U = 0 on boundary
%
%   F is input as an M x N x P matrix of Chebyshev coefficients. The
%   equation is solved using an M x N x P discretization.
%
%   POISSON_CUBE( F, TOL ) solves to the specified error tolerance.

% DEVELOPER'S NOTE:
%
% METHOD: Spectral method (in coefficient space). We use a C^{(3/2)} basis
% to discretize the equation, resulting in a discretization of the form
% (kron(A,A,I) + kron(A,I,A) + kron(I,A,A)) X(:) = F(:), where A is a
% symmetric tridiagonal matrix.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation is solved by a
% three-level nested alternating direction implicit (ADI) method.
%
% SOLVE COMPLEXITY:  O(M*N*P*log(MAX(M,N,P))^3*log(1/eps)) with M*N*P =
% total degrees of freedom. This is not theoretically justified.
% 
% AUTHORS: Dan Fortunato (dan.fortunato@gmail.com)
%          Alex Townsend (townsend@cornell.edu)
%
% The fast Poisson solver is based on:
%
% D. Fortunato and A. Townsend, Fast Poisson solvers for spectral methods,
% in preparation, 2017.

if ( nargin < 2 )
    tol = 1e-13;
end

[m, n, p] = size( F );

% Convert the RHS to C^{(3/2)} coefficients
F = cheb2ultra_tensor( F );

% Construct Tm, Tn, Tp
jjm = (0:m-1)';
dsub = -1./(2*(jjm+3/2)).*(jjm+1).*(jjm+2)*1/2./(1/2+jjm+2);
dsup = -1./(2*(jjm+3/2)).*(jjm+1).*(jjm+2)*1/2./(1/2+jjm);
d = -dsub - dsup;
Mm = spdiags([dsub d dsup], [-2 0 2], m, m);
invDm = spdiags(-1./(jjm.*(jjm+3)+2), 0, m, m);
Tm = invDm * Mm;

jjn = (0:n-1)';
dsub = -1./(2*(jjn+3/2)).*(jjn+1).*(jjn+2)*1/2./(1/2+jjn+2);
dsup = -1./(2*(jjn+3/2)).*(jjn+1).*(jjn+2)*1/2./(1/2+jjn);
d = -dsub - dsup;
Mn = spdiags([dsub d dsup], [-2 0 2], n, n);
invDn = spdiags(-1./(jjn.*(jjn+3)+2), 0, n, n);
Tn = invDn * Mn;

jjp = (0:p-1)';
dsub = -1./(2*(jjp+3/2)).*(jjp+1).*(jjp+2)*1/2./(1/2+jjp+2);
dsup = -1./(2*(jjp+3/2)).*(jjp+1).*(jjp+2)*1/2./(1/2+jjp);
d = -dsub - dsup;
Mp = spdiags([dsub d dsup], [-2 0 2], p, p);
invDp = spdiags(-1./(jjp.*(jjp+3)+2), 0, p, p);
Tp = invDp * Mp;

% Diagonally scale the RHS, as in the 2D Poisson solver
[ii, jj, kk] = ndgrid( jjm, jjn, jjp );
F = -F  ./ (ii.*(ii+3)+2) ./ (jj.*(jj+3)+2) ./ (kk.*(kk+3)+2);

Im = eye( m );
In = eye( n );
Ip = eye( p );

% Let U = kron(Tp,Tn,Im), V = kron(Ip,Tn,Tm), W = kron(Tp,In,Tm).
% Then Poisson can be written as (U + V + W) x = f with x = X(:), f = F(:).
% We will solve this using a three-level nested ADI iteration.

X = zeros( m, n, p );

% Calculate ADI shifts based on bounds on the eigenvalues of Tm, Tn, Tp
innertol = 1e-16;
a = @(n) -39*n^-4;
b = @(n) -4/pi^2;
% TODO: Calculate shifts based on m, n, and p
[p1, ~] = ADIshifts(b(n)^2, a(n)^2, -a(n)^2, -b(n)^2, tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADI LEVEL 1: Solve ((U+V) + W) x = f %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(p1)
    fprintf('Outer iteration: %g / %g\n', i, numel(p1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ADI LEVEL 2: Solve ((U+p1(i)I/2) + (V+p1(i)I/2)) x = f - (W-p1(i)I) x_i %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RHSi = F + p1(i)*X;
    for j = 1:n
        RHSi(:,j,:) = squeeze(RHSi(:,j,:)) - Tm*squeeze(X(:,j,:))*Tp';
    end
    [p2, ~] = ADIshifts(b(n)^2+p1(i)/2, a(n)^2+p1(i)/2, -a(n)^2+p1(i)/2, -b(n)^2+p1(i)/2, innertol);
    for j = 1:numel(p2)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ADI LEVEL 3: Solve ((U+p1(i)I/2) + p2(j)I) x = f - (W-p1(i)I) x_i - ((V+p1(i)I/2) - p2(j)I) x_j %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This decouples in one dimension, yielding: T*X_k*T' + s*X_k = F_k
        % We can rewrite this as: s*T\X_k + X_k*T' = T\F_k
        % and solve it fast using ADI.
        RHSj = RHSi - p1(i)/2*X + p2(j)*X;
        for k = 1:p
            RHSj(:,:,k) = squeeze(RHSj(:,:,k)) - Tm*squeeze(X(:,:,k))*Tn';
        end
        s = p1(i)/2 + p2(j);
        [p3, q3] = ADIshifts(-5/a(n)*s, -1/b(n)*s, a(n), b(n), innertol);
        for k = 1:m
            RHSk = squeeze(RHSj(k,:,:));
            Xk = zeros(n, p);
            FF = Tn \ RHSk;
            for l = 1:numel(p3)
                Xk = ( FF - (s*(Tn\Xk) + p3(l)*Xk) ) / (p3(l)*Ip - Tp');
                Xk = (s*In+q3(l)*Tn) \ ( RHSk - Tn*Xk*(q3(l)*Ip-Tp') );
            end
            X(k,:,:) = Xk;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ADI LEVEL 3: Solve ((V+p1(i)I/2) + p2(j)I) x = f - (W-p1(i)I) x_i - ((U+p1(i)I/2) - p2(j)I) x_j %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This decouples in one dimension, yielding: T*X_k*T' + s*X_k = F_k
        % We can rewrite this as: s*T\X_k + X_k*T' = T\F_k
        % and solve it fast using ADI.
        RHSj = RHSi - p1(i)/2*X + p2(j)*X;
        for k = 1:m
            RHSj(k,:,:) = squeeze(RHSj(k,:,:)) - Tn*squeeze(X(k,:,:))*Tp';
        end
        s = p1(i)/2 + p2(j);
        [p3, q3] = ADIshifts(-5/a(n)*s, -1/b(n)*s, a(n), b(n), innertol);
        for k = 1:p
            RHSk = squeeze(RHSj(:,:,k));
            Xk = zeros(m, n);
            FF = Tm \ RHSk;
            for l = 1:numel(p3)
                Xk = ( FF - (s*(Tm\Xk)+p3(l)*Xk) ) / (p3(l)*In-Tn');
                Xk = (s*Im+q3(l)*Tm) \ ( RHSk - Tm*Xk*(q3(l)*In-Tn') );
            end
            X(:,:,k) = Xk;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ADI LEVEL 2: Solve (W + p1(i)I) x = f - ((U+V)-p1(i)I) x_i %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This decouples in one dimension, yielding: T*X_k*T' + s*X_k = F_k
    % We can rewrite this as: s*T\X_k + X_k*T' = T\F_k
    % and solve it fast using ADI.
    RHSi = F + p1(i)*X;
    for j = 1:p
        RHSi(:,:,j) = RHSi(:,:,j) - Tm*X(:,:,j)*Tn';
    end
    for j = 1:m
        RHSi(j,:,:) = squeeze(RHSi(j,:,:)) - Tn*squeeze(X(j,:,:))*Tp';
    end
    s = p1(i);
    [p3, q3] = ADIshifts(-5/a(n)*s, -1/b(n)*s, a(n), b(n), innertol);
    for j = 1:n
        RHSj = squeeze(RHSi(:,j,:));
        Xk = zeros(m, p);
        FF = Tm \ RHSj;
        for l = 1:numel(p3)
            Xk = ( FF - (s*(Tm\Xk)+p3(l)*Xk) ) / (p3(l)*Ip-Tp');
            Xk = (s*Im+q3(l)*Tm) \ ( RHSj - Tm*Xk*(q3(l)*Ip-Tp') );
        end
        X(:,j,:) = Xk;
    end
end

% Convert back to Chebyshev
X = ultra1mx2cheb_tensor( X );

end

function X = cheb2ultra_tensor( X )
%CHEB2ULTRA_TENSOR Convert tensor Chebyshev coefficients to C^{(3/2)}
[m, n, p] = size( X );
for k = 1:p
    X(:,:,k) = cheb2ultra( cheb2ultra( X(:,:,k) ).' ).';
end
S = leg2ultra_mat( p );
for i = 1:m
    for j = 1:n
        X(i,j,:) = S * cheb2leg( cheb2leg( permute(X(i,j,:),[3 2 1]) ).' ).';
    end
end
end

function X = ultra1mx2cheb_tensor( X )
%ULTRA1MX2CHEB_TENSOR Convert tensor (1-x.^2)T_k coefficients to T_k
[m, n, p] = size( X );
for k = 1:p
    X(:,:,k) = ultra1mx2cheb( ultra1mx2cheb( X(:,:,k) ).' ).';
end
S = ultra1mx2leg_mat( p );
for i = 1:m
    for j = 1:n
        X(i,j,:) = leg2cheb( leg2cheb( S * permute(X(i,j,:),[3 2 1]) ).' ).';
    end
end
end
