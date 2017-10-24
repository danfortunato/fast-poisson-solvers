function M = leg2cheb_mat( N )
% Construct the leg2cheb conversion matrix.

% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

M = zeros(N, N);
for j = 0:N-1
    for k = j:2:N-1
        M(j+1, k+1) = 2/pi*vals((k-j)+1).*vals((k+j)+1);
    end
end
M(1,:) = .5*M(1,:);

end
