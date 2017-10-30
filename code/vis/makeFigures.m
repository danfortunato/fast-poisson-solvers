function makeFigures( name, writeToDisk )
%MAKEFIGURES   Make figures for the paper.
%
%   name: The name of the figure to generate.
%         Options are 'all' (default) or one of the following:
%         - 'FiniteDifferenceTimings'
%         - 'SquareSolution'
%         - 'SquareTimings'
%         - 'CylinderSolution'
%         - 'CylinderTimings'
%         - 'SphereSolution'
%         - 'CubeSolution'
%         - 'CubeTimings'
%
%   writeToDisk: Logical flag indicating whether to write the figure to
%                disk (1, default) or just display the figure (0).

close all
if ( nargin < 1 ), name = 'all'; end
if ( nargin < 2 ), writeToDisk = 1; end

figs = { 'FiniteDifferenceTimings' @FiniteDifferenceTimings 'vector' ;
         'SquareSolution'          @SquareSolutionFigure    'raster' ;
         'SquareTimings'           @SquareTimingsFigure     'vector' ;
         'CylinderSolution'        @CylinderSolutionFigure  'raster' ;
         'CylinderTimings'         @CylinderTimingsFigure   'vector' ;
         'SphereSolution'          @SphereSolutionFigure    'raster' ;
         'CubeSolution'            @CubeSolutionFigure      'raster' ;
         'CubeTimings'             @CubeTimingsFigure       'vector' };

for i = 1:length(figs)
    [figname, figfun, figformat] = figs{i,:};
    if ( strcmp(name, figname) || strcmp(name, 'all') )
        feval( figfun );
        writefig( figname, writeToDisk, figformat );
    end
end

end

function c = colorpalette()
%COLORPALETTE   Define the colors to be used for all figures.
    c = magma();
end

function writefig( filename, writeToDisk, format )
%WRITEFIG   Write a figure to disk.
shg
if ( writeToDisk )
    savefig(strcat(filename, '.fig'));
    if ( strcmp(format, 'vector') )
        print('-depsc', strcat(filename, '.eps'));
    elseif ( strcmp(format, 'raster') )
        print('-dpng', '-r500', strcat(filename, '.png'));
    end
end
end

function FiniteDifferenceTimings()
    nn = floor(logspace(1,3.7,50));
    t_adi3  = zeros(size(nn));
    t_adi6  = zeros(size(nn));
    t_adi13 = zeros(size(nn));
    t_fft   = zeros(size(nn));

    j = 1;
    for n = nn
        % ADI
        F = ones(n,n);

        s = tic;
        X_ADI3 = FiniteDifference_ADI(F, 1e-3);
        t_adi3(j) = toc(s);

        s = tic;
        X_ADI6 = FiniteDifference_ADI(F, 1e-6);
        t_adi6(j) = toc(s);

        s = tic;
        X_ADI13 = FiniteDifference_ADI(F, 1e-13);
        t_adi13(j) = toc(s);

        % FFT
        f = @(x,y) 1+0*x;
        s = tic;
        X_FFT = FiniteDifference_FFT(f, n);
        t_fft(j) = toc(s);

        fprintf('n = %g\n', n)
        j = j + 1;
    end

    loglog(nn, t_fft,   'LineWidth', 2), hold on
    loglog(nn, t_adi3,  'LineWidth', 2)
    loglog(nn, t_adi6,  'LineWidth', 2)
    loglog(nn, t_adi13, 'LineWidth', 2)
    loglog(nn, 2e-7*nn.^2.*log(nn), 'k--', 'LineWidth', 2), hold off
    legend('FFT',                               ...
           ['ADI, ' char(1013) ' = 10^{ -3}' ], ...
           ['ADI, ' char(1013) ' = 10^{ -6}' ], ...
           ['ADI, ' char(1013) ' = 10^{ -13}'], ...
           'Location', 'NorthWest');
    set(gca, 'FontSize', 16)
    xlim([min(nn) max(nn)])
    ylim([1e-4 10])
end

function SquareSolutionFigure()
    m = 200; n = 200;
    f = chebfun2( @(x,y) -100*x.*sin(20*pi*x.^2.*y).*cos(4*pi*(x+y)) );
    F = coeffs2( f, m, n );
    X = poisson_rectangle( F );
    u = chebfun2( X, 'coeffs' );
    plot(u)
    view(2)
    colormap(colorpalette())
    colorbar('FontSize', 16)
    axis square
    axis off
end

function SquareTimingsFigure()
    nn_lyap = floor(logspace(1,3.7,40));
    t_lyap = zeros(size(nn_lyap));
    j = 1;
    for n = nn_lyap
        fprintf('n = %g\n', n);
        F = ones(n, n);
        fprintf('  lyap: ');
        s = tic;
        X_LYAP = fastPoisson2D_lyap( F, n );
        t_lyap(j) = toc(s);
        fprintf('%g s\n', t_lyap(j));
        j = j + 1;
    end

    nn_adi  = floor(logspace(1,4,40));
    t_adi3  = zeros(size(nn_adi));
    t_adi6  = zeros(size(nn_adi));
    t_adi13 = zeros(size(nn_adi));
    j = 1;
    for n = nn_adi
        fprintf('n = %g\n', n);
        F = ones(n, n);

        s = tic;
        X_ADI3 = poisson_rectangle( F, 1e-3 );
        t_adi3(j) = toc(s);
        fprintf('  ADI 1e-3: %g s\n', t_adi3(j));

        s = tic;
        X_ADI6 = poisson_rectangle( F, 1e-6 );
        t_adi6(j) = toc(s);
        fprintf('  ADI 1e-6: %g s\n', t_adi6(j));

        s = tic;
        X_ADI13 = poisson_rectangle( F, 1e-13 );
        t_adi13(j) = toc(s);
        fprintf('  ADI 1e-13: %g s\n', t_adi13(j));

        j = j + 1;
    end

    loglog(nn_adi,  t_adi3,  'LineWidth', 2), hold on,
    loglog(nn_adi,  t_adi6,  'LineWidth', 2)
    loglog(nn_adi,  t_adi13, 'LineWidth', 2)
    loglog(nn_lyap, t_lyap,  'LineWidth', 2)
    nn = nn_adi(nn_adi > 500);
    loglog(nn, 2.2e-8*nn.^2.*log(nn).^2, 'k--', 'LineWidth', 2)
    loglog(nn, 6.1e-9*nn.^3, 'k--', 'LineWidth', 2), hold off
    legend(['ADI, ' char(1013) ' = 10^{ -3}' ], ...
           ['ADI, ' char(1013) ' = 10^{ -6}' ], ...
           ['ADI, ' char(1013) ' = 10^{ -13}'], ...
           'Bartels-Stewart', 'Location', 'NorthWest')
    set(gca, 'FontSize', 16)
    xlim([min(nn_adi) max(nn_adi)])
    ylim([2e-4 400])
end

function CylinderSolutionFigure()

    function CFS = vals2coeffs( VALS )
    % Convert to Chebyshev--Fourier--Chebyshev coefficients
    [n1, n2, n3] = size( VALS );
    CFS = VALS;
    for k = 1:n3
        CFS(:,:,k) = chebtech2.vals2coeffs( CFS(:,:,k) );
        CFS(:,:,k) = trigtech.vals2coeffs( CFS(:,:,k).' ).';
    end
    for j = 1:n2
        vj = reshape( CFS(:,j,:), n1, n3 );
        vj = chebtech2.vals2coeffs( vj.' ).';
        CFS(:,j,:) = reshape( vj, n1, 1, n3 );
    end
    end

    n1 = 55; n2 = 84; n3 = 42;
    r = chebpts( n1 );
    th = pi*trigpts( n2 );
    z = chebpts( n3 );
    [rr, tt, zz] = ndgrid( r, th, z );
    xx = rr.*cos(tt);
    yy = rr.*sin(tt);
    v = @(x,y,z) (1-x.^2-y.^2).*(1-z.^2).*(z.*cos(4*pi*(x.^2))+cos(4*pi*y.*z));
    f = lap(chebfun3(v));
    F = vals2coeffs(f(xx,yy,zz));
    X = poisson_cylinder( F );
    colormap(colorpalette())
    cylinderplot(X, 'coeffs')
end

function CylinderTimingsFigure()
    nn = floor(logspace(1,2.6,50));
    t  = zeros(size(nn));
    j = 1;
    for n = nn
        fprintf('n = %g\n', n);
        F = ones(n, n, n);
        s = tic;
        X = poisson_cylinder( F );
        t(j) = toc(s);
        fprintf('%g s\n', t(j));
        j = j + 1;
    end

    loglog(nn, t, 'LineWidth', 2), hold on
    nn2 = nn(nn > 110);
    loglog(nn2, 1e-7*nn2.^3.*log(nn2).^2, 'k--', 'LineWidth', 2), hold off
    set(gca, 'FontSize', 16)
    xlim([min(nn) max(nn)])
    ylim([7e-3 500])
end

function SphereSolutionFigure()

    function CFS = vals2coeffs( VALS )
    % Convert to Chebyshev--Fourier--Fourier coefficients
    [n1, n2, n3] = size( VALS );
    CFS = VALS;
    for k = 1:n3
        CFS(:,:,k) = chebtech2.vals2coeffs( CFS(:,:,k) );
        CFS(:,:,k) = trigtech.vals2coeffs( CFS(:,:,k).' ).';
    end
    for j = 1:n2
        vj = reshape( CFS(:,j,:), n1, n3 );
        vj = trigtech.vals2coeffs( vj.' ).';
        CFS(:,j,:) = reshape( vj, n1, 1, n3 );
    end
    end

    n1 = 31; n2 = 32; n3 = 32;
    r = chebpts( n1 );
    lam = pi*trigpts( n2 );
    th = pi*trigpts( n3 );
    [rr, tt, ll] = ndgrid( r, th, lam );
    k = 2;
    rhs = @(r, th, lam) -2*(2*abs(k)+3).*r.^abs(k).*sin(lam).^abs(k).*exp(1i*k*th);
    F = rhs(rr, tt, ll);
    F = vals2coeffs( F );
    X = poisson_solid_sphere( F );
    colormap(colorpalette())
    sphereplot( X )
end

function CubeSolutionFigure()
    m = 10; n = 10; p = 10;
    u = chebfun3( @(x,y,z) (1-x.^2).*(1-y.^2).*(1-z.^2).*cos(x.*y.*z.^2) );
    f = lap( u );
    F = coeffs3( f, m, n, p );
    X = poisson_cube( F );
    colormap(colorpalette())
    cubeplot( X )
end

function CubeTimingsFigure()
    nn = floor(logspace(1,1.7,17));
    t  = zeros(size(nn));
    j = 1;
    for n = nn
        fprintf('n = %g\n', n);
        F = ones(n, n, n);
        s = tic;
        X = poisson_cube( F );
        t(j) = toc(s);
        fprintf('%g s\n', t(j));
        j = j + 1;
    end

    loglog(nn, t, 'LineWidth', 2), hold on
    nn2 = nn(nn > 25);
    loglog(nn2, 7.6e-4*nn2.^3.*log(nn2).^3, 'k--', 'LineWidth', 2), hold off
    set(gca, 'FontSize', 16)
    xlim([min(nn) max(nn)])
    ylim([10 10000])
end
