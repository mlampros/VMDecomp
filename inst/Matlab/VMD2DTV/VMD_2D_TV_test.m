% 2D Compact Variational Mode Decomposition -- TEST SCRIPT 
% that reproduces figures in [1]
%
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% dominique.zosso@montana.edu
% http://www.math.montana.edu/dzosso
% Initial release 2018-05-10 (c) 2014--
%
% THIS CODE IS PROVIDED AS IS.
%
% Input and Parameters:
% ---------------------
% signal     - the space domain signal (2D) to be decomposed
% alpha      - narrowbandedness of subsignals coefficient (scalar)
% beta       - L1 penalty coefficient of spatial mode support (scalar)
% gamma      - Spatial support TV-term: heat diffusion coefficient
% delta      - Artifact classification threshold (inf for no artifacts)
% rho        - data fidelity coefficient
% rho_k      - u-v splitting coefficient
% tau        - time-step of dual ascent for data ( pick 0 for noise-slack )
% tau_k      - time-step of dual ascent for u-v splitting
% t          - spatial support TV-term: time-step of ODE/PDE
% K          - the number of modes to be recovered
% DC         - true, if the first mode is put and kept at DC (0-freq)
% init       - 0 = all omegas start initialized radially uniformly
%              1 = all omegas start initialized randomly on half plane
%              2 = all omegas initalized by user graphical input
% X_tol      - tolerances of convergence criterion; typically around 1e-7
% N          - maximum number of iterations
% M          - number of submodes
% A_phase    - 2D VMD - 2D-TV-VMD - 2D-TV-VMD-Seg scheduling
%
% When using this code, please do cite our papers:
% -----------------------------------------------
%
% [1] D. Zosso, K. Dragomiretskiy, A.L. Bertozzi, P.S. Weiss, Two-Dimensional 
% Compact Variational Mode Decomposition, Journal of Mathematical Imaging 
% and Vision, 58(2):294–320, 2017.
% DOI:10.1007/s10851-017-0710-z
%
% [2] K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing, 62(3):531-544, 2014.
% DOI:10.1109/TSP.2013.2288675
%
% [3] K. Dragomiretskiy, D. Zosso, Two-Dimensional Variational Mode
% Decomposition, EMMCVPR 2015, Hong Kong, LNCS 8932:197-208, 2015.
% DOI:10.1007/978-3-319-14612-6_15
%

clear all;
close all;
clc;

mkdir('output');


parfor problem = 1:10
    
    
    switch problem
        case 1
            % jerome's texture data
            texture = load('texture.mat');
            
            signal = texture.f;
            alpha = 1000;       % bandwidth constraint
            beta = 0.5;         % L1 penalty on A (compact support)
            gamma = 500;         % TV weight
            delta = inf;        % threshold value for artifact detection
            rho = 10;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 2.5;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 2.5;        % gradient ascent step size for u/v lagrange multipliers
            t = 1.5;             % gradient descent step size for A heat propogation
            
            K = 5;              % number of modes
            M = 1;              % number of submodes (per mode)
            DC = 1;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 1e-10;      % tolerance for u convergence
            A_tol = 1e-4;    % tolerance for A convergence
            omega_tol = 1e-10;  % tolerance for omega convergence
            
            N = 130;             % max number of iterations
            
            A_phase = [100,Inf]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
            
        case 2
            % peptide data
            F = imread('AB_pH_no_scale.png');
            f = double(mean(F,3));
            f = f - mean(mean(f(:)));
            f = conv2(f, fspecial('log', 7,1.25), 'same');
            f = conv2(f, fspecial('log', 7,1.5), 'same');
            
            %             f = downsample(f',2);
            %             f = downsample(f',2);
            signal = f;
            
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = inf;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 6;              % number of modes
            M = 1;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 250;             % max number of iterations
            
            A_phase = [100,150]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
            
        case 3
            % peptide data
            F = imread('H10A.png');
            f = double(mean(F,3));
            f = f - mean(mean(f(:)));
            f = conv2(f, fspecial('log', 7,1.25), 'same');
            f = conv2(f, fspecial('log', 7,1.5), 'same');
            
            %             f = downsample(f',2);
            %             f = downsample(f',2);
            signal = f;
            
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 75;         % TV weight
            delta = 3.5;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 3;              % number of modes
            M = 1;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 200;             % max number of iterations
            
            A_phase = [100,150]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
        case 4
            % cage molecule
            
            f_1 = 40;
            H = 256;
            [X,Y] = meshgrid((1:H)/H, (1:H)/H);
            
            M1 = zeros(H,H);
            M2 = zeros(H,H);
            
            M1(1:H/2, :) = 1;
            M2(H/2+1:H, :) = 1;
            
            A = 0.1;
            B = 1;
            v_1 = (cos(2*pi*f_1/sqrt(1+A^2)*(X+A*Y)) + cos(2*pi*f_1/sqrt(1+A^2)*(A*X-Y)) ).*M1;
            v_2 = (cos(2*pi*f_1/sqrt(1+B^2)*(X+B*Y)) + cos(2*pi*f_1/sqrt(1+B^2)*(B*X-Y)) ).*M2;
            
            
            f = v_1 + v_2;
            f = double(f);
            f = f - mean(f(:));
            
            signal = f;
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = Inf;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 2;              % number of modes
            M = 2;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 200;             % max number of iterations
            
            A_phase = [100,150]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
            
        case 5
            % cage molecule
            %f = imread('synthetic.png');
            f_1 = 50;
            H = 256;
            [X,Y] = meshgrid( (1:H)/H, (1:H)/H );
            
            M1 = zeros(H,H);
            M2 = zeros(H,H);
            
            M1(1:H/2, :) = 1;
            M2(H/2+1:H, :) = 1;
            
            cx=H/2;cy=H/2;ix=H;iy=H;r=80;
            [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
            M3 = double((x.^2+y.^2)<=r^2);
            M3_OP = ones(H,H) - M3;
            
            
            
            theta_1 = 0;
            theta_2 = pi/12;
            theta_3 = pi/4;
            
            v_1 = (cos(2*pi*f_1*(cos(theta_1)*X + sin(theta_1)*Y)) + cos(2*pi*f_1*(cos(theta_1+pi/3)*X + sin(theta_1+pi/3)*Y)) + cos(2*pi*f_1*(cos(theta_1+2*pi/3)*X + sin(theta_1+2*pi/3)*Y)) ).*M1.*M3_OP;
            v_2 = (cos(2*pi*f_1*(cos(theta_2)*X + sin(theta_2)*Y)) + cos(2*pi*f_1*(cos(theta_2+pi/3)*X + sin(theta_2+pi/3)*Y)) + cos(2*pi*f_1*(cos(theta_2+2*pi/3)*X + sin(theta_2+2*pi/3)*Y)) ).*M2.*M3_OP;
            v_3 = (cos(2*pi*f_1*(cos(theta_3)*X + sin(theta_3)*Y)) + cos(2*pi*f_1*(cos(theta_3+pi/3)*X + sin(theta_3+pi/3)*Y)) + cos(2*pi*f_1*(cos(theta_3+2*pi/3)*X + sin(theta_3+2*pi/3)*Y)) ).*M3;
            
            f = v_1 + v_2 + v_3;
            f = double(f);
            f = f - mean(f(:));
            
            signal = f;
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = Inf;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 3;              % number of modes
            M = 3;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 200;             % max number of iterations
            
            A_phase = [100,150]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
        case 6
            % cage molecule
            F = imread('synthetic_eq.png');
            f = double(F);
            f = f - mean(f(:));
            f = conv2(f, fspecial('log', 7,1.25), 'same');
            f = conv2(f, fspecial('log', 7,1.5), 'same');
            
            signal = f;
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = inf;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 1;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 1;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 5;              % number of modes
            M = 3;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 500;             % max number of iterations
            
            A_phase = [10,20]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
        case 7
            % cage molecule
            F = imread('synthetic_eq.png');
            f = double(F);
            f = f - mean(f(:));
            f = conv2(f, fspecial('log', 7,1.25), 'same');
            f = conv2(f, fspecial('log', 7,1.5), 'same');
            
            signal = f;
            alpha = 20000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 500;         % TV weight
            delta = 20;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 5;              % number of modes
            M = 3;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = cat( 3, [0.1908    0.1818    0.1666    0.1425    0.1096; ...
                0.0196    0.0585    0.0951    0.1282    0.1545], ...
                [0.0763    0.0390    0.0001   -0.0391   -0.0761; ...
                0.1739    0.1871    0.1910    0.1872    0.1740], ...
                [-0.1102   -0.1403   -0.1671   -0.1831   -0.1890; ...
                0.1545    0.1276    0.0965    0.0588    0.0202]);           % initialize omegas
            
            u_tol = 10^-8;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 10^-8;  % tolerance for omega convergence
            
            N = 100;             % max number of iterations
            
            A_phase = [10,20]; % A propogation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
        case 8
            
            % Chirps
            
            f_1 = 10;
            f_2 = 5;
            f_3 = 8;
            
            H = 256;
            [X,Y] = meshgrid((1:H)/H, (1:H)/H);
            
            F1 = ( cos(f_1*2*pi*((X+2).^2-(Y+3).^2)) );
            F2 = ( cos(f_2*2*pi*(2/3*(X+0.5).^2-1/3*(Y+1).^2)) );
            F3 = ( cos(f_3*2*pi*(3/3*(X+1).^2+3/3*(Y+2).^2)) );
            
            M1 = zeros(H,H);
            M2 = zeros(H,H);
            M3 = zeros(H,H);
            
            M1(1:H/2, :) = 1;
            M2(H/2+1:H, :) = 1;
            M3(:,H/4:3*H/4) = 1;
            
            f = F1.*M1 + F2.*M2 + F3.*M3;
            
            f = double(f);
            f = f - mean(f(:));
            
            signal = f;
            
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 1000;         % TV weight
            delta = inf;        % threshold value for artifact detection
            rho = 7;           % fidelity weight
            rho_k = 10;        % u/v coupling weight
            tau = 1;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 1;        % gradient ascent step size for u/v lagrange multipliers
            t = 1;             % gradient descent step size for A heat propogation
            
            K = 3;              % number of modes
            M = 1;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 1e-10;      % tolerance for u convergence
            A_tol = 1e-4;    % tolerance for A convergence
            omega_tol = 1e-8;  % tolerance for omega convergence
            
            N = 200;             % max number of iterations
            
            A_phase = [100,Inf]; % A propagation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
            
        case 9
            
            f = imread('Normal.png');
            f = mean(f(:,:,:),3);
            
            f = f - mean(f(:));
            
            signal = f;
            
            alpha = 1500;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = 30;        % threshold value for artifact detection
            rho = 150;           % fidelity weight
            rho_k = 20;        % u/v coupling weight
            tau = 0;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 2;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 2;              % number of modes
            M = 1;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;%cat(3, [-0.11; 0.11], [0.11; 0.11]);           % initialize omegas
            
            
            
            u_tol = 1e-10;      % tolerance for u convergence
            A_tol = 1e-4;    % tolerance for A convergence
            omega_tol = 1e-8;  % tolerance for omega convergence
            
            N = 500;             % max number of iterations
            
            A_phase = [inf,Inf]; % A propagation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
            
        case 10
            F = imread('PFC.jpg');
            f = mean(F(:,:,:),3);
            
            f = f - mean(f(:));
            
            %f = conv2(f, fspecial('log', 7,1), 'same');
            f = conv2(f, fspecial('log', 7,5.0), 'same');
            
            f = downsample(f',2);
            f = downsample(f',2);
            
            signal = f;
            alpha = 2000;       % bandwidth constraint
            beta = 1;         % L1 penalty on A (compact support)
            gamma = 250;         % TV weight
            delta = 2.5;        % threshold value for artifact detection
            rho = 10;           % fidelity weight
            rho_k = 50;        % u/v coupling weight
            tau = 0.1;          % gradient ascent step size for fid. lagrange multipliers
            tau_k = 0.1;        % gradient ascent step size for u/v lagrange multipliers
            t = 2.5;             % gradient descent step size for A heat propogation
            
            K = 4;              % number of modes
            M = 3;              % number of submodes (per mode)
            DC = 0;             % includes DC part (first mode at DC)
            init = 0;           % initialize omegas
            
            u_tol = 1e-9;      % tolerance for u convergence
            A_tol = K*10^-4;    % tolerance for A convergence
            omega_tol = 1e-8;  % tolerance for omega convergence
            
            N = 500;             % max number of iterations
            
            A_phase = [20,35]; % A propagation phases [a,b]
            % iterations 1:a --> 2D VMD
            % iterations a:b --> 2D TV VMD
            % iterations b:end --> 2D TV VMD Segmentation
    end
    
    
    for typenum = 1:4
        
        aphase = [inf inf];
        delta2 = inf;
        
        switch typenum
            case 1
                typename = '2D'; % only run 2D VMD without artifacts
            case 2
                typename = 'TV'; % switch on TV part at some point
                if ~isfinite(A_phase(1)) % TV was never intended
                    continue;
                else
                    aphase(1) = A_phase(1);
                end
            case 3
                typename = 'SEG'; % switch on TV part at some point
                if ~all(isfinite(A_phase)) % TV&Seg was never intended
                    continue;
                else
                    aphase = A_phase;
                end
            case 4
                typename = 'X';
                if ~isfinite(delta) % artifacts was never intended
                    continue;
                else
                    aphase = A_phase; % whatever model
                    delta2 = delta; % but this time with artifacts
                end
        end
        
        
        fname = [ 'output/' num2str(problem) '_' typename ];
        
        
        
        [u, v, omega, A, X] = VMD_2D_TV(signal, alpha, beta, gamma, delta2, rho, rho_k, tau, tau_k, t, K, DC, init, u_tol, omega_tol, A_tol, N, M, aphase );
        
        
        if numel(init) == 1
            inits = init;
        else 
            inits = 3;
        end
        
        fid = fopen([fname '.txt'], 'w');
        fprintf( fid, 'alpha = %g\nbeta = %g\ngamma = %g\ndelta = %g\nrho = %g\nrho_k = %g\ntau=%g\ntau_k = %g\nt = %g\nK = %u\nDC = %u\ninit = %u\nu_tol = %g\nomega_tol = %g\nA_tol = %g\nN = %u\nM = %u\nA_phase = [%u,%u]\n',...
            alpha, beta, gamma, delta2, rho, rho_k, tau, tau_k, t, K, DC, inits, u_tol, omega_tol, A_tol, N, M, aphase(1), aphase(2)   );
        fclose(fid);
        
        
        cnt = '';
        oldsignal = signal; % save this one for later
        
        for c = 1:2
            
            if c==2 % improve contrast by cutting negative image values
                
                signal( signal < 0 ) = 0;
                u( u < 0 ) = 0;
                v( v < 0 ) = 0;
                cnt = '_HI';
            end
            
            borders = zeros( size(signal) );
            
            for k = 1:K
                for m = 1:M
                    tmp = squeeze(u(:,:,k,m));
                    tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                    imwrite( tmp, [fname '_u_' num2str(k) '_' num2str(m) cnt '.png'] );
                    tmp = squeeze(v(:,:,k,m));
                    tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                    imwrite( tmp, [fname '_v_' num2str(k) '_' num2str(m) cnt '.png'] );
                    tmp = squeeze(A(:,:,k).*v(:,:,k,m));
                    tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                    imwrite( tmp, [fname '_Av_' num2str(k) '_' num2str(m) cnt '.png'] );
                    tmp = squeeze(A(:,:,k).*u(:,:,k,m));
                    tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                    imwrite( tmp, [fname '_Au_' num2str(k) '_' num2str(m) cnt '.png'] );
                    tmp = abs(fftshift(fft2(squeeze(u(:,:,k,m)))));
                    tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                    imwrite( tmp, [fname '_uhat_' num2str(k) '_' num2str(m) cnt '.png'] );
                end
                tmp = squeeze(A(:,:,k));
                tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                imwrite( tmp, [fname '_A_' num2str(k) cnt '.png'] );
                
                tmp = signal; % normalize to same scale as input signal
                tmp2 =  A(:,:,k).*sum(v(:,:,k,:),4);
                tmp2 = (tmp2 - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                imwrite( tmp2, [fname '_Av_' num2str(k) cnt '.png'] );
                tmp2 = A(:,:,k).*sum(u(:,:,k,:),4);
                tmp2 = (tmp2 - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                imwrite( tmp2, [fname '_Au_' num2str(k) cnt '.png'] );
                
                borders = borders | bwperim(squeeze(A(:,:,k)));
            end
            tmp = squeeze(X);
            tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp, [fname '_X' cnt '.png'] );
            
            tmp = squeeze(abs(fftshift(fft2(signal))));
            tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp, [fname '_fhat' cnt '.png'] );
            
            tmp = squeeze(signal);
            tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp, [fname '_f' cnt '.png'] );
            
            V = tmp;
            S = X;
            H = 0.5+zeros(size(signal));
            fX = hsv2rgb( cat(3, H,S,V) );
            imwrite( fX, [fname '_fX' cnt '.png'] );
            
            S = borders;
            H = zeros(size(signal));
            fA = hsv2rgb( cat(3, H,S,V) );
            imwrite( fA, [fname '_fA' cnt '.png'] );
            
            S = borders | X;
            H = 0.5*X;
            fAX = hsv2rgb( cat(3, H,S,V) );
            imwrite( fAX, [fname '_fAX' cnt '.png'] );
            
            tmp = signal; % normalize to same scale as input signal
            tmp2 = sum( A.*sum(v,4), 3 );
            tmp2 = (tmp2 - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp2, [fname '_Av' cnt '.png'] );
            tmp2 = sum( A.*sum(u,4), 3 );
            tmp2 = (tmp2 - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp2, [fname '_Au' cnt '.png'] );
            
            tmp = signal - sum( A.*sum(v,4), 3 );
            tmp = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
            imwrite( tmp, [fname '_fAv' cnt '.png'] );
            
            
        end
        
        signal = oldsignal; % restore signal for further computations
    end
    
    disp( [ num2str(problem) ' done.' ] );
end