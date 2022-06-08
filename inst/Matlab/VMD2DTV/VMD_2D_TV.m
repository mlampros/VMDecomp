function [u, v, omega, A, X] = VMD_2D_TV(signal, alpha, beta, gamma, delta, rho, rho_k, tau, tau_k, t, K, DC, init, u_tol, omega_tol, A_tol, N, M, A_phase )
% 2D Compact Variational Mode Decomposition -- core function
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


    


%% Variable Declaration
% Resolution of image
[Hy,Hx] = size(signal);
[X,Y] = meshgrid((1:Hx)/Hx, (1:Hy)/Hy);

% Spectral Domain discretization
fx = 1/Hx;
fy = 1/Hy;
freqs_1 = X - 0.5 - fx;
freqs_2 = Y - 0.5 - fy;

% N iterations at most, 2 spatial coordinates, K modes, M submodes
omega = zeros(N, 2, K, M);

% Storage matrices for (Fourier) modes. Iterations are not recorded.
u_hat = zeros(Hy,Hx,K,M);

u = u_hat;
u_old = u;
v = u;

% Augmented Lagrangian variables

lambda_k = u;   %linking variable u/v   ~rho_k
lambda = zeros(Hy,Hx);      %data fidelity          ~rho

% Spatial support variables
A = ones(Hy,Hx,K);
A_old = ones(Hy,Hx,K);

% Artifact map
X = zeros(Hy,Hx);

%% Initialization of omega_k

switch numel(init)
    case 1   % use built-in omega initialization
        switch init
            case 0    % Case 0: Radially Uniform
                % if DC, keep first mode at 0,0
                if DC
                    maxK = K-1;
                else
                    maxK = K;
                end
                radius = 0.3;
                for k = DC + (1:maxK)
                    for m = 1:M
                        omega(1,1,k,m) = radius*cos(pi*(k-1+(m-1)*maxK)/maxK/M);
                        omega(1,2,k,m) = radius*sin(pi*(k-1+(m-1)*maxK)/maxK/M);
                    end
                end
            case 1   % Case 1: Random on Half-Plane
                for k=1:K
                    for m = 1:M
                        omega(1,1,k,m) = rand()-1/2;
                        omega(1,2,k,m) = rand()/2;
                    end
                end
                % DC component (if expected)
                if DC
                    omega(1,1,1,:) = 0;
                    omega(1,2,1,:) = 0;
                end
            case 2      % Case 2: User Input Graphically
                figure;
                imagesc(abs(fftshift(fft2(signal))));
                [x,y] = ginput(K*M);
                x = x/Hx - 0.5;
                y = y/Hy - 0.5;
                
                % If DC = 1, first selection is overwritten to origin
                if DC
                    maxK = K-1;
                    omega(1,1,1,:) = 0;
                    omega(1,2,1,:) = 0;
                else
                    maxK = K;
                end
                
                for k = DC + (1:maxK)
                    for m=1:M
                        omega(1,1,k,m) = x((k-1)*M + m);
                        omega(1,2,k,m) = y((k-1)*M + m);
                    end
                end
            otherwise
                error('init parameter has inappropriate value');
                
        end
    case 2*K*M % use given omega list for initialization; should be 2xKxM
        if size(init,1) ~= 2 || size(init,2) ~= K
            error('init parameter has inappropriate size');
        end
        omega(1, :,:,:) = init;
    otherwise
        error('init parameter has inappropriate size');
end

%% Main loop for iterative updates

% Stopping criteria tolerances
uDiff=inf;
ADiff=inf;
omegaDiff = inf;

% Stores the sum of A_j*v_j for all j not equal to k
sum_Avk = 0;


% Main loop - run until convergence or max number of iterations

n = 1;
while((uDiff > u_tol || ADiff > A_tol || omegaDiff > omega_tol) && n < N) || n <= max(isfinite(A_phase).*A_phase)
    
    % Modes
    for k=1:K
        % Submodes
        for m=1:M
            
            % Compute the halfplane spectral mask for the 2D "analytic signal"
            HilbertMask = ( sign(freqs_1*omega(n,1,k,m) + freqs_2*omega(n,2,k,m) ) + 1 );
            
            % Update accumulator
            if m == 1
                if k == 1
                    sum_Avk = sum_Avk + A(:,:,end).*v(:,:,end,end) - A(:,:,1).*v(:,:,1,1);
                else
                    sum_Avk = sum_Avk + A(:,:,k-1).*v(:,:,k-1,end) - A(:,:,k).*v(:,:,k,1);
                end
            else
                sum_Avk = sum_Avk + A(:,:,k).*v(:,:,k,m-1) - A(:,:,k).*v(:,:,k,m);
            end
            
            
            % Update v (time domain averaging)
            v(:,:,k,m) = ( rho_k*u(:,:,k,m) + lambda_k(:,:,k,m) + rho*A(:,:,k).*( signal - sum_Avk + lambda/rho ).*(1-X) ) ./ ( rho_k + rho*(1-X).*A(:,:,k).^2 );
            
            % Update u_hat (analytic signal spectrum via Wiener filter)
            u_hat(:,:,k,m) = (fftshift(fft2(rho_k*v(:,:,k,m) - lambda_k(:,:,k,m))).*HilbertMask) ./ ( rho_k + 2*alpha*( (freqs_1 - omega(n,1,k,m)).^2 + (freqs_2 - omega(n,2,k,m) ).^2) );
            
            % Update center frequencies (first mode is kept at omega = 0 if DC = 1)
            
            if ~DC || k > 1
                
                % update signal frequencies as center of mass of power spectrum
                omega(n+1,1,k,m) = sum(sum(freqs_1.*(abs(u_hat(:,:,k,m)).^2)))/sum(sum(abs(u_hat(:,:,k,m)).^2));
                omega(n+1,2,k,m) = sum(sum(freqs_2.*(abs(u_hat(:,:,k,m)).^2)))/sum(sum(abs(u_hat(:,:,k,m)).^2));
                
                % keep omegas on same halfplane (top half)
                if omega(n+1,2,k,m) < 0
                    omega(n+1,:,k,m) = -omega(n+1,:,k,m);
                end
            end
            
            % recover full spectrum (and signal) from analytic signal spectrum
            u(:,:,k,m) = real(ifft2(ifftshift(squeeze(u_hat(:,:,k,m)))));
            
        end
        
        % No MBO/TV-term propagation in phase I (2D VMD, only)
        
        % Individual MBO for TV-term Propagation in phase II (2D TV VMD)
        if(n >= A_phase(1)  && n < A_phase(2))
            % Reconstruction Fidelity + Area Penalty + Segmentation Penalty
            A(:,:,k) = A(:,:,k) + t*( -beta + 2*rho*sum(v(:,:,k,:),4).*( signal - sum(A.*sum(v,4),3) + A(:,:,k).*sum(v(:,:,k,:),4) + lambda/rho ).*(1-X) );
            A(:,:,k) = A(:,:,k)./(1 + t*2*rho*(1-X).*sum(v(:,:,k,:),4).^2);
            
            % Project to characteristic range
            A( A > 1 ) = 1;
            A( A < 0 ) = 0;
            
            % Propagate by heat equation
            A(:,:,k) = ifft2(  fft2( A(:,:,k) )./ (1 + t*gamma*ifftshift(freqs_1.^2 + freqs_2.^2)) );
            
            % individual MBO thresholding [0,1] (no segmentation constraint)
            
            A(:,:,k) = (A(:,:,k) >= 0.5);
        end
        
    end
    
    % Joint MBO prop. with segmentation, "Winner Takes All", phase III
    
    if(n >= A_phase(2))
        sumAv = sum(A.*sum(v,4),3);
        for k = 1:K
            A(:,:,k) = A(:,:,k) + t*( -beta + 2*rho*sum(v(:,:,k,:),4).*( signal - sumAv + A(:,:,k).*sum(v(:,:,k,:),4) + lambda/rho ) );
            A(:,:,k) = A(:,:,k)./(1 + t*2*rho*sum(v(:,:,k,:),4).^2);
            A(:,:,k) = ifft2(  fft2( A(:,:,k) )./ (1 + t*gamma*ifftshift(freqs_1.^2 + freqs_2.^2)) );
        end
        A = reshape(A, Hx*Hy, K);
        [~,I] = max(A,[],2);
        A = reshape(full(sparse(1:Hx*Hy, I, 1, Hx*Hy, K)), Hx, Hy, K);
        
    end
    
    % Artifact thresholding
    DF = (signal - sum(A.*sum(v,4),3));
    X = DF.^2 >= delta;
    
    
    % data fidelity dual ascent
    lambda = lambda + tau*DF;
    % update Lagrangian multiplier variables via gradient ascent
    lambda_k = lambda_k + tau_k*( u - v );
    
    
    % Update counter
    n = n+1;
    
    % Tolerance calculation for stopping criteria
    uDiff = norm(u(:)-u_old(:)).^2 / norm(u(:)).^2 / (Hx*Hy);
    ADiff = norm(A(:)-A_old(:),1) / (Hx*Hy);
    
    omegaDiff = norm(omega(n,:) - omega(n-1,:)).^2;
    
    % Storage of n-th iteration
    u_old = u;
    A_old = A;
    
    
    % Debug: Display subsequent relative differences + Graphs
%     
%     disp(uDiff);
%     disp(ADiff);
%     disp(omegaDiff);
%     disp(n);
%     
%     if (mod(n,10)==0)
%         
%         for k = 1:K
%             
%             
%             
%             for m = 1:M
%                 
%                 % Plot modes u or v
%                 subplot(M+2,K,k + (m-1)*K);
%                 imagesc(u(:,:,k,m));
%                 colormap jet;
%                 axis equal;
%                 axis off;
%                 
%                 % Plot supports A
%                 subplot(M+2,K,k + (M)*K);
%                 imagesc(A(:,:,k));
%                 colormap jet;
%                 axis equal;
%                 axis off;
%                 
%                 subplot(M+2,K,k + (M+1)*K);
%                 imagesc(X);
%                 
%                 
%                 colormap jet;
%                 axis equal;
%                 axis off;
%                 
%                 drawnow;
%                 
%             end
%         end
%     end
    
end

omega = omega(n, :,:,:);