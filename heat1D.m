close all; clear

% 1D transient heat problem with an imposed temperature at both boundaries
% comparison to analytical solution (a series until infinity)

% problem parameters
ne = 50;           % number of elements
alpha = 1;         % thermal diffusivity coefcient
ox = 0;            % origin of the domain
L = 1;             % length of the domain
deltat = 1e-3;     % time step size
nstep = 100;       % number of time steps
plotevery = 10;    % plot every ... time step
nplot = 200;       % number of points for plotting analytical solution
nseries = 100;     % number of terms in analytical solution
ana = 1;           % plot analytical sol 1: Cartesian 2: axisymm.
keepK = 1;         % 1: keep elem matrix between time steps  0: do not keep
coorsys = 0;       % 0: Cartesion, 1: axisymmetric
BCs = [0 0];       % BCs(1) for left and BCs(2) right side, with values:
                   % 0: Dirichlet, 1: Neumann, 2: Robin
T0    = 1;         % initial temperature
Twall = [0.1 0.1]; % wall temperature in case BCs = 0
flux = [0 5];      % heat flux / (rho * cp) in case BCs = 1
hheat = [10 5];    % heat transfer coeff. / (rho * cp) in case BCs = 2
Tinf = [0 0];      % temperature at inf. in case BCs = 2

% derived parameters
h = L/ne;       % element size
xnod = ox:h:(ox+L);   % nodal coordinates
time = 0;       % initial time

% the location of the integration points (parametric domain is [-1:1])
xi = [-1/sqrt(3) 1/sqrt(3)];

% the weights in the integration points
w = [1 1];

% the shape functions for each node in each integration point
% Ne(ndf,ninti)
Ne = [1/2*(1-xi); 1/2*(1+xi)];

% derivaties of the shape function wrt to xi in each integration point
% gradxiNe(ndf,ninti)
gradxiNe = [-1/2 -1/2; 1/2 1/2];

% define the initial solution of the temperature
sol = zeros(ne+1,1)+T0;

% plot the solution
figure(1); hold on
plot(ox:h:(ox+L),sol,'-','color',[0 0 1],'LineWidth',2)
xplot = ox:L/nplot:(ox+L);
if any ( ana == [1,2] )
    plot(xplot,0.*xplot+T0,'-k','LineWidth',2)
end

% load the Bessel roots (these were found with chebfun)
if ana == 2
    br = load('31831_Bessel_roots.txt');
end
    
% store the old solution
solold = sol;

% time stepping 
for step = 1:nstep
    
    % update the current time
    time = time + deltat;

    disp(['step = ',num2str(step),' time = ',num2str(time)])
    
    % initialize the element matrix and rhs 
    if ( step == 1 && keepK ) || ~keepK
        K = zeros(ne+1,ne+1);
    end
    f = zeros(ne+1,1);

    % loop over the elements and assemble the system matrix
    for ie = 1:ne
        
        % get the current nodal coordinates
        re = [xnod(ie) xnod(ie+1)]';
        
        % determine the Jacobian in each integration point
        J = transpose(gradxiNe'*re);

        % derivaties of the shape function wrt to x in each integration 
        % point gradNe(ndf,ninti)
        gradNe = gradxiNe./J;
        
        % the coordinates in the integration points
        r = transpose(Ne'*re);

        % integration in a cylindrical coordinate system
        if coorsys == 1 
            J = J*2*pi.*r;
        end
        
        % get the previous solution in the nodes
        Told = [solold(ie) solold(ie+1)];
        
        % system matrix only has to be assembled once (does not change)
        if ( step == 1 && keepK ) || ~keepK
            
            % initialize local system matrix
            Ke = 0; 

            % loop over the integration points
            for k = 1:2
              Ke = Ke + (1/deltat)*w(k)*Ne(:,k)*Ne(:,k)'*J(k) + ...
                           + alpha*w(k)*gradNe(:,k)*gradNe(:,k)'*J(k);
            end

            % add to the global system matrix
            K(ie:ie+1,ie:ie+1) = K(ie:ie+1,ie:ie+1) + Ke;
        
        end
        
        % initialize local RHS
        fe = 0;
        
        % loop over the integration points
        for k = 1:2
          fe = fe + (1/deltat)*w(k)*Ne(:,k)*Ne(:,k)'*Told'*J(k);
        end       
        
        % add to the total RHS vector
        f(ie:ie+1) = f(ie:ie+1) + fe;  
        
    end
    
    % get factor depending on coordinate system
    facL = 1; facR = 1;
    if coorsys == 1
        facL = 2*pi*ox; % left boundary
        facR = 2*pi*(ox+L); % right boundary
    end
        
    % add Neumann boundary conditions
    if ( BCs(1) == 1 )            
        f(1) = f(1) - facL*flux(1);
    end

    if ( BCs(2) == 1 )    
        f(end) = f(end) - facR*flux(2);
    end
    
    % add Robin boundary conditions
    if ( BCs(1) == 2 )   
        if ( step == 1 && keepK ) || ~keepK
            K(1,1) = K(1,1) + facL*hheat(1);
        end
        f(1) = f(1) + facL*hheat(1)*Tinf(1);
    end

    if ( BCs(2) == 2 ) 
        if ( step == 1 && keepK ) || ~keepK
            K(end,end) = K(end,end) + facR*hheat(2);
        end        
        f(end) = f(end) + facR*hheat(2)*Tinf(2);
    end    
            
    % add essential boundary conditions
    % NOTE: see userguide of eztfem for details on this approach
    if any ( BCs == 0 )
        
        % make a system vector with the essential BCs
        uess = zeros(ne+1,1);    
        
        % indices of essential BCs and indices of "free" dofs
        if BCs(1) == 0 && BCs(2) ~= 0
          iess = [1]; 
          iu = 2:ne+1;
          uess(iess(1)) = Twall(1);
        elseif BCs(1) ~= 0 && BCs(2) == 0
          iess = [ne+1]; 
          iu = 1:ne;
          uess(iess(1)) = Twall(2);
        else
          iess = [1 ne+1]; 
          iu = 2:ne;
          uess([iess(1),iess(2)]) = Twall;
        end        


        if ( step == 1 && keepK ) || ~keepK

          % extract Kup
          Kup = K(iu,iess);

          % modify K
          K(iu,iess) = 0;                     % Kup = 0
          K(iess,iu) = 0;                     % Kpu = 0
          K(iess,iess) = eye(length(iess));   % Kpp = I

        end

        % modify f
        f(iu) = f(iu) - Kup * uess(iess) ; % fu = fu - Kup * up
        f(iess) = uess(iess) ;             % fp = up
    
    end

    % solve the linear system
    sol = K\f;
    
    % store the old solution
    solold = sol;

    % postprocessing
    if mod(step,plotevery) == 0
    
        % the numerical solution
        plot(ox:h:(ox+L),sol,'-','color',[step/nstep 0 1],'LineWidth',2)
        
        if ana == 1

            % calculate coefficients of the anlytical solution once
            % NOTE: assume Twall = 0 and T0 = 1 here, and compensate for the 
            % difference later
            if exist('Bn','var') == 0
               Bn = zeros(nseries,1);
               for n = 1:nseries
                  Bn(n) = -2*(-1+(-1)^n)/(n*pi);
               end
            end

            % calculate the analytical solution
            solex = zeros(nplot+1,1);
            for n = 1:nseries
               solex = solex + (Bn(n)*sin(n*pi*xplot/L)*exp(-(n*pi/L)^2*alpha*time))';
            end

            % plot the analytical solution using real values for T0 and Twall
            plot(xplot,Twall+solex*(T0-Twall),'-k','LineWidth',2)
            
        elseif ana == 2
            
            % only calculate Bessel functions once
            if exist('rr','var') == 0
               rr = zeros(nseries,nplot+1);
               brr = zeros(nseries,1);
               for n = 1:nseries
                   rr(n,:) = besselj(0,br(n)*xplot);
                   brr(n) = br(n)*besselj(1,br(n));
               end
            end

            % calculate the analytical solution
            solex = zeros(nplot+1,1);
            for n = 1:nseries
               solex = solex + 2*exp(-br(n)^2*time)*rr(n,:)' ./ brr(n);
            end
            
            % plot the analytical solution using real values for T0 and Twall
            plot(xplot,Twall+solex*(T0-Twall),'-k','LineWidth',2)
            
        end
        
    end
    
end
   
% make the plot a bit nicer
xlabel('$x$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
ax = gca; 
ax.FontSize = 24;
xlim([ox ox+L])
%ylim([0 1]);

disp(['outward flux at left boundary = ',num2str((sol(2) - sol(1)) / h)])
if BCs(1) == 1
    disp(['which should be = ',num2str(flux(1))])
elseif BCs(1) == 2
    disp(['which should be = ',num2str(hheat(1)*(sol(1)-Tinf(1)))])
end

disp(['outward flux at right boundary = ',num2str(-(sol(end) - sol(end-1)) / h)])
if BCs(2) == 1
    disp(['which should be = ',num2str(flux(2))])
elseif BCs(2) == 2
    disp(['which should be = ',num2str(hheat(2)*(sol(end)-Tinf(2)))])
end
