close all; clear

% 1D transient heat problem with an imposed temperature at both boundaries
% comparison to analytical solution (a series until infinity)

% problem parameters
ne = 50;          % number of elements
alpha = 1;        % diffusion coefcient
L = 0.1;          % length of the domain
deltat = 1e-5;    % time step size
nstep = 100;      % number of time steps
plotevery = 10;   % plot every ... time step
Twall = 0.1;      % wall temperature
T0    = 1;        % initial temperature
nplot = 200;      % number of points for plotting analytical solution
nseries = 100;    % number of terms in analytical solution
keepK = 1;        % 1: keep elem matrix between time steps  0: do not keep

% derived parameters
h = L/ne;       % element size
xnod = 0:h:L;   % nodal coordinates
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
sol(1) = Twall;
sol(end) = Twall;

% plot the solution
figure(1); hold on
plot(0:h:L,sol,'-','color',[0 0 1],'LineWidth',2)
xplot = 0:L/nplot:L;
plot(xplot,0.*xplot+T0,'-k','LineWidth',2)

% store the old solution
solold = sol;

% time stepping 
for step = 1:nstep
    
    % update the current time
    time = time + deltat;

    disp(['step = ',num2str(step),' time = ',num2str(time)])
    
    % initialize the element matrix and rhs 
    if step == 1 && keepK
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

        % get the previous solution in the nodes
        Told = [solold(ie) solold(ie+1)];
        
        % system matrix only has to be assembled once (does not change)
        if step == 1 && keepK
            
            % initialize local system matrix
            Ke = 0; 

            % loop over the integration points
            for k = 1:2
              Ke = Ke + w(k)*Ne(:,k)*Ne(:,k)'*J(k) + ...
                           + deltat*alpha*w(k)*gradNe(:,k)*gradNe(:,k)'*J(k);
            end

            % add to the global system matrix
            K(ie:ie+1,ie:ie+1) = K(ie:ie+1,ie:ie+1) + Ke;
        
        end
        
        % initialize local RHS
        fe = 0;
        
        % loop over the integration points
        for k = 1:2
          fe = fe + w(k)*Ne(:,k)*Ne(:,k)'*Told'*J(k);
        end       
        
        f(ie:ie+1) = f(ie:ie+1) + fe;  
        
    end

    % add boundary condition on the left and right boundary
    % NOTE: see userguide of eztfem for details on this approach

    % indices of essential BCs
    iess = [1 ne+1];

    % indices of "free" dofs
    iu = 2:ne;

    % make a system vector with the essential BCs
    uess = zeros(ne+1,1);
    uess(iess) = Twall;

    if step == 1 && keepK

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

    % solve the linear system
    sol = K\f;
    
    % store the old solution
    solold = sol;

    % postprocessing
    if mod(step,plotevery) == 0
    
        % the numerical solution
        plot(0:h:L,sol,'-','color',[step/nstep 0 1],'LineWidth',2)
        
        % calculate coefficients of the anlytical solution
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

    end
    
end
   
% make the plot a bit nicer
xlabel('$x$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
ax = gca; 
ax.FontSize = 24;
ylim([min([Twall T0]) max([Twall T0])]);