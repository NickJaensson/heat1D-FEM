close all; clear

% problem parameters
ne = 40;       % number of elements
alpha = 1;      % diffusion coefcient
R = 1;          % length of the domain
deltat = 1e-4;  % time step size
nstep = 50;      % number of time steps
Twall = 1;      % wall temperature
T0    = 0;      % initial temperature

% derived parameters
h = R/ne;       % element size
xnod = 0:h:R;   % nodal coordinates
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

% initialize the element matrix and rhs 
K = zeros(ne+1,ne+1);
f = zeros(ne+1,1);

% define the initial solution of the temperature
sol = zeros(ne+1,1);
sol(1) = Twall;

% plot the solution
figure(1); hold on
plot(0:h:R,sol,'-','color',[0 0 1],'LineWidth',2)

% store the old solution
solold = sol;

% time stepping 
for step = 1:nstep
    
    % update the current time
    time = time + deltat;

    disp(['step = ',num2str(step),' time = ',num2str(time)])
    
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
        
        % initialize local system matrix and RHS
        Ke = 0; fe = 0;
        
        % loop over the integration points
        for k = 1:2
          Ke = Ke + w(k)*Ne(:,k)*Ne(:,k)'*J(k) + ...
                       + deltat*alpha*w(k)*gradNe(:,k)*gradNe(:,k)'*J(k);
          fe = fe + w(k)*Ne(k)*Ne(:,k)'*Told'*J(k);
        end

        % add to the global system matrix
        K(ie:ie+1,ie:ie+1) = K(ie:ie+1,ie:ie+1) + Ke;
        f(ie:ie+1) = f(ie:ie+1) + fe;  
        
    end

    % add boundary condition on the left boundary
    K(1,:) = 0; K(1,1) = 1; f(1) = Twall;
    K(end,:) = 0; K(end,end) = 1; f(end) = 0;

    % solve the linear system
    sol = K\f;
    
    % store the old solution
    solold = sol;

    % plot the solution
    if step == nstep
        plot(0:h:R,sol,'-','color',[step/nstep 0 1],'LineWidth',2)
        plot(0:h:R,erfc((0:h:R)/sqrt(4*alpha*time)),'-k','LineWidth',2)
    end
    
end
    
% make the plot a bit nicer
xlabel('$x$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
ax = gca; 
ax.FontSize = 24;

