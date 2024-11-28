% sch_2d_adi: Solves 2D Schrödinger equation using ADI scheme.
% 
% Inputs:
%
%   tmax:   Maximum integration time
%   level:  Discretization level
%   lambda: dt/dx
%   idtype: Selects initial data type
%   idpar:  Vector of initial data parameters
%   vtype:  Selects potential type
%   vpar:   Vector of potential parameters
%
% Outputs:
%
%   x:      Column vector of x coordinates [nx]
%   y:      Column vector of y coordinates [ny]
%   t:      Column vector of t coordinates [nt]
%   psi:    Array of computed psi values [nt x nx x ny]
%   psire:  Array of computed psi_re values [nt x nx x ny]
%   psiim:  Array of computed psi_im values [nt x nx x ny]
%   psimod: Array of computed sqrt(psi psi*) values [nt x nx x ny]
%   v:      Array of potential values [nx x ny]
function [x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)

    % Define mesh and derived parameters
    nx = 2^level + 1;               ny = nx;
    x  = linspace(0.0, 1.0, nx);    y  = x;
    dx = x(2) - x(1);               dy = dx;
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    t  = (0 : nt-1) * dt;

    % Define meshgrid for populating psi(x,y,0) and V(x,y)
    [X, Y] = meshgrid(x, y);

    % Initialize solution, and set initial data
    psi = zeros(nt, nx, ny);
    if idtype == 0
        % Exact family 
        psi(1, :, :) = sin(idpar(1)*pi*x)' * sin(idpar(2)*pi*y);
    elseif idtype == 1
        % Boosted Gaussian
        % Create variable names for function parameters 
        x0      = idpar(1);      y0 = idpar(2);    
        delta_x = idpar(3); delta_y = idpar(4); 
        p_x     = idpar(5);     p_y = idpar(6);   

        % Calculate psi(x, y, 0)
        psi_0 = exp(1i*p_x*X) .* exp(1i*p_y*Y) ...
             .* exp(-(((X - x0).^2)/delta_x^2 + ((Y - y0).^2)/delta_y^2));
        psi(1, :, :) = reshape(psi_0, [1, nx, ny]);
    else
        fprintf('sch_2d_adi: Invalid idtype=%d\n', idtype);
        return
    end
    % Set boundary conditions of initial data to zero
    % t = 0:   ψ(0,y,t) = ψ(1,y,t) = ψ(x,0,t) = ψ(x,1,t) = 0
    psi(1, 1, :)  = 0; 
    psi(1, :, 1)  = 0; 
    psi(1, nx, :) = 0; 
    psi(1, :, ny) = 0; 

    % Initialize time-independent potential
    v = zeros(nx,ny);
    if vtype == 0
        % No potential - leave unchanged
    elseif vtype == 1
        % Rectangular barrier or well 
        % Create variable names for function parameters 
        x_min = vpar(1);   x_max = vpar(2);    
        y_min = vpar(3);   y_max = vpar(4); 
        Vc    = vpar(5);
        
        % Calculate V(x, y)
        v((X >= x_min & X <= x_max) & (Y >= y_min & Y <= y_max)) = Vc;
    elseif vtype == 2
        % Double slit
        % Create variable names for function parameters 
        x1 = vpar(1);   x2 = vpar(2);    
        x3 = vpar(3);   x4 = vpar(4); 
        Vc = vpar(5); 
        j_prime = (ny - 1)/4 + 1;

        % Calculate V(x, y)
        Vc_indices = (X <= x1) | (X >= x2 & X <= x3) | (X >= x4);
        v(Vc_indices, j_prime) = Vc;
        v(Vc_indices, j_prime + 1) = Vc;
    else
        fprintf('sch_2d_adi: Invalid vtype=%d\n', vtype);
        return
    end

    % Define sparse matrix diagonals for first ADI eqn
    dl = (-1i*dt/(2*dx^2)) * ones(nx, 1);
    d  = (1 + 1i*dt/(dx^2)) * ones(nx, 1);
    du = dl;
    % Impose boundary conditions
    d(1)     = 1.0;
    du(2)    = 0.0;
    dl(nx-1) = 0.0;
    d(nx)    = 1.0;
    % Compute sparse matrix for first ADI eqn
    A_half = spdiags([dl d du], -1:1, nx, nx);

    % Loop that iterates each time step 
    for n = 1:nt-1
        % reshape ψ to create a 2d matrix at this timestep
        psi_n = reshape(psi(n,:,:), nx, ny);
        % Create matrix for ψ^(n+1/2)
        psi_half = zeros(nx,ny);

        % Solve tridiagonal system for each j (row)
        for j = 2:ny-1
            % Array for holding the RHS of the first ADI eqn
            f = zeros(nx,1);

            % Compute RHS of first ADI eqn in stages. 
            % Ensure boundary conditions are maintained by using 2:nx-1
            f(2:nx-1) = (1i*dt/(2*dy^2)) * (psi_n(2:nx-1, j+1) + psi_n(2:nx-1, j-1)) ...
                      + (1 - 1i*dt*(1/dy^2 + v(2:nx-1,j)/2)) .* psi_n(2:nx-1, j);
            f(2:nx-1) = (1i*dt/(2*dx^2)) * (f(1:nx-2) + f(3:nx)) + ...
                        (1 - 1i*dt/dx^2) * f(2:nx-1);

            % Solve first ADI system
            psi_half(:,j) = A_half \ f;
            % Impose boundary conditions
            psi_half(1,:)  = 0.0;
            psi_half(:,1)  = 0.0;
            psi_half(nx,:) = 0.0;
            psi_half(:,ny) = 0.0;
        end

        % Define upper and lower sparse matrix diagonals for second ADI eqn
        dl = (-1i*dt/(2*dy^2)) * ones(ny, 1);
        du = dl;
        % Impose boundary conditions
        du(2)    = 0.0;
        dl(nx-1) = 0.0;

        % Solve tridiagonal system for each i (column)
        for i = 2:nx-1
            % Define middle sparse matrix diagonal for second ADI eqn
            v_i = reshape(v(i,:), ny, 1);
            d = 1 + 1i*dt/dy^2 + (1i*dt/2)*v_i;
            % Impose boundary conditions 
            d(1)  = 1.0;
            d(nx) = 1.0;

            % Compute sparse matrix for second ADI eqn
            A_full = spdiags([dl d du], -1:1, ny, ny);

            % Compute RHS of second ADI eqn 
            f = reshape(psi_half(i,:), ny, 1);

            % Solve second ADI system
            psi(n+1, i, :) = A_full \ f;
            % Impose boundary conditions 
            psi(n+1, 1, :)  = 0.0;
            psi(n+1, :, 1)  = 0.0;
            psi(n+1, nx, :) = 0.0;
            psi(n+1, :, ny) = 0.0;
        end
    end 

    % Compute real, imaginary, and modulus of each entry in psi
    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);

    % Convert to column vectors
    x = x.';
    y = y.';
    t = t.';
end