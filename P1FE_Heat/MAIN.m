%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conforming FE P1 code for  u_t = Delta u + f with  Neuman BC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long;
tic
%% Final times
T=1;
itermax=1000;

123

%% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_4.mat'};%'mesh1_05.mat'};%'mesh1_6.mat'};%'mesh1_7.mat'};
nbmeshes=size(meshes,1);

%% Errors for each meash initiation
MAXL2error = zeros(nbmeshes, 1);
L1W11error = zeros(nbmeshes, 1);

%% Order of convergence initiation
ocMAXL2error = zeros(nbmeshes - 1, 1);
ocL1W11error = zeros(nbmeshes - 1, 1);

% To see the results printed in file
fid = fopen('results.txt','w');

%% Loop over each mesh in the sequence
for imesh=1:nbmeshes

    %% Load mesh here!
    loadmesh=strcat('load matlab_meshes/', meshes{imesh}); 
    str = sprintf('\n %s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    disp('mesh loaded');

    %% Compute real centers of mass and mesh size
    cg=gravity_centers(ncell, cell_v, vertex, area);
    h(imesh)=max(abs(diam));
    %% Time steps
    Ndt(imesh) = ceil(T/h(imesh)^2);
    dt = T/Ndt(imesh);

    str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
    forkprint(fid,str);
 
    %% Assemble Mass and Stifness (diffusion) Matrices
    A = assemble_diffusion_system(cell_v, ncell, nvert, vertex);
    M = assemble_mass_system(area, ncell, nvert, cell_v);

    %% Initial condition, interpolation of initial condition, and error
    U_pre = test_cases(0,vertex)';
    [L2error_ini, W11error_ini] = compute_norms(U_pre - test_cases(0,vertex)', area, vertex, cell_v);

    %% Print the scheme solution and exact solution at initial condition to view in Paraview
    write_solution_vtk(U_pre, strcat('VTKout/solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    write_solution_vtk(U_pre, strcat('VTKout/exact_solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    
    %% Left hand side matrices of the scheme
    G = M + dt*A; 

    %% Error norms initiation
    L2error = zeros(Ndt(imesh), 1);
    W11error = zeros(Ndt(imesh), 1);

    %% Time stepping starts here!
    for idt = 1 : Ndt(imesh);
        b = assemble_source(idt * dt, cell_v, ncell, nvert, area, cg);
        rhs = M * U_pre + dt * b;
        U = G\rhs; 
        U_pre = U; % Current time scheme solution is retiring

        %% Errors
        [L2error(idt), W11error(idt)] = compute_norms(U - test_cases(idt * dt,vertex)', area, vertex, cell_v);

        %% Create files for visualising the approximate and exact solution 
        write_solution_vtk(U,strcat('VTKout/p1_c_heat_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the scheme solution at current time in Paraview
        write_solution_vtk(test_cases(idt*dt, vertex)',strcat('VTKout/exact_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the exact solution at current time in Paraview

    end 
    
    %% Time norms: L infinity, L1 norm, and L2 norm
    MAXL2error(imesh) = max(L2error);
    L1W11error(imesh) = dt * sum(W11error);  

    str = sprintf('Mesh %i. MAXL2error=%4.2e, L1W11error:%4.2e',imesh, MAXL2error(imesh), L1W11error(imesh));
    forkprint(fid,str);

    Time(imesh)=toc;
end 
str = sprintf(['\n Elapsed time is ' num2str(toc) 'seconds']);
forkprint(fid,str); 

%% Orders of convergence
ocL2error=zeros(nbmeshes-1,1);
ocH10error=zeros(nbmeshes-1,1);
for imesh = 1:nbmeshes-1
    ocMAXL2error(imesh)=log(MAXL2error(imesh)/MAXL2error(imesh+1)) / log(h(imesh)/h(imesh+1));
    ocL1W11error(imesh)=log(L1W11error(imesh)/L1W11error(imesh+1)) / log(h(imesh)/h(imesh+1));
end 

str = sprintf('\n MaxL1 convergence rate:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',MAXL2error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',MAXL2error(imesh),ocMAXL2error(imesh-1));
        forkprint(fid,str);
    end
end

str = sprintf('\n L1W11 convergence rate:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
    if (imesh==1)
        str = sprintf('\t%4.2e\n',L1W11error(imesh));
        forkprint(fid,str);
    else
        str = sprintf('\t%4.2e \t %4.2e\n',L1W11error(imesh), ocL1W11error(imesh-1));
        forkprint(fid,str);
    end
end

fclose(fid);

% Write data file
fid = fopen('data_rates.dat','w');
fprintf(fid,'meshsize timestep MaxL2error L1W11error Time \n');
for i=1:nbmeshes
    fprintf(fid,'%f %f %f %f %f\n',h(i),T/Ndt(i),MAXL2error(i),L1W11error(i),Time(i));
end;
fclose(fid);

%% Plot of the MAXL2 error of u and L1W11 error Gradu with respect of h in a logarithmic scale
subplot(2,1,1);
loglog(h,MAXL2error,'*-');
axis([0.01 1 1e-5 1]);
hold on
loglog(h,h.^2); % I compute the vector in order to have a slope equal to 2
title('Errors in $L^{\infty}L^2$ norm on the functions','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
ylabel('$\left\| u- \overline{u} \right\|_{L^{\infty}(L^2)}$','Interpreter','latex')
legend('L^{\infty}L^2 error','Slope equal to 2','Location','southeast');
grid on
subplot(2,1,2);
loglog(h,L1W11error,'*-');
axis([0.01 1 0.0001 1]);
hold on
loglog(h,h);
title('Errors in $L^1L^{1}$ norm on the gradients','Interpreter','latex')
xlabel('$h$','Interpreter','latex')
ylabel('$\left\| \nabla u-\nabla \overline{u} \right\|_{L^1(L^1)}$','Interpreter','latex')
legend('L^1W^{1,1} error','Slope equal to 1','Location','southeast')
grid on

