%Lid Driven Cavity Flow Problem
%MAE 435
%Velocity-Pressure was written: by Victor Ong
%5/19/2020
% clear workspace
clear
clc
close all

% define variables
nx = 10;
ny = nx;

sf=zeros(nx,ny);                % stream function
vort=zeros(nx,ny);              % vorticity
w=zeros(nx,ny);                 

u = zeros(nx,ny+1);             % u-velocity
un = zeros(nx,ny+1);            % new u-velocity
ut = zeros(nx,ny+1);
uu = zeros(nx,ny);

v = zeros(nx+1,ny);             % v-velocity
vn = zeros(nx+1,ny);            % new v-velocity
vt = zeros(nx+1,ny);
vv = zeros(nx,ny);

p = zeros(nx+1,ny+1);
pp = zeros(nx+1,ny+1);

h = 1/(nx-1);
dt = 0.001;

Re = 10;                        % Reynolds number
alpha = 1/Re;
beta = 1.5;

t=0;

% Initial condition
UwallTop = 1;
UwallBot = 0;
error=1;
maxError = 0.03;
maxSteps = 100000;
maxedOut = 0;

% Convergence test selection:
%1=Pressure based, 2=Velocity based, 3=stream function based
testSelect = 1;

it = 0;
step = 0;
while error>maxError && step<maxSteps
    error = 0;
    it = it + 1;
    step = step +1;
    disp(['Iteration: ', int2str(it)]);
    
    %Vorticity-Stream
    w=sf;
    for i=2:nx-1
        for j=2:ny-1
            sf(i,j)=0.25*beta*(sf(i+1,j)+sf(i-1,j)...
                +sf(i,j+1)+sf(i,j-1)+h*h*vort(i,j))+(1.0-beta)*sf(i,j);
        end
    end
    vort(2:nx-1,1)=-2.0*sf(2:nx-1,2)/(h*h);
    vort(2:nx-1,ny)=-2.0*sf(2:nx-1,ny-1)/(h*h)-2.0/h;
    vort(1,2:ny-1)=-2.0*sf(2,2:ny-1)/(h*h);
    vort(nx,2:ny-1)=-2.0*sf(nx-1,2:ny-1)/(h*h);
    for i=2:nx-1
        for j=2:ny-1
            w(i,j)=-0.25*((sf(i,j+1)-sf(i,j-1))*(vort(i+1,j)-vort(i-1,j))...
                -(sf(i+1,j)-sf(i-1,j))*(vort(i,j+1)-vort(i,j-1)))/(h*h)...
                +1/Re*(vort(i+1,j)+vort(i-1,j)+vort(i,j+1)+vort(i,j-1)-4.0*vort(i,j))/(h*h);
        end
    end
    vort(2:nx-1,2:ny-1)=vort(2:nx-1,2:ny-1)+dt*w(2:nx-1,2:ny-1);
    

    % Velocity-Pressure 
        % u-v Boundary conditions
        u(:,1) = 2*UwallBot -u(:,2);            % u bottom wall
        u(:,ny+1) = 2*UwallTop-u(:,ny);         % u top wall
        u(1,:) = 0;                             % u left wall
        u(nx,:) = 0;                            % u right wall
        
        v(:,1) = 0;                             % v bottom wall
        v(:,ny) = 0;                            % v top wall
        v(1,:) = -v(2,:);                       % v left wall
        v(nx+1,:) = -v(nx,:);                   % v right wall

        % solve (advect) x-momentum
        for i=2:nx-1
            for j=2:ny
                u_ip1j = ( 0.5*(u(i+1,j) + u(i,j) ))^2;
                u_ij = ( 0.5*( u(i,j) + u(i-1,j)  ))^2;
                uv_ij = (0.5*(u(i,j+1) + u(i,j) )) * (0.5*(v(i,j) + v(i+1,j)));                
                uv_ijm1 = (0.5*(u(i,j) + u(i,j-1) )) * (0.5*(v(i,j-1) + v(i+1,j-1)));
                
                Ax =  (1/h) * ( u_ip1j -  u_ij + uv_ij - uv_ijm1);
                Dx = (alpha/h^2) * (u(i+1,j) - 2*u(i,j) + u(i-1,j) + u(i,j+1) - 2*u(i,j) + u(i,j-1));
                
                ut(i,j) = u(i,j) + dt*(-Ax + Dx);      
            end
        end
        
        % solve (advect) y-momentum
        for i=2:nx
            for j=2:ny-1
                v_ijp1 = ( 0.5*(v(i,j+1) + v(i,j) ))^2;
                v_ij = ( 0.5*( v(i,j) + v(i,j-1)  ))^2;
                uv_ij = (0.5*(u(i,j+1) + u(i,j) )) * (0.5*(v(i,j) + v(i+1,j)));
                uv_im1j = (0.5*(u(i-1,j+1) + u(i-1,j) )) * (0.5*(v(i,j) + v(i-1,j)));
                
                Ay =  (1/h) * ( uv_ij - uv_im1j + v_ijp1 - v_ij);
                Dy = (alpha/h^2) * (v(i+1,j) - 2*v(i,j) + v(i-1,j) + v(i,j+1) - 2*v(i,j) + v(i,j-1));
                
                vt(i,j) = v(i,j) + dt*(-Ay + Dy);
            end
        end
                        
        % Pressure Boundary treatment
        uuOld = uu;
        vvOld = vv;
        pOld = pp;
        for i=2:nx
           for j=2:ny
               
               if (i>2 && i<nx) && (j>2 && j<ny)
                   mult = 1/4; 
               elseif (j==2 || j==ny) && (i>2 && i<nx)
                   mult = 1/3;
               elseif (i==2 || i==nx) && (j>2 && j<ny)
                   mult = 1/3;
               elseif (i==2 && j==2) || (i==nx && j==2) || (i==2 && j==ny) || (i==nx && j==ny)
                   mult = 1/2;                                                              
               end
          
               p(i,j) = beta*(mult*(p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1))) ...
                   - (mult*h/dt)*(ut(i,j) - ut(i-1,j) + vt(i,j) - vt(i,j-1)) ...
                   + (1-beta)*p(i,j);
               
               % Projection
                un(i,j) = ut(i,j) - (dt/h)*(p(i+1,j) - p(i,j));
                vn(i,j) = vt(i,j) - (dt/h)*(p(i,j+1) - p(i,j));

                % Interpolate
                uu(i,j) = 0.5*(un(i,j)+un(i,j-1));
                vv(i,j) = 0.5*(vn(i,j)+vn(i-1,j));
                pp(i,j) = 0.25*(p(i,j)+p(i-1,j)+p(i,j-1)+p(i-1,j-1));
                %pp(i,j) = 0.5*(p(i,j)+p(i-1,j-1));
                
                % Convergence test
                if testSelect == 1
                    % Convergence test based on pressure
                    test = abs(pp(i,j) - pOld(i,j));
                elseif testSelect == 2
                    % Convergence test based on velocity
                    test = abs(uu(i,j) - uuOld(i,j)) + abs(vv(i,j) - vvOld(i,j));
                elseif testSelect == 3
                    % convergence test based on stream function
                    test = abs(w(i,j)-sf(i,j));
                end               
                error = error + test;                               
           end
        end
  
        % update u and v
        u = un;
        v = vn;
        t=t+dt;
        
        if(step==maxSteps)
            maxedOut = 1;
        end
         
end

if maxedOut==1
    disp ('Iteration maxed out! Either the solution diverges or more iterations are required!')
else

    fh = figure();
    fh.WindowState = 'maximized';
    fontSize = 15;
    
    X = linspace(0,1,nx);
    Y = linspace(0,1,ny);    
    figure(1)
    quiver(X,Y,uu,vv);
    
    subplot(231), contourf(X,Y,rot90(fliplr(uu)),'LineColor','none'), axis('square');
    title('X Velocity Contour')
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar
    set(gca,'FontSize',fontSize)
    
    
    subplot(232), contourf(X,Y,rot90(fliplr(vv)),'LineColor','none'), axis('square');
    title('Y Velocity Contour')
    xlabel('X (m)');
    ylabel('Y (m)');    
    colorbar
    set(gca,'FontSize',fontSize)
    
    Xp = linspace(0,1,nx+1);
    Yp = linspace(0,1,ny+1);
    subplot(233), contourf(Xp,Yp,rot90(fliplr(pp)),'LineColor','none'), axis('square');
    title('Pressure Contour') 
    xlabel('X (m)');
    ylabel('Y (m)');    
    colorbar
    set(gca,'FontSize',fontSize)

    subplot(234), contourf(X,Y,rot90(fliplr(vort)), 10,'LineColor','none'), axis('square'), colorbar;
    title('Vorticity Contour')
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar
    set(gca,'FontSize',fontSize)
    
    subplot(235), contourf(X,Y,rot90(fliplr(sf)), 10,'LineColor','none' ), axis('square'), colorbar;
    title('Stream Function Contour Plot')
    xlabel('X (m)');
    ylabel('Y (m)');
    colorbar
    set(gca,'FontSize',fontSize)
    
    subplot(236)    
    mt=uu(ceil(nx/2),:);
    plot(mt,Y), axis('square');
    xlabel('U velocity (m/s)');
    ylabel('Y (m)');
    title('U vs Y plot');
    colorbar
    set(gca,'FontSize',fontSize)
end

    
    
%     subplot(235)    
%     nt=vv(ceil(ny/2),:);
%     plot(X,-nt), axis('square');
%     xlabel('V velocity (m/s)');
%     ylabel('X (m)');
%     title('V vs X plot');
%     colorbar
    
