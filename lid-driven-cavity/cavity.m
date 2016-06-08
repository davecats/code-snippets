clear all
close all

nx=60;
ny=60;

t=0;
dt=0.001;
tmax=10;
Re=100;
alpha=0.6;
tol=1e-4;
maxiter=500;

ni=1/Re;
Un = 1; Us = 0; Uw = 0; Ue=0;
Vn = 0; Vs = 0; Vw = 0; Ve=0;

x = linspace(0,1,nx);
y = linspace(0,1,ny);
[X,Y] = meshgrid(x,y);
dx=x(2)-x(1); dy=y(2)-y(1);

U=zeros(ny,nx+1);
V=zeros(ny+1,nx);
P=zeros(ny,nx);   b=P;

% FD Matrix for solution of pressure correction
A=zeros(nx*ny,nx*ny);
eq=@(i,j) j+(i-1)*ny;
for i=1:nx; for j=1:ny
            A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j))+2*dt*(dy/dx+dx/dy);
            if i<nx; A(eq(i,j),eq(i+1,j)) = A(eq(i,j),eq(i+1,j))-dt*dy/dx; end
            if i>1;  A(eq(i,j),eq(i-1,j)) = A(eq(i,j),eq(i-1,j))-dt*dy/dx; end
            if j<ny; A(eq(i,j),eq(i,j+1)) = A(eq(i,j),eq(i,j+1))-dt*dx/dy; end
            if j>1;  A(eq(i,j),eq(i,j-1)) = A(eq(i,j),eq(i,j-1))-dt*dx/dy; end
end; end
iXm=round(nx/2); iYm=round(ny/2);
A(eq(iXm,iYm),:)=A(eq(iXm,iYm),:)*0; A(eq(iXm,iYm),eq(iXm,iYm)) = 1; % XXX is this needed?
B=inv(A);

% Apply boundary condition
V(:,1)=Vw; V(:,end)=Ve; V(1,:)=2*Vs-V(2,:); V(end,:)=2*Vn-V(end-1,:);  % XXX is U=0 or U=Ubar in corners?
U(1,:)=Us; U(end,:)=Un; U(:,1)=2*Ue-U(:,2); U(:,end)=2*Uw-U(:,end-1);  % XXX is V=0 or V=Ubar in corners?
% Time loo
while t<tmax
    t
    uplot = (U(:,1:nx)+U(:,2:nx+1))/2;
    vplot = (V(1:ny,:)+V(2:ny+1,:))/2;
    Len = sqrt(uplot.^2+vplot.^2+eps);
    uplot=uplot./Len; vplot=vplot./Len;
    pcolor(X,Y,P);
        axis equal
%         axis([0 dimx 0 dimy])
        hold on 
        colormap(jet)
        colorbar
        shading interp
        quiver(X,Y,uplot,vplot,0.6,'k-');
        %title({['2D Cavity Flow with Re = ',num2str(UN*dimx/nu)];['CFL = ',num2str(UN/dx*dt),', \itt = ',num2str(t)]})
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('Spatial co-ordinate (y) \rightarrow')
        drawnow;
        hold off 
        pause(0.1)
    
    res=1e10; iter=0;
    P=P.*0; u=U; v=V;
    while res>tol && iter<=maxiter
     % Solve the velocity 
    uuW=0.25*(U(2:end-1,1:end-2)+U(2:end-1,2:end-1)).^2;
    uuE=0.25*(U(2:end-1,3:end)+U(2:end-1,2:end-1)).^2;
    uvS=0.25*(U(1:end-2,2:end-1)+U(2:end-1,2:end-1)).*(V(2:end-2,1:end-1)+V(2:end-2,2:end));
    uvN=0.25*(U(3:end,2:end-1)+U(2:end-1,2:end-1)).*(V(3:end-1,1:end-1)+V(3:end-1,2:end));
    u(2:end-1,2:end-1) = U(2:end-1,2:end-1) + dt/(dx*dy)*(    ...
          ni*( ...
               dy*(U(2:end-1,3:end  )-U(2:end-1,2:end-1))/dx  ... 
              -dy*(U(2:end-1,2:end-1)-U(2:end-1,1:end-2))/dx  ...
              +dx*(U(3:end  ,2:end-1)-U(2:end-1,2:end-1))/dy  ...
              -dx*(U(2:end-1,2:end-1)-U(1:end-2,2:end-1))/dy   ...
             ) ...
          +dy*(uuW-uuE)+dx*(uvS-uvN)  ... check here
          -dx*dy*( P(2:end-1,2:end)-P(2:end-1,1:end-1) )/dx  ... 
                                                         );
    uvW=0.25*(V(2:end-1,2:end-1)+V(2:end-1,1:end-2)).*(U(1:end-1,2:end-2)+U(2:end,2:end-2));
    uvE=0.25*(V(2:end-1,2:end-1)+V(2:end-1,3:end)).*(U(1:end-1,3:end-1)+U(2:end,3:end-1));
    vvS=0.25*(V(2:end-1,2:end-1)+V(1:end-2,2:end-1)).^2;
    vvN=0.25*(V(2:end-1,2:end-1)+V(3:end,2:end-1)).^2;
    v(2:end-1,2:end-1) = V(2:end-1,2:end-1) + dt/(dx*dy)*(    ...
          ni*( ...
               dy*(V(2:end-1,3:end  )-V(2:end-1,2:end-1))/dx  ... 
              -dy*(V(2:end-1,2:end-1)-V(2:end-1,1:end-2))/dx  ...
              +dx*(V(3:end  ,2:end-1)-V(2:end-1,2:end-1))/dy  ...
              -dx*(V(2:end-1,2:end-1)-V(1:end-2,2:end-1))/dy   ...
             ) ...
          +dy*(uvW-uvE)+dx*(vvS-vvN)  ... check here
          -dx*dy*( P(2:end,2:end-1)-P(1:end-1,2:end-1) )/dy  ... 
                                                         );  
    % Apply boundary condition
    v(1,:)=2*Vs-v(2,:); v(end,:)=2*Vn-v(end-1,:); v(:,1)=Vw; v(:,end)=Ve;   % XXX is U=0 or U=Ubar in corners?
    u(:,1)=2*Ue-u(:,2); u(:,end)=2*Uw-u(:,end-1); u(1,:)=Us; u(end,:)=Un;    % XXX is V=0 or V=Ubar in corners?
    % Continuity equation residual
    b=zeros(ny,nx);
    b(2:end-1,2:end-1)=-dy*(u(2:end-1,3:end-1)-u(2:end-1,2:end-2))-dx*(v(3:end-1,2:end-1)-v(2:end-2,2:end-1));
    b(2:end-1,1)=-dy*(u(2:end-1,2)-u(2:end-1,1))-dx*(v(3:end-1,1)-v(2:end-2,1));
    b(2:end-1,nx)=-dy*(u(2:end-1,nx+1)-u(2:end-1,nx))-dx*(v(3:end-1,nx)-v(2:end-2,nx));
    b(1,2:end-1)=-dy*(u(1,3:end-1)-u(1,2:end-2))-dx*(v(2,2:end-1)-v(1,2:end-1));
    b(ny,2:end-1)=-dy*(u(ny,3:end-1)-u(ny,2:end-2))-dx*(v(ny+1,2:end-1)-v(1,2:end-1));
    % Boundary condition for Pressure correction
    b(iYm,iXm)=0; %P(iYm,iXm); % Here set zero pressure
    % Laplacian of p'
    P(:) = alpha*B*b(:)+P(:);
    % Compute residual
    res = max(abs(b(:))); iter = iter+1;
    end
    iter
    U=u; V=v;
    % Apply boundary condition XXX probably not needed anymore
    V(1,:)=2*Vs-V(2,:); V(end,:)=2*Vn-V(end-1,:); V(:,1)=Vw; V(:,end)=Ve;   % XXX is U=0 or U=Ubar in corners?
    U(:,1)=2*Ue-U(:,2); U(:,end)=2*Uw-U(:,end-1); U(1,:)=Us; U(end,:)=Un;    % XXX is V=0 or V=Ubar in corners?
    % Update time
    t=t+dt;
end
