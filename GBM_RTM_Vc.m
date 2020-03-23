% clear;
% load('Seis');
dt=0.0005;    % 
% tt=-0.04:dt:0.04;
fs=1/dt;
% fm=50;
% A=1;
% wave=A*(1-2*(pi*fm*tt).^2).*exp(-(pi*fm*tt).^2);
T=0.5;        % 
% wave(round(T/dt))=0;    % 
% wave=fliplr(Seis(:,100));
% w2=fft(wave);
% Vv=load('V1.mat');Vv=Vv.V1;       % �����ٶ�,m/s
% [Nz Nx] = size(Vv);
% pml=50;
% m = Nx + 2*pml;
% n = Nz + 2*pml;
v0=2000;
Vv=v0*ones(300,300);
dz=4;                              % ���������С����λm
dx=4;                              % ���������С����λm




% %% Ray tracing
% 
% x0=20;      % x coordinate of takeoff point
% z0=1;       % z coordinate of takeoff point
% a0=pi/18;          % takeoff angle (radians) -pi/2~pi/2
% 
% c=cos(a0);
% s=sin(a0);
% 
% ray=zeros(nt,3);
% ray(1,1)=x0;
% ray(1,2)=z0;
% ray(1,3)=a0;
% 
% xx=x0*dx;
% zz=z0*dz;
% xt=x0;
% zt=z0;
% at=a0;
% for i=2:nt
%     ddx=Vv(zt,xt)*s*dt;
%     ddz=Vv(zt,xt)*c*dt;
%     xx=xx+ddx;
%     zz=zz+ddz;
%     dvx=(Vv(zt,xt+1)-Vv(zt,xt))/dx;
%     dvz=(Vv(zt+1,xt)-Vv(zt,xt))/dz;
%     da=-dvx*c*dt+dvz*s*dt;
%     at=at+da;
%     xt=floor(xx/dx);
%     zt=floor(zz/dz);
%     c=cos(at);
%     s=sin(at);
% %     scatter(xt,zt);
%     ray(i,1)=xt;
%     ray(i,2)=zt;
%     ray(i,3)=at;
% end
Uf=zeros(300,300,T*fs);
Ut=zeros(300,300,T*fs);
for x0=40:30:251
z0=0;
% a=0*pi/6;
% fn=20;
wave=fliplr(Seis(:,x0));
w2=fft(wave);
for fn=0:1/T:fs-1
w=fn*pi*2;
u0=w2(fn*T+1);
U=zeros(300,300);
% a=atan(-(150-x0)/150);
for a1=-pi/3:pi/3:pi/3
for x=1:300
    for z=1:300
%       a0=atan((x-x0)*dx/(z-z0)/dz);
%       a1=a+a0;
      s=cos(a1)*sqrt((x-x0)^2*dx^2+(z-z0)^2*dz^2);
      n=sin(a1)*sqrt((x-x0)^2*dx^2+(z-z0)^2*dz^2);
      U(z,x)=U(z,x)+u0*(v0^2/(v0^2/4+i*v0*s))^(1/2)*exp(i*w*s/v0-w/2/(v0^2/4+i*v0*s)*n^2);
    end
end
end

% U1=U1+U;
% imagesc(real(U));
Uf(:,:,fn/2+1)=U+Uf(:,:,fn/2+1);
end
x0
end
tic;
for x=1:300
    for z=1:300
        Ut(z,x,:)=ifft(Uf(z,x,:));
    end
end
toc;
% 
% for i=1:1000
%     imagesc(real(Ut(:,:,i)));
%     shading interp;      %Ϊʹͼ��ƽ�����в�ֵ
%     axis square;         %��������������ϵ
%     colormap('gray');
%     pause(1e-3);
%     
% end



% Figure;
tic;
M1=max(real(Ut),[],3);
M2=max(imag(Ut),[],3);
imagesc(M1+M2);
toc;