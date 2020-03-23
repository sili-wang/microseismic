%clear 
dt=0.0005;     % time sampling interval
T=0.5;        % time length
nt=T/dt;      % number of time samples

tt=-0.04:dt:0.04;
fs=1/dt;
fm=50;
A=1;
% wave=A*(1-2*(pi*fm*tt).^2).*exp(-(pi*fm*tt).^2);
% T=0.5;        % 
% wave(round(T/dt))=0;    % 
% w=fft(wave);


% Vv=load('marmvel.mat');Vv=Vv.V1;       % m/s
% Vv=imgaussfilt(Vv,4);
Vv=load('marmsmooth.mat');Vv=Vv.V1;
% v0=2000;
% Vv=v0*ones(300,300);


[nz nx] = size(Vv);               % nx number of x samples
                                  % nz number of z samples
wr=5;
Va=mean(mean(Vv));
L0=2*pi*Va/wr;

                                  
dz=5;                              % z sampling interval m
dx=5;                              % x sampling interval m

d=2000*ones(nz,nx);                    % density kg/m^3
Us=zeros(300,300,1000);
imagesc(Vv);
hold on;
%% Ray tracing
for x0=10:30:290
% x0=100;      % x coordinate of takeoff point
z0=2;       % z coordinate of takeoff point
a0=atan((150-x0)/150);         % takeoff angle (radians) -pi/2~pi/2
wave=fliplr(Seis(:,x0))';
w=fft(wave);

P=j/Vv(z0,x0);
Q=wr*L0^2/Vv(z0,x0);


c=cos(a0);
s=sin(a0);

xx=x0*dx;
zz=z0*dz;
xt=x0;
zt=z0;
at=a0;
ss=0;
Uf=zeros(300,300,1000);
ray=zeros(nt,6);
ray(1,1)=xx;
ray(1,2)=zz;
ray(1,3)=a0;
ray(1,4)=P;
ray(1,5)=Q;
ray(1,6)=s;
for i=2:nt
    ddx=Vv(zt,xt)*s*dt;
    ddz=Vv(zt,xt)*c*dt;
    xx=xx+ddx;
    zz=zz+ddz;
    ss=ss+sqrt(ddx^2+ddz^2);
    Q=Q+dt*Vv(zt,xt)^2*P;
    P=P-Q/Vv(zt,xt)*dt*(((Vv(zt,xt+1)-2*Vv(zt,xt)+Vv(zt,xt-1))/(dx^2))*c^2-2*(Vv(zt+1,xt+1)-Vv(zt,xt+1)-Vv(zt+1,xt)+Vv(zt,xt))/(dx*dz)*c*s+((Vv(zt+1,xt)-2*Vv(zt,xt)+Vv(zt-1,xt))/(dz^2))*s^2);
    
    dvx=(Vv(zt,xt+1)-Vv(zt,xt))/dx;
    dvz=(Vv(zt+1,xt)-Vv(zt,xt))/dz;
    da=-dvx*c*dt+dvz*s*dt;
    at=at+da;
    xt=floor(xx/dx);
    zt=floor(zz/dz);
    c=cos(at);
    s=sin(at);
    scatter(xt,zt);
    
    
    ray(i,1)=xx;
    ray(i,2)=zz;
    ray(i,3)=at;
    ray(i,4)=P;
    ray(i,5)=Q;
    ray(1,6)=ss;
    
    for n=0:4:(ss/4+60)
        Usn=w*(Vv(zt,xt)/Q)^(1/2).*exp(j*2*pi*(0:2:1999)*(i-1)*dt+j*2*pi*(0:2:1999)/2*P/Q*n^2);
        x=floor((xx-n*c)/dx);
        z=floor((zz+n*s)/dz);
        if x>0 && x<300 && z>0 && z<300
           Uf(z,x,:)=Usn; 
        end
        x=floor((xx+n*c)/dx);
        z=floor((zz-n*s)/dz);
        if x>0 && x<300 && z>0 && z<300
           Uf(z,x,:)=Usn; 
        end
    end

    
    
end

Us=Us+Uf;
end


for x=1:300
    for z=1:300
        Ut(z,x,:)=ifft(Us(z,x,:));
    end
end


for i=1:500
    imagesc(real(Ut(:,:,i)));
    shading interp;      %Ϊʹͼ��ƽ�����в�ֵ
    axis square;         %��������������ϵ
    colormap('gray');
    pause(1e-3);
    
end


