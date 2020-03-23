%%***********************************************
function M=RTM(channelnum,Seis,Vv)
%     h = figure;

% filename = 'testAnimated1.gif';

dt=0.0005;    % time step, s
tt=-0.04:dt:0.04;
fm=50;
A=1;
wave=fliplr(Seis');
wave=wave';
T=0.4;        % transmit time, s


% Vv=load('marmsmooth.mat');Vv=Vv.V1;       % velocity map, m/s



[Nz Nx] = size(Vv);
pml=50;
m = Nx + 2*pml;
n = Nz + 2*pml;

dz=4;                              % grid step size on z, m
dx=4;                              % grid step size on z, m

d=2000*ones(n,m);                    % density,kg/m^3

V=zeros(n,m);
V(pml+1:pml+Nz,pml+1:pml+Nx) = Vv;
for i=1:pml
    V(pml+1:pml+Nz,i) = Vv(:,1);
    V(pml+1:pml+Nz,i+Nx+pml) = Vv(:,Nx);
end
for i=1:pml
    V(i,:) = V(pml+1,:);
    V(i+pml+Nz,:) = V(pml+Nz,:);
end
clear Vv;
Vmax=0;         % max velocity of P wave
for i=1:n
    for k=1:m
        if V(i,k) > Vmax
            Vmax=V(i,k);
        end
    end
end
R=1e-6;
ddx=zeros(n,m); ddz=zeros(n,m);
plx=pml*dx;      % pml size
plz=pml*dz;
for i=1:n
    for k=1:m
        %% zone 1
        if i>=1 && i<=pml && k>=1 && k<=pml
            x=pml-k;z=pml-i;
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
            ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        elseif i>=1 && i<=pml && k>m-pml && k<=m
            x=k-(m-pml);z=pml-i;
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
            ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        elseif i>n-pml && i<=n && k>=1 && k<=pml
            x=pml-k;z=i-(n-pml);
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
            ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        elseif i>n-pml && i<=n && k>m-pml && k<=m
            x=k-(m-pml);z=i-(n-pml);
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);
            ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        %% zone 2
        elseif i<=pml && k>pml && k<m-pml+1
            x=0;z=pml-i;
            ddx(i,k)=0;ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        elseif  i>n-pml && i<=n && k>pml && k<=m-pml
            x=0;z=i-(n-pml);
            ddx(i,k)=0;ddz(i,k)=-log(R)*3*Vmax*z^2/(2*plz^2);
        %% zone 3
        elseif i>pml && i<=n-pml && k<=pml
            x=pml-k;z=0;
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);ddz(i,k)=0;
        elseif i>pml && i<=n-pml && k>m-pml && k<=m
            x=k-(m-pml);z=0;
            ddx(i,k)=-log(R)*3*Vmax*x^2/(2*plx^2);ddz(i,k)=0;
        end
    end
end

K=zeros(n,m);   
for i=3:n-2
        for k=3:m-2
            K(i,k)=d(i,k)*V(i,k)^2;
        end
end
clear V;
shotnum = 1;
%for shot = 1:shotnum
z0=pml;                            % source location on z
x0=pml+channelnum;                       % source location on x


p=zeros(n,m);    
Pp=zeros(Nz,Nx);  
px=zeros(n,m);   
pz=zeros(n,m);   
Vx=zeros(n,m); 
Vz=zeros(n,m); 

Wavefeild = zeros(Nz,Nx,ceil(T/dt));

left = 0;right = 0;down = 0;radius = 0;

ttt = 0;
for t=dt:dt:T
    ttt = ttt + 1;
     p(z0,pml+1:pml+80)=wave(ttt,channelnum-39:channelnum+40);
%      p(z0,pml+40:30:pml+251)=wave(ttt,40:30:251);  
%    p(z0,pml+1:pml+300)=wave(ttt,:);
    

i=4:m-4;
k=4:m-4;
            Vz(i,k)=((1-0.5.*dt.*ddz(i,k)).*Vz(i,k)-dt.*((p(i+1,k)-p(i,k)).*1.171875-(p(i+2,k)-p(i-1,k)).*0.065104166666667+(p(i+3,k)-p(i-2,k)).*0.0046875)./d(i,k)./dz)./(1+0.5.*dt.*ddz(i,k));
            Vx(i,k)=((1-0.5.*dt.*ddx(i,k)).*Vx(i,k)-dt.*((p(i,k+1)-p(i,k)).*1.171875-(p(i,k+2)-p(i,k-1)).*0.065104166666667+(p(i,k+3)-p(i,k-2)).*0.0046875)./d(i,k)./dx)./(1+0.5.*dt.*ddx(i,k));

            pz(i,k)=((1-0.5.*dt.*ddz(i,k)).*pz(i,k)-K(i,k).*dt.*((Vz(i,k)-Vz(i-1,k)).*1.171875-(Vz(i+1,k)-Vz(i-2,k)).*0.065104166666667+(Vz(i+2,k)-Vz(i-3,k)).*0.0046875)./dz)./(1+0.5.*dt.*ddz(i,k));
            px(i,k)=((1-0.5.*dt.*ddx(i,k)).*px(i,k)-K(i,k).*dt.*((Vx(i,k)-Vx(i,k-1)).*1.171875-(Vx(i,k+1)-Vx(i,k-2)).*0.065104166666667+(Vx(i,k+2)-Vx(i,k-3)).*0.0046875)./dx)./(1+0.5.*dt.*ddx(i,k)); 
            p(i,k) = px(i,k) + pz(i,k);

    
    for i=1:n-2*pml
        for k=1:m-2*pml
            Pp(i,k)=p(i+pml,k+pml);
        end
    end
%     imagesc(Pp*10e9);
%     shading interp;      % interp to make image smooth
%     axis square;         % square coordinate system
%     colormap('gray');
%     pause(1e-4);
    Wavefeild(:,:,ttt) = p(pml+1:pml+Nz,pml+1:pml+Nx);
%         drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if t == dt 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
%     if mod(ttt,50)==0
%         time = ['Current time = ' num2str(t) ' of total ' num2str(T) 's'];disp(time);
%     end

end
    M = max(Wavefeild,[],3);
    
    
% M1=zeros(Nz,Nx);
% for channelnum=30:30:271
%     M1=M1+M(:,:,channelnum);
% end
% imagesc(M1);
%name = ['Wavefeild' num2str(channelnum)];
%save(name,'Wavefeild');




