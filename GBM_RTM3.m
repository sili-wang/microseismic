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
wr=40;
Va=mean(mean(Vv));
L0=2*pi*Va/wr;

                                  
dz=4;                              % z sampling interval m
dx=4;                              % x sampling interval m

d=2000*ones(nz,nx);                    % density kg/m^3
Us=zeros(300,300,1000);
imagesc(Vv);
hold on;
for i=1:60
    
end


for x=1:300
    for z=1:300
        Usn=zeros(1,1000);
        for i=1:2
        s=fliplr(sqrt(((z-1)*dz)^2+((x-i)*dx)^2));
        Va=0;
        for j=1:z
            x0=floor((x-i)*(j-1)/(z)+i);
            Va=Va+Vv(j,x0);
        end
        Va=Va/z;
        u0=fft(Seis(:,i));
        w=0:2:1999;
        Usn=Usn+u0'.*((Va^2/(Va^2/4+i*Va*s))^(1/2)*exp(i*w*s/Va));
        end
        Uf(z,x,:)=Usn;
    end
end


for x=1:300
    for z=1:300
        Ut(z,x,:)=ifft(Us(z,x,:));
    end
end


for i=1:1000
    imagesc(real(Ut(:,:,i)));
    shading interp;      %Ϊʹͼ��ƽ�����в�ֵ
    axis square;         %��������������ϵ
    colormap('gray');
    pause(1e-3);
    
end


