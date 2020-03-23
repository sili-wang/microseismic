% A=zeros(300,300);
% B1=A;
% B2=A;
% for i=2:299
% 
%     name1 = ['Wavefeild' num2str(i-1)];
%     name2 = ['Wavefeild' num2str(i)];
%     name3 = ['Wavefeild' num2str(i+1)];
%     W1=load(name1);
%     W1=W1.Wavefeild;
%     W2=load(name2);
%     W2=W2.Wavefeild;
%     W3=load(name3);
%     W3=W3.Wavefeild;
% for j=1:1000
%     A1=W2(:,:,j).*W1(:,:,j);
%     B1=B1+A1;
%     A2=W2(:,:,j).*W3(:,:,j);
%     B2=B2+A2;
% end
%     A=A+B1+B2;
% end




% Stacked method
Stack=zeros(300,300,1000);
for i=1:300
    name = ['Wavefeild' num2str(i)]; 
    W1=load(name);
    Stack=Stack+W1.Wavefeild;
end

M = max(Stack,[],3);
imagesc(M)

% W1=load(name);
% Cross=W1.Wavefeild;
% for i=2:300
%     name = ['Wavefeild' num2str(i)]; 
%     W1=load(name);
%     Cross=Cross.*W1.Wavefeild/max(max(max(W1.Wavefeild)));
%     i
% end
M = max(Cross,[],3);
imagesc(M)

