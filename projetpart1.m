clear all, close all, clc
dt=0.001;
n=10000;
t=0:dt:(n-1)*dt;

%CI 1
yc1(1,:)=[0.5,0];
%methode d'Euler
y1e(1,:)=yc1(1,:);
for i=1:n-1
    y1e(i+1,:)=y1e(i,:)+dt*fctepart1(t(i),y1e(i,:));
end

%methode ordre 2 (rk2)
y12(1,:)=yc1(1,:);
for i=1:n-1
    k=fctepart1(t(i),y12(i,:));
    y12(i+1,:)=y12(i,:)+dt*fctepart1(t(i)+dt/2,y12(i,:)+dt/2*k);
end

% ordre 4
y14=zeros(n,2);
y14(1,:)=yc1(1,:);
for i=1:n-1
    k1=fctepart1(t(i),y14(i,:));
    k2=fctepart1(t(i)+dt/2,y14(i,:)+dt/2*k1);
    k3=fctepart1(t(i)+dt/2,y14(i,:)+dt/2*k2);
    k4=fctepart1(t(i)+dt,y14(i,:)+dt*k3);
    y14(i+1,:)=y14(i,:)+dt*(k1+2*k2+2*k3+k4)/6;
end

%CI 2 
yc2(1,:)=[1,-2.5];
%methode d'Euler
y2e(1,:)=yc2(1,:);
for i=1:n-1
    y2e(i+1,:)=y2e(i,:)+dt*fctepart1(t(i),y2e(i,:));
end
%methode ordre 2 (rk2)
y22(1,:)=yc2(1,:);
for i=1:n-1
    k=fctepart1(t(i),y22(i,:));
    y22(i+1,:)=y22(i,:)+dt*fctepart1(t(i)+dt/2,y22(i,:)+dt/2*k);
end
% ordre 4
y24=zeros(n,2);
y24(1,:)=yc2(1,:);
for i=1:n-1
    k1=fctepart1(t(i),y24(i,:));
    k2=fctepart1(t(i)+dt/2,y24(i,:)+dt/2*k1);
    k3=fctepart1(t(i)+dt/2,y24(i,:)+dt/2*k2);
    k4=fctepart1(t(i)+dt,y24(i,:)+dt*k3);
    y24(i+1,:)=y24(i,:)+dt*(k1+2*k2+2*k3+k4)/6;
end

%trace
subplot(1,3,1)
plot(y1e(:,1),y1e(:,2),'o','MarkerIndices',1:30:length(y1e(:,2)))
hold on ;
plot(y12(:,1),y12(:,2),'x','MarkerIndices',1:30:length(y12(:,2)))
plot(y14(:,1),y14(:,2))
legend('Euler','rk2','rk4')
title('Portrait de phase avec CI1')
xlabel('θ')
ylabel('dθ/dt')

subplot(1,3,2)
plot(y2e(:,1),y2e(:,2),'o','MarkerIndices',1:30:length(y2e(:,2)))
hold on ;
plot(y22(:,1),y22(:,2),'x','MarkerIndices',1:30:length(y22(:,2)))
plot(y24(:,1),y24(:,2))
title('Portrait de phase avec CI2')
xlabel('θ')
ylabel('dθ/dt')
legend('Euler','rk2','rk4')


%erreur absolue
err=abs(y1e(:,1)-y12(:,1));
subplot(1,3,3)
plot(t,err)
title('Ecart absolu Euler-Rk2')
legend('Ecart absolu')
xlabel('t')
ylabel('écart absolu')

