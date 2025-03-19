clear all, close all, clc
dt=10;
n=25000;
t=0:dt:(n-1)*dt;
mt=5.9736e24;
G=6.67e-11;

%CI 1
yc1(1,:)=[4.223e7,0,0,3071];
%methode ode23
[t,y11]=ode23(@fctepart2,t,yc1(1,:));
%methode ordre 2 (rk2)
y12(1,:)=yc1(1,:);
for i=1:n-1
    k1=fctepart2(1,y12(i,:));
    k2 = fctepart2(1, (y12(i,:)' + dt/2 * k1) )';
    y12(i+1,:)=y12(i,:)+dt*k2;
end

%CI 2
yc2(1,:)=[4.223e7,0,0,1.2*3071];
% methode ode23
[t,y21]=ode23(@fctepart2,t,yc2(1,:));
%methode ordre 2 (rk2)
y22(1,:)=yc2(1,:);
for i=1:n-1
    k1=fctepart2(1,y22(i,:));
    k2 = fctepart2(1, (y22(i,:)' + dt/2 * k1) )';
    y22(i+1,:)=y22(i,:)+dt*k2;
end

%part 3, calcul des energies 
Emode231=1/2.*100.*(y11(:,4).^2+y11(:,2).^2)-(mt*G*100)./((y11(:,1).^2+y11(:,3).^2).^(1/2));
Emode232=1/2.*100.*(y21(:,4).^2+y21(:,2).^2)-(mt*G*100)./((y21(:,1).^2+y21(:,3).^2).^(1/2));
EPod231=-(mt*G*100)./((y11(:,1).^2+y11(:,3).^2).^(1/2));
ECode231=1/2.*100.*(y11(:,4).^2+y11(:,2).^2);
EPode232=-(mt*G*100)./((y21(:,1).^2+y21(:,3).^2).^(1/2));
ECode232=1/2.*100.*(y21(:,4).^2+y21(:,2).^2);

%plot des orbites pour la methode ODE23 avec les 2 CI
subplot(3,2,1)
plot(y11(:,1),y11(:,3),y21(:,1),y21(:,3))
legend('CI1','CI2' )
title('ODE23')
xlabel('x')
ylabel('y')
axis equal

%plot des orbites pour la methode RK2 avec les 2CI
subplot(3,2,2)
plot(y12(:,1),y12(:,3),y22(:,1),y22(:,3))
legend('CI1','CI2' )
title('RK2/heun')
xlabel('x')
ylabel('y')
axis equal

%Plot des energies pour les CI1
subplot(3,2,3)
plot(t,Emode231,t,EPod231,t,ECode231)
title('Ode23 Energie CI1')
legend('EM','EP','EC')
xlabel('t')
ylabel('E')

%plot des energies pour les CI2
subplot(3,2,4)
plot(t,Emode232,t,EPode232,t,ECode232)
title('Ode23 Energie CI2')
legend('EM','EP','EC')
xlabel('t')
ylabel('E')

% err absolu
err=abs(y11(:,1)-y12(:,1));
subplot(3,2,5.5)
plot(t,err)
title('Ecart absolu Heun-Ode23')
xlabel('t')
ylabel('écart absolu')


