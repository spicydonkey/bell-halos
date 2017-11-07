az_idx=85;       %[1,99]
el_idx=32;       %[1,50]

a=amp_m{2};
p=squeeze(P_rabi_m{2}(az_idx,el_idx,:));

feq_rabi=@(p,a) cos(abs(p(1)+p(2)*a+p(3)*a.^2+p(4)*a.^3)).^2;
p0=[0.1,0.1,0.1,0.1];

fopts=statset('TolFun',1e-12,'TolX',1e-12,'MaxIter',1e3);

f=fitnlm(a,p,feq_rabi,p0,'Options',fopts);

x=linspace(0,0.55);
y=feval(f,x);

figure();
hold on;
plot(a,p,'o');
plot(x,y);
ylim([0,1]);