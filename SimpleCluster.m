% SNN Classification
% 2 groups with gaussian distribution with mean equals to zero
%
% each group has 12 points.
% supervised adaptation method 
% axonal delays equals to 1.
%
% Classes (output): 2
% Codification: lineal with one neuron as reference 
% Receptive fields: 0
% Synapses type : simples
%
%-------------------------------------------------------------------------
% Prepara Entrada

a = 0;
b = 100;

aux1 = [(b-a).*rand(1,1) + a, (b-a).*rand(1,1) + a ;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a];

%-------------------------------------------------------------------------
tam=size(aux1,1); 
in_neu=size(aux1,2)+1; 
out_neu=2; 
rho=0.1; 
tau=1.84; 
maxd=10; 
tmax=30; 
teta=2.04; 
%-------------------------------------------------------------------------
% Training
%tra=50; trb=(tam/2)+tra; 
tra=6; trb=(tam/2)+tra; 
aus=kodieren_ln(aux1,maxd,rho,0); 
%ca=[2.0 2.5];%ca=[7.8 6.5];
ca=round((sum(aus(1:tra,:))./tra)./rho).*rho; % class 1 mean
%cb=[8.2 6.4];%cb=[3 3.2];
%cb=round((sum(aus(51:trb,:))./tra)./rho).*rho; %
cb=round((sum(aus(7:trb,:))./tra)./rho).*rho; %  class 2 mean
d=maxd-[ca(1) cb(1); ca(2) cb(2); maxd maxd]; % delays of centers
w=[1 1;1 1;1 1];%[1 1.42;1 1.38;1 1];
%-------------------------------------------------------------------------
% Tesing the net
class=zeros(1,tam); 
aus=kodieren_ln(aux1,maxd,rho,0); 
aus=[aus ones(tam,1).*maxd]; 
for k=1:tam
    t=0;
    neu=0;
    while neu==0 && t<=tmax
        out=zeros(1,out_neu);
        for j=1:out_neu
            for i=1:in_neu
                dt=t-aus(k,i);
                sai=w(i,j)*((dt-d(i,j))/tau)*exp(1-((dt-d(i,j))/tau));
                if sai<0
                    sai=0;
                end
                out(j)=out(j)+sai;
            end
            subplot(out_neu,1,j); hold on;
            plot(t,out(j),'.',t,teta,'r.'); drawnow;
        end
        if max(out)>=teta
            neu=find(out==max(out));
            inst=t;
        end
        t=t+rho;
    end
    if size(neu,2)>1 
        neu=0; % class = 0
    end
    if neu==0 % 
        class(k)=0; % class = 0
    end
    class(k)=neu; 
    %pause
    %h=gcf; % pega hadle da figura
    %close(h); % fecha figura
end
%-------------------------------------------------------------------------

ptos={'+k' 'ok' 'or' '+b' 'xk' 'sr' 'db'...
    '.k' '*r' 'ob' '+k' 'xr' 'sb' 'dk'...
    '.r' '*b' 'ok' '+r' 'xb' 'sk' 'dr'};
figure; hold on;
for i=1:out_neu
    a=find(class==i);
    for n=1:size(a,2)
        plot(aux1(a(n),1),aux1(a(n),2),char(ptos(i)));
    end
end
a=find(class==0);
for n=1:size(a,2)
    plot(aux1(a(n),1),aux1(a(n),2),char(ptos(out_neu+1)));
end
%--------------------------------------------------------------------------
%function plota(p1, p2, px, py, teta)
% subplot(p1,1,p2); hold on;
% plot(px,py,'.',px,teta,'r.'); drawnow;
