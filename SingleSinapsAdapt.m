%A.3.1 First method - Time constant adaptation
% single_syn_5a - classification with SNN
% It's equivalent to simple synapses. It transforms multiple synapses
% in to simples.
% Usa o resultado de cl_spk_5a.m com 16
% neuronios de entrada e cinco de saida.
%

% Each group has 12 points.
%
% Classes (output): 5
% Codification: Fixed receptive fields
% Number of receptive fields : 1
% Type of synapses : simples
% Delays: Constant of time EPSP
%
%-------------------------------------------------------------------------
% Input
clear all;
%load('single_syn.mat'); % carrega vetores de pesos w
% Em single_syn.mat estao:
% aux1 = Analog inputs
% aus = Coded Analog inputs
% dx = Values of tau ( time EPSP )
% wx = Weights ( max time of EPSP )
% epsp = .
%-------------------------------------------------------------------------
% Initial parameters

a = 0;
b = 100;

dx = [(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;]

aus = [(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;]  

 aux1 = [(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;]    
   
   
wx = [(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;
       (b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a,(b-a).*rand(1,1) + a;]      
   
[in_neu out_neu]=size(dx); % number of inputs and outputs
conj=size(aus,1); % number of groups to classify
rho=0.1; % time resolution
teta=12;
tmax=30;
%-------------------------------------------------------------------------
% Classification starts
class=zeros(1,conj); % new table of classes
for n=1:conj
    in=aus(n,:);
    ndx=find(in~=-1); % fired neurons
    tam=size(ndx,2);
    t=0; % 
    neu=0;
    % h=figure; hold on;
    while neu==0 & t<=tmax
        out=zeros(1,out_neu);
        for j=1:out_neu % for each output
            for i=1:tam
                if wx(ndx(i),j)==0 | dx(ndx(i),j)==0
                    sai=0;
                else
                    dt=(t-in(ndx(i)))/dx(ndx(i),j);
                    sai=wx(ndx(i),j)*dt*exp(1-dt);
                    if sai<0
                        sai=0;
                    end
                end
                out(j)=out(j)+sai;
            end
            % plot(t,out(j),t,teta);
        end
        if max(out)>=teta
            neu=find(out==max(out));
            inst=t;
        end
        t=t+rho;
    end
    % pause; close(h);
    if size(neu,2)>1 % equals...
        neu=0; % class = 0
    end
    class(n)=neu; % saves a class
end
%-------------------------------------------------------------------------
% Plotting result and testing
ptos={'.b' '*k' 'or' '+b' 'xk' 'sr' 'db'...
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