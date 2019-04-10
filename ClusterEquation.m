%A.3.3 Using equation (5.9)
% Classification using SNN
% It uses 5 groups with gaussian distribution with zero mean.
% Each group has 12 points.
%
% Classes (output): 5
% Codification: Fixed Receptive fields
% Number of receptive fields: 1
% Type of synapses: simples
% Training : weight and time adaptation 
%
%-------------------------------------------------------------------------
% Input
clear all;

a = 0;
b = 100;

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

%-------------------------------------------------------------------------

% Initial parameters
out_neu=5; % number of output neurons (classes)
rho=0.1; % time
%conj=125; %
conj=12; % number of training samples
%tau=3; % time EPSP
rf=[8 8 8 8 8]; % number of input neurons
sigma=[1/(1.5*(rf(1)-2)) 1/(1.5*(rf(2)-2)) 1/(1.5*(rf(1)-2)) 1/(1.5*(rf(1)-2)) 1/(1.5*(rf(1)-2))]; % receptive fields amplitude
f_cut=0.9;
maxd=10; % codification interval
%-------------------------------------------------------------------------
% learning parameters
tau_max=9.8; tau_min=0.1; % min and max time
w_max=6.98; w_min=0; % min and max weights
max_epoch=10; % max number of epochs
t_max=30; % training max time
eta_w=0.035; % learing rate  b/ weights
beta=0.26; % learning window of weights
nu=3; % neihborhood 
dx=2.8; % shifting of W(t) / left (+) or rigth
%-------------------------------------------------------------------------
% parameters calculation
kapa=1-(nu^2)/(2*log(beta/(beta+1))); % number of learning window
in_neu=sum(rf); % number of input neurons
teta=12; % boundary
%-------------------------------------------------------------------------
% weight initializing and time 
w=w_min+rand(in_neu,out_neu).*(w_max-w_min);
tau=tau_min+rand(in_neu,out_neu).*(tau_max-tau_min);
%tau=tau_max.*(0.2452.*exp(-0.07922.*(w./w_max))...
% +0.09115.*exp(2.149.*(w./w_max)));
%-------------------------------------------------------------------------
% training starts..
tic;
aus=kodieren_rf(aux1,rf,maxd,sigma,f_cut,rho,0,0); % codifica entrada
ctrl_1=1; % 
h=waitbar(ctrl_1/max_epoch,sprintf('Executing  %i of %i iteration',...
    ctrl_1, max_epoch));
f_times=zeros(2,1);
while ctrl_1<=max_epoch
    for z=1:conj
        tr=aus(z,:);
        [neu, inst]=wer_spuckt(tr,t_max,w,tau,rho,teta,0);
        f_times=[f_times [neu;inst]];
        if neu~=0 & size(neu,2)==1
            for i=1:in_neu
                if tr(i)==-1
                    w(i,neu)=w(i,neu)-eta_w;
                    tau(i,neu)=tau_max*(0.2452*exp(-0.07922*(w(i,neu)/...
                        w_max))+0.09115*exp(2.149*(w(i,neu)/w_max)));
                else % weight fixing and time
                    delta_t=tr(i)-inst; % 
                    w(i,neu)=w(i,neu)+eta_w.*((1+beta).*...
                        exp(-(((delta_t+dx).^2)./(2*kapa-2)))-beta);
                    tau(i,neu)=tau_max*(0.2452*exp(-0.07922*(w(i,neu)/...
                        w_max))+0.09115*exp(2.149*(w(i,neu)/w_max)));
                end
                % if w(i,neu)<=0
                % tau(i,neu)=tau_min;
                % end
                if w(i,neu)>w_max
                    w(i,neu)=w_max;
                elseif w(i,neu)<w_min
                    w(i,neu)=w_min;
                end
                % if tau(i,neu)>tau_max
                % tau(i,neu)=tau_max;
                % elseif tau(i,neu)<tau_min
                % tau(i,neu)=tau_min;
                % end
            end
        end
    end
    ctrl_1=ctrl_1+1; 
    waitbar(ctrl_1/max_epoch,h,sprintf('Executing %i of %i iteration',...
        ctrl_1, max_epoch));
    teta=teta+(0.3*teta)/max_epoch;
    %show_w(w);
end
close(h);
toc;
%-------------------------------------------------------------------------
% test
[conj,in_neu]=size(aus);
class=zeros(1,conj); % new classes table
for i=1:conj
    tr=aus(i,:); % conjunto de teste
    [neu,inst]=wer_spuckt(tr,t_max,w,tau,rho,teta,0); % who fired?
    if size(neu,2)>1 % is it equals...
        neu=0; % class = 0
    end
    class(i)=neu; % it saves the class
end
% Plotting result
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