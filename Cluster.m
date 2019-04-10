% SNN Clustering
%
% it has 5 groups with gaussian distribution with mean equals to zero.
% each group has 12 points.
%
% output classes: 5
% Codification: fixed receptive fields
% Number of receptive fields: 1
% Type of synapses: multiples
%
%-------------------------------------------------------------------------

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
% Parameters initialization
out_neu=5; 
rho=0.1; 
%conj=125; 
conj=6; % training data
tau=3; % 
rf=[8 8 8 8 8]; % number of neurons by input
sigma=[1/(1.5*(rf(1)-2)) 1/(1.5*(rf(2)-2)) 1/(1.5*(rf(2)-2)) 1/(1.5*(rf(2)-2)) 1/(1.5*(rf(2)-2))]; 
f_cut=0.9;
maxd=10; 
%-------------------------------------------------------------------------
% Learning parameters
d=[0 1 2 3 4 5 6 7 8 9 10 11 12]; % sub-sinapses delays
w_max=1; w_min=0; % max and min weight
max_epoch=3; % number of epochs
t_max=30; % training max time 
eta=0.35; % learning rate
beta=0.2; % learning window
nu=5.0; % neighborhood
dx=2.3; % learning window motion
% p/ left (+) or right (-)
%-------------------------------------------------------------------------
% Parameters calculation
kapa=1-(nu^2)/(2*log(beta/(beta+1))); % number of learning window
in_neu=sum(rf); % number of input neurons
ssin=size(d,2); % number of sub synapses
%teta=1.5*(in_neu/out_neu)*(w_max-w_min)*ssin/2;
teta=12; % boundary
%-------------------------------------------------------------------------
% Weight initialize
w=zeros(in_neu,out_neu,ssin);
for i=1:in_neu
    for j=1:out_neu
        w(i,j,:)=w_min+rand(1,ssin).*(w_max-w_min);
    end
end
%-------------------------------------------------------------------------
% Training
tic;
aus=kodieren_rf(aux1,rf,maxd,sigma,f_cut,rho,0,0); % input codification
ctrl_1=1; % 
delta_t=zeros(1,1,ssin);
h=waitbar(ctrl_1/max_epoch,sprintf('Executando %i de %i iteracoes',ctrl_1, max_epoch));
%show_w(w);
f_times=zeros(2,1);
while ctrl_1<=max_epoch
    for z=1:conj
        tr=aus(z,:);
        [neu, inst]=wer_schiesst(tr,t_max,w,d,tau,rho,teta,0);
        f_times=[f_times [neu;inst]];
        if neu~=0 && size(neu,2)==1
            for i=1:in_neu
                if tr(i)==-1
                    w(i,neu,:)=w(i,neu,:)-eta;
                else
                    delta_t(1,1,:)=tr(i)+d-inst; % delta_t
                    w(i,neu,:)=w(i,neu,:)+eta.*((1+beta).*...
                        exp(-(((delta_t+2.3).^2)./(2*kapa-2)))-beta);
                end
                ndx=find(w(i,neu,:)<w_min); % boundary
                for u=1:size(ndx,1)
                    w(i,neu,ndx(u))=w_min;
                end
                ndx=find(w(i,neu,:)>w_max); % boundary
                for u=1:size(ndx,1)
                    w(i,neu,ndx(u))=w_max;
                end
            end
        end
    end
    ctrl_1=ctrl_1+1; % 
    waitbar(ctrl_1/max_epoch,h,sprintf('Executing %i de %i iteration',...
    ctrl_1, max_epoch));
    teta=teta+(0.3*teta)/max_epoch;
   % show_w(w);
end
close(h);
toc;
%-------------------------------------------------------------------------
% Testing net
[conj,in_neu]=size(aus);
class=zeros(1,conj); % classes table
for i=1:conj
    tr=aus(i,:); % number of testing data
    [neu,inst]=wer_schiesst(tr,t_max,w,d,tau,rho,teta,0); % who fired?
    if size(neu,2)>1 % equal...
        neu=0; % class = 0
    end
    class(i)=neu; % it saves the class
end
%-------------------------------------------------------------------------
% Grap the result of the test
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