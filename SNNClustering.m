% cl_img_snn - Image classification with SNN
%
% Name of the file 'OriginalImage.jpg'
% Classes (output): 3
% Codification: fixed receptive fields
% Number of Receptive Fields: 1
% Type of synapses: multiples
%--------------------------------------------------------------------------
clear all
%
% file transformation to double to be able to process it
aux=double(imread('OriginalImage.jpg')); 
[x1 y1 z1]=size(aux);
%--------------------------------------------------------------------------
% 3D transformation to 2D  x1 x y1 x z1 --> (x1 * y1) x z1
aux1=zeros((x1*y1),z1); conta=1; % aux1 = analog Inputs
for i=1:x1
    for j=1:y1
        for k=1:z1
            aux1(conta,k)=aux(i,j,k);
        end
        conta=conta+1;
    end
end
%--------------------------------------------------------------------------
% Initial Parameters 
out_neu=4; % output neurons (classes)
rho=0.1; % time
conj=32768; % !!!! % number of training groups 
tau=3; % time EPSP
rf=[6 6 6]; % number of neurons by input
% receptive fields amplitud
sigma=[1/(1.5*(rf(1)-2)) 1/(1.5*(rf(2)-2)) 1/(1.5*(rf(3)-2))];
f_cut=0.9;
maxd=20; % interval of encoding
%-------------------------------------------------------------------------
% Learning parameters
d=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; % sub synapses delays
w_max=1; w_min=0; % max and min weights
max_epoch=3; % max epoch
t_max=30; % max time of training
eta=0.35; % learning rate
beta=0.2; % learning window
nu=5.0; % neighborhood
%-------------------------------------------------------------------------
% Parameters calculation
kapa=1-(nu^2)/(2*log(beta/(beta+1))); % number of learning window
in_neu=sum(rf); % number of imput neurons
ssin=size(d,2); % number of sub synapses
teta=9.0; % boundary
%-------------------------------------------------------------------------
% Weight initialization
w=zeros(in_neu,out_neu,ssin);
for i=1:in_neu
    for j=1:out_neu
        w(i,j,:)=w_min+rand(1,ssin).*(w_max-w_min);
    end
end
%-------------------------------------------------------------------------
% Training starts..
aus=kodieren_rf(aux1,rf,maxd,sigma,f_cut,rho,0,0); % input encoding
ctrl_1=1; 
delta_t=zeros(1,1,ssin);
h=waitbar(ctrl_1/max_epoch,sprintf('Executing %i de %i iterations',...
    ctrl_1, max_epoch));
while ctrl_1<=max_epoch
    for z=1:conj
        tr=aus(z,:);
        [neu, inst]=wer_schiesst(tr,t_max,w,d,tau,rho,teta,0);
        if neu~=0 & size(neu,2)==1
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
    ctrl_1=ctrl_1+1; % incrementa contador de epocas
    waitbar(ctrl_1/max_epoch,h,sprintf('Executing %i de %i iterations',...
        ctrl_1, max_epoch));
    teta=teta+(0.3*teta)/max_epoch;
end
close(h);
%-------------------------------------------------------------------------
%Define classes from a  SNN trained
[conj,in_neu]=size(aus);
class=zeros(1,conj); % new table of classes
for i=1:conj
    tr=aus(i,:); % samples for testing
    [neu,inst]=wer_schiesst(tr,t_max,w,d,tau,rho,teta,0); % who fired?
    if size(neu,2)>1 % iquals...
        neu=out_neu+1; % class = n de classes+1
    end
    
    if neu == 0
        neu = 1;
    end
    
    class(i)=neu; % saves a class
end
% calculation of center of a class
c=zeros(out_neu+1,size(rf,2)); % 
for i=1:out_neu
    a=find(class==i);
    for n=1:size(a,2)
        c(i,1)=c(i,1)+aux1(a(n),1);
        c(i,2)=c(i,2)+aux1(a(n),2);
        c(i,3)=c(i,3)+aux1(a(n),3);
    end
    c(i,:)=c(i,:)./size(a,2);
end
%a=find(class==out_neu+1);
%for n=1:size(a,2)
% c(out_neu+1,:)=[256 0 0]; % 
%end
c_int=round(c); % centers must be integer
%imagy=uint8(zeros((x1*y1),z1)); % new output area unsigned int
imagy= uint8(zeros((x1*y1),z1)); % new output area unsigned int
for i=1:(x1*y1) % substitution of the vector of center classes
        imagy(i,:)=c_int(class(i),:);
end
imagz=uint8(zeros(x1,y1,z1)); conta=1; %gets 3D array
for i=1:x1
    for j=1:y1
        for k=1:z1
            imagz(i,j,k)=imagy(conta,k);
        end
        conta=conta+1;
    end
end
imwrite(imagz,'output_snn.jpg','jpg','Quality',100); % saves image
% results plotting
r=aux1(:,1);
g=aux1(:,2);
b=aux1(:,3);
figure(1); plot3(r,g,b,'.'); grid on;
box on; axis tight;
title('Points Distribution','FontSize',14);
xlabel('Red','FontSize',14);
ylabel('Green','FontSize',14);
zlabel('Blue','FontSize',14);
% centers plotting
cr=c(:,1);
cg=c(:,2);
cb=c(:,3);
figure(2); hold on; grid on;
box on; axis([0 256 0 256 0 256]);
plot3(cr,cg,cb,'r*');
plot3([0; 255],[0; 255],[0; 255],'r');
title('SNN Center Location','FontSize',14);
xlabel('Red','FontSize',14);
ylabel('Green','FontSize',14);
zlabel('Blue','FontSize',14);