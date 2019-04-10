function [neu, inst] = wer_spuckt(in, tmax, w, tau, rho, teta, plt)
%
%  It returns if the neuron fires or not over a time instant.

%
%--------------------------------------------------------------------------
[ein,aus]=size(w); % ein = inputs (neurons), aus = output
ndx=find(in~=-1); % it has fired inputs
tam=size(ndx,2); % how many neurons fired?
if tam==0 % if nothing fired
    neu=0; inst=-1;
    return
end
t=0; % 
neu=0;
while neu==0 & t<=tmax
    out=zeros(1,aus);
    for j=1:aus
        for i=1:tam
            dt=(t-in(ndx(i)))/tau(ndx(i),j);
            sai=w(ndx(i),j)*dt*exp(1-dt);
            if sai<0
                sai=0;
            end
            out(j)=out(j)+sai;
        end
        if plt==1
            plota(aus,j,t,out(j),teta);
        end
    end
    if max(out)>=teta
        neu=find(out==max(out));
        inst=t;
    end
    t=t+rho;
end
if neu==0 
    inst=-1;
end
%--------------------------------------------------------------------------
function plota(p1, p2, px, py, teta)
subplot(p1,1,p2); hold on;
plot(px,py,'.',px,teta,'r.'); drawnow;