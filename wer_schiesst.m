function [neu, inst] = wer_schiesst(in, tmax, w, d, tau, rho, teta, plt)
%
% It returns if the neuron fires or not over a time instant.
%--------------------------------------------------------------------------
[ein,aus,ssin]=size(w); % ein = inputs (neurons), aus = outputs
ndx=find(in~=-1); % input fired neurons
tam=size(ndx,2);
if tam==0 % if no fires
    neu=0; inst=-1;
    return
end
t=0; % 
neu=0;
while neu==0 & t<=tmax
    out=zeros(1,aus);
    for j=1:aus
        for i=1:tam
            dt=t-in(ndx(i));
            for k=1:ssin
                dtt=(dt-d(k))/tau;
                sai=w(ndx(i),j,k)*dtt*exp(1-dtt);
                if sai<0
                    sai=0;
                end
                out(j)=out(j)+sai;
            end
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
if neu==0 %
    inst=-1;
end
%--------------------------------------------------------------------------
