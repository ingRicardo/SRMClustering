function aus = kodieren_rf(ein, rf, maxd, sigma, f_cut, rho, typ, plt)

% It encode analog values by fixed receptive fields 
%
t_neur=sum(rf); % total of neurons
[t_cj,t_in]=size(ein); % groups (t_cj) and analog inputs. (t_in)
aus=zeros(t_cj,t_neur); % it crates output matrix
%sigma=0.233./rf; %  max < 0.9*dmax
%
% Inputs
n=0; % 
for i=1:t_in
    ct=zeros(1,rf(i)); % center of curves
    %
    % It validates the number of neurons
    if rf(i)<3
        error(' Number of receptive fields must be greater than 3 ');
    end
    % Normalizes inputs no interval (0,1)
    range=max(ein(:,i))-min(ein(:,i)); % 
    ein(:,i)=(ein(:,i)-min(ein(:,i)))./range; %  normalized variables from 0 to 1
    % Center of curves
    if typ==1
        cdis=1/rf(i); % distance between centers
        for j=1:rf(i)
            ct(j)=(j-0.5)*cdis;
        end
    elseif typ==0
        cdis=1/(rf(i)-2); %  distance between centers
        for j=1:rf(i)
            ct(j)=(j-1.5)*cdis;
        end
    else
        error('Invalida parameter it mas be  0 or 1');
    end
    % Delays
    for j=1:t_cj
        for k=1:rf(i)
            aux=maxd-maxd.*exp(-((ein(j,i)-ct(k))^2)./(2.*sigma(i)^2));
            aus(j,k+n)=round(aux*(1/rho))*rho; % 
            if aus(j,k+n)>maxd*f_cut % the neuron fires
                aus(j,k+n)=-1; % or delay of -1
            end
        end
    end
    n=n+rf(i); % 
end
%
% Plot results
if plt==1
    figure;
    for i=1:t_cj
        for j=1:t_neur
            subplot(t_neur,1,j); hold on;
            H_line=plot([aus(i,j) aus(i,j)],[0 1],'b');
            set(H_line,'LineWidth',1);
            axis([0 maxd 0 2]);
            Ha_ax=gca;
            ylabel(sprintf('%1g',j),'FontSize',8);
            set(Ha_ax,'YTick',[]);
            set(Ha_ax,'XTick',[0 maxd]);
            %set(Ha_ax,'XTick',[0 maxd*0.2 maxd*0.4 maxd*0.6 maxd*0.8 maxd]);
            set(Ha_ax,'XGrid','on','XTickLabel',' ');
            if j==t_neur
                set(Ha_ax,'XTickLabelMode','auto');
                xlabel('t (ms)','FontSize',12);
            end
            if j==1
                title(sprintf('Classes: %1g - Inputs: %1g - ...Samples: %1g', t_in,t_neur,t_cj),'FontSize',14);
            end
        end
    end
    figure; hold on;
    x=0:0.01:1;
    for i=1:rf(end)
        y=maxd.*exp(-(x-ct(i)).^2./(2*sigma(end)^2));
        plot(x,y);
    end
    plot([0 1],[maxd*(1-f_cut) maxd*(1-f_cut)],'r');
elseif plt~=0
    error('Invalid plt parameter : it must be 0 or 1');
end