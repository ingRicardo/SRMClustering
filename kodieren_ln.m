function aus = kodieren_ln(ein, maxd, rho, plt)

%
[t_cj,t_in]=size(ein);
%

for i=1:t_in
    range=max(ein(:,i))-min(ein(:,i)); 
    ein(:,i)=(ein(:,i)-min(ein(:,i)))./range; 
end

aus=round((ein.*maxd)./rho).*rho;

aus
if plt==1
    figure;
    for i=1:t_cj
        for j=1:t_in
            subplot(t_in,1,j); hold on;
            H_line=plot([aus(i,j) aus(i,j)],[0 1],'b');
            set(H_line,'LineWidth',1);
            axis([0 maxd 0 2]);
            Ha_ax=gca;
            ylabel(sprintf('%1g',j),'FontSize',8);
            set(Ha_ax,'YTick',[]);
            
            set(Ha_ax,'XTick',[0 maxd]);
            %set(Ha_ax,'XTick',[0 maxd*0.2 maxd*0.4 maxd*0.6 maxd*0.8 maxd]);
            set(Ha_ax,'XGrid','on','XTickLabel',' ');
            if j==t_in
                set(Ha_ax,'XTickLabelMode','auto');
                xlabel('t (ms)','FontSize',12);
            end
            if j==1
                title(sprintf('Classes: %1g - Inputs: %1g - ...Samples: %1g', t_in,t_in,t_cj),'FontSize',14);
            end
        end
    end
elseif plt~=0
    error('Invalid plt parameter: 0,1');
end