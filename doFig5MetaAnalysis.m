%do PL range meta analysis
clear
ncmap=jet(100);
bincmap=[1 0 0; 0 0 1; 1 1 0.5];
plotflag=0;
temp=dir('*\*\*.csv');
% temp(1:2)=[];
nf=length(temp);
q=0;
temp1=[];
temp2=[];
for f=1:nf
    fname=temp(f).name;
    %load data
    dat=load([temp(f).folder,'\',fname]);
    
   
    
    %define plot symbol
    %   circle = spikes, x = LFP, up triangle = BOLD, * = MEG, down
    %   triangle = EEG, square = widefield imaging
    %   blue = awake, red = non-awake, green = in vitro
    psshape='.'; pscol='b';
    if any(strfind(temp(f).folder,'spike')); psshape='o';  end
    if any(strfind(temp(f).folder,'LFP')); psshape='x'; end
    if any(strfind(temp(f).folder,'EEG')); psshape='v'; end
    if any(strfind(temp(f).folder,'ECOG')); psshape='+'; end
    if any(strfind(temp(f).folder,'MEG')); psshape='*'; end
    if any(strfind(temp(f).folder,'Widefield')); psshape='s'; end
    if any(strfind(temp(f).folder,'vitro')); pscol='g'; end
    if any(strfind(temp(f).folder,'BOLD')); psshape='^'; end
    if any(strfind(temp(f).folder,'non-awake')); pscol='r'; end
%      if any(strfind(fname,'Ponce')); pscol='m';  end
    ps=[psshape,pscol];
    
    %remove any non-unique x values
    [~,IA,~] = unique(dat(:,1));
    dat=dat(IA,:);
    
    %interpolate to constant number of points per decade
%     rsiz=min(dat(:,1)):0.03:max(dat(:,1));
    rsiz=min(dat(:,1)):0.1:max(dat(:,1));    
    ns=length(rsiz);
    [rdat]=interp1(dat(:,1),dat(:,2),rsiz,'pchip');
    
    %extract DT from filename
    tind=strfind(fname,'DT');
    DT=str2num(fname(tind+2:tind+6))*0.0001;
    if ~any(DT); DT=0.00001; end
    
    
    %find linear fit with largest range that meets gof
    gofthresh=0.99;
    leeway=0.3;%*(max(dat(:,2))-min(dat(:,2)));
    bestplr=0;
    expon=0;
    for i=1:ns-1
        for j=2:ns
            plr = rsiz(j)-rsiz(i);
            if plr>bestplr
                p = polyfit(rsiz(i:j),rdat(i:j),1);
                fitdat = polyval(p,rsiz(i:j));
                gof = sum(abs(rdat(i:j)-fitdat)<leeway)/length(fitdat);
                if gof>gofthresh
                    bestplr=plr; %power law range
                    expon=p(1); %exponent
                    
                    if plotflag
                        figure(1)
                        plot(dat(:,1),dat(:,2),'.')
                        hold on;
                        plot(rsiz ,rdat,'.')
                        plot(rsiz ,rdat+leeway,'g')
                        plot(rsiz ,rdat-leeway,'g')
                        plot(rsiz(i:j),fitdat,'m')
                        hold off
                        title([fname,"..PLR=",num2str(plr)])
                        
                    end
                end
            end
        end
    end
%     pause

    figure(4)
    if bestplr>0
    semilogx(DT+randn*DT*0.05,bestplr,'.','Color',[1 1 1]*0.7)
    hold on
    q=q+1;
    end

    figure(3)
    subplot(321)
    semilogx(DT,bestplr,ps)
    hold on
    
    temp1=[temp1 DT];
    temp2=[temp2 bestplr]; 

    subplot(322)
    semilogx(DT,expon,ps)
    hold on
    
    
    if any(strfind(temp(f).folder,'spike')) && ~any(strfind(temp(f).folder,'non-awake')) && ~any(strfind(temp(f).folder,'vitro'))
        
        if(any(strfind(fname,'linbin'))); bintype=1; else; bintype=2; end
        NN=str2num(fname(tind-5:tind-1));
%         cind=round(100*(log10(NN)+0)/3);  if cind>100; cind=100; end
        if NN<60; bintype=3; end
        subplot(323)
%         semilogx(DT,bestplr,'+','color',ncmap(cind,:))
        semilogx(DT,bestplr,'o','color',bincmap(bintype,:))
        hold on
    
        
        subplot(324)
%         semilogx(DT,expon,'+','color',ncmap(cind,:))
        semilogx(DT,expon,'o','color',bincmap(bintype,:))
        hold on
    else
        subplot(323)
        semilogx(DT,bestplr,'.','color',[1 1 1]*0.7)
        hold on
        
        subplot(324)
        semilogx(DT,expon,'.','color',[1 1 1]*0.7)
        hold on
    end


    if any(strfind(temp(f).folder,'LFP')) || any(strfind(temp(f).folder,'EEG')) || any(strfind(temp(f).folder,'MEG'))
        
        pcol=[0 0 1];
        if(any(strfind(temp(f).folder,'vitro'))); ps='.'; pcol=[1 1 1]*0.7; end
        if(any(strfind(temp(f).folder,'non-awake'))); ps='.'; pcol=[1 1 1]*0.7; end
       
        subplot(325)
        semilogx(DT,bestplr,ps,'color',pcol)
        hold on
        
        
    else
        subplot(325)
        semilogx(DT,bestplr,'.','color',[1 1 1]*0.7)
        hold on
        
    end
    
    if any(strfind(temp(f).folder,'spikes')) && any(strfind(temp(f).folder,'vitro')) 
        
%         if(any(strfind(temp(f).folder,'vitro'))); bintype=1; else; bintype=2; end
        bintype=1;
       
        subplot(326)
        semilogx(DT,bestplr,ps,'color',bincmap(bintype,:))
        hold on
        
        temp1=[temp1 DT];
        temp2=[temp2 bestplr]; 
        
    else
        subplot(326)
        semilogx(DT,bestplr,'.','color',[1 1 1]*0.7)
        hold on
       
    end
    
end

figure(3)
subplot(321)
hold off
xlabel('\Delta T')
ylabel('PLR')

subplot(322)
hold off
xlabel('\Delta T')
ylabel('exponent \tau')

subplot(323)
hold off
xlabel('\Delta T')
ylabel('PLR')
title('awake spikes')

subplot(324)
hold off
xlabel('\Delta T')
ylabel('exponent \tau')

subplot(325)
hold off
xlabel('\Delta T')
ylabel('PLR')
title('awake LFP, EEG, MEG')

subplot(326)
hold off
xlabel('\Delta T')
ylabel('PLR')
title('in vitro spikes')

figure(5)
semilogx(temp1,temp2,'.')
[rho,pval]=corr(temp1',temp2','type','spearman')