%this code creates example time series that demonstrate temporal scale
%invariance
clear

T=360000; %one hour
dt=0.01; %model time step assumption in seconds
tvec=(1:T)*dt;

col=[linspace(177/255,0,7)' linspace(45/255,188/255,7)' linspace(143/255,161/255,7)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model 1 - critical        
x=zeros(1,T);
phi1=0.5;
phi2=0.5;
for i=3:T
    x(i) = phi1*x(i-1) + phi2*x(i-2) + randn;   
end

%model space
figure(20)
plot(phi1,phi2,'.k','MarkerSize',4)
hold on

%time series stack
figure(21)
subx=x;
subx=(subx-mean(subx))/max(subx);
lsx=length(subx);
q=1;
while lsx>360
    tvec=(1:lsx)/lsx;
    mask1=tvec<1/3;
    mask2=tvec>2/3;
    if q<7
        plot(tvec(mask1),subx(mask1)-q*2,'Color',col(q,:))
        hold on
        plot(tvec(mask2),subx(mask2)-q*2,'Color',col(q,:))
        plot(tvec(~mask1 & ~mask2),subx(~mask1 & ~mask2)-q*2,'Color',col(q+1,:))
    else
        plot(tvec,subx-q*2,'Color',col(q,:))
        hold on
    end
    
    subx=subx(round(lsx/3):round(2*lsx/3));
    subx=(subx-mean(subx));
    subx=subx/max(subx);
    lsx=length(subx);
    q=q+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model 2 - stable, near critical        
x=zeros(1,T);
phi1=0.49;
phi2=0.49;
for i=3:T
    x(i) = phi1*x(i-1) + phi2*x(i-2) + randn;   
end

%model space
figure(20)
plot(phi1,phi2,'.k','MarkerSize',4)

%time series stack
figure(21)
subx=x;
subx=(subx-mean(subx))/max(subx);
lsx=length(subx);
q=1;
while lsx>360
    tvec=(1:lsx)/lsx;
    mask1=tvec<1/3;
    mask2=tvec>2/3;
    if q<7
        plot(tvec(mask1)+1.1,subx(mask1)-q*2,'Color',col(q,:))
        hold on
        plot(tvec(mask2)+1.1,subx(mask2)-q*2,'Color',col(q,:))
        plot(tvec(~mask1 & ~mask2)+1.1,subx(~mask1 & ~mask2)-q*2,'Color',col(q+1,:))
    else
        plot(tvec+1.1,subx-q*2,'Color',col(q,:))
        hold on
    end
    subx=subx(round(lsx/3):round(2*lsx/3));
    subx=(subx-mean(subx));
    subx=subx/max(subx);
    lsx=length(subx);
    q=q+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model 3 - near white noise        
x=zeros(1,T);
phi1=0.1;
phi2=0.1;
for i=3:T
    x(i) = phi1*x(i-1) + phi2*x(i-2) + randn;   
end

%model space
figure(20)
plot(phi1,phi2,'.k','MarkerSize',4)
hold off

%time series stack
figure(21)
subx=x;
subx=(subx-mean(subx))/max(subx);
lsx=length(subx);
q=1;
while lsx>360
    tvec=(1:lsx)/lsx;
    mask1=tvec<1/3;
    mask2=tvec>2/3;
    if q<7
        plot(tvec(mask1)+2.2,subx(mask1)-q*2,'Color',col(q,:))
        hold on
        plot(tvec(mask2)+2.2,subx(mask2)-q*2,'Color',col(q,:))
        plot(tvec(~mask1 & ~mask2)+2.2,subx(~mask1 & ~mask2)-q*2,'Color',col(q+1,:))
    else
        plot(tvec+2.2,subx-q*2,'Color',col(q,:))
        hold on
    end
    subx=subx(round(lsx/3):round(2*lsx/3));
    subx=(subx-mean(subx));
    subx=subx/max(subx);
    lsx=length(subx);
    q=q+1;
end
hold off