clear all; close all;

%% Create some pseudo-data for demo purposes 



dt = 1/24; %sampling frequency 1/hr
t = dt:dt:30;


y1 = 5*cos(t*2*pi*2.0); % Semi-diurnal (2/day) signal
y2 = 2*cos(t*2*pi*1.0); % small Diurnal (1/day) signal
y3 = 10*cos(t*2*pi*0.2); % large Low-frequency (0.2/day) signal

Y = y1+y2+y3+randn(size(t))*3; %superimpose signals and and some moderate noise

figure(1);clf
plot(t,Y,'k','linewidth',3)
grid on
xlabel('Time [days]','fontweight','bold')
ylabel('Signal [m]','fontweight','bold')
set(gca,'fontsize',14)
title('Sample time-series','fontsize',20)

%% Start calculating the power-spectrum

%remove trendline for powerspectrum
coef = polyfit(t,Y,1); % 1 = linear fit
Y_fit = coef(1)*t+coef(2);
Y_detrend = Y-Y_fit;


L = length(t);
k = (2*pi/L)*[0:L/2-1 -L/2:-1]; ks = fftshift(k); %set up frequency grid
    f1 = ks*24/(2*pi);

    
    
   Sgt = (1/24)*fft(Y_detrend);       
          powspec = fftshift(Sgt.*conj(Sgt))/(1/24*L);

    idpos = find(f1>0);
    f_pos = f1(idpos);
    S_pos = powspec(idpos)*2; %Since the powerspectrum is symmetric for a real Y, we simplify by looking at just the posti
 
figure(2); clf
loglog(f_pos,S_pos,'k','linewidth',2)    
    


ind = find(f_pos > 2.0); ind2 = find(f_pos <=2);
f0 = f_pos(ind); y0 = S_pos(ind);

clear f2 y2
for i = 1:length(f0)/3
    f2(i) = 1/3*(f0(i*3-2) + f0(i*3-1) + f0(i*3));
    y2(i) = 1/3*(y0(i*3-2) + y0(i*3-1) + y0(i*3));
end

ff = [f_pos(ind2) f2]; Syyf = [S_pos(ind2) y2];

hold on
loglog(ff,Syyf,'r')

a = 1;

            dfnew = 3*2*a;

            lo1 = chi2inv(0.975,2*a); hi1 = chi2inv(0.025,2*a);
            hi1 = 2*a/hi1; lo1 = 2*a/lo1;
            
            lo2 = chi2inv(0.975,3*2*a); hi2 = chi2inv(0.025,3*2*a);
            hi2 = dfnew/hi2; lo2 = dfnew/lo2;
 
            hold on
      line1 = .35*ones(size(f_pos));
        line2 = zeros(size(f_pos));
        line2(ind2) = line1(1)*hi1/lo1;
        line2(ind) = line1(1)*hi2/lo2;
        
        loglog(f_pos,line1,'k',f_pos,line2,'k','linewidth',2)
        pos1 = (line1(1)/lo1);
        loglog(.11,pos1,'k.','markersize',20)
        loglog([.11 .11],[line1(1) line2(1)],'k','linewidth',2)
        text(.12,pos1*.95,'95% CI','fontsize',12,'fontweight','bold')


 
 
 %% Average over 3 frequency bins for high frequencies
            ind = find(fpos > 4.1); ind2 = find(fpos<=4.1); % Frequencies in [cpd]
            f0 = fpos(ind); y0 = Syy(ind);

            for i = 1:length(f0)/3
                  f2(i) = 1/3*(f0(i*3-2) + f0(i*3-1) + f0(i*3));
                  y2(i) = 1/3*(y0(i*3-2) + y0(i*3-1) + y0(i*3));
            end

            ff = [fpos(ind2) f2]; Syyf = [Syy(ind2) y2];
        
        loglog(ff,Syyf,'k','linewidth',2)
        %loglog(ff,Syyf,'r*','markersize',5)
    

        
        [a,b] = size(Syy_ind);
        cd /Users/thennon/Argo/Testing
            dfnew = 3*2*a;

            lo1 = chi2inv(0.975,2*a); hi1 = chi2inv(0.025,2*a);
            hi1 = 2*a/hi1; lo1 = 2*a/lo1;
            
            lo2 = chi2inv(0.975,3*2*a); hi2 = chi2inv(0.025,3*2*a);
            hi2 = dfnew/hi2; lo2 = dfnew/lo2;
        cd /Users/thennon/Argo/Testing/NewSpec
        
        line1 = .35*ones(size(fpos));
        line2 = zeros(size(fpos));
        line2(ind2) = line1(1)*hi1/lo1;
        line2(ind) = line1(1)*hi2/lo2;
        
        loglog(fpos,line1,'k',fpos,line2,'k','linewidth',2)
        pos1 = (line1(1)/lo1);
        loglog(.11,pos1,'k.','markersize',20)
        loglog([.11 .11],[line1(1) line2(1)],'k','linewidth',2)
        text(.12,pos1*.95,'95% CI','fontsize',12,'fontweight','bold')
        %text(1,line1(1)*.8,['df = ',num2str(2*a)],'fontsize',12,'fontweight','bold')
        %text(8,line1(1)*.8,['df = ',num2str(dfnew)],'fontsize',9,'fontweight','bold')
        
        text(9,10^1.7,['Lat: ',num2str(avlatR),' +/- ' ,num2str(slatR)],'fontsize',12,'fontweight','bold')    
        text(9,10^1.5,['Lon: ',num2str(avlonR),' +/- ', num2str(slonR)],'fontsize',12,'fontweight','bold')
        text(9,10^1.3,['Record: ',num2str(record),' Days'],'fontsize',12,'fontweight','bold')
        text(9,10^1.1,['Segments: ',num2str(del/24),' Days'],'fontsize',12,'fontweight','bold')
        
        legend('S_{av}','GM','location','southeast')




