%% power_spec_demo.m
% Use pseudo data to demonstrate technique used to create a power spectrum
% with error bars

clear all; close all;

%% Create some pseudo-data for demo purposes 

dt = 1/24; %sampling frequency 1/hr
t = dt:dt:30;


y1 = 5*cos(t*2*pi*2.0); % Semi-diurnal (2/day) signal
y2 = 2*cos(t*2*pi*1.0); % small Diurnal (1/day) signal
y3 = 10*cos(t*2*pi*0.2); % large Low-frequency (0.2/day) signal

Y = y1+y2+y3+randn(size(t))*3; %superimpose signals and and some moderate noise
% We will assume that Y has units of meters.

hf = figure(1);clf
set(hf,'position',[50 50 1300 750])

subplot(121)
plot(t,Y,'k','linewidth',3)
grid on
xlabel('Time [days]','fontweight','bold')
ylabel('Signal [m]','fontweight','bold')
set(gca,'fontsize',20)
title('Sample time-series','fontsize',24)


%% Calculate the power spectrum

%remove trendline for powerspectrum
coef = polyfit(t,Y,1); % 1 = linear fit
Y_fit = coef(1)*t+coef(2);
Y_detrend = Y-Y_fit;


L = length(t);
k = (2*pi/L)*[0:L/2-1 -L/2:-1]; ks = fftshift(k); %set up frequency grid
f = ks*24/(2*pi); % Convert to cycles/day

S = (1/24)*fft(Y_detrend);       
SPEC = fftshift(S.*conj(S))/(1/24*L);

id_pos = find(f>0);
f_pos = f(id_pos);
S_pos = SPEC(id_pos)*2; %Since the powers pectrum is symmetric if all 'Y' is real, 
                           %we simplify by looking at just the postive frequency values
 
%% Check parseval's theorem... variance in TS should be equal to variance in spectrum

var_ts = var(Y);
var_spec = sum(SPEC)*(f(2)-f(1));

var_diffs = var_ts-var_spec; %this should be a  pretty small number 
    
%% Let's average over some higher frequency to increase sample size and decrease error...


ind = find(f_pos > 2.0); ind2 = find(f_pos <=2); %split data into >2 and <2 cycles/day
f0 = f_pos(ind); y0 = S_pos(ind);

clear f2 y2
for i = 1:length(f0)/3
    f2(i) = 1/3*(f0(i*3-2) + f0(i*3-1) + f0(i*3)); %take average frequency
    y2(i) = 1/3*(y0(i*3-2) + y0(i*3-1) + y0(i*3)); %take average spectum amplitude
end

% cat the averaged data (>2cpd) with the original data <2cpd.
ff = [f_pos(ind2) f2]; Syyf = [S_pos(ind2) y2];

%% Plot the data
subplot(122)

    loglog(ff,Syyf,'k','linewidth',2); hold on %plot power spectrum in log-scale

    % Now plot 95% CI ranges, one for the data <2cpd with one df, and one for
    % the averaged data (>2cpd) with 3df. a = 1; 

    a = 1; %sample size = 1, for one 30-day segment. Using multiple 30-day segments would increase this number.
    df_0 = 2*a; 
    lo_CI_0 = chi2inv(0.975,2*a); hi_CI = chi2inv(0.025,2*a); %
    hi_CI = 2*a/hi_CI; lo_CI_0 = 2*a/lo_CI_0;

    df_2 = 3*2*a; %degrees of freedom when averaging over frequency bins. 
    lo_CI_2 = chi2inv(0.975,3*2*a); hi_CI_2 = chi2inv(0.025,3*2*a);
    hi_CI_2 = df_2/hi_CI_2; lo_CI_2 = df_2/lo_CI_2;
 
            hold on
            
    % Draw 95% CI range
    setpoint = 1e-4;
    line_bottom = setpoint*ones(size(f_pos));
    line_top = zeros(size(f_pos));
    line_top(ind2) = line_bottom(1)*hi_CI/lo_CI_0; % for <2cpd
    line_top(ind) = line_bottom(1)*hi_CI_2/lo_CI_2; % for >2cpd
        
    loglog(f_pos,line_bottom,'k',f_pos,line_top,'k','linewidth',2)
    pos1 = (line_bottom(1)/lo_CI_0);
    loglog(.11,pos1,'ko','markersize',10,'markerfacecolor',[0 0 0])
    loglog([.11 .11],[line_bottom(1) line_top(1)],'k','linewidth',2)
    text(.12,pos1*.95,'95% CI','fontsize',12,'fontweight','bold')

    set(gca,'xlim',[10^(-1.8) 10^(1.5)],'ylim',[1e-5 1e5],'fontsize',14)
    grid on
    xlabel('Frequency [cpd]','fontsize',20,'fontweight','bold')
    ylabel('Amplitude [m^2 cpd^{-1}]','fontsize',20,'fontweight','bold')
    title('Demo Power Spectrum','fontsize',24)
    
 



