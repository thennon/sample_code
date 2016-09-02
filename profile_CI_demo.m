%%%% A script to calculate the average profile of respiration
%%%% with 95% error bars, and plot up the results.
%%%%
%%%% this script requires MATLAB statistics toolbox and m_map (free online)



clear all; close all;

%% Load Data / generate pseudo-data for demo purposes

% load(pathname/respiration_data.mat) %data with N profiles of respiration (p = pressure, ... 
% ... R = respiration rate) and longitude/latitude where data is collected

%Generate pseudo-data for demo
N = 20; %number of profiles used
dz = 25; % spacing of 

p_vec = 0:dz:2000; %vector of pressures
p_mat = repmat(p_vec,N,1)'; %matrix of pressures for multiple profiles

R = 5*exp(-p_mat/150)+randn(size(p_mat)); %make artificial respiration profiles that ~ match real world patterns

lon = -150+randn(N,1)*10; % Assume data come from off of Hawaii
lat = 20 + randn(N,1)*4;

%% Calculate average, and 95% CI, assuming Gaussian distribution


for i = 1:length(p_vec) % for loop at each pressure level
    R_p = R(i,:); % respiration rates at pressure (i). I define R_p for clarity, but is not needed explicitly in next few lines.
    R_av(i) = nanmean(R_p); %Calculate the average, removing invalid results.
    
    idn  = find(isnan(R_p)==0); % find indices of valid data
    [mu,sig,muci,sigci] = normfit(R_p(idn)); % <-- REQUIRES STATISTICS TOOLBOX FROM MATLAB
 
    R_CI(i,:) = muci;
end

%%
figure(1); clf 
    fill([R_CI(:,1); flipud(R_CI(:,2))],[p_vec fliplr(p_vec)],[1 1 1]*0.5,'edgecolor','none') %CI intervals = shaded gray
    hold on 
    plot(R_av,p_vec,'k','linewidth',3) % mean = thick black
    set(gca,'ydir','reverse','ylim',[0 500],'fontsize',16) %reverse y-axis (so ocean surface is on top), and zoom to top 500 dbar
    xlabel('Respiration [mol C kg^{-3}]','fontweight','bold')
    ylabel('Pressure [dbar]','fontweight','bold')
    title('Respiration near Hawaii','fontsize',24)
    grid on
    
    % Add in small globe for reference on spatial location <--- REQUIRES 'm_map' ROUTINES
     axes('position',[.7 .03 .30 .30]) %Define size of subplot
     m_proj('ortho','lat',20,'long',-150); %make globe centered on hawaii
     m_coast('patch','k'); % plot black coastline
     m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',[-80:40:80]);
     hold on
     m_plot(lon,lat,'k.','markersize',20,'color','m') %plot profile locations in magenta
    
%% Save data and figure
% save('PATHNAME/respiration_hawaii','p_vec','R_av','R_CI','lon','lat')
% saveas('PATHNAME/respiration_hawaii','gcf','jpg')
