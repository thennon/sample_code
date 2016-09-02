# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:05:13 2016

@author: thennon
"""
#%% LOAD DATA

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import date

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

game_csv = pd.read_csv('team-game-statistics.csv')
team_csv = pd.read_csv('team.csv')
conf_csv = pd.read_csv('conference.csv')

#%% DEFINE PARAMETERS 
year     = 2013.;
scorefcs = -12.;    
weeks    = 20;
num = range(0,30)
win_bonus = 1.;
fcs_penalty = 0.;
iterations = 10;

#%% GAME STATISTICS

teamcode = game_csv['Team Code'].values

rushtd = game_csv['Rush TD'].values
passtd = game_csv['Pass TD'].values
korttd = game_csv['Kickoff Ret TD'].values
ptrttd = game_csv['Punt Ret TD'].values
fmrttd = game_csv['Fum Ret TD'].values
inrttd = game_csv['Int Ret TD'].values
msrttd = game_csv['Misc Ret TD'].values
fgm = game_csv['Field Goal Made'].values
xp1 = game_csv['Off XP Kick Made'].values
xp2 = game_csv['Off 2XP Made'].values
dxp2= game_csv['Def 2XP Made'].values
safe = game_csv['Safety'].values

points = rushtd*6+passtd*6+korttd*6+ptrttd*6+fmrttd*6+inrttd*6+msrttd*6+fgm*3+xp1*1+xp2*2+dxp2*2+safe*2;

#%% FIND DATE OF GAME
gamecode = game_csv['Game Code'].values

gameday = [0]*len(gamecode)

for i in range(0,len(gamecode)):
    gamecodestr = str(gamecode[i])
    gameyear = int(gamecodestr[-8:-4])
    gamemonth = int(gamecodestr[-4:-2])
    gd1 = gamecodestr[-2]
    gd2 = gamecodestr[-1]
    gameday1 = int(gd1+gd2)
    
    gameday[i] = date.toordinal(datetime.date(gameyear,gamemonth,gameday1))

startdate = date.toordinal(datetime.date(2013,8,26))
enddate = date.toordinal(datetime.date(2013,9,2))
bowlseason = date.toordinal(datetime.date(2013,12,9))
L = 125.

#%% Define weeks of CFB (Tuesday-Monday, generally)

frame = np.zeros((weeks,2))

for i in range(0,weeks):
    a = i
    frame[i,0] = startdate+7*(i)
    frame[i,1] = enddate+7*(i)
    if startdate+7*(i-1)>= bowlseason:
        frame[i,0] = bowlseason
        frame[i,1] = bowlseason+100
        weeks = i+1
        break

#%% TEAMS & Initialize team score

teams = team_csv['Team Code'].values
team_name = team_csv['Name'].values
team_name[L] = 'FCS'
teams = teams[0:L]
#teams = np.concatenate((teams,[-999]))

teamscore = np.zeros([L+1,1])
opponents = np.ones([L+1,weeks])*np.NaN
weekscore = np.ones([L+1,weeks])*np.NaN


#%% Come up with score weighting 
#here I weight scores so that the first 10 points are most important, points 10-20 are half as important, 20-30 are a quarter, etc. 

frac = float(2.0/4.0);           

a0 = np.arange(.1,10.1,.1)
a1 = np.arange(.1,10.1,.1)*frac**1+a0[-1]
a2 = np.arange(.1,10.1,.1)*frac**2+a1[-1]
a3 = np.arange(.1,10.1,.1)*frac**3+a2[-1]
a4 = np.arange(.1,10.1,.1)*frac**4+a3[-1]
a5 = np.arange(.1,10.1,.1)*frac**5+a4[-1]
a6 = np.arange(.1,10.1,.1)*frac**6+a5[-1]
a7 = np.arange(.1,10.1,.1)*frac**7+a6[-1]
a8 = np.arange(.1,10.1,.1)*frac**8+a7[-1]
a9 = np.arange(.1,10.1,.1)*frac**9+a8[-1]

scoreweight1 = np.concatenate((a0,a1,a2,a3,a4,a5,a6,a7,a8,a9))
scoreweight2 = -np.flipud(scoreweight1)

scoreweight = np.concatenate((scoreweight2,[0],scoreweight1))
sc = np.arange(-100,100.1,.1)


#%%
teamscore_weeks = np.NaN*np.ones((L+1,weeks*iterations))

temp_weekscore = np.zeros((L+1,weeks))
temp_teamav = np.zeros((L+1,1))

#%% Meat of the code
# - Essentially, start the teams at zero values, and record score differences for entire season. Weight scores according to above, and save the average as 'teamscore'. 
# - Then, run the season again, using new team values (now non-zero) to see how actual game score difference compares to expected game score difference (from above teamscores).
# - Repeat this process X times until teamscores converge to steady values. 

for I in range(0,iterations):
    for i in range(0,weeks):
        week = np.arange(frame[i,0],frame[i,1]+1,1)
        idw = np.where((week[0]<=gameday) & (gameday <= week[-1]))
        
        code = gamecode[idw[0]]
        tcode = teamcode[idw[0]]
        pts = points[idw[0]]
        
        lookup = np.unique(code)
        
        for jj in range(0,len(lookup)):
            id = np.where(code==lookup[jj])
                        
            team0 = tcode[id[0][0]]
            tid0 = np.where(team0 == teams)
            team1 = tcode[id[0][1]]
            tid1 = np.where(team1 == teams)
            
            opponents[tid0,i] = team1
            opponents[tid1,i] = team0
            
            s0 = float(shape(tid0)[1])
            s1 = float(shape(tid1)[1])
            
            if s0 < 1.:
                tid0 = len(teams)
                teamscore[tid0] = scorefcs
                opponents[tid1,i] = -999
                temp_teamav[tid0] = scorefcs
            if s1 < 1.:
                tid1 = len(teams)
                teamscore[tid1] = scorefcs
                opponents[tid0,i] = -999
                temp_teamav[tid1] = scorefcs

            pts0 = pts[id[0][0]]
            pts1 = pts[id[0][1]]           
            
            act_team0 = pts0-pts1
            act_team1 = pts1-pts0
                
            exp_team0 = temp_teamav[tid0]                
            exp_team1 = temp_teamav[tid1]
            
            diff0 = act_team0 + exp_team1
            diff1 = act_team1 + exp_team0
           
            tscore0 = scoreweight[find_nearest(sc,act_team0)]+exp_team1
            tscore1 = scoreweight[find_nearest(sc,act_team1)]+exp_team0
            
            if tid0 == L+1:
                tscore1 = tscore1 - fcs_penalty
            elif tid1 == L+1:
                tscore0 = tscore1 - fcs_penalty
            
            if pts0 > pts1:
                tscore0 = tscore0+win_bonus
                tscore1 = tscore1-win_bonus
            elif pts1 > pts0:
                tscore0 = tscore0-win_bonus
                tscore1 = tscore1+win_bonus
            
            temp_weekscore[tid0,i] = tscore0
            temp_weekscore[tid1,i] = tscore1
            
            weekscore[tid0,i] = tscore0
            weekscore[tid1,i] = tscore1
            
            teamscore_weeks[tid0,i+I*weeks] = tscore0
            teamscore_weeks[tid1,i+I*weeks] = tscore1
        
    id1 = np.where(temp_weekscore == 0)
    temp_weekscore[id1] = np.NaN
    
    for jj in range(0,len(temp_teamav)):
        tw = temp_weekscore[jj,:]
        idgood = np.where(np.isnan(tw)==False)
        temp_teamav[jj] = np.mean(tw[idgood[0]])
        
teamscore = temp_teamav

#%%

teamscore[-1] = scorefcs#np.NaN

x = np.flipud(np.sort(teamscore,axis=0))
teamorder = np.flipud(np.argsort(teamscore,axis=0))




#%% Calculate SOS
oppscore = np.NaN*np.zeros([len(teams),weeks])
SOS_av = np.NaN*np.zeros([len(teams)+1,1])

for ii in range(0,len(teams)):
    opp = opponents[ii,:]
    for jj in range(0,len(opp)):
        if (np.isnan(opp[jj]) == False) & (opp[jj] > 0):
            oppscore[ii,jj] = teamscore[np.where(teams==opp[jj])[0][0]]
        elif np.isnan(opp[jj]) == True:
            oppscore[ii,jj] = np.NaN
        elif opp[jj] < 0 :
            oppscore[ii,jj] = scorefcs
    idgood = np.where(np.isnan(oppscore[ii,:])==False)
    SOS_av[ii] = np.mean(oppscore[ii,idgood[0]])
  
          
#%%
SOS_order = np.squeeze(SOS_av[teamorder])
abc = ["%.2f" % v for v in SOS_av]          


#%% Make a bar plot of top 25 teams
n25 = 25; 

figure(1)
barh(pos,scorelist[0:n25], align='center',height=0.6,color = [0.7,0.7,0.7],edgecolor='none')
barh(pos,SOS_order[0:n25],align='center',height=0.35,color = [0.9,0,0],edgecolor = 'none')

yname = [None]*n25
for ii in range(0,n25):
    yname[ii] = np.squeeze(team_name[teamorder[ii]])+' ('+str(ii+1)+')'


yticks(pos, (yname),fontweight='bold',fontsize=14)
xlabel('TEAM SCORE',fontsize=18)
title('TOP 25  [2013]',fontsize=24)
#grid(True)
plt.gca().invert_yaxis()
xlim([-6,21])
ylim([25.5,-0.5])
matplotlib.rcParams.update({'font.size': 16})

plt.plot([0,0],[-1, 26],color='black',linewidth = 5)

show()

#%% MAKE A LIST FOR VIEWING

list_teams = team_name[teamorder]
scorelist = squeeze(teamscore[teamorder])
rank = np.arange(1,L+2,1)
SOSorder = SOS_av[teamorder]

list_tscore = ["%.2f" % v for v in scorelist]
list_rank = ["%.0f" % v for v in rank]
list_SOS = ["%.2f" % v for v in SOSorder]
  
for ii in range(0,len(twodecimals)):
    str0 = list_tscore[ii]
    str1 = list_rank[ii]
    str2 = list_teams[ii]
    str3 = list_SOS[ii]
    
    list_tscore[ii] = str0.rjust(6)
    list_rank[ii] = str1.ljust(5)
    list_teams[ii] = str2[0].ljust(22)
    list_SOS[ii] = str3.rjust(8)

full_list = transpose([list_rank,squeeze(list_teams),list_tscore,list_SOS])

print('Full List:')
print(full_list)

N = 10
topN = full_list[0:N,:]

print('Top '+str(N)+':')
print(topN)