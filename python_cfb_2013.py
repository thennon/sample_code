# This is a script to assemble data from the raw format and calculate team rankings based on
# game score differential and expected differential. 


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
fgm    = game_csv['Field Goal Made'].values
xp1    = game_csv['Off XP Kick Made'].values
xp2    = game_csv['Off 2XP Made'].values
dxp2   = game_csv['Def 2XP Made'].values
safe   = game_csv['Safety'].values

# add all the points together for game score (oddly not included in original data)
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

#%% Define weeks of CFB (Tuesday-Monday, generally)

startdate = date.toordinal(datetime.date(2013,8,26)) #first day of first week
enddate = date.toordinal(datetime.date(2013,9,2)) #last day of first week
bowlseason = date.toordinal(datetime.date(2013,12,9)) # start of bowl season
L = 125. #total # of teams


frame = np.zeros((weeks,2)) #frame will have start date and end date for each week

for i in range(0,weeks):
    a = i
    frame[i,0] = startdate+7*(i)
    frame[i,1] = enddate+7*(i)
    if startdate+7*(i-1)>= bowlseason: # if bowl season, use all dates afterward (+100)
        frame[i,0] = bowlseason
        frame[i,1] = bowlseason+100
        weeks = i+1
        break

#%% TEAMS & Initialize team score

teams = team_csv['Team Code'].values
team_name = team_csv['Name'].values
team_name[L] = 'FCS' #add in generic FCS team (I'm only doing FBS here)
teams = teams[0:L]

#initialize vectors
teamscore = np.zeros([L+1,1])
opponents = np.ones([L+1,weeks])*np.NaN
weekscore = np.ones([L+1,weeks])*np.NaN


#%% Come up with score weighting 
#here I weight scores so that the first 10 points are most important, points 10-20 are half as important, 20-30 are a quarter, etc. 

frac = (1.0/2.0); #each 10 points the weight diminshes by half (tuneable)       

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
scoreweight2 = -np.flipud(scoreweight1) # add in negative scores as well for when team loses

sc = np.arange(-100,100.1,.1) #What the score is
scoreweight = np.concatenate((scoreweight2,[0],scoreweight1)) #the weight of the score


#%% Meat of the code
# - Essentially, start the teams at zero values, and record score differences for entire season. Weight scores according to above, and save the average as 'teamscore'. 
# - Then, run the season again, using new team values (now non-zero) to see how actual game score difference compares to expected game score difference (from above teamscores).
# - Repeat this process X times until teamscores converge to steady values. 

#initialize more matricies
teamscore_weeks = np.NaN*np.ones((L+1,weeks*iterations))
temp_teamav = np.zeros((L+1,1))


for I in range(0,iterations): #for each season iteriation (we'll repeat 10x)
    for i in range(0,weeks): #for each week of the season (1-16ish)
        week = np.arange(frame[i,0],frame[i,1]+1,1)
        idw = np.where((week[0]<=gameday) & (gameday <= week[-1]))
        
        code = gamecode[idw[0]]
        tcode = teamcode[idw[0]]
        pts = points[idw[0]]
        
        lookup = np.unique(code) #each game played that week
        
        for jj in range(0,len(lookup)): # for game in the week
            id = np.where(code==lookup[jj]) #game index
                        
            team0 = tcode[id[0][0]] 
            tid0 = np.where(team0 == teams) #team 0 ID
            team1 = tcode[id[0][1]]
            tid1 = np.where(team1 == teams) #team 1 ID
            
            opponents[tid0,i] = team1 # who is the opponent?
            opponents[tid1,i] = team0
            
            s0 = float(shape(tid0)[1]) #is there a defined opponent?
            s1 = float(shape(tid1)[1])
            
            if s0 < 1.: #if opponent is not defined, it is an FCS opponent
                tid0 = len(teams)
                teamscore[tid0] = scorefcs #give generic FCS a team score
                opponents[tid1,i] = -999 # and ID
                temp_teamav[tid0] = scorefcs
            if s1 < 1.:
                tid1 = len(teams)
                teamscore[tid1] = scorefcs
                opponents[tid0,i] = -999
                temp_teamav[tid1] = scorefcs

            pts0 = pts[id[0][0]] #how many points did team 0 score?
            pts1 = pts[id[0][1]]           
            
            act_team0 = pts0-pts1 #what is the difference in the scores?
            act_team1 = pts1-pts0
                
            exp_team0 = temp_teamav[tid0] # what did we expect to happen?         
            exp_team1 = temp_teamav[tid1]
            
            diff0 = act_team0 + exp_team1 # add actual results to expected results 
            diff1 = act_team1 + exp_team0
            #the above lines say that if a team beats a really bad team by a lot, they won't score incredibly well
            # but if a team beats a great team by a lot, they'll score really well.
            # and if a team loses to a great team by a field goal, they'll still score pretty well. 
                  
            tscore0 = scoreweight[find_nearest(sc,act_team0)]+exp_team1 #weight the actual/expected difference by score weight. 
            tscore1 = scoreweight[find_nearest(sc,act_team1)]+exp_team0
            
            if tid0 == L+1:
                tscore1 = tscore1 - fcs_penalty # can choose to penalize for playing fcs teams (currently set to 0)
            elif tid1 == L+1:
                tscore0 = tscore1 - fcs_penalty
            
            if pts0 > pts1:
                tscore0 = tscore0+win_bonus #give a slight bonus for actually winning the game
                tscore1 = tscore1-win_bonus
            elif pts1 > pts0:
                tscore0 = tscore0-win_bonus
                tscore1 = tscore1+win_bonus
            
            
            weekscore[tid0,i] = tscore0 #store each score 
            weekscore[tid1,i] = tscore1
            
            teamscore_weeks[tid0,i+I*weeks] = tscore0 #store each score for every season iteration (see how scores evolve with iterations)
            teamscore_weeks[tid1,i+I*weeks] = tscore1
        
    id1 = np.where(weekscore == 0)
    weekscore[id1] = np.NaN
    
    for jj in range(0,len(temp_teamav)): #take average score for each team
        tw = weekscore[jj,:]
        idgood = np.where(np.isnan(tw)==False) #don't use nans
        temp_teamav[jj] = np.mean(tw[idgood[0]])
        
teamscore = temp_teamav

#%%

teamscore[-1] = scorefcs #add in FCS to list (just for comparison)


teamorder = np.flipud(np.argsort(teamscore,axis=0)) #create an ordered vector (top=best)


#%% Calculate SOS
oppscore = np.NaN*np.zeros([len(teams),weeks]) 
SOS_av = np.NaN*np.zeros([len(teams)+1,1])

for ii in range(0,len(teams)): # for each team
    opp = opponents[ii,:] # the team IDs for opponents
    for jj in range(0,len(opp)):
        if (np.isnan(opp[jj]) == False) & (opp[jj] > 0): #if non-nan and non-FCS, use the score
            oppscore[ii,jj] = teamscore[np.where(teams==opp[jj])[0][0]]
        elif np.isnan(opp[jj]) == True: # nans beget nans
            oppscore[ii,jj] = np.NaN
        elif opp[jj] < 0 : #if fcs use fcs score
            oppscore[ii,jj] = scorefcs
    idgood = np.where(np.isnan(oppscore[ii,:])==False) #don't use nans in average
    SOS_av[ii] = np.mean(oppscore[ii,idgood[0]]) #average opponent ranking (strength of schedule rating)
               
SOS_order = np.squeeze(SOS_av[teamorder])
#%% make a clean list for viewing data (team rank, teamscore, and SOS)

list_teams = team_name[teamorder]
scorelist = squeeze(teamscore[teamorder])
rank = np.arange(1,L+2,1)
SOSorder = SOS_av[teamorder]

list_tscore = ["%.2f" % v for v in scorelist]
list_rank = ["%.0f" % v for v in rank]
list_SOS = ["%.2f" % v for v in SOSorder]
  
for ii in range(0,len(list_tscore)):
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


#%% Make a bar plot of top 25 teams, and include SOS
n25 = 25; 
pos = arange(n25)+.5

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
