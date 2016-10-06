import pickle
import urllib2
try:    from bs4 import BeautifulSoup
except:    from BeautifulSoup import BeautifulSoup
from pandas import DataFrame, date_range, Series, datetools, tseries, rolling_mean, rolling_std, ols, concat, HDFStore, merge, read_csv, read_pickle
from datetime import datetime


ALPHA = .1
SPORTS_YEAR = 2012

DATE_RANGES = {2011:  date_range(datetime(2011, 1, 1), datetime(2011, 12,31)),
               2012: date_range(datetime(2012, 1, 1), datetime(2012, 12,31))}

analysis_range = DATE_RANGES[SPORTS_YEAR]

TEAM_URLS = ['http://espn.go.com/LEAGUE/team/schedule/_/name/TEAM/year/YEAR/seasontype/2/half/1/',
            'http://espn.go.com/LEAGUE/team/schedule/_/name/TEAM/year/YEAR/seasontype/2/half/2/',
            'http://espn.go.com/LEAGUE/team/schedule/_/name/TEAM/year/YEAR/']

def loadSingleYearSeason(abbrev, league, season_year, dataframe, ignore_fall):
    dataframe['temp'] = zeros(shape(dataframe)[0])
    season = Series(index = dataframe.index)
    for url in TEAM_URLS:
        data_url = url.replace('YEAR', str(season_year)).replace('TEAM', abbrev).replace('LEAGUE',league)
        print data_url
        soup = BeautifulSoup(urllib2.urlopen(data_url).read())
        dates = []
        for row in soup.findAll("tr", {"class" : lambda L: L and (L.startswith('evenrow') or L.startswith('oddrow'))}):
            tds = row('td')
            if(league == 'nfl'):
                try:
                    [date, win, loss, tie] = [tds[1].text, 
                                              {'W':1, 'L':0, 'T':0}[tds[3].li()[0].contents[0]], 
                                              {'W':0, 'L':1, 'T':0}[tds[3].li()[0].contents[0]], 
                                              {'W':0, 'L':0, 'T':1}[tds[3].li()[0].contents[0]]]
                except (IndexError, TypeError):
                    continue
            else:
                try:
                    try:
                        [date, win, loss, tie] = [tds[0].text, 
                                                  {'W':1, 'L':0, 'T':0}[tds[2].li()[0].contents[0]], 
                                                  {'W':0, 'L':1, 'T':0}[tds[2].li()[0].contents[0]], 
                                                  {'W':0, 'L':0, 'T':1}[tds[2].li()[0].contents[0]]]
                    except TypeError:
                        continue
                except IndexError:
                    [date, win, loss, tie] = [tds[0].contents[0], 0, 0, 1]

            try:
                date1 = datetools.to_datetime(date + ' ' +str({True:season_year, False:season_year-1}[ignore_fall]),errors='raise')
            except ValueError:
                print season_year, date, ignore_fall
                continue

            if(ignore_fall):
                if(league == 'nba'):
                    if(date1.month > 7): continue
                if(league == 'nhl'):
                    if(date1.month > 7): continue
            else:
                if(league == 'nba'):
                    if(date1.month < 11): continue
                if(league == 'nhl'):
                    assert(season_year!=2013)
                    if(season_year == 2012):
                        if(date1 < datetime(2011, 10, 6)): continue

            date = datetime({True:season_year, False:season_year-1}[ignore_fall], date1.month, date1.day)
            if(league == 'nfl'):
                season_start = {2011: datetime(2011,9,8), 2012:datetime(2012,9,5)}[season_year]
                if(date < season_start): 
                    continue

            if(dataframe.ix[date]['temp'] == 0):
                if(isnan(dataframe.ix[date][team+'_wins'])):
                    dataframe.ix[date][team+'_wins'] = win
                    dataframe.ix[date][team+'_losses'] = loss
                    dataframe.ix[date][team+'_ties'] = tie
                    dataframe.ix[date][team+'_num_events'] = 1
                else:
                    dataframe.ix[date][team+'_wins'] += win
                    dataframe.ix[date][team+'_losses'] += loss
                    dataframe.ix[date][team+'_ties'] += tie
                    dataframe.ix[date][team+'_num_events'] += 1

            season.ix[datetools.to_datetime(date,errors='raise')] = 1
            dates.append(date)

        try:
            [first_date, last_date] =  [min(dates),max(dates)]
            if((first_date != None) and (last_date != None)):  dataframe['temp'].ix[date_range(first_date, last_date)] = 1
        except ValueError:
            pass

    dataframe = dataframe.drop('temp', 1)
    return {'start':season.index[isnan(season)==False][0],
            'end':season.index[isnan(season)==False][-1]}


TEAMS FROM 5 BOROUGH ANALYSIS, 2011 and 2012
league = {'mets':'mlb', 'yankees':'mlb', 'knicks':'nba', 'nets':'nba',
          'giants':'nfl', 'jets':'nfl', 'rangers':'nhl'}
abbrev = {'mets':'nym', 'yankees':'nyy', 'knicks':'ny', 'nets':'bkn',
          'giants':'nyg', 'jets':'nyj', 'rangers':'nyr'}
if(SPORTS_YEAR == 2011): teams = ['rangers', 'knicks', 'mets', 'yankees', 'giants', 'jets']
else:  teams = ['rangers','knicks', 'mets', 'yankees', 'giants', 'jets', 'nets']



sports = DataFrame(index=analysis_range)
start_end_dates = {}
for team in teams:
    sports[team+'_wins'] = array(nan)
    sports[team+'_losses'] = array(nan)
    sports[team+'_num_events'] = array(nan)
    sports[team+'_ties'] = array(nan)
    sports[team+'_mask'] = zeros(len(analysis_range))
    
    if(league[team] in ['nhl', 'nba']):
        if(team != 'nets'):
            team_start_end_dates = loadSingleYearSeason( abbrev[team], league[team], SPORTS_YEAR ,sports, ignore_fall=True)
            start_end_dates[team] = [team_start_end_dates]
            sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1
        else:
            start_end_dates[team] = []

        if( not ((league[team]  == 'nhl') and ( SPORTS_YEAR == 2012))):
            team_start_end_dates = loadSingleYearSeason(abbrev[team], league[team], SPORTS_YEAR+1,sports, ignore_fall=False)
            start_end_dates[team].append(team_start_end_dates)
            sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1
    else:
        team_start_end_dates = loadSingleYearSeason(abbrev[team], league[team], SPORTS_YEAR ,sports, ignore_fall=True)
        start_end_dates[team] = team_start_end_dates
        sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1

win_cols = filter(lambda x:'wins' in x, sports.columns)
sports['num_events'] = map(lambda x:sum(1- isnan(sports.ix[analysis_range[x]][win_cols].values)), 
                           range(len(analysis_range)))
sports['sum_wins'] = zeros(len(analysis_range))
sports['sum_losses'] = zeros(len(analysis_range))
sports['sports_composite'] = zeros(len(zeros(len(analysis_range))))
sports['composite_mask'] = zeros(len(zeros(len(analysis_range))))

for team in teams:
    sports['sum_wins'] += sports[team+'_wins'].fillna(0)
    sports['sum_losses'] += sports[team+'_losses'].fillna(0)
    sports['sports_composite'] += sports[team+'_wins'].fillna(0) - sports[team+'_losses'].fillna(0)

sports['p_win'] = sports['sum_wins'] / (sports['sum_wins'] + sports['sum_losses']) #sports['num_events']
sports['p_loss'] =sports['sum_losses'] /  (sports['sum_wins'] + sports['sum_losses']) #sports['num_events']

sports['p_win'] = map(lambda x:{True:nan, False:x}[x==inf],sports.p_win)
sports['p_loss'] = map(lambda x:{True:nan, False:x}[x==inf],sports.p_loss)

figure()
for team in teams:
    sports[team+'_pe'] = zeros(len(sports.index))  
    sports[team+'_pe_with_na'] = map(lambda x:nan, range(len(sports.index)))
    avg_rec = Series(index=sports.index)
    avg_rec_with_na = Series(map(lambda x:nan, range(len(sports.index))), index=sports.index)

    if(league[team] in ['nhl', 'nba']):
        for season_num in range(2):
            avg = .5
            if((league[team]  == 'nhl') and ( SPORTS_YEAR == 2012) and (season_num==1)): continue
            if((team  == 'nets') and ( SPORTS_YEAR == 2012) and (season_num==1)): continue

            for date in date_range(start_end_dates[team][season_num]['start'],start_end_dates[team][season_num]['end']):
                if(isnan(sports.ix[date][team + '_wins'])): 
                    avg_rec.set_value(date, avg)
                    sports.set_value(date,team+'_pe', 0)
                    continue

                pe = sports.ix[date][team + '_wins'] - avg
                pe = (sports.ix[date][team + '_wins'] / sports.ix[date][team + '_num_events']) - avg
                sports.set_value(date,team+'_pe', pe)
                sports.set_value(date,team+'_pe_with_na', pe)
                avg += ALPHA*pe
                avg_rec.set_value(date, avg)
                avg_rec_with_na.set_value(date, avg)

        subplot(211)
        avg_rec.plot(label=team)

    else:
        avg_rec = Series(index = sports.index)
        avg = .5

        for date in date_range(start_end_dates[team]['start'],start_end_dates[team]['end']):
            if(isnan(sports.ix[date][team + '_wins'])): 
                avg_rec.set_value(date, avg)
                sports.set_value(date,team+'_pe', 0)
                continue

            pe = (sports.ix[date][team + '_wins'] / sports.ix[date][team + '_num_events']) - avg
            sports.set_value(date,team+'_pe', pe)
            sports.set_value(date,team+'_pe_with_na', pe)
            avg += ALPHA*pe
            avg_rec.set_value(date, avg)
            avg_rec_with_na.set_value(date, avg)

        subplot(211)
        avg_rec.plot(label=team) 

    subplot(212)
    sports[team+'_pe'].plot(label=team)

subplot(211)
legend(prop={'size':11},ncol=2)
subplot(212)
legend(prop={'size':11},ncol=2)
show()


sports['sum_pe_with_na'] = zeros(len(analysis_range))
sports['sum_pe'] = zeros(len(analysis_range))
for team in teams:
    sports['sum_pe'] += sports[team+'_pe'].fillna(0)
    sports['sum_pe_with_na'] += sports[team+'_pe'].fillna(0)

for date in sports.index:
    if(sports.ix[date]['num_events'] == 0): sports.set_value(date, 'sum_pe_with_na', nan)

sports['mean_pe'] = (sports['sum_pe']/sports['num_events']).replace(inf, nan)

pickle.dump(sports, open('data/sports/sports_metro_area_'+str(SPORTS_YEAR)+'_alpha_'+str(ALPHA)+'.dat', 'wb'), True)

other_league = {'miami': {'panthers':'nhl', 'marlins':'mlb', 'dolphins':'nfl', 'heat':'nba'},
                'phoenix': {'coyotes':'nhl', 'diamondbacks':'mlb', 'cardinals':'nfl', 'suns':'nba'},
                'minneapolis':{'wild':'nhl', 'twins':'mlb', 'vikings':'nfl', 'timberwolves':'nba'},
                'chicago':{'bears':'nfl', 'cubs':'mlb', 'white-sox':'mlb', 'bulls':'nba', 'blackhawks':'nhl'},
                'boston': {'bruins':'nhl', 'red-sox':'mlb', 'patriots':'nfl', 'celtics':'nba'},
                'la': {'dodgers':'mlb', 'clippers':'nba','lakers':'nba', 'kings':'nhl'},
                'dallas':{'rangers':'mlb', 'mavericks':'nba', 'cowboys':'nfl', 'stars':'nhl'},
                'denver':{'rockies':'mlb', 'nuggets':'nba', 'broncos':'nfl', 'avalanche':'nhl'},
                'philadelphia':{'phillies':'mlb', '76ers':'nba', 'eagles':'nfl', 'flyers':'nhl'},
                'sfbay':{'49ers':'nfl','raiders':'nfl', 'giants':'mlb', 'athletics':'mlb', 'warriors':'nba', 'sharks':'nhl'},
                'dc':{'redskins':'nfl','nationals':'mlb',  'wizards':'nba', 'capitals':'nhl'},
                'detroit':{'lions':'nfl','tigers':'mlb',  'pistons':'nba', 'redwings':'nhl'}}
                
                
other_abbrev = {'miami': {'panthers':'fla', 'marlins':'mia', 'dolphins':'mia', 'heat':'mia'},
                'phoenix': {'coyotes':'ari', 'diamondbacks':'ari', 'cardinals':'ari', 'suns':'phx'},
                'minneapolis':{'wild':'min', 'twins':'min', 'vikings':'min', 'timberwolves':'min'},
                'chicago': {'bears':'chi', 'cubs':'chc', 'white-sox':'chw', 'bulls':'chi', 'blackhawks':'chi'},
                'boston': {'bruins':'bos', 'red-sox':'bos', 'patriots':'ne', 'celtics':'bos'},
                'la': {'dodgers':'lad', 'clippers':'lac','lakers':'lal', 'kings':'la'},
                'dallas':{'rangers':'tex', 'mavericks':'dal', 'cowboys':'dal', 'stars':'dal'},
                'denver':{'rockies':'col', 'nuggets':'den', 'broncos':'den', 'avalanche':'col'},
                'philadelphia':{'phillies':'phi', '76ers':'phi', 'eagles':'phi', 'flyers':'phi'},
                'sfbay':{'49ers':'sf','raiders':'oak', 'giants':'sf', 'athletics':'oak', 'warriors':'gs', 'sharks':'sj'},
                'dc':{'redskins':'wsh','nationals':'wsh',  'wizards':'wsh', 'capitals':'wsh'},
                'detroit':{'lions':'det','tigers':'det',  'pistons':'det', 'redwings':'det'}}


for city  in ['dallas', 'chicago', 'la',  'philadelphia', 'boston', 'sfbay']:
    other_sports = DataFrame(index=analysis_range)
    erostart_end_dates = {}
    for team in other_abbrev[city].keys():
        other_sports[team+'_wins'] = array(nan)
        other_sports[team+'_losses'] = array(nan)
        other_sports[team+'_ties'] = array(nan)
        other_sports[team+'_num_events'] = array(nan)
        other_sports[team+'_mask'] = zeros(len(analysis_range))

        if(other_league[city][team] in ['nhl', 'nba']):
            team_start_end_dates = loadSingleYearSeason(other_abbrev[city][team], other_league[city][team], SPORTS_YEAR ,other_sports, ignore_fall=True)
            start_end_dates[team] = [team_start_end_dates]
            other_sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1

            if( not ((other_league[city][team]  == 'nhl') and ( SPORTS_YEAR == 2012))):
                team_start_end_dates = loadSingleYearSeason(other_abbrev[city][team], other_league[city][team], SPORTS_YEAR+1,other_sports, ignore_fall=False)
                start_end_dates[team].append(team_start_end_dates)
                other_sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1

        else:
            team_start_end_dates = loadSingleYearSeason(other_abbrev[city][team], other_league[city][team], SPORTS_YEAR ,other_sports, ignore_fall=True)
            start_end_dates[team] = team_start_end_dates
            other_sports[team+'_mask'].ix[team_start_end_dates['start']:team_start_end_dates['end']] = 1

    win_cols = filter(lambda x:'wins' in x, other_sports.columns)

    other_sports['num_events'] = map(lambda x:sum(1- isnan(other_sports.ix[analysis_range[x]][win_cols].values)), 
                                     range(len(analysis_range)))

    other_sports['sum_wins'] = zeros(len(zeros(len(analysis_range))))
    other_sports['sum_losses'] = zeros(len(zeros(len(analysis_range))))
    other_sports['sports_composite'] = zeros(len(zeros(len(analysis_range))))

    for team in other_abbrev[city].keys():
        other_sports['sum_wins'] += other_sports[team+'_wins'].fillna(0)
        other_sports['sum_losses'] += other_sports[team+'_losses'].fillna(0)
        other_sports['sports_composite'] += other_sports[team+'_wins'].fillna(0) - other_sports[team+'_losses'].fillna(0)

    other_sports['p_win'] = other_sports['sum_wins'] / (other_sports['sum_wins'] + other_sports['sum_losses']) 
    other_sports['p_loss'] = other_sports['sum_losses'] /  (other_sports['sum_wins'] + other_sports['sum_losses']) 

    other_sports['p_win'] = map(lambda x:{True:nan, False:x}[x==inf],other_sports.p_win)
    other_sports['p_loss'] = map(lambda x:{True:nan, False:x}[x==inf],other_sports.p_loss)

    other_sports['sum_streak'] = zeros(shape(other_sports)[0])
    other_sports['sum_loss_streak'] = zeros(shape(other_sports)[0])

    for team in other_abbrev[city].keys():
        if(other_league[city][team] in ['nba', 'nhl']): 
            streak_length =0 
            streak_lengths = Series(index= analysis_range,data =zeros(len(analysis_range))-1)
            loss_streak_length =0 
            loss_streak_lengths = Series(index= analysis_range,data =zeros(len(analysis_range))-1)

            for season_num in range(2):
                if((other_league[city][team]  == 'nhl') and ( SPORTS_YEAR == 2012) and (season_num==1)): continue

                for date in date_range(start_end_dates[team][season_num]['start'],start_end_dates[team][season_num]['end']):
                    if(other_sports.ix[date][team+'_wins'] ==1):        streak_length += 1
                    elif(other_sports.ix[date][team+'_losses'] == 1):        streak_length = 0

                    streak_lengths.ix[date] = streak_length

                for date in date_range(start_end_dates[team][season_num]['start'],start_end_dates[team][season_num]['end']):
                    if(other_sports.ix[date][team+'_losses'] ==1):        loss_streak_length += 1
                    elif(other_sports.ix[date][team+'_wins'] == 1):         loss_streak_length = 0
                    loss_streak_lengths.ix[date] = loss_streak_length

        else:
            streak_length =0 
            streak_lengths = Series(index= analysis_range,data =zeros(len(analysis_range))-1)
            loss_streak_length =0 
            loss_streak_lengths = Series(index= analysis_range,data =zeros(len(analysis_range))-1)

            for date in date_range(start_end_dates[team]['start'],start_end_dates[team]['end']):
                if(other_sports.ix[date][team+'_wins'] ==1):      streak_length += 1
                elif(other_sports.ix[date][team+'_losses'] == 1):     streak_length = 0
                streak_lengths.ix[date] = streak_length

            for date in date_range(start_end_dates[team]['start'],start_end_dates[team]['end']):
                if(other_sports.ix[date][team+'_losses'] ==1):        loss_streak_length += 1
                elif(other_sports.ix[date][team+'_wins'] == 1):         loss_streak_length = 0
                loss_streak_lengths.ix[date] = loss_streak_length

        other_sports[team+'_streak'] = map(lambda x:{True:0, False:x}[x==-1],streak_lengths) 
        other_sports[team+'_loss_streak'] = map(lambda x:{True:0, False:x}[x==-1],loss_streak_lengths) 
        other_sports['sum_streak']+= other_sports[team+'_streak']
        other_sports['sum_loss_streak']+= other_sports[team+'_loss_streak']

    other_sports['date'] = other_sports.index

    figure()
    for team in other_abbrev[city].keys():
        other_sports[team+'_pe'] = zeros(len(sports.index)) #array(nan) #zeros(len(sports.index))-1
        if(other_league[city][team] in ['nhl', 'nba']):
            for season_num in range(2):
                avg_rec = []
                avg = .5

                if((other_league[city][team]  == 'nhl') and ( SPORTS_YEAR == 2012) and (season_num==1)): continue
                if((team  == 'nets') and ( SPORTS_YEAR == 2012) and (season_num==1)): continue

                for date in date_range(start_end_dates[team][season_num]['start'],start_end_dates[team][season_num]['end']):
                    if(isnan(other_sports.ix[date][team + '_wins'])): 
                        other_sports.set_value(date,team+'_pe', 0)
                        continue

                    pe = other_sports.ix[date][team + '_wins'] - avg
                    pe = (other_sports.ix[date][team + '_wins'] / other_sports.ix[date][team + '_num_events']) - avg
                    other_sports.set_value(date,team+'_pe', pe)
                    avg += ALPHA*pe
                    avg_rec.append(avg)
                plot(avg_rec, label=team+str(season_num+1))

        else:
            avg_rec = []
            avg = .5

            for date in date_range(start_end_dates[team]['start'],start_end_dates[team]['end']):
                if(isnan(other_sports.ix[date][team + '_wins'])): 
                    other_sports.set_value(date,team+'_pe', 0)
                    continue

                pe = (other_sports.ix[date][team + '_wins'] / other_sports.ix[date][team + '_num_events']) - avg
                other_sports.set_value(date,team+'_pe', pe)
                avg += ALPHA*pe
                avg_rec.append(avg)
            plot(avg_rec,label=team)

    legend()
    title(str(ALPHA))
    show()

    other_sports['sum_pe'] = zeros(len(analysis_range))
    other_sports['sum_pe_with_na'] = zeros(len(analysis_range))

    for team in other_abbrev[city].keys():
        other_sports['sum_pe'] += other_sports[team+'_pe'].fillna(0)
        other_sports['sum_pe_with_na'] += other_sports[team+'_pe'].fillna(0)

    for date in sports.index:
        if(other_sports.ix[date]['num_events'] == 0): other_sports.set_value(date, 'sum_pe_with_na', nan)

    other_sports['mean_pe'] = (other_sports['sum_pe']/other_sports['num_events']).replace(inf, nan)

    pickle.dump(other_sports, open('data/'+city+'_sports_'+str(SPORTS_YEAR)+ '_alpha_'+str(ALPHA)+'.dat', 'wb'), True)


