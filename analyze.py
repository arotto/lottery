import gzip, os
from numpy import *
from pylab import *
from scipy.stats import zscore, pearsonr, linregress, scoreatpercentile
from pandas import DataFrame, date_range, Series, datetools, tseries, rolling_mean, rolling_std, ols, concat, HDFStore, merge, read_csv, read_pickle
from pandas.io.data import get_data_yahoo
#import statsmodels.api as sm
import urllib2, csv, re, cPickle, calendar, time
from datetime import datetime
import statsmodels.api as sm
import statsmodels.formula.api as smf

def pd_zscore(x): return (x - mean(x)) / std(x)

DATE_RANGES = {2011:  date_range(datetime(2011, 1, 1), datetime(2011, 12,31)),
               2012: date_range(datetime(2012, 1, 1), datetime(2012, 12,31)),
               'BOTH': date_range(datetime(2011, 1, 1), datetime(2012, 12,31))}
ANALYSIS_YEAR = 2012
analysis_range = DATE_RANGES[ANALYSIS_YEAR]
SPORTS_ALPHA = 0.1
WEATHER_ALPHA = 0.1
OTHER_CITY_LIST = OTHER_CITY_LIST = ['dallas', 'chicago', 'la',  'philadelphia', 'boston', 'sfbay']

lott = DataFrame.from_csv('data/'+str(ANALYSIS_YEAR)+'_lottery_sales.csv', index_col='SalesDate')

all_zips = lott.ZIP.unique()

demographics = read_pickle('data/demographics.dat')
demographics['z_ses'] = zscore(demographics.ses)

zoning = read_pickle('data/zoning.dat')
sports = read_pickle('data/sports/sports_'+str(ANALYSIS_YEAR)+'_alpha_'+str(SPORTS_ALPHA)+'.dat')

other_sports = {}
for other_city in OTHER_CITY_LIST:
    other_sports[other_city] = read_pickle('data/sports/'+other_city+'_sports_'+str(ANALYSIS_YEAR)+'_alpha_'+str(SPORTS_ALPHA)+'.dat')

[all_other_sports_sum_pe, all_other_sports_sum_pe_with_na] = [zeros(len(analysis_range)),zeros(len(analysis_range))]
for other_city in OTHER_CITY_LIST:
    all_other_sports_sum_pe += other_sports[other_city]['sum_pe']
    all_other_sports_sum_pe_with_na += other_sports[other_city]['sum_pe_with_na']

other_sports['all'] = DataFrame({'sum_pe':Series( all_other_sports_sum_pe, index=analysis_range),
                                 'sum_pe_with_na':Series( all_other_sports_sum_pe_with_na, index=analysis_range)})

weather = DataFrame.from_csv('data/weather/weather.csv', index_col='EST')
weather['PrecipitationIn'] = weather.PrecipitationIn.replace('T', '0')
weather['PrecipitationIn'] = weather.PrecipitationIn.astype(float)
weather['temp_exp_avg'] = Series(index=weather.index) 
weather['temp_pe'] = Series(index=weather.index) 
avg_temp = weather.ix[datetime(2010,12,31)]['MeanTemperatureF']
for date in weather.index[1:]:
    pe = (weather.ix[date]['MeanTemperatureF'] - avg_temp)
    avg_temp += WEATHER_ALPHA*pe
    weather.set_value(date, 'temp_exp_avg', avg_temp)
    weather.set_value(date, 'temp_pe', pe)

solar = DataFrame.from_csv('data/weather/solar_irradiance_' +str(ANALYSIS_YEAR) + '.csv')
solar['date'] = solar.index
solar['time'] = map(lambda x:time.strptime(x, "%H:%M").tm_hour, solar['Time (HH:MM)'])
solar['ghi'] = solar['GHI (W/m^2)']
solar['dni'] = solar['DNI (W/m^2))'] 

ghi = solar[solar.ghi != 0].groupby('date').mean()['ghi'].dropna()
dni = solar[logical_and(~isnan(solar.dni), solar.dni != 0)].groupby('date').mean()['dni']

solar_df = DataFrame({'ghi':ghi[logical_and(ghi.index >= analysis_range[0], ghi.index <= analysis_range[-1])],
                      'dni':dni[logical_and(dni.index >= analysis_range[0], dni.index <= analysis_range[-1])], 
                      'temp':weather.MeanTemperatureF[logical_and(weather.index >= analysis_range[0],
                                                                  weather.index <= analysis_range[-1])]})

for month in range(12):
    solar_df[datetools.MONTHS[month]] = map(int, (solar_df.index.month-1) == month)

irradiance = DataFrame(index = analysis_range)
irradiance['dni'] = solar_df.dni
irradiance['ghi'] = solar_df.ghi

irradiance['ghi_exp_avg'] = Series(index=irradiance.index) 
irradiance['ghi_pe'] = Series(index=irradiance.index) 
irradiance['dni_exp_avg'] = Series(index=irradiance.index) 
irradiance['dni_pe'] = Series(index=irradiance.index) 

avg_ghi = mean(irradiance['ghi'])
avg_dni = mean(irradiance['dni'])
for date in irradiance.index[1:]:
    pe = (irradiance.ix[date]['ghi'] - avg_ghi)
    avg_ghi += WEATHER_ALPHA*pe
    irradiance.set_value(date, 'ghi_exp_avg', avg_ghi)
    irradiance.set_value(date, 'ghi_pe', pe)

    if(not isnan(irradiance.ix[date]['dni'])):
        dni_pe = (irradiance.ix[date]['dni'] - avg_dni)
        avg_dni += WEATHER_ALPHA*dni_pe
        irradiance.set_value(date, 'dni_exp_avg', avg_dni)
        irradiance.set_value(date, 'dni_pe', dni_pe)
    else:
        irradiance.set_value(date, 'dni_exp_avg', avg_dni)
        irradiance.set_value(date, 'dni_pe', nan)

irradiance.to_pickle('data/weather/irradiance_'+str(ANALYSIS_YEAR)+'_alpha_'+str(WEATHER_ALPHA)+'.dat')

other_irradiance = {}
for other_city in OTHER_CITY_LIST:
    if(other_city in ['boston', 'philadelphia']): continue
    other_weather = DataFrame.from_csv('data/weather/'+other_city+'_weather.csv', index_col='EST')
    other_solar = DataFrame.from_csv('data/weather/'+other_city+'_solar_irradiance_'+str(ANALYSIS_YEAR)+'.csv')
    other_solar['date'] = other_solar.index
    other_solar['time'] = map(lambda x:time.strptime(x, "%H:%M").tm_hour, other_solar['Time (HH:MM)'])
    other_solar['ghi'] = other_solar['GHI (W/m^2)']
    other_solar['dni'] = other_solar['DNI (W/m^2))'] 

    other_solar = other_solar[other_solar.index.year==ANALYSIS_YEAR]
    
    other_ghi = other_solar[other_solar.ghi != 0].groupby('date').mean()['ghi']
    other_dni = other_solar[other_solar.dni != 0].groupby('date').mean()['dni']
    other_solar_df = DataFrame({'ghi': other_ghi[logical_and(other_ghi.index >= analysis_range[0], other_ghi.index <= analysis_range[-1])],
                                'dni':other_dni[logical_and(other_dni.index >= analysis_range[0], other_dni.index <= analysis_range[-1])],
                                'temp':other_weather.MeanTemperatureF[logical_and(other_weather.index >= analysis_range[0],
                                                                                  other_weather.index <= analysis_range[-1])]})

    other_ghi_solar_regress_obj = smf.ols(formula = 'ghi ~ 1 + temp', data=other_solar_df) 
    other_ghi_solar_res = other_ghi_solar_regress_obj.fit()
    other_dni_solar_regress_obj = smf.ols(formula = 'dni ~ 1 + temp', data=other_solar_df) 
    other_dni_solar_res = other_dni_solar_regress_obj.fit()

    other_irradiance[other_city]['dni_exp_avg'] = Series(index=other_irradiance[other_city].index) 
    other_irradiance[other_city]['dni_pe'] = Series(index=other_irradiance[other_city].index) 
    avg_dni = mean(other_irradiance[other_city]['dni'])
    for date in other_irradiance[other_city].index[1:]:    
        if(not isnan(other_irradiance[other_city].ix[date]['dni'])):
            dni_pe = (other_irradiance[other_city].ix[date]['dni'] - avg_dni)
            avg_dni += WEATHER_ALPHA*dni_pe
            other_irradiance[other_city].set_value(date, 'dni_exp_avg', avg_dni)
            other_irradiance[other_city].set_value(date, 'dni_pe', dni_pe)
        else:
            other_irradiance[other_city].set_value(date, 'dni_exp_avg', avg_dni)
            other_irradiance[other_city].set_value(date, 'dni_pe', nan)

all_other_ghi = Series(index=analysis_range, data= zeros(len(analysis_range)))
all_other_dni = Series(index=analysis_range, data= zeros(len(analysis_range)))
for other_city in OTHER_CITY_LIST:
    if(other_city in ['boston', 'philadelphia']): continue
    all_other_ghi += other_irradiance[other_city]['ghi']
    all_other_dni += other_irradiance[other_city]['dni']

all_other_ghi = all_other_ghi / (len(OTHER_CITY_LIST)-2)
all_other_dni = all_other_dni / (len(OTHER_CITY_LIST)-2)



other_irradiance['all'] = DataFrame({'dni':all_other_dni, 'ghi':all_other_ghi})

demographics = merge(demographics, zoning, on='ZIP')
demographics.index = demographics.ZIP.astype(int)

hurricanes = {  'irene':[datetime(2011,8,27), datetime(2011,8,28), datetime(2011,8,29)],
                'sandy': date_range(datetime(2012,10,29),datetime(2012,11,1)) }
storm_regressor = Series(data=zeros(len(DATE_RANGES['BOTH'])),index=DATE_RANGES['BOTH'])
hurricane_regressor = Series(data=zeros(len(DATE_RANGES['BOTH'])),index=DATE_RANGES['BOTH'])

for name, dates in hurricanes.iteritems():
    for date in dates:
        hurricane_regressor.ix[date] = 1

temp_weather = weather[weather.index.year != 2010]
storm_regressor = Series(data=zeros(len(DATE_RANGES['BOTH'])),index=DATE_RANGES['BOTH'])
storm_regressor.ix[ temp_weather[logical_and(map(lambda x:'Snow' in str(temp_weather.Events.ix[x]), temp_weather.index), 
                                        temp_weather.MeanVisibilityMiles < 5)].index.tolist() ] = 1

days_of_week_regressors = DataFrame(index=analysis_range)
for day_of_week in range(7):
    days_of_week_regressors[datetools.DAYS[day_of_week]] = map(int, analysis_range.dayofweek==day_of_week)

month_of_year_regressors = DataFrame(index=analysis_range)
for month_of_year in range(1,13):
    month_of_year_regressors[datetools.MONTHS[month_of_year-1]] = map(int, analysis_range.month==month_of_year)

week_of_year_regressors = DataFrame(index=analysis_range)
for week_of_year in range(1,53):
    week_of_year_regressors['WEEK_'+str(week_of_year)] = map(int, analysis_range.week==week_of_year)


csv_file = open('data/holidays.csv', 'rU')
reader = csv.DictReader(csv_file)
holidays = {}
for row in reader:
    holidays[datetools.to_datetime(row['date'], errors='raise')] = row['holiday']+'_'+str(datetools.to_datetime(row['date'], errors='raise').year)[-2:]

holiday_regressor = Series(index=DATE_RANGES['BOTH'], data=zeros(len(DATE_RANGES['BOTH'])))
for date in holiday_regressor.index:
    if(date in holidays.keys()): holiday_regressor[date] = 1

ind_holiday_regressors = DataFrame(index=analysis_range)
for holiday_date, name in holidays.iteritems():
    if(holiday_date in analysis_range): 
        if(not '*' in name):
            holiday_name = name.replace(' ', '').replace('\'','').replace('*','2').replace('.', '').upper()
            ind_holiday_regressors[holiday_name] = Series( map(int, analysis_range == holiday_date),index=analysis_range)

for holiday_date, name in holidays.iteritems():
    if(holiday_date in analysis_range): 
        if('*' in name):
            holiday_name = name.replace(' ', '').replace('\'','').replace('*','').replace('.', '').upper()
            ind_holiday_regressors[holiday_name].ix[holiday_date] = 1

payday_regressors = DataFrame(index=analysis_range)
payday_regressors['FIRST_OF_MONTH'] = zeros(len(analysis_range))
payday_regressors['FIFTEENTH_OF_MONTH'] = zeros(len(analysis_range))

for date in analysis_range:
    if(date.day == 1):
        if((not datetools.isBusinessDay(date)) or sum(ind_holiday_regressors.ix[date]) > 0):
            for i in range(5):
                if((date-i)  not in analysis_range):   break
                if(datetools.isBusinessDay(date-i)):
                    payday_regressors['FIRST_OF_MONTH'].ix[date-i]  = 1
                    break
        else:
            payday_regressors['FIRST_OF_MONTH'].ix[date]  = 1
    elif(date.day == 15):
        if((not datetools.isBusinessDay(date)) or sum(ind_holiday_regressors.ix[date]) > 0):
            for i in range(5):
                if((date-i)  not in analysis_range): break
                if(datetools.isBusinessDay(date-i)):
                    payday_regressors['FIFTEENTH_OF_MONTH'].ix[date-i]  = 1
                    break
        else:
            payday_regressors['FIFTEENTH_OF_MONTH'].ix[date]  = 1

lott['composite'] = (lott['Take5'] + lott['Win4Day'] + lott['Win4Eve'] +lott['QuickDraw'] + \
                         lott['Pick10'] +lott['NumbersDay'] + lott['NumbersEve'] ) #/ lott['RetailerCount']

lott['composite_per_capita'] = zeros(lott.shape[0])
grouped = lott.groupby('ZIP')
for zip, zip_lott in grouped:
    lott.ix[lott.ZIP==zip,'composite_per_capita'] = lott[lott.ZIP==zip]['composite'] / demographics.ix[zip]['population_over_18']


low_residential_zips = demographics[demographics.prop_residential < .1].index
for exclude_zip in low_residential_zips:
    lott = lott[(lott.ZIP != exclude_zip)]
print 'excluded ', len(low_residential_zips), 'low density'


if(ANALYSIS_YEAR == 2012):
    for zip in lott.ZIP.unique():
        zip_df = lott[lott.ZIP == zip]
        pre_storm_mean = mean(zip_df[:datetime(2012, 10,27)].composite_per_capita)

        storm_mean = mean(zip_df[datetime(2012, 11,1):datetime(2012, 12,31)].composite_per_capita)
        if(storm_mean/pre_storm_mean < .75):
            print 'DROPPING FOR SANDY:', zip, demographics.ix[zip]['desc'], demographics.ix[zip]['population'], storm_mean/pre_storm_mean
            lott = lott[lott.ZIP != zip]

        

print len(unique(lott.ZIP)), 'ZIPs total'

regress_df = DataFrame() 

lott['composite_per_capita'] = zeros(lott.shape[0])
grouped = lott.groupby('ZIP')
for zip, zip_lott in grouped:
    print '***', zip, demographics.ix[zip]['desc'], demographics.ix[zip]['population']

    lott.ix[lott.ZIP==zip,'composite_per_capita'] = lott[lott.ZIP==zip]['composite'] / demographics.ix[zip]['population_over_18']
    zip_df = DataFrame(columns = [], index=zip_lott.index) 
    zip_df['purchase'] = zip_lott['composite'] / demographics.ix[zip]['population_over_18'] 
    zip_df['log_purchase'] = log(zip_lott['composite'] / demographics.ix[zip]['population_over_18'] )
    zip_df['purchase_per_income'] = (zip_lott['composite'] / demographics.ix[zip]['population_over_18'] ) / demographics.ix[zip]['income'] 
    zip_df['POWERBALL'] = pd_zscore(log(jackpot.powerball))
    zip_df['MEGAMILLION'] = pd_zscore(log(jackpot.megamillion))
    zip_df['NUM_EVENTS'] = pd_zscore(sports.num_events.shift(1))
    zip_df['temp_pe'] = pd_zscore(weather.temp_pe) #pd_zscore( (rolling_mean(weather.MeanTemperatureF,1) - rolling_mean(weather.MeanTemperatureF,20))/ rolling_mean(weather.MeanTemperatureF,20) ) 
    zip_df['z_ses'] = demographics.ix[zip]['z_ses']
    zip_df['ZIP'] = repeat(zip, len(zip_lott)) #len(zip_lott)-20)
    zip_df['sports_composite'] = pd_zscore(sports.sports_composite.shift(1))#sports.sum_wins#(rolling_mean(sports.sum_wins, 5))
    zip_df['sports_wins'] = pd_zscore(sports.sum_wins.shift(1))#sports.sum_wins#(rolling_mean(sports.sum_wins, 5))
    zip_df['sports_losses'] = pd_zscore(sports.sum_losses.shift(1))#sports.sum_wins#(rolling_mean(sports.sum_wins, 5))
    zip_df['sports_p_win'] = (sports.p_win.shift(1))
    zip_df['sports_p_loss'] = (sports.p_loss.shift(1))
    zip_df['streak_composite'] = (sports.sum_streak.shift(1)) #(rolling_mean(sports.sum_streak, 5))
    zip_df['log_streak_composite'] = log(sports.sum_streak.shift(1)+1) #(rolling_mean(sports.sum_streak, 5))
    zip_df['loss_streak_composite'] = pd_zscore(sports.sum_loss_streak.shift(1)) #(rolling_mean(sports.sum_streak, 5))
    zip_df['streak_composite_lag_0'] = sports.sum_streak #(rolling_mean(sports.sum_streak, 5))
    zip_df['sports_sum_pe'] = pd_zscore(sports.sum_pe.shift(1))
    zip_df['sports_sum_pe_with_na'] = pd_zscore(sports.sum_pe_with_na.shift(1))
    zip_df['ghi_pe'] = pd_zscore(irradiance.ghi_pe)
    zip_df['dni_pe'] = pd_zscore(irradiance.dni_pe)

    if(ANALYSIS_YEAR == 2011): teams = ['rangers', 'rangers','knicks', 'mets', 'yankees', 'giants', 'jets']
    else:  teams = ['rangers', 'rangers','knicks', 'mets', 'yankees', 'giants', 'jets', 'nets']

    for team in teams:
        zip_df[team+'_wins'] = sports[team+'_wins']
        zip_df[team+'_streak'] = sports[team+'_streak']

    for other_city in OTHER_CITY_LIST + ['all']:
        zip_df[ {True: 'all_other', 
                 False:other_city}[other_city=='all'] +'_sum_pe'] = pd_zscore(other_sports[other_city].sum_pe.shift(1)).copy()
        zip_df[ {True: 'all_other', 
                 False:other_city}[other_city=='all'] +'_sum_pe_with_na'] = pd_zscore(other_sports[other_city].sum_pe_with_na.shift(1)).copy()
        if(other_city in ['boston', 'philadelphia']): continue
        if(other_city != 'all'):
            zip_df[other_city+'_dni_pe'] = pd_zscore(other_irradiance[other_city].dni_pe).copy()
        

    zip_df['STORM'] = storm_regressor
    zip_df['HURRICANE'] = hurricane_regressor
    zip_df['HOLIDAYS'] = holiday_regressor    
    zip_df['BAD_WEATHER'] = storm_regressor + hurricane_regressor

    for col in ind_holiday_regressors.columns: zip_df[col.split('_')[0]] = ind_holiday_regressors[col]
    for col in payday_regressors.columns:    zip_df[col] = payday_regressors[col]
    for col in days_of_week_regressors.columns:    zip_df[col] = days_of_week_regressors[col]
    for col in month_of_year_regressors.columns:    zip_df[col] = month_of_year_regressors[col]
#    for col in week_of_year_regressors.columns:    zip_df[col] = week_of_year_regressors[col]

    if((sum(isinf(zip_df['log_purchase'])) > 0) or (sum(isnan(zip_df['log_purchase'])) > 0)):
        print '\t missing' ,sum(isinf(zip_df['log_purchase'])) + sum(isnan(zip_df['log_purchase'])), 'observations'

    if(ANALYSIS_YEAR == 2012):
        zip_df = concat([ zip_df.ix[date_range(datetime(2012,1,1),datetime(2012,10,27))],
                          zip_df.ix[date_range(datetime(2012,11,4),datetime(2012,12,31))]])

    regress_df = concat([regress_df, zip_df])  


regress_df = regress_df[~isinf(regress_df['log_purchase'])]
regress_df = regress_df[~isnan(regress_df['log_purchase'])]

regress_df.to_csv('regress/regress_input_'+str(ANALYSIS_YEAR)+'.csv')


f_in = open('regress/regress_input_'+str(ANALYSIS_YEAR)+'.csv', 'rb')
f_out = gzip.open('regress/regress_input_'+str(ANALYSIS_YEAR)+'.csv.gz', 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.system('rm regress/regress_input_'+str(ANALYSIS_YEAR)+'.csv')

os.system('scp -P 757 regress/regress_input_'+str(ANALYSIS_YEAR)+'.csv.gz rotto@impatience.cns.nyu.edu:~/dawlab/lottery/regress/')


mean_lott = lott.groupby(lambda x:x.date).aggregate(mean)
mean_lott.index = analysis_range
'''
mean_lott['composite'] = (mean_lott['Take5'] + mean_lott['Win4Day'] + mean_lott['Win4Eve'] +mean_lott['QuickDraw'] + mean_lott['Pick10'] +  mean_lott['NumbersDay'] + mean_lott['NumbersEve'] )
#mean_lott = mean_lott.drop(['Lotto', 'MegaMillions', 'Megaplier', 'SweetMillion', 'PowerBall', 'PowerPlay'],1)
mean_lott['sports_composite'] = sports.sum_wins
for day, regressor in days_of_week_regressors.iteritems():   mean_lott[day] = regressor
for day, regressor in payday_regressors.iteritems():    mean_lott[day] = regressor
mean_lott['HOLIDAYS'] = holiday_regressor
'''

