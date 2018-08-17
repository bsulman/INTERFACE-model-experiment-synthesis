import pandas
import os.path
import numpy
from pylab import *

datadir = '../../Data/Model-output'
obsdata_dir = '../../Data/Meta-analysis'

if 'CORPSE_loaded' not in globals():
    CORPSE_loaded={}
    MIMICS_loaded={}
    Daycent_loaded={}
    LBL_loaded={}
    MEND_loaded={}

# Functions to load data from each model. Data will be put into a pandas dataframe
# Each call will load an individual experiment based on clay, litter quality, and manipulations specified
def load_CORPSE(clay='medium',litterquality='high',experiment='control',verbose=False):
    clay_strs = {'high':'higherclay','low':'lowclay','medium':'highclay'}
    quality_strs = {'high':'highquality','low':'lowquality'}

    expt_strs = {
                    'control':'control',
                    'warming2':'warming_2',
                    'warming5':'warming_5',
                    'litter_addition_30':'total_addition_30',
                    'litter_addition_100':'total_addition_100',
                    'labile_addition':'labile_addition',
                    'warming_x_litter':'warming_2_total_addition_30',
                    'warming_x_priming':'warming_2_labile_addition',
                    'litter_removal':'litter_removal'
    }

    if '*'.join((experiment,clay,litterquality)) in CORPSE_loaded.keys():
        return CORPSE_loaded['*'.join((experiment,clay,litterquality))]

    filename = 'CORPSE-output-%s_%s-%s.csv'%(clay_strs[clay],quality_strs[litterquality],expt_strs[experiment])
    if verbose:
        print filename
    data = pandas.read_csv(os.path.join(datadir,'CORPSE/2017-03',filename))

    data['total_protectedC']=data['protectedC']
    data['total_unprotectedC']=data['unprotectedC_fast']+data['unprotectedC_slow']+data['unprotectedC_deadmic']
    data['total_C']=data['total_protectedC']+data['total_unprotectedC']+data['microbeC']
    data['total_microbeC']=data['microbeC']
    data['total_litterC']=0.0

    # Resample to monthly resolution
    # Not doing any averaging, just grabbing first value. Probably ok since they don't change fast

    data['month']=data.index/30
    data['year']=data['month']/12.0
    data=data[::30].set_index('month')
    data['CO2flux']=data['CO2'].diff() #kgC/month

    CORPSE_loaded['*'.join((experiment,clay,litterquality))]=data
    return data


def load_MIMICS(clay='medium',litterquality='high',experiment='control',verbose=False):
    clay_strs={'medium':'CLAY20','low':'CLAY5','high':'CLAY70'}
    quality_strs={'high':'LIG16.6','low':'LIG24.4'}

    expt_strs = {
                    'control':'Control',
                    'warming2':'Warm2',
                    'warming5':'Warm5',
                    'litter_addition_30':'1.3xLitter',
                    'litter_addition_100':'2.0xLitter',
                    'labile_addition':'Prime_1.3',
                    'warming_x_litter':'Warm2_1.3xLitter',
                    'warming_x_priming':'Warm2_1.3xprime',
                    'litter_removal':'0xLitter',
                    'litter_removal2':'0xLitter_lowMicTurn'
    }

    if '*'.join((experiment,clay,litterquality)) in MIMICS_loaded.keys():
        return MIMICS_loaded['*'.join((experiment,clay,litterquality))]

    filename='MIMICS_%s_%s_%s.csv'%(expt_strs[experiment],clay_strs[clay],quality_strs[litterquality])
    if verbose:
        print filename

    if expt_strs[experiment] is None:
        print  "Experiment %s for MIMICS does not exist"%experiment
        data = load_MIMICS(clay,litterquality,'control')
        data['CO2flux']=numpy.nan
        data['total_protectedC']=numpy.nan
        data['total_unprotectedC']=numpy.nan
        data['total_microbeC']=numpy.nan
        data['total_C']=numpy.nan
        data['total_litterC']=numpy.nan
        return data

    data = pandas.read_csv(os.path.join(datadir,'MIMICS/data',filename))
    data['month']=data.index/30
    data['year']=data['month']/12.0
    data['CO2flux']=data['CO2']*1e-6*100**2*20  # mgC/cm3 -> kgC/m2

    # Based on conversation with Will, SOMc (chemically protected SOM) is meant to represent particulate organic C
    # This is consistent with how "unprotected C" is defined in this analysis (i.e. not physically inaccessible to microbes)
    data['total_protectedC']=(data['SOMp'])*1e-6*100**2*20  # mgC/cm3 -> kgC/m2
    data['total_unprotectedC']=(data['SOMa']+data['SOMc'])*1e-6*100**2*20  # mgC/cm3 -> kgC/m2
    data['total_litterC']=(data['LITm']+data['LITs'])*1e-6*100**2*20  # mgC/cm3 -> kgC/m2
    data['total_microbeC']=(data['MICr']+data['MICk'])*1e-6*100**2*20  # mgC/cm3 -> kgC/m2
    data['total_C']=data['total_protectedC']+data['total_unprotectedC']+data['total_microbeC']

    MIMICS_loaded['*'.join((experiment,clay,litterquality))]=data[::30]

    return data[::30]

# Daycent has a directory with monthly results. Just use that instead of downsampling.
# Note that Melannie said "resp(1) is an annual accumulator for heterotrophic respiration rather than a monthly value"

def load_Daycent(clay='medium',litterquality='high',experiment='control',verbose=False):
    clay_strs = {'high':'high_protect','medium':'med_protect','low':'low_protect'}
    quality_strs = {'high':'HQ','low':'LQ'}

    if '*'.join((experiment,clay,litterquality)) in Daycent_loaded.keys():
        return Daycent_loaded['*'.join((experiment,clay,litterquality))]

    # These file names need to be filled with clay and litter quality for each experiment
    T_strs = {
                    'control':'20C',
                    'warming2':'22C',
                    'warming5':'25C',
                    'litter_addition_30':'20C',
                    'litter_addition_100':'20C',
                    'labile_addition':'20C',
                    'warming_x_litter':'22C',
                    'warming_x_priming':'22C',
                    'litter_removal':'20C'
    }
    litteraddition_strs = {
                    'control':'500',
                    'warming2':'500',
                    'warming5':'500',
                    'litter_addition_30':'650',
                    'litter_addition_100':'X2',
                    'labile_addition':'L30',
                    'warming_x_litter':'650',
                    'warming_x_priming':'L30',
                    'litter_removal':'0'
    }

    if litteraddition_strs[experiment] is None:
        print  "Experiment %s for Daycent does not exist"%experiment
        data = load_Daycent(clay,litterquality,'control')
        data['CO2flux']=numpy.nan
        data['total_protectedC']=numpy.nan
        data['total_unprotectedC']=numpy.nan
        data['total_microbeC']=numpy.nan
        data['total_C']=numpy.nan
        data['total_litterC']=numpy.nan
        return data

    filename = 'exp_%s%s.%s.%s.lis.csv'%(quality_strs[litterquality],litteraddition_strs[experiment],clay_strs[clay],T_strs[experiment])
    if verbose:
        print filename


    data=pandas.read_csv(os.path.join(datadir,'DAYCENT','Monthly_2017-03-20',filename))

    data.drop('Unnamed: 13',axis=1,inplace=True)
    data.drop(120,inplace=True,axis=0)
    data.drop(121,inplace=True,axis=0)
    data['total_protectedC']=data['som3c'] # Passive pool
    data['total_unprotectedC']=data['somsc']-data['som3c'] # Everything except passive pool
    data['total_litterC']=data['somtc']-data['somsc']
    data['total_C']=data['somtc']
    data['total_microbeC']=numpy.nan

    data['month']=data.index.values
    data.set_index('month',inplace=True)
    data.rename(columns={'time':'year'},inplace=True)

    CO2flux=data['resp(1)'].diff()
    CO2flux[CO2flux<0]=data['resp(1)'][CO2flux<0]
    data['CO2flux']=CO2flux

    data['year'][data['year']<2000]=data['year'][data['year']<2000]+9
    data['year'][data['year']>2000]=data['year'][data['year']>2000]-floor(data['year'].iloc[0])

    Daycent_loaded['*'.join((experiment,clay,litterquality))]=data.loc[0:731]

    return data.loc[0:731]


# LBL (RESOM) is all in one csv file. Just return cross section for this experiment

def load_LBL(clay='medium',litterquality='high',experiment='control',verbose=False):
    clay_strs = {'high':'high clay','low':'low clay','medium':'high clay'}
    quality_strs = {'high':'high quality litter','low':'low quality litter'}

    if '*'.join((experiment,clay,litterquality)) in LBL_loaded.keys():
        return LBL_loaded['*'.join((experiment,clay,litterquality))]

    expt_strs = {
                    'control':'Control',
                    'warming2':'+2C',
                    'warming5':'+5C',
                    'litter_addition_30':'+30% total input',
                    'litter_addition_100':'+100% total input',
                    'labile_addition':'+30% labile inputs',
                    'warming_x_litter':'+2C and +30% total input',
                    'warming_x_priming':'+2C and +30% labile input',
                    'litter_removal':'No input'
    }

    if expt_strs[experiment] is None:
        print  "Experiment %s for LBL does not exist"%experiment
        data = load_LBL(clay,litterquality,'control')
        data['CO2flux']=numpy.nan
        data['total_protectedC']=numpy.nan
        data['total_unprotectedC']=numpy.nan
        data['total_microbeC']=numpy.nan
        data['total_C']=numpy.nan
        data['total_litterC']=numpy.nan
        return data

    # data=pandas.read_csv(os.path.join(datadir,'LBL Output','LBL_model_interface_expanded.csv'))
    data=pandas.read_csv(os.path.join(datadir,'RESOM','LBL_model_Interface_extra_runs_sameProtectedPools.csv'))
    data.drop('Unnamed: 0',axis=1,inplace=True)
    data['month']=data['Year']*12
    data.rename(columns={'Year':'year'},inplace=True)
    data.set_index(['month','Clay content','Litter quality','Experiment'],inplace=True)

    data=data.xs((clay_strs[clay],quality_strs[litterquality],expt_strs[experiment]),level=('Clay content','Litter quality','Experiment'))

    # Medium and high clay content both have name "high clay"
    # Medium comes first, higher clay comes second
    if clay=='medium':
        data=data.iloc[:548]
    if clay=='high':
        data=data.iloc[548:]

    data['total_protectedC']=data['Adsorbed monomers (gC m-2)']+data['Adsorbed enzymes (gC m-2)']
    data['total_unprotectedC']=data['Free monomers  (gC m-2)']+data['Polymers (gC m-2)']+data['Free enzymes (gC m-2)']
    data['total_C']=data['SOC']
    data['total_microbeC']=data['Microbial structural biomass (gC m-2)']+data['Microbial reserve biomass (gC m-2)']
    data['CO2flux']=data['CO2 production (gC m-2 d-1)']*30 # gC/m2/month
    data['total_litterC']=0.0

    LBL_loaded['*'.join((experiment,clay,litterquality))]=data

    return data



def load_MEND(clay='medium',litterquality='high',experiment='control',verbose=False):
    clay_strs = {'high':'sand10%','low':'sand50%','medium':'sand10%'}
    quality_strs = {'high':'lignin17%','low':'lignin25%'}

    expt_strs = {
                    # 'control': no control experiment: Need to just fill with first line value
                    'warming2':'Exp1_%s_%s_2C',
                    'warming5':'Exp2_%s_%s_5C',
                    'litter_addition_30':'Exp3_%s_%s_litter30%%',
                    'litter_addition_100':'Exp4_%s_%s_litter100%%',
                    'labile_addition':'Exp5_%s_%s_labile30%%',
                    'warming_x_litter':'Exp6_%s_%s_2C_litter30%%',
                    'warming_x_priming':'Exp7_%s_%s_2C_labile30%%',
                    'litter_removal':'Exp8_%s_%s_litter0'
    }

    if '*'.join((experiment,clay,litterquality)) in MEND_loaded.keys():
        return MEND_loaded['*'.join((experiment,clay,litterquality))]

    if experiment is not 'control':
        if expt_strs[experiment] is None:
            print  "Experiment %s for MEND does not exist"%experiment
            data = load_MEND(clay,litterquality,'control')
            data['CO2flux']=numpy.nan
            data['total_protectedC']=numpy.nan
            data['total_unprotectedC']=numpy.nan
            data['total_microbeC']=numpy.nan
            data['total_C']=numpy.nan
            data['total_litterC']=numpy.nan
            return data
        filename = 'MEND_' + expt_strs[experiment]%(quality_strs[litterquality],clay_strs[clay]) + '.out'
        if verbose:
            print filename
        data=pandas.read_table(os.path.join(datadir,'MEND',filename),delim_whitespace=True,skiprows=1)

        data['total_protectedC']=data['MOM_C']+data['QOM_C']
        data['total_unprotectedC']=data['POM1_C']+data['POM2_C']+data['DOM_C']+data['ENZ_C']
        data['total_C']=data['TOM_C']
        data['total_microbeC']=data['MB_C']
        ndays=pandas.date_range(datetime.date(2000,12,1),datetime.date(2050,12,31),freq='M').days_in_month
        data['CO2flux']=data['CO2'].diff()/ndays
        data['total_litterC']=0.0

        # Convert months and years into starting index
        # Note that first 10 years of control are not included and assumed to equal the Mon=0 line
        yr=numpy.floor(data['Mon']/100)
        month=data['Mon']-yr*100
        data['year']=yr+10-yr[1]+month/12.0
        data['month']=(yr+10-yr[1]+month/12.0)*12.0

        # Add first 10 years back
        first10=pandas.DataFrame(index=numpy.arange(120),columns=data.columns,dtype=float)
        first10[data.columns]=data.iloc[0].values
        yr=numpy.arange(10).repeat(12)
        month=numpy.arange(len(yr))%12+1
        first10['year']=yr+month/12.0
        first10['month']=(yr+month/12.0)*12.0
        # Gangsheng says baseline CO2 flux is  2.8539E-04 mgC/cm3/hour
        CO2flux_control=2.8539e-04*24
        first10['CO2flux']=CO2flux_control
        data['CO2flux'].loc[1]=CO2flux_control


        MEND_loaded['*'.join((experiment,clay,litterquality))]=pandas.concat((first10,data.iloc[1:])).iloc[0:731]
        return pandas.concat((first10,data.iloc[1:])).iloc[0:731]
    else:
        # Need to read a file with correct clay and litter quality but fill with first val
        data=load_MEND(clay,litterquality,'warming2')
        data[data.columns.drop(['year','month','CO2'])]=data.iloc[0][data.columns.drop(['year','month','CO2'])].values
        data['CO2flux']=2.8539e-04*24
        MEND_loaded['*'.join((experiment,clay,litterquality))]=data.iloc[0:731]
        return data.iloc[0:731]

expt_titles = {
                'control':'Control',
                'warming2':'2 $^\circ$C warming',
                'warming5':'5 $^\circ$C warming',
                'litter_addition_30':'+30% litter',
                'litter_addition_100':'+100% litter',
                'labile_addition':'+30% labile',
                'warming_x_litter':'2 $^\circ$C warming x +30% litter',
                'warming_x_priming':'2 $^\circ$C warming x +30% labile',
                'litter_removal':'No litter inputs'
}

def plot_models(clay='medium',litterquality='high',expt='warming2',variable='total_C',do_log=False,**kwargs):


    if expt not in expt_titles.keys():
        raise ValueError(
        '%s is not a defined experiment. '%expt +
        'Available experiments: %s' %str(expt_titles.keys())
        )

    corpse_control=load_CORPSE(clay,litterquality,'control')
    corpse_expt=load_CORPSE(clay,litterquality,experiment=expt)
    daycent_control=load_Daycent(clay,litterquality,'control')
    daycent_expt=load_Daycent(clay,litterquality,experiment=expt)
    lbl_control=load_LBL(clay,litterquality,'control')
    lbl_expt=load_LBL(clay,litterquality,experiment=expt)
    mend_control=load_MEND(clay,litterquality,'control')
    mend_expt=load_MEND(clay,litterquality,experiment=expt)
    mimics_control=load_MIMICS(clay,litterquality,'control')
    mimics_expt=load_MIMICS(clay,litterquality,experiment=expt)
    if expt=='litter_removal':
        mimics_expt2=load_MIMICS(clay,litterquality,experiment='litter_removal2')


    from pylab import figure,clf,plot,title,xlabel,ylabel,draw,xlim
    if do_log:
        h_CORPSE=plot(corpse_control['year']-10,log(pandas.rolling_mean(corpse_expt[variable],12,center=True)/pandas.rolling_mean(corpse_control[variable],12,center=True)),label='CORPSE',**kwargs)
        h_Daycent=plot(daycent_control['year']-10,log(pandas.rolling_mean(daycent_expt[variable],12,center=True)/pandas.rolling_mean(daycent_control[variable],12,center=True)),label='DAYCENT',**kwargs)
        h_LBL=plot(lbl_control['year']-10,log(pandas.rolling_mean(lbl_expt[variable],12,center=True)/pandas.rolling_mean(lbl_control[variable],12,center=True)),label='RESOM',**kwargs)
        h_MEND=plot(mend_control['year'][6::12]-10,log(mend_expt[variable][6::12]/mend_control[variable][6::12]),label='MEND',**kwargs)
        h_MIMICS=plot(mimics_control['year']-10,log(pandas.rolling_mean(mimics_expt[variable],12,center=True)/pandas.rolling_mean(mimics_control[variable],12,center=True)),label='MIMICS',**kwargs)
        if expt=='litter_removal':
             h_MIMICS2=plot(mimics_control['year']-10,log(pandas.rolling_mean(mimics_expt2[variable],12,center=True)/pandas.rolling_mean(mimics_control[variable],12,center=True)),label='MIMICS (DDM)',ls='--',c=h_MIMICS[0].get_color(),**kwargs)
        ylabel('Log response ratio')
    else:
        h_CORPSE=plot(corpse_control['year']-10,pandas.rolling_mean(corpse_expt[variable],12,center=True)/pandas.rolling_mean(corpse_control[variable],12,center=True)*100-100,label='CORPSE',**kwargs)
        h_Daycent=plot(daycent_control['year']-10,pandas.rolling_mean(daycent_expt[variable],12,center=True)/pandas.rolling_mean(daycent_control[variable],12,center=True)*100-100,label='DAYCENT',**kwargs)
        h_LBL=plot(lbl_control['year']-10,pandas.rolling_mean(lbl_expt[variable],12,center=True)/pandas.rolling_mean(lbl_control[variable],12,center=True)*100-100,label='RESOM',**kwargs)
        h_MEND=plot(mend_control['year'][6::12]-10,mend_expt[variable][6::12]/mend_control[variable][6::12]*100-100,label='MEND',**kwargs)
        h_MIMICS=plot(mimics_control['year']-10,pandas.rolling_mean(mimics_expt[variable],12,center=True)/pandas.rolling_mean(mimics_control[variable],12,center=True)*100-100,label='MIMICS',**kwargs)
        if expt=='litter_removal':
             h_MIMICS2=plot(mimics_control['year']-10,pandas.rolling_mean(mimics_expt2[variable],12,center=True)/pandas.rolling_mean(mimics_control[variable],12,center=True)*100-100,label='MIMICS (DDM)',ls='--',c=h_MIMICS[0].get_color(),**kwargs)
        ylabel('Difference from control (%)')
    plot([-5,55],[0,0],'k--')


    xlabel('Year')
    xlim(-5,50)
    draw()

    return [h_CORPSE[0],h_Daycent[0],h_LBL[0],h_MEND[0],h_MIMICS[0]]


def plot_models_range(expt=['warming2'],variable='total_C',do_log=False,**kwargs):

    if isinstance(expt,str):
        expts=[expt]
    else:
        expts=expt

    corpse_control=[]
    corpse_expt=[]
    daycent_control=[]
    daycent_expt=[]
    lbl_control=[]
    lbl_expt=[]
    mend_control=[]
    mend_expt=[]
    mimics_control=[]
    mimics_expt=[]

    for expt in expts:
        if expt not in expt_titles.keys():
            raise ValueError(
            '%s is not a defined experiment. '%expt +
            'Available experiments: %s' %str(expt_titles.keys())
            )
        if expt=='litter_removal':
            mimics_expt2=[]
        for clay in ['low','medium','high']:
            for litterquality in ['low','high']:
                corpse_control.append(load_CORPSE(clay,litterquality,'control').set_index('year')[variable])
                corpse_expt.append(load_CORPSE(clay,litterquality,experiment=expt).set_index('year')[variable])
                daycent_control.append(load_Daycent(clay,litterquality,'control').set_index('year')[variable])
                daycent_expt.append(load_Daycent(clay,litterquality,experiment=expt).set_index('year')[variable])
                lbl_control.append(load_LBL(clay,litterquality,'control').set_index('year')[variable])
                lbl_expt.append(load_LBL(clay,litterquality,experiment=expt).set_index('year')[variable])
                mend_control.append(load_MEND(clay,litterquality,'control').set_index('year')[variable])
                mend_expt.append(load_MEND(clay,litterquality,experiment=expt).set_index('year')[variable])
                mimics_control.append(load_MIMICS(clay,litterquality,'control').set_index('year')[variable])
                mimics_expt.append(load_MIMICS(clay,litterquality,experiment=expt).set_index('year')[variable])
                if expt=='litter_removal':
                    mimics_expt2.append(load_MIMICS(clay,litterquality,experiment='litter_removal2').set_index('year')[variable])

    t_corpse=corpse_control[0].index.values-10
    corpse_control=column_stack(corpse_control)
    corpse_expt=column_stack(corpse_expt)
    t_daycent=daycent_control[0].index.values-10
    daycent_control=column_stack(daycent_control)
    daycent_expt=column_stack(daycent_expt)
    t_lbl=lbl_control[0].index.values-10
    lbl_control=column_stack(lbl_control)
    lbl_expt=column_stack(lbl_expt)
    t_mend=mend_control[0].index.values[6::12]-10
    mend_control=column_stack(mend_control)
    mend_expt=column_stack(mend_expt)
    t_mimics=mimics_control[0].index.values-10
    mimics_control=column_stack(mimics_control)
    mimics_expt=column_stack(mimics_expt)
    if expt=='litter_removal':
        mimics_expt2=column_stack(mimics_expt2)


    def fillmodelrange(t,control,expt,color,smoothing=12,label=None,**kwargs):
        if do_log:
            h_lines=plot(t,log(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)),c=color,alpha=0.25,lw=0.5,label=None,**kwargs)
            h_meanline=plot(t,log(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).mean(axis=1),c=color,alpha=0.9,lw=2.0,label=label,**kwargs)
            h_fill=fill_between(t,log(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).min(axis=1),
                            log(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).max(axis=1),alpha=0.1,facecolor=color,label=None)
            ylabel('Log response ratio')
        else:
            h_lines=plot(t,(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True))*100-100,c=color,alpha=0.25,lw=0.5,label=None,**kwargs)
            h_lines=plot(t,(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).mean(axis=0)*100-100,c=color,alpha=0.9,lw=2.0,label=label,**kwargs)
            h_fill=fill_between(t,(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).min(axis=1)*100-100,
                            log(pandas.rolling_mean(expt,smoothing,center=True)/pandas.rolling_mean(control,smoothing,center=True)).max(axis=1)*100-100,alpha=0.1,facecolor=color,label=None)
            ylabel('Difference from control (%)')
        return h_meanline[0]

    from pylab import figure,clf,plot,title,xlabel,ylabel,draw,xlim
    h_CORPSE=fillmodelrange(t_corpse,corpse_control,corpse_expt,color='C0',label='CORPSE',**kwargs)
    h_Daycent=fillmodelrange(t_daycent,daycent_control,daycent_expt,color='C1',label='DAYCENT',**kwargs)
    h_LBL=fillmodelrange(t_lbl,lbl_control,lbl_expt,label='RESOM',color='C2',**kwargs)
    h_MEND=fillmodelrange(t_mend,mend_control[6::12],mend_expt[6::12],smoothing=1,color='C3',label='MEND',**kwargs)
    h_MIMICS=fillmodelrange(t_mimics,mimics_control,mimics_expt,color='C4',label='MIMICS',**kwargs)
    if expt=='litter_removal':
         h_MIMICS2=fillmodelrange(t_mimics,mimics_control,mimics_expt2,color='C4',label='MIMICS (DDM)',ls='--',**kwargs)
    ylabel('Log response ratio')

    plot([-5,55],[0,0],'k--',lw=0.5)


    xlabel('Year')
    xlim(-5,50)
    draw()

    return [h_CORPSE,h_Daycent,h_LBL,h_MEND,h_MIMICS]




if __name__=='__main__':
    import sys
    from pylab import show,figure,clf,legend,Line2D,subplot,draw

    save_figs=False
    config_leg_ax=5

    def avgdata(data,var,maxyears=30,justoneyear=False):
        avg=zeros(maxyears)
        if justoneyear:
            for yr in xrange(maxyears):
                xx=(data['year']>=yr)&(data['year']<yr+1)
                avg[yr]=data[var][xx].mean()
        else:
            for yr in xrange(maxyears):
                xx=(data['year']<yr+1)
                avg[yr]=data[var][xx].mean()
        return avg

    def stddata(data,var,maxyears=30,justoneyear=False):
        avg=zeros(maxyears)
        if justoneyear:
            for yr in xrange(maxyears):
                xx=(data['year']>=yr)&(data['year']<yr+1)
                avg[yr]=data[var][xx].std()
        else:
            for yr in xrange(maxyears):
                xx=(data['year']<yr+1)
                avg[yr]=data[var][xx].std()
        return avg


    # log response ratio of total SOC and CO2 flux under warming
    litterquality='high'
    for clay in ['high','medium','low']:
        for var in ['total_C','CO2flux']:
            corpse_control=load_CORPSE(clay,litterquality,'control')
            daycent_control=load_Daycent(clay,litterquality,'control')
            lbl_control=load_LBL(clay,litterquality,'control')
            mend_control=load_MEND(clay,litterquality,'control')
            mimics_control=load_MIMICS(clay,litterquality,'control')

            corpse_warming2=load_CORPSE(clay,litterquality,'warming2')
            daycent_warming2=load_Daycent(clay,litterquality,'warming2')
            lbl_warming2=load_LBL(clay,litterquality,'warming2')
            mend_warming2=load_MEND(clay,litterquality,'warming2')
            mimics_warming2=load_MIMICS(clay,litterquality,'warming2')

            corpse_warming5=load_CORPSE(clay,litterquality,'warming5')
            daycent_warming5=load_Daycent(clay,litterquality,'warming5')
            lbl_warming5=load_LBL(clay,litterquality,'warming5')
            mend_warming5=load_MEND(clay,litterquality,'warming5')
            mimics_warming5=load_MIMICS(clay,litterquality,'warming5')

            corpse_doublelitter=load_CORPSE(clay,litterquality,'litter_addition_100')
            daycent_doublelitter=load_Daycent(clay,litterquality,'litter_addition_100')
            lbl_doublelitter=load_LBL(clay,litterquality,'litter_addition_100')
            mend_doublelitter=load_MEND(clay,litterquality,'litter_addition_100')
            mimics_doublelitter=load_MIMICS(clay,litterquality,'litter_addition_100')

            corpse_nolitter=load_CORPSE(clay,litterquality,'litter_removal')
            daycent_nolitter=load_Daycent(clay,litterquality,'litter_removal')
            lbl_nolitter=load_LBL(clay,litterquality,'litter_removal')
            mend_nolitter=load_MEND(clay,litterquality,'litter_removal')
            mimics_nolitter=load_MIMICS(clay,litterquality,'litter_removal')


            # Save mean and standard deviation of model output to csv files for statistical analysis
            pandas.DataFrame({
                             'corpse_control':avgdata(corpse_control,var),
                             'corpse_2C':avgdata(corpse_warming2,var),
                             'corpse_5C':avgdata(corpse_warming5,var),
                             'corpse_doublelitter':avgdata(corpse_doublelitter,var,justoneyear=True),
                             'corpse_nolitter':avgdata(corpse_nolitter,var,justoneyear=True),
                             'daycent_control':avgdata(daycent_control,var),
                             'daycent_2C':avgdata(daycent_warming2,var),
                             'daycent_5C':avgdata(daycent_warming5,var),
                             'daycent_doublelitter':avgdata(daycent_doublelitter,var,justoneyear=True),
                             'daycent_nolitter':avgdata(daycent_nolitter,var,justoneyear=True),
                             'lbl_control':avgdata(lbl_control,var),
                             'lbl_2C':avgdata(lbl_warming2,var),
                             'lbl_5C':avgdata(lbl_warming5,var),
                             'lbl_doublelitter':avgdata(lbl_doublelitter,var,justoneyear=True),
                             'lbl_nolitter':avgdata(lbl_nolitter,var,justoneyear=True),
                             'mend_control':avgdata(mend_control,var),
                             'mend_2C':avgdata(mend_warming2,var),
                             'mend_5C':avgdata(mend_warming5,var),
                             'mend_doublelitter':avgdata(mend_doublelitter,var,justoneyear=True),
                             'mend_nolitter':avgdata(mend_nolitter,var,justoneyear=True),
                             'mimics_control':avgdata(mimics_control,var),
                             'mimics_2C':avgdata(mimics_warming2,var),
                             'mimics_5C':avgdata(mimics_warming5,var),
                             'mimics_doublelitter':avgdata(mimics_doublelitter,var,justoneyear=True),
                             'mimics_nolitter':avgdata(mimics_nolitter,var,justoneyear=True),
                             }).to_csv('%s_%sclay_mean.csv'%(var,clay))

            pandas.DataFrame({
                             'corpse_control':stddata(corpse_control,var),
                             'corpse_2C':stddata(corpse_warming2,var),
                             'corpse_5C':stddata(corpse_warming5,var),
                             'corpse_doublelitter':stddata(corpse_doublelitter,var,justoneyear=True),
                             'corpse_nolitter':stddata(corpse_nolitter,var,justoneyear=True),
                             'daycent_control':stddata(daycent_control,var),
                             'daycent_2C':stddata(daycent_warming2,var),
                             'daycent_5C':stddata(daycent_warming5,var),
                             'daycent_doublelitter':stddata(daycent_doublelitter,var,justoneyear=True),
                             'daycent_nolitter':stddata(daycent_nolitter,var,justoneyear=True),
                             'lbl_control':stddata(lbl_control,var),
                             'lbl_2C':stddata(lbl_warming2,var),
                             'lbl_5C':stddata(lbl_warming5,var),
                             'lbl_doublelitter':stddata(lbl_doublelitter,var,justoneyear=True),
                             'lbl_nolitter':stddata(lbl_nolitter,var,justoneyear=True),
                             'mend_control':stddata(mend_control,var),
                             'mend_2C':stddata(mend_warming2,var),
                             'mend_5C':stddata(mend_warming5,var),
                             'mend_doublelitter':stddata(mend_doublelitter,var,justoneyear=True),
                             'mend_nolitter':stddata(mend_nolitter,var,justoneyear=True),
                             'mimics_control':stddata(mimics_control,var),
                             'mimics_2C':stddata(mimics_warming2,var),
                             'mimics_5C':stddata(mimics_warming5,var),
                             'mimics_doublelitter':stddata(mimics_doublelitter,var,justoneyear=True),
                             'mimics_nolitter':stddata(mimics_nolitter,var,justoneyear=True),
                             }).to_csv('%s_%sclay_std.csv'%(var,clay))




    # Two 4-panel figures. First is total and CO2flux with data. Second breaks out unprotected and protected
    from pylab import arange,bar,xticks,ylabel,title,legend,text

    def letter_label(plotnum=None,ax=None,xpos=0.05,ypos=1.05):
        if ax is None:
            ax=gca()
        from string import ascii_lowercase
        if plotnum is None:
            plotnum=ax.get_subplotspec().num1
        return text(xpos,ypos,'('+ascii_lowercase[plotnum]+')',transform=ax.transAxes)



    do_log=True
    claycontent='low'
    litqual='high'
    ms=3
    elw=1.0

    dirt_soc_data=pandas.read_csv(os.path.join(obsdata_dir,'new DIRT soc.csv'))
    treatment=dirt_soc_data['treat']
    yr=dirt_soc_data['Study.duration..y.']
    if not do_log:
        cntrl=dirt_soc_data['control.SOC.mean']
        treat=dirt_soc_data['treat.SOC.mean']
        cntrl_sd=dirt_soc_data['control.SOC.sd']/sqrt(dirt_soc_data['control.N'])
        treat_sd=dirt_soc_data['treat.SOC.SD']/sqrt(dirt_soc_data['treat.SOC.N'])
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=dirt_soc_data['mean LRR']
        err=sqrt(dirt_soc_data['study variance'])


    # litteraddfig=figure(9,figsize=(12,4));clf()
    # litteradd_subplots2=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,left=0.6,wspace=0.4,right=0.98,bottom=0.2)
    # litteradd_subplots1=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,right=0.42,left=0.05,wspace=0.4,bottom=0.2)
    litteraddfig=figure(9,figsize=(8,7));clf()
    # litteradd_subplots2=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,left=0.05,wspace=0.4,right=0.98,bottom=0.6)
    # litteradd_subplots1=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,right=0.98,left=0.05,wspace=0.4,top=0.5)

    # warmingfig=figure(8,figsize=(13,3.56));clf()
    # warming_subplots2=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,left=0.6,wspace=0.5,right=0.98,bottom=0.2)
    # warming_subplots1=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,right=0.42,left=0.06,wspace=0.5,bottom=0.2)
    warmingfig=figure(8,figsize=(8,7));clf()
    # warming_subplots2=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,left=0.05,wspace=0.5,right=0.98,bottom=0.6)
    # warming_subplots1=matplotlib.gridspec.GridSpec(ncols=2,nrows=1,right=0.98,left=0.06,wspace=0.5,top=0.5)


    # ax=litteraddfig.add_subplot(litteradd_subplots1[1])
    ax=litteraddfig.add_subplot(222)
    sca(ax)

    var='total_C'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    title('Total soil carbon')
    if do_log:
        ylabel('Log response ratio')
    else:
        ylabel('Difference from control (%)')

    DL=treatment=='DL'
    NI=treatment=='NI'
    errorbar(yr[DL],ratio[DL],yerr=err[DL],ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms,elinewidth=elw,label='Experiments')
    # errorbar(yr[NI],ratio[NI]*100-100,yerr=err[NI]*100,ls='None',marker='o',markerfacecolor='w',markeredgecolor='gray',ecolor='gray')
    xlim(-5,52)
    # text(0.05,0.9,'c',transform=gca().transAxes)
    letter_label(xpos=0.03)


    dirt_resp_data=pandas.read_csv(os.path.join(obsdata_dir,'new DIRT resp.csv'))
    treatment=dirt_resp_data['treat']
    yr=dirt_resp_data['Study.duration..y.']
    if not do_log:
        cntrl=dirt_resp_data['mean.control']
        treat=dirt_resp_data['mean.treat']
        cntrl_sd=dirt_resp_data['sd.control']/sqrt(dirt_resp_data['N.control'])
        treat_sd=dirt_resp_data['sd.treat']/sqrt(dirt_resp_data['N.treat'])
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=dirt_resp_data['LRR']
        err=sqrt(dirt_resp_data['variance'])

    # sca(litteraddfig.add_subplot(litteradd_subplots1[0]))
    sca(litteraddfig.add_subplot(221))
    var='CO2flux'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    title('CO$_2$ flux rate')

    DL=treatment=='DL'
    NI=treatment=='NI'
    errorbar(yr[DL],ratio[DL],yerr=err[DL],ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',label='Experiments',ms=ms,elinewidth=elw)
    # errorbar(yr[NI],ratio[NI]*100-100,yerr=err[NI]*100,ls='None',marker='o',markerfacecolor='w',markeredgecolor='gray',ecolor='gray')
    xlim(-5,52)
    # text(0.05,0.9,'d',transform=gca().transAxes)
    letter_label()
    legend(loc='upper left',bbox_to_anchor=(-.05,-0.2),ncol=6,handletextpad=0.5)

    # leg=legend(loc='best',fontsize='small')
    # leg.get_frame().set_alpha(1.0)

    warm_soc_data=pandas.read_csv(os.path.join(obsdata_dir,'new warm soc.csv'))
    yr=warm_soc_data['Duration.of.Study..y.']
    if not do_log:
        cntrl=warm_soc_data['SOC.mean.control']
        treat=warm_soc_data['SOC.mean.warmed']
        cntrl_sd=warm_soc_data['SOC.st.dev.control']/sqrt(warm_soc_data['SOC.control.N'])
        treat_sd=warm_soc_data['SOC.st.dev.warmed']/sqrt(warm_soc_data['SOC.warm.N'])
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=warm_soc_data['LRR mean']
        err=sqrt(warm_soc_data['Study variance']) #Fix this


    # sca(warmingfig.add_subplot(warming_subplots1[1]))
    sca(warmingfig.add_subplot(222))
    ax=gca()

    var='total_C'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    title('Total soil carbon')
    if do_log:
        ylabel('Log response ratio')
    else:
        ylabel('Difference from control (%)')

    errorbar(yr,ratio,yerr=err,ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms,elinewidth=elw,label='Experiments')
    xlim(-5,50)
    # ylim(-1.1,1.1)
    # ylim(-0.3,0.35)
    # text(0.05,0.9,'a',transform=gca().transAxes)
    letter_label(xpos=0.03)

    #
    from mpl_toolkits.axes_grid import inset_locator
    # iax=inset_locator.inset_axes(gca(),width='35%',height='40%',loc='upper right')
    # # iax=axes([0,0,1,1],label='ax1')
    # # ip1=inset_locator.InsetPosition(ax, [0.55,0.55,0.5,0.5])
    # # iax.set_axes_locator(ip1)
    # iax.tick_params(length=2.0)
    # # plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    # plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    # errorbar(yr,ratio,yerr=err,ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms*0.7,elinewidth=elw*0.7)
    # if do_log:
    #     # ylim(-0.2,0.01)
    #     ylim(-1.1,0.6)
    # else:
    #     ylim(-10,1)
    # xlim(-5,50)
    # xlabel('');ylabel('')
    # xticks([0,25,50],fontsize='small')
    # yticks(fontsize='small')


    # warm_resp_data=pandas.read_csv('new warming resp.csv')
    warm_resp_data=pandas.read_csv(os.path.join(obsdata_dir,'warm resp effects with Tres.csv'))
    warm_resp_data=warm_resp_data[warm_resp_data['Model.Exp']=='Exp']
    yr=warm_resp_data['Duration.of.Study..y.']
    if not do_log:
        cntrl=warm_resp_data['resp.mean.control..umol.m.2.s.1.']
        treat=warm_resp_data['resp.mean.warmed']
        cntrl_sd=warm_resp_data['resp.st.error']
        treat_sd=warm_resp_data['resp.st.error.warmed']
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=warm_resp_data['yi']
        err=sqrt(warm_resp_data['vi'])


    # sca(warmingfig.add_subplot(warming_subplots1[0]))
    sca(warmingfig.add_subplot(221))
    ax1=gca()
    var='CO2flux'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    title('CO$_2$ flux rate')

    errorbar(yr,ratio,yerr=err,ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms,elinewidth=elw,label='Experiments')
    xlim(-5,50)
    # ylim(-.6,0.6)
    # text(0.05,0.9,'b',transform=gca().transAxes)
    letter_label()
    legend(loc='upper left',bbox_to_anchor=(-.05,-0.2),ncol=6,handletextpad=0.5)
    # legend(fontsize='small',loc='upper right')

    # iax1=axes([0,0,1,1],label='ax1')
    # ip1=inset_locator.InsetPosition(ax1, [0.45,0.5,0.25,0.3])
    iax1=inset_locator.inset_axes(gca(),width='48%',height='52%',loc=1)
    # iax1.set_axes_locator(ip1)
    iax1.tick_params(length=2.0)
    # inax=inset_locator.inset_axes(gca(),width='20%',height='50%',loc=9)
    # inax.tick_params(length=2.0)
    # plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    errorbar(yr,ratio,yerr=err,ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms*0.7,elinewidth=elw*0.7,label='Experiments')
    if do_log:
        ylim(-0.35,0.8)
    else:
        ylim(-5,25)
    xlim(-1,5.5)
    xticks([0,2,4])
    xlabel('');ylabel('')
    xticks(fontsize='small')
    yticks(fontsize='small')


######################################
# Protected and unprotected fractions

    do_log=True
    if do_log:
        ytext='Log response ratio'
    else:
        ytext='Difference from control (%)'

    # claycontent='medium'
    # litqual='low'

    # f9=figure(9,figsize=(9,7));clf()
    # sca(litteraddfig.add_subplot(litteradd_subplots2[0]))
    sca(litteraddfig.add_subplot(223))
    var='total_unprotectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    title('Unprotected C')
    ylabel(ytext)
    # text(0.05,1.05,'c',transform=gca().transAxes)
    letter_label(2)
    xlim(-5,52)


    # sca(litteraddfig.add_subplot(litteradd_subplots2[1]))
    sca(litteraddfig.add_subplot(224))
    var='total_protectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_addition_100',variable=var,zorder=10,do_log=do_log)
    title('Protected C')
    # leg=legend(loc='best',fontsize='small')
    # leg.set_zorder(12)
    # text(0.05,1.05,'d',transform=gca().transAxes)
    letter_label(3,xpos=0.03)
    xlim(-5,52)

    # sca(warmingfig.add_subplot(warming_subplots2[0]))
    sca(warmingfig.add_subplot(223))
    var='total_unprotectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    title('Unprotected C')
    ylabel(ytext)
    # text(0.05,1.05,'a',transform=gca().transAxes)
    letter_label(2)

    # ax2=warmingfig.add_subplot(warming_subplots2[1]);sca(ax2)
    ax2=warmingfig.add_subplot(224);sca(ax2)
    var='total_protectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log)
    title('Protected C')
    # text(0.05,1.05,'b',transform=gca().transAxes)
    letter_label(3,xpos=0.03)
    # ylim(-0.2,0.02)
    xlim(-1,50)

    iax2=axes([0,0,1,1],label='ax2')
    ip2=inset_locator.InsetPosition(ax2, [0.5,0.50,0.4,0.3])
    iax2.set_axes_locator(ip2)
    iax2.tick_params(length=2.0)
    # inset_locator.inset_axes(gca(),width='40%',height='40%',bbox_to_anchor=(0.7,0.15,0.4,0.4))
    # plot_models(clay=claycontent,litterquality=litqual,expt='warming2',variable=var,zorder=10,do_log=do_log,axes=iax2)
    plot_models_range(expt=['warming2','warming5'],variable=var,zorder=10,do_log=do_log,axes=iax2)
    if do_log:
        ylim(-0.2,0.02)
    else:
        ylim(-25,2)
    xlim(-1,50)
    xlabel('');ylabel('')
    xticks([0,25,50],fontsize='small');yticks(fontsize='small')

    warmingfig.tight_layout(h_pad=3.0)
    litteraddfig.tight_layout(h_pad=3.0)


    figure(11);clf()
    corpse_control_prot=[]
    daycent_control_prot=[]
    lbl_control_prot=[]
    mend_control_prot=[]
    mimics_control_prot=[]
    corpse_control_unprot=[]
    daycent_control_unprot=[]
    lbl_control_unprot=[]
    mend_control_unprot=[]
    mimics_control_unprot=[]
    corpse_control_litter=[]
    daycent_control_litter=[]
    lbl_control_litter=[]
    mend_control_litter=[]
    mimics_control_litter=[]
    for clay in ['low','medium','high']:
        for litterquality in ['low','high']:
            corpse=load_CORPSE(clay,litterquality,'control')
            corpse_control_prot.append(corpse['total_protectedC'].mean())
            corpse_control_unprot.append(corpse['total_unprotectedC'].mean())
            corpse_control_litter.append(corpse['total_litterC'].mean())
            daycent=load_Daycent(clay,litterquality,'control')
            daycent_control_prot.append(daycent['total_protectedC'].mean())
            daycent_control_unprot.append(daycent['total_unprotectedC'].mean())
            daycent_control_litter.append(daycent['total_litterC'].mean())
            lbl=load_CORPSE(clay,litterquality,'control')
            lbl_control_prot.append(lbl['total_protectedC'].mean())
            lbl_control_unprot.append(lbl['total_unprotectedC'].mean())
            lbl_control_litter.append(lbl['total_litterC'].mean())
            mend=load_MEND(clay,litterquality,'control')
            mend_control_prot.append(mend['total_protectedC'].mean())
            mend_control_unprot.append(mend['total_unprotectedC'].mean())
            mend_control_litter.append(mend['total_litterC'].mean())
            mimics=load_MIMICS(clay,litterquality,'control')
            mimics_control_prot.append(mimics['total_protectedC'].mean())
            mimics_control_unprot.append(mimics['total_unprotectedC'].mean())
            mimics_control_litter.append(mimics['total_litterC'].mean())

    corpse_control_prot=array(corpse_control_prot)
    daycent_control_prot=array(daycent_control_prot)
    lbl_control_prot=array(lbl_control_prot)
    mend_control_prot=array(mend_control_prot)
    mimics_control_prot=array(mimics_control_prot)
    corpse_control_unprot=array(corpse_control_unprot)
    daycent_control_unprot=array(daycent_control_unprot)
    lbl_control_unprot=array(lbl_control_unprot)
    mend_control_unprot=array(mend_control_unprot)
    mimics_control_unprot=array(mimics_control_unprot)
    corpse_control_litter=array(corpse_control_litter)
    daycent_control_litter=array(daycent_control_litter)
    lbl_control_litter=array(lbl_control_litter)
    mend_control_litter=array(mend_control_litter)
    mimics_control_litter=array(mimics_control_litter)

    x=arange(5)
    unprotectedfrac=numpy.array([mean((corpse_control_unprot+corpse_control_litter)/(corpse_control_unprot+corpse_control_litter+corpse_control_prot)),
                                 mean((daycent_control_unprot+daycent_control_litter)/(daycent_control_unprot+daycent_control_litter+daycent_control_prot)),#/1000,
                                mean((lbl_control_unprot+lbl_control_litter)/(lbl_control_unprot+lbl_control_litter+lbl_control_prot)),#/1000,
                                mean((mend_control_unprot+mend_control_litter)/(mend_control_unprot+mend_control_litter+mend_control_prot)),#*20*1e-6*100**2,
                                mean((mimics_control_unprot+mimics_control_litter)/(mimics_control_unprot+mimics_control_litter+mimics_control_prot))])

    unprotectedfrac_std=numpy.array([std((corpse_control_unprot+corpse_control_litter)/(corpse_control_unprot+corpse_control_litter+corpse_control_prot)),
                                 std((daycent_control_unprot+daycent_control_litter)/(daycent_control_unprot+daycent_control_litter+daycent_control_prot)),#/1000,
                                std((lbl_control_unprot+lbl_control_litter)/(lbl_control_unprot+lbl_control_litter+lbl_control_prot)),#/1000,
                                std((mend_control_unprot+mend_control_litter)/(mend_control_unprot+mend_control_litter+mend_control_prot)),#*20*1e-6*100**2,
                                std((mimics_control_unprot+mimics_control_litter)/(mimics_control_unprot+mimics_control_litter+mimics_control_prot))])



    prot_h=bar(x,1-unprotectedfrac,yerr=unprotectedfrac_std,color='k',edgecolor='k',width=0.8,label='Protected',ecolor=(0.3,0.3,0.3),capsize=5.0)
    unprot_h=bar(x,unprotectedfrac,bottom=1-unprotectedfrac,color='w',edgecolor='k',width=0.8,label='Unprotected')

    # text(1,unprotected[1]+protected[1]+litter[1]+0.1,'%s clay, %s litter quality'%(clay.capitalize(),litterquality.capitalize()),rotation=90,ha='center',va='bottom')
    xticks(arange(5),['CORPSE','Daycent','RESOM','MEND','MIMICS'],ha='center',rotation=30)
    yticks(arange(0,1.1,0.2))
    ylabel('Total carbon fraction')
    title('Protected SOC fraction in control simulations')
    legend(handles=[prot_h,unprot_h],loc='best',fontsize='medium')
    ylim(0,1.2)
    # text(0.05,0.95,'e',transform=gca().transAxes)
    # letter_label()

    tight_layout()



    # Litter removal
    f10=figure(10,figsize=(9,7.5));clf()
    f10.add_subplot(221)


    dirt_soc_data=pandas.read_csv(os.path.join(obsdata_dir,'new DIRT soc.csv'))
    treatment=dirt_soc_data['treat']
    yr=dirt_soc_data['Study.duration..y.']
    if not do_log:
        cntrl=dirt_soc_data['control.SOC.mean']
        treat=dirt_soc_data['treat.SOC.mean']
        cntrl_sd=dirt_soc_data['control.SOC.sd']/sqrt(dirt_soc_data['control.N'])
        treat_sd=dirt_soc_data['treat.SOC.SD']/sqrt(dirt_soc_data['treat.SOC.N'])
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=dirt_soc_data['mean LRR']
        err=sqrt(dirt_soc_data['study variance'])


    var='total_C'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    title('Total soil carbon')
    if do_log:
        ylabel('Litter removal\nLog response ratio')
    else:
        ylabel('Litter removal\nDifference from control (%)')

    DL=treatment=='DL'
    NI=treatment=='NI'
    errorbar(yr[NI],ratio[NI],yerr=err[NI],ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',ms=ms,elinewidth=elw)
    # errorbar(yr[NI],ratio[NI]*100-100,yerr=err[NI]*100,ls='None',marker='o',markerfacecolor='w',markeredgecolor='gray',ecolor='gray')
    xlim(-5,55)
    text(0.05,1.05,'a',transform=gca().transAxes)

    dirt_resp_data=pandas.read_csv(os.path.join(obsdata_dir,'new DIRT resp.csv'))
    treatment=dirt_resp_data['treat']
    yr=dirt_resp_data['Study.duration..y.']
    if not do_log:
        cntrl=dirt_resp_data['mean.control']
        treat=dirt_resp_data['mean.treat']
        cntrl_sd=dirt_resp_data['sd.control']/sqrt(dirt_resp_data['N.control'])
        treat_sd=dirt_resp_data['sd.treat']/sqrt(dirt_resp_data['N.treat'])
        ratio=treat/cntrl
        err=sqrt((cntrl_sd/cntrl)**2+(treat_sd/treat)**2)*ratio*100
        ratio=ratio*100-100
    else:
        ratio=dirt_resp_data['LRR']
        err=sqrt(dirt_resp_data['variance'])

    f10.add_subplot(222)
    var='CO2flux'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    title('CO$_2$ flux rate')
    legend(fontsize='small')
    text(0.05,1.05,'b',transform=gca().transAxes)

    DL=treatment=='DL'
    NI=treatment=='NI'
    errorbar(yr[NI],ratio[NI],yerr=err[NI],ls='None',marker='o',markerfacecolor='gray',markeredgecolor='gray',ecolor='gray',label='Experiments',ms=ms,elinewidth=elw)
    # errorbar(yr[NI],ratio[NI]*100-100,yerr=err[NI]*100,ls='None',marker='o',markerfacecolor='w',markeredgecolor='gray',ecolor='gray')
    xlim(-5,55)


    f10.add_subplot(223)
    var='total_unprotectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    title('Unprotected C')
    ylabel('Litter removal\n%s'%ytext)
    text(0.05,1.05,'c',transform=gca().transAxes)

    ax8=f10.add_subplot(224)
    var='total_protectedC'
    # h1=plot_models(clay=claycontent,litterquality=litqual,expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    h1=plot_models_range(expt='litter_removal',variable=var,zorder=10,do_log=do_log)
    title('Protected C')
    text(0.05,1.05,'d',transform=gca().transAxes)

    iax8=axes([0,0,1,1],label='ax8')
    # ip8=inset_locator.InsetPosition(ax8, [0.55,0.25,0.4,0.4])
    ip8=inset_locator.InsetPosition(ax8, [0.1,0.1,0.3,0.3])
    iax8.set_axes_locator(ip8)
    iax8.tick_params(length=2.0)
    # inset_locator.inset_axes(gca(),width='40%',height='40%',bbox_to_anchor=(0.7,0.15,0.4,0.4))
    # plot_models(clay=claycontent,litterquality=litqual,expt='litter_removal',variable=var,zorder=10,do_log=do_log,axes=iax8)
    plot_models_range(expt='litter_removal',variable=var,zorder=10,do_log=do_log,axes=iax8)
    if do_log:
        ylim(-0.2,0.02)
    else:
        ylim(-25,2)
    xlim(-1,50)
    xlabel('');ylabel('')
    xticks([0,25,50],fontsize='small');yticks(fontsize='small')

    tight_layout()
    draw()


    show()
