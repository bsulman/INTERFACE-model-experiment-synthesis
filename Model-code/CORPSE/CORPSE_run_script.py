import CORPSE
from pylab import *
import pandas

# 5% clay
params_lowclay={
    'vmaxref':[1500,50,600], #Relative maximum enzymatic decomp rates
    'Ea':[37e3,54e3,50e3],    # Activation energy
    'kC':[0.01,0.01,0.01],    # Michaelis-Menton parameter
    'gas_diffusion_exp':2.5,  # Determines suppression of decomp at high soil moisture
    'minMicrobeC':1e-3,       #Minimum microbial biomass (fraction of total C)
    'Tmic':0.15,       # Microbial lifetime
    'et':0.5,          # Fraction of turnover not converted to CO2
    'eup':[0.6,0.05,0.6], # Carbon uptake efficiency
    'tProtected':75.0,    # Protected C turnover time (years)
    'protection_rate':[1.0,0.00005,1.0], # Protected carbon formation rate (year-1)
}

# 20% clay
params_highclay=params_lowclay.copy()
params_highclay['protection_rate']=array(params_lowclay['protection_rate'])*CORPSE.prot_clay(20.0)/CORPSE.prot_clay(5.0)
params_higherclay=params_lowclay.copy()
params_higherclay['protection_rate']=array(params_lowclay['protection_rate'])*CORPSE.prot_clay(70.0)/CORPSE.prot_clay(5.0)

# c_lowclay_highquality=CORPSE.soil_carbon_cohort(litterC=[0.0018,0.5535,0.0073],livingMicrobeC=0.025,protectedC=[0.215,0.0,1.14],params=params_lowclay)
# c_highclay_highquality=CORPSE.soil_carbon_cohort(litterC=[0.0018,0.5535,0.0073],livingMicrobeC=0.025,protectedC=[0.4,0.0,2.5],params=params_highclay)
# c_lowclay_lowquality=CORPSE.soil_carbon_cohort(litterC=[0.0022587,0.605026,0.005075],livingMicrobeC=0.012595,protectedC=[0.215,0.0,0.778],params=params_lowclay)
# c_highclay_lowquality=CORPSE.soil_carbon_cohort(litterC=[0.002195,0.5946,0.005646],livingMicrobeC=0.013072,protectedC=[0.408,0.02922,1.8762],params=params_highclay)

c_lowclay_highquality = CORPSE.soil_carbon_cohort(litterC=[0.004849,0.3673,0.006732], protectedC=[0.3637,0.001377,0.5049], livingMicrobeC=0.02342)
c_lowclay_lowquality = CORPSE.soil_carbon_cohort(litterC=[0.002171,0.622,0.004328], protectedC=[0.1628,0.002333,0.3246], livingMicrobeC=0.01188)
c_highclay_highquality = CORPSE.soil_carbon_cohort(litterC=[0.004849,0.3673,0.006732], protectedC=[0.7108,0.002691,0.9867], livingMicrobeC=0.02342)
c_highclay_lowquality = CORPSE.soil_carbon_cohort(litterC=[0.002171,0.622,0.004327], protectedC=[0.5831,0.008352,1.162], livingMicrobeC=0.01188)
c_higherclay_highquality = CORPSE.soil_carbon_cohort(litterC=[0.004849,0.3673,0.006732], protectedC=[0.7108,0.002691,0.9867], livingMicrobeC=0.02342)
c_higherclay_lowquality = CORPSE.soil_carbon_cohort(litterC=[0.002171,0.622,0.004327], protectedC=[0.5831,0.008352,1.162], livingMicrobeC=0.01188)

c_lowclay_highquality.set_params(params_lowclay)
c_lowclay_lowquality.set_params(params_lowclay)
c_highclay_highquality.set_params(params_highclay)
c_highclay_lowquality.set_params(params_highclay)
c_higherclay_highquality.set_params(params_higherclay)
c_higherclay_lowquality.set_params(params_higherclay)

do_spinup=False

configs=[
'lowclay_highquality',
'lowclay_lowquality',
'highclay_highquality',
'highclay_lowquality',
'higherclay_highquality',
'higherclay_lowquality'
]

if do_spinup:
    # Just do control runs to save time
    sims=['control']
else:
    sims=[
        'control',
        'warming_2',
        'warming_5',
        'labile_addition',
        'total_addition_30',
        'total_addition_100',
        'warming_2_total_addition_30',
        'warming_2_labile_addition',
        'litter_removal'
        ]

cohorts=dict()
cohorts['lowclay_highquality']=dict()
cohorts['highclay_highquality']=dict()
cohorts['higherclay_highquality']=dict()
cohorts['lowclay_lowquality']=dict()
cohorts['highclay_lowquality']=dict()
cohorts['higherclay_lowquality']=dict()

for sim in sims:
    cohorts['lowclay_highquality'][sim]=c_lowclay_highquality.copy()
    cohorts['highclay_highquality'][sim]=c_highclay_highquality.copy()
    cohorts['higherclay_highquality'][sim]=c_highclay_highquality.copy()
    cohorts['lowclay_lowquality'][sim]=c_lowclay_lowquality.copy()
    cohorts['highclay_lowquality'][sim]=c_higherclay_lowquality.copy()
    cohorts['higherclay_lowquality'][sim]=c_higherclay_lowquality.copy()


nsteps=365*60


def make_outputs(nsteps):
    out=dict()
    for sim in sims:
        out[sim]=dict()
        out[sim]['unprotectedC']=zeros((nsteps,3))
        out[sim]['protectedC']=zeros((nsteps,3))
        out[sim]['microbeC']=zeros(nsteps)
        out[sim]['CO2']=zeros(nsteps)

    return out

outputs=dict()
outputs['lowclay_highquality']=make_outputs(nsteps)
outputs['lowclay_lowquality']=make_outputs(nsteps)
outputs['highclay_highquality']=make_outputs(nsteps)
outputs['highclay_lowquality']=make_outputs(nsteps)
outputs['higherclay_highquality']=make_outputs(nsteps)
outputs['higherclay_lowquality']=make_outputs(nsteps)


dt=1.0/365.0
start=365*10

temperature=zeros(nsteps)+273.15+20
theta=zeros(nsteps)+0.5

T_warming_2=temperature.copy()
T_warming_2[start:]=T_warming_2[start:]+2.0
T_warming_5=temperature.copy()
T_warming_5[start:]=T_warming_5[start:]+5.0

T={
    'control':temperature,
    'warming_2':T_warming_2,
    'warming_5':T_warming_5,
    'labile_addition':temperature,
    'total_addition_30':temperature,
    'total_addition_100':temperature,
    'warming_2_total_addition_30':T_warming_2,
    'warming_2_labile_addition':T_warming_2,
    'litter_removal':temperature
}

inputs_highquality=column_stack([zeros((nsteps,1))+150.0,zeros((nsteps,1))+350.0,zeros((nsteps,1))+0.0])/1000  # kgC Per year
inputs_lowquality=column_stack([zeros((nsteps,1))+50.0,zeros((nsteps,1))+450.0,zeros((nsteps,1))+0.0])/1000  # kgC Per year

inputs_highquality_labile_addition=inputs_highquality.copy()
inputs_highquality_labile_addition[start:,0]=inputs_highquality_labile_addition[start:,0]*1.3

inputs_lowquality_labile_addition=inputs_lowquality.copy()
inputs_lowquality_labile_addition[start:,0]=inputs_lowquality_labile_addition[start:,0]*1.3

inputs_highquality_total_addition_30=inputs_highquality.copy()
inputs_highquality_total_addition_30[start:,:]=inputs_highquality_total_addition_30[start:,:]*1.3
inputs_highquality_total_addition_100=inputs_highquality.copy()
inputs_highquality_total_addition_100[start:,:]=inputs_highquality_total_addition_100[start:,:]*2.0

inputs_lowquality_total_addition_30=inputs_lowquality.copy()
inputs_lowquality_total_addition_30[start:,:]=inputs_lowquality_total_addition_30[start:,:]*1.3
inputs_lowquality_total_addition_100=inputs_lowquality.copy()
inputs_lowquality_total_addition_100[start:,:]=inputs_lowquality_total_addition_100[start:,:]*2.0

inputs_highquality_litter_removal=inputs_highquality.copy()
inputs_highquality_litter_removal[start:,:]=0.0
inputs_lowquality_litter_removal=inputs_lowquality.copy()
inputs_lowquality_litter_removal[start:,:]=0.0

inputs={'lowclay_highquality':{},'highclay_highquality':{},'lowclay_lowquality':{},'highclay_lowquality':{},'higherclay_highquality':{},'higherclay_lowquality':{}}
inputs['lowclay_highquality']['control']=inputs_highquality
inputs['lowclay_highquality']['warming_2']=inputs_highquality
inputs['lowclay_highquality']['warming_5']=inputs_highquality
inputs['lowclay_highquality']['labile_addition']=inputs_highquality_labile_addition
inputs['lowclay_highquality']['total_addition_30']=inputs_highquality_total_addition_30
inputs['lowclay_highquality']['total_addition_100']=inputs_highquality_total_addition_100
inputs['lowclay_highquality']['warming_2_total_addition_30']=inputs_highquality_total_addition_30
inputs['lowclay_highquality']['warming_2_labile_addition']=inputs_highquality_labile_addition
inputs['lowclay_highquality']['litter_removal']=inputs_highquality_litter_removal

inputs['highclay_highquality']['control']=inputs_highquality
inputs['highclay_highquality']['warming_2']=inputs_highquality
inputs['highclay_highquality']['warming_5']=inputs_highquality
inputs['highclay_highquality']['labile_addition']=inputs_highquality_labile_addition
inputs['highclay_highquality']['total_addition_30']=inputs_highquality_total_addition_30
inputs['highclay_highquality']['total_addition_100']=inputs_highquality_total_addition_100
inputs['highclay_highquality']['warming_2_total_addition_30']=inputs_highquality_total_addition_30
inputs['highclay_highquality']['warming_2_labile_addition']=inputs_highquality_labile_addition
inputs['highclay_highquality']['litter_removal']=inputs_highquality_litter_removal

inputs['higherclay_highquality']['control']=inputs_highquality
inputs['higherclay_highquality']['warming_2']=inputs_highquality
inputs['higherclay_highquality']['warming_5']=inputs_highquality
inputs['higherclay_highquality']['labile_addition']=inputs_highquality_labile_addition
inputs['higherclay_highquality']['total_addition_30']=inputs_highquality_total_addition_30
inputs['higherclay_highquality']['total_addition_100']=inputs_highquality_total_addition_100
inputs['higherclay_highquality']['warming_2_total_addition_30']=inputs_highquality_total_addition_30
inputs['higherclay_highquality']['warming_2_labile_addition']=inputs_highquality_labile_addition
inputs['higherclay_highquality']['litter_removal']=inputs_highquality_litter_removal

inputs['lowclay_lowquality']['control']=inputs_lowquality
inputs['lowclay_lowquality']['warming_2']=inputs_lowquality
inputs['lowclay_lowquality']['warming_5']=inputs_lowquality
inputs['lowclay_lowquality']['labile_addition']=inputs_lowquality_labile_addition
inputs['lowclay_lowquality']['total_addition_30']=inputs_lowquality_total_addition_30
inputs['lowclay_lowquality']['total_addition_100']=inputs_lowquality_total_addition_100
inputs['lowclay_lowquality']['warming_2_total_addition_30']=inputs_lowquality_total_addition_30
inputs['lowclay_lowquality']['warming_2_labile_addition']=inputs_lowquality_labile_addition
inputs['lowclay_lowquality']['litter_removal']=inputs_lowquality_litter_removal

inputs['highclay_lowquality']['control']=inputs_lowquality
inputs['highclay_lowquality']['warming_2']=inputs_lowquality
inputs['highclay_lowquality']['warming_5']=inputs_lowquality
inputs['highclay_lowquality']['labile_addition']=inputs_lowquality_labile_addition
inputs['highclay_lowquality']['total_addition_30']=inputs_lowquality_total_addition_30
inputs['highclay_lowquality']['total_addition_100']=inputs_lowquality_total_addition_100
inputs['highclay_lowquality']['warming_2_total_addition_30']=inputs_lowquality_total_addition_30
inputs['highclay_lowquality']['warming_2_labile_addition']=inputs_lowquality_labile_addition
inputs['highclay_lowquality']['litter_removal']=inputs_lowquality_litter_removal

inputs['higherclay_lowquality']['control']=inputs_lowquality
inputs['higherclay_lowquality']['warming_2']=inputs_lowquality
inputs['higherclay_lowquality']['warming_5']=inputs_lowquality
inputs['higherclay_lowquality']['labile_addition']=inputs_lowquality_labile_addition
inputs['higherclay_lowquality']['total_addition_30']=inputs_lowquality_total_addition_30
inputs['higherclay_lowquality']['total_addition_100']=inputs_lowquality_total_addition_100
inputs['higherclay_lowquality']['warming_2_total_addition_30']=inputs_lowquality_total_addition_30
inputs['higherclay_lowquality']['warming_2_labile_addition']=inputs_lowquality_labile_addition
inputs['higherclay_lowquality']['litter_removal']=inputs_lowquality_litter_removal


# Run model for all configs and sims
for conf in configs:
    for sim in sims:
        print 'Starting ',conf,sim
        for step in xrange(nsteps):
            if (step*dt)%10 ==0:
                print 'Year: %d'% (step*dt)


            output=cohorts[conf][sim].update(T[sim][step],theta[step],dt)
            cohorts[conf][sim].check_validity()
            outputs[conf][sim]['unprotectedC'][step,:]=cohorts[conf][sim].litterC
            outputs[conf][sim]['protectedC'][step]=cohorts[conf][sim].protectedC.sum()
            outputs[conf][sim]['microbeC'][step]=cohorts[conf][sim].livingMicrobeC
            outputs[conf][sim]['CO2'][step]=cohorts[conf][sim].CO2

            cohorts[conf][sim].add_carbon(inputs[conf][sim][step,:]*dt)


t=arange(nsteps)/365.0

plotstyles={
    'control':'-',
    'warming_2':'--',
    'warming_5':'--',
    'labile_addition':':',
    'total_addition_30':'-.',
    'total_addition_100':'-.',
    'warming_2_total_addition_30':'-.',
    'warming_2_labile_addition':':',
    'litter_removal':'-.'
}

plotcolors={
    'control':'k',
    'warming_2':[0.8,0.5,0.0],
    'warming_5':[1.0,0.0,0.0],
    'labile_addition':'g',
    'total_addition_30':[0.0,0.8,0.0],
    'total_addition_100':[0.0,0.7,0.6],
    'warming_2_total_addition_30':'m',
    'warming_2_labile_addition':'m',
    'litter_removal':'b'
}

plotlegends={
    'control':'Control',
    'warming_2':'+2$^\circ$ C',
    'warming_5':'+5$^\circ$ C',
    'labile_addition':'+30% Labile',
    'total_addition_30':'+30% Total',
    'total_addition_100':'+100% Total',
    'warming_2_total_addition_30':r'+30% total, +2$^\circ$ C',
    'warming_2_labile_addition':r'+30% labile, +2$^\circ$ C',
    'litter_removal':'No litter'
}

def plot_results(mode='total'):
    for plotnum,conf in enumerate(configs):
        subplot(2,3,plotnum+1)
        for sim in sims:
            if mode=='total':
                y=outputs[conf][sim]['microbeC']+\
                    outputs[conf][sim]['unprotectedC'].sum(axis=1)+\
                    outputs[conf][sim]['protectedC'].sum(axis=1)
            elif mode=='protected':
                y=outputs[conf][sim]['protectedC'].sum(axis=1)
            elif mode=='unprotected':
                y=outputs[conf][sim]['unprotectedC'].sum(axis=1)
            else:
                raise ValueError('Invalid plotting mode')

            h=plot(t,y,c=plotcolors[sim],ls=plotstyles[sim],label=plotlegends[sim],lw=2.0)

        title(conf.replace('_',', ').title()+' litter')
        ylabel('Carbon pools (kgC/m$^2$)')
        xlabel('Time (years)')

        if plotnum==2:
            leg=legend(fontsize='small',loc='best')
            leg.get_frame().set_alpha(0.5)

    subplots_adjust(hspace=0.25)
    figtext(0.5,0.95,mode.capitalize()+' C',fontsize=20,ha='center')

    draw()

figure(1);clf()
plot_results('total')
figure(2);clf()
plot_results('protected')
figure(3);clf()
plot_results('unprotected')

show()

def flatten_output(output):
    t=arange(len(output['CO2']))
    df=pandas.DataFrame(index=t)
    df['unprotectedC_fast']=output['unprotectedC'][:,0]
    df['unprotectedC_slow']=output['unprotectedC'][:,1]
    df['unprotectedC_deadmic']=output['unprotectedC'][:,2]
    df['protectedC']=output['protectedC'].sum(axis=1)
    df['CO2']=output['CO2']
    df['microbeC']=output['microbeC']

    return df

def flatten_output_all(output):
    df_out=pandas.DataFrame()
    for conf in output.keys():
        for expt in output[conf].keys():
            df=flatten_output(output[conf][expt])
            clay,quality=conf.split('_')
            df['clay']=clay
            df['quality']=quality
            df['Experiment']=expt

            df['Day']=arange(len(df))

            df_out=pandas.concat([df_out,df])

    df_out.set_index(['Day','clay','quality','Experiment'],inplace=True)
    df_out.sort_index(inplace=True)
    return df_out

def save_output_csv(output,directory=None):
    for conf in output.keys():
        for expt in output[conf].keys():
            df=flatten_output(output[conf][expt])
            fname='CORPSE-output-%s-%s.csv'%(conf,expt)
            if directory is not None:
                import os.path
                fname=os.path.join(directory,fname)
            df.to_csv(fname,index_label='Day')

def save_output_pickle(output,directory=None):
    df=flatten_output_all(output)
    fname='CORPSE-output-python-dataframe.pik'
    if directory is not None:
        import os.path
        fname=os.path.join(directory,fname)
    df.to_pickle(fname)


def steady_state_protected(protected,output):
    return (output['protected_produced']/dt)/(output['protected_turnover_rate']/protected)


if do_spinup:
    for conf in configs:
        sim='control'
        total=outputs[conf][sim]['microbeC']+outputs[conf][sim]['unprotectedC'].sum(axis=1)+outputs[conf][sim]['protectedC'].sum(axis=1)
        print conf
        print 'Total C change: %1.3f (%1.2f%%)'%(total[-1]-total[0],(total[-1]-total[0])/total[0]*100)
        unprot=outputs[conf][sim]['unprotectedC'].sum(axis=1)
        print 'Unprotected C change: %1.3f (%1.2f%%)'%(unprot[-1]-unprot[0],(unprot[-1]-unprot[0])/unprot[0]*100)
        prot=outputs[conf][sim]['protectedC'].sum(axis=1)
        print 'Protected C change: %1.3f (%1.2f%%)'%(prot[-1]-prot[0],(prot[-1]-prot[0])/prot[0]*100)


    print 'To make spun up cohorts:'
    for conf in configs:
        sim='control'
        c=cohorts[conf][sim].copy()
        o=c.update(temperature[-1],theta[-1],dt)
        p=steady_state_protected(c.protectedC,o)
        s='soil_carbon_cohort(litterC=[%1.4g,%1.4g,%1.4g], protectedC=[%1.4g,%1.4g,%1.4g], livingMicrobeC=%1.4g)'\
                    %(c.litterC[0],c.litterC[1],c.litterC[2],
                      p[0],p[1],p[2],
                      c.livingMicrobeC)
        print 'c_%s = CORPSE.%s'%(conf,s)
