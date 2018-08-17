# This should be a complete list of parameters with descriptions
# The set_params function will check that all are present and there are no extras
expected_params={	'vmaxref': 'Relative maximum enzymatic decomp rates (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
        	'minMicrobeC':	   'Minimum microbial biomass (fraction of total C)',
        	'Tmic': 'Microbial lifetime at 20C (years)',
        	'et':		  'Fraction of turnover not converted to CO2',
        	'eup': 'Carbon uptake efficiency (length 3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            }

#All soils: slope=0.4833,intercept=2.3282
#Alfisols: slope=0.5945, intercept=2.2788
def prot_clay(claypercent,slope=0.4833,intercept=2.3282,BD=1.15,porosity=0.4):
    ''' Calculate protection rate as a function of clay content, based on sorption isotherms from Mayes et al (2012) Table 3
    Calculates Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using bulk density
    Typically used as relative value for calculating protection_rate parameter.
    claypercent: Soil % clay (out of 100)
    slope: Either for all soils or a soil order, from Mayes paper
    intercept: Either for all soils or a soil order, from Mayes paper
    BD: Soil bulk density in g/cm3
    '''
    from numpy import log10
    prot=1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6)
    return prot

class soil_carbon_cohort:
    def __init__(self,
            litterC=[0.0,0.0,0.0],
            protectedC=[0.0,0.0,0.0],
            livingMicrobeC=0.0,
            CO2=0.0,
            params=None):
        '''Initialize a soil carbon cohort.
           litterC: Unprotected C, available for decomposition. (kgC/m2)
           protectedC: Physically or chemically protected, inaccessible to microbes. (kgC/m2)
           livingMicrobeC: Microbes that decompose the carbon. (kgC/m2)
           CO2: Cumulative CO2 production (kgC/m2)
           originalC: Cumulative carbon added to cohort, used for conservation checks. (kgC/m2)
           params: Parameter set (see set_params method)

           Vectors (for litterC and protectedC) are [Labile, Resistant, Microbial necromass]

           '''

        from numpy import array

        self.litterC=array(litterC)
        self.protectedC=array(protectedC)
        self.livingMicrobeC=livingMicrobeC
        self.CO2=CO2
        self.originalC=sum(litterC+protectedC)+livingMicrobeC+CO2

        if params is not None:
            self.set_params(params)
        else:
            self.params=None


    def set_params(self,params):
        '''params: dictionary containing parameter values. Should contain these fields (showing reasonable default values):
                 vmaxref=[2500,600,2000]; Relative maximum enzymatic decomp rates
                 Ea=[37e3,54e3,50e3];     Activation energy
                 kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
                 gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
                 minMicrobeC=1e-3;       Minimum microbial biomass (fraction of total C)
                 Tmic=0.15;        Microbial turnover rate
                 et=0.5;           Fraction of turnover not converted to CO2
                 eup=[0.6,0.05,0.6];  Carbon uptake efficiency
                 tProtected=75.0;     Protected C turnover time (years)
                 protection_rate=[1.0,0.0,1.0];  Protected carbon formation rate (year-1)'''

        from numpy import iterable,array
        self.params=params
        unused_params=expected_params.copy()
        for k in self.params.keys():
            if k not in expected_params:
                raise ValueError('Parameter set contains unexpected parameter %s'%k)
            unused_params.pop(k)
            if iterable(self.params[k]):
                self.params[k]=array(self.params[k])
        if len(unused_params)>0:
            for k in unused_params.keys():
                print ('Missing parameter: %s [%s]'%(k,unused_params[k]))
            raise ValueError('Missing parameters: %s'%unused_params.keys())


    def add_carbon(self,litterC=[0.0,0.0,0.0],protectedC=[0.0,0.0,0.0],
                        livingMicrobeC=0.0, CO2=0.0):
        'Add carbon to the cohort, keeping track of total for conservation checks.'

        from numpy import array

        if(any(array(litterC)<0.0)):
            raise ValueError('litterC must be >= 0')
        if(any(array(protectedC)<0.0)):
            raise ValueError('protectedC must be >= 0')
        if livingMicrobeC<0.0:
            raise ValueError('livingMicrobeC must be >=0.0')
        if CO2<0.0:
            raise ValueError('CO2 must be >=0.0')

        self.litterC=self.litterC+array(litterC)
        self.protectedC=self.protectedC+array(protectedC)
        self.livingMicrobeC=self.livingMicrobeC+livingMicrobeC
        self.CO2=self.CO2+CO2
        self.originalC=self.originalC+sum(litterC+protectedC)+livingMicrobeC+CO2

    def update(self,T,theta,dt):
        '''Update the cohort, with decomposition, microbial growth, etc.
           T: Temperature (K)
           theta: Soil water content (fraction of saturation)
           dt: Time step (years)

           Returns a dictionary of outputs.'''

        if T<0.0:
            raise ValueError('T must be >=0')

        if theta<0.0:
            theta=0.0
        elif theta>1.0:
            theta=1.0

        totalResp=0.0;

        et=self.params['et']
        eup=self.params['eup']

        # Calculate maximum potential C decomposition rate
        tempResp=self.Resp(self.litterC,self.livingMicrobeC,T,theta)

        # Carbon loss cannot exceed size of pool
        tempResp[dt*tempResp > self.litterC] = self.litterC[dt*tempResp > self.litterC]/dt


        # Microbial turnover
        microbeTurnover=max(0.0,(self.livingMicrobeC-self.params['minMicrobeC']*sum(self.litterC))/self.params['Tmic']);   # kg/m2/yr
        maintenance_resp=microbeTurnover*(1.0-et)

        deadmic_C_produced=dt*microbeTurnover*et   # actual fraction of microbial turnover
        self.litterC[2]=self.litterC[2]+deadmic_C_produced  # kg/m2

        # CO2 production and cumulative CO2 produced by cohort
        CO2prod=dt*(sum(tempResp*(1.0-eup))+maintenance_resp) # kg/m2
        self.CO2=self.CO2+CO2prod  # kg/m2

        microbeGrowth=sum(tempResp*eup);

        self.livingMicrobeC=self.livingMicrobeC + dt*(microbeGrowth-microbeTurnover);

        # Update the amount of organic C and N in the cohort after the decomposition process

        self.litterC=self.litterC-dt*tempResp;     # kg/m2
        totalResp=totalResp+tempResp;             # kg/m2/yr


        # Update protected carbon
        protectedCturnover = self.protectedC/self.params['tProtected'] ;

        if(sum(self.litterC)>0.0):

            newProtectedC = self.params['protection_rate']*self.litterC*dt;
            #  kg/m2      =   yr-1             kg/m2        yr

        else:
            newProtectedC = [0.0,0.0,0.0];

        newProtectedC[newProtectedC>self.litterC] = self.litterC[newProtectedC>self.litterC];
        self.protectedC = self.protectedC + newProtectedC - dt*protectedCturnover;
        self.litterC = self.litterC - newProtectedC + dt*protectedCturnover;

        protected_produced=newProtectedC;  # kg/m2
        protected_turnover_rate=protectedCturnover;  # kg/m2/dt

        outputs={
            'decomp':totalResp,
            'protected_produced':protected_produced,
            'protected_turnover_rate':protected_turnover_rate,
            'CO2prod':CO2prod
            }

        return outputs



    # Decomposition rate
    def Resp(self,Ctotal,Chet,T,theta):
        '''Chet        heterotrophic (microbial) C, living microbial biomass in the cohort
        T,theta     temperature (k), theta (fraction of 1.0)
        Ctotal      Substrate C (3-value vector)'''


        if(sum(Ctotal)==0.0 or theta==0.0 or Chet==0.0):
            return 0.0


        Resp=self.Vmax(T)*theta**3.0*(Ctotal)*Chet/(sum(Ctotal)*self.params['kC']+Chet)*(1.0-theta)**self.params['gas_diffusion_exp'];

        return Resp



    def Vmax(self,T):
        Tref=293.15;
        Rugas=8.314472;

        from numpy import exp

        # Normalization value
        alpha=self.params['vmaxref']/exp(-self.params['Ea']/(Rugas*Tref));
        Vmax=alpha*exp(-self.params['Ea']/(Rugas*T));
        return Vmax


    def totalC(self,includeCO2=True):
        if includeCO2:
            return sum(self.litterC+self.protectedC)+self.livingMicrobeC+self.CO2;
        else:
            return sum(self.litterC+self.protectedC)+self.livingMicrobeC

    def check_validity(self,tolerance=1e-6):
        'Check for valid values and carbon conservation'

        from numpy import isfinite

        bad=False;
        if any(self.litterC<0) or any(self.protectedC<0):
            bad=True;

        if self.livingMicrobeC<0 or self.CO2<0 or self.originalC<0:
            bad=True;

        if any(~isfinite(self.litterC)) or any(~isfinite(self.protectedC)) or ~isfinite(self.livingMicrobeC) or \
                ~isfinite(self.CO2) or ~isfinite(self.originalC):
            bad=True;

        if bad:
            print(self)
            raise RuntimeError('Cohort bad')

        # Check carbon conservation
        if abs(self.totalC()-self.originalC) > tolerance:
            print self
            raise RuntimeError('Cohort C conservation violated. Sum = %1.2d, Original = %1.2d'%(self.totalC(),self.originalC))



    def __mul__(self,value):
        # if not is_numlike(value):
        #     raise ValueError('Must multiply by a number')

        return soil_carbon_cohort(
            litterC=self.litterC*value,
            protectedC=self.protectedC*value,
            livingMicrobeC=self.livingMicrobeC*value,
            CO2=self.CO2*value,params=self.params)

    def __rmul__(self,value):
        return(self*value)

    def __div__(self,value):
        return self*(1.0/value)

    def __add__(self,cohort2):
        if not isinstance(cohort2,soil_carbon_cohort):
            raise ValueError('Cohorts can only be added to other cohorts')

        # Note: Resulting cohort will have parameters of self, not cohort2 (if they differ)
        return soil_carbon_cohort(
            litterC=self.litterC+cohort2.litterC,
            protectedC=self.protectedC+cohort2.protectedC,
            livingMicrobeC=self.livingMicrobeC+cohort2.livingMicrobeC,
            CO2=self.CO2+cohort2.CO2,params=self.params)

    def __neg__(self):
        return self*-1

    def __sub__(self,cohort2):
        return self + -cohort2

    def copy(self):
        out=self*1.0
        return out


    def __str__(self):
        s='soil_carbon_cohort(litterC=[%1.2g,%1.2g,%1.2g], protectedC=[%1.2g,%1.2g,%1.2g], livingMicrobeC=%1.2g, CO2=%1.2g, originalC=%1.2g)'\
                    %(self.litterC[0],self.litterC[1],self.litterC[2],
                      self.protectedC[0],self.protectedC[1],self.protectedC[2],
                      self.livingMicrobeC,self.CO2,self.originalC)
        return s

    def __repr__(self):
        return self.__str__()
