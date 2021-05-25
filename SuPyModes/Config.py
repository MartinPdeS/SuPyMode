PROPERTIES = ['Index', 'ITR', 'Field', 'Beta', 'Axes']

AdiabaticDict = { 'name' : r'Adiabatic criterion $\frac{1}{\rho} \frac{d\rho}{dz} \,\,$',
                  'unit' : r' $[ \mu m^{-1} ]$'}

CouplingDict = { 'name' : r'Coupling',
                 'unit' : r' '}

IndexDict = { 'name' : r'Effective Index',
              'unit' : r' [1]'}

BetaDict = { 'name' : r'Beta value $\beta_{i,j}$',
             'unit' : r' [1]'}

BasePlotKwarg = {'Index'    : {'xlim'  : None,
                               'ylim'  : None,
                               'xscale' : 'lin',
                               'yscale' : 'lin'},

                 'Coupling' : {'xlim'  : None,
                               'ylim'  : None,
                               'xscale' : 'lin',
                               'yscale' : 'lin'},

                 'Adiabatic': {'xlim'  : None,
                               'ylim'  : [1e-7, 1],
                               'xscale' : 'lin',
                               'yscale' : 'log'}
                               }
