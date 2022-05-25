PROPERTIES = ['Index', 'Field', 'Beta', 'Axes']

BasePlotKwarg = {'Index'    : { 'name' : r'Effective Index',
                                'unit' : r' [1]',
                                'xlim'  : None,
                                'ylim'  : None,
                                'xscale' : 'lin',
                                'yscale' : 'lin'},

                 'Beta'    : { 'name' : r'Beta',
                                         'unit' : r' [1]',
                                         'xlim'  : None,
                                         'ylim'  : None,
                                         'xscale' : 'lin',
                                         'yscale' : 'lin'},

                 'Coupling' : { 'name' : r'Coupling',
                                'unit' : r' ',
                                'xlim'  : None,
                                'ylim'  : None,
                                'xscale' : 'lin',
                                'yscale' : 'lin'},

                 'Adiabatic': { 'name'   : r'Adiabatic criterion $\frac{1}{\rho} \frac{d\rho}{dz} \,\,$',
                                'unit'   : r' $[ \mu m^{-1} ]$',
                                'xlim'   : None,
                                'ylim'   : [None, 1],
                                'xscale' : 'lin',
                                'yscale' : 'log'}
                               }
