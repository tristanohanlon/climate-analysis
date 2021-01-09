
import numpy as np
import os
import math
import constants
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib import rc


location = constants.network
labels = ["$\Delta cl$", "$\Delta lwp$", "$\Delta iwp$"]

# GLobal

SWcre_CLerror = [ [9.3] , [3.3] ]
SWcre_LWPerror = [ [3.8] , [8.1] ]
SWcre_IWPerror = [ [4.4] , [2.1] ]

LWcre_CLerror = [ [3.2] , [4.5] ] 
LWcre_LWPerror = [ [4.7] , [1.9] ] 
LWcre_IWPerror = [ [1.6] , [2.9] ]

fig, ( ax1, ax2 ) = plt.subplots(2, 1, figsize=(6, 6))
fig.suptitle('Sensitivity Spread in Cloud Radiative Effect from Ensemble Mean')

ax1.errorbar(0, 0, yerr = SWcre_CLerror, capsize=5)
ax1.errorbar(1, 0, yerr = SWcre_LWPerror, capsize=5)
ax1.errorbar(2, 0, yerr = SWcre_IWPerror, capsize=5)

ax2.errorbar(0, 0, yerr = LWcre_CLerror, capsize=5)
ax2.errorbar(1, 0, yerr = LWcre_LWPerror, capsize=5)
ax2.errorbar(2, 0, yerr = LWcre_IWPerror, capsize=5)

ax1.axhline(y=0, label = 'Ensemble Mean', color = 'black', linestyle='--')
ax2.axhline(y=0, color = 'black', linestyle='--')
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels)
ax1.set_yticks(np.arange(-12,13,2))
ax2.set_yticks(np.arange(-12,13,2))
ax1.set_ylim([-12,12])
ax2.set_ylim([-12,12])
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels)
ax1.legend(bbox_to_anchor=(1, 1))
ax1.set_ylabel('$\Delta$SW CRE $W m^{-2}$')
ax2.set_ylabel('$\Delta$LW CRE $W m^{-2}$')
plt.savefig( location + 'Images/RRTGM/' + "global_enemble_cre_bias_error_bars.pdf", format="pdf", bbox_inches='tight')
plt.show()




# SO


SWcre_CLerror = [ [6.1] , [3.5] ]
SWcre_LWPerror = [ [6.6] , [9.1] ]
SWcre_IWPerror = [ [4.1] , [3.4] ]

LWcre_CLerror = [ [1.5] , [3.8] ] 
LWcre_LWPerror = [ [5.8] , [3.9] ] 
LWcre_IWPerror = [ [2.8] , [2.0] ]

fig, ( ax1, ax2 ) = plt.subplots(2, 1, figsize=(6, 6))
fig.suptitle('SO Sensitivity Spread in Cloud Radiative Effect from Ensemble Mean')

ax1.errorbar(0, 0, yerr = SWcre_CLerror, capsize=5)
ax1.errorbar(1, 0, yerr = SWcre_LWPerror, capsize=5)
ax1.errorbar(2, 0, yerr = SWcre_IWPerror, capsize=5)

ax2.errorbar(0, 0, yerr = LWcre_CLerror, capsize=5)
ax2.errorbar(1, 0, yerr = LWcre_LWPerror, capsize=5)
ax2.errorbar(2, 0, yerr = LWcre_IWPerror, capsize=5)

ax1.axhline(y=0, label = 'Ensemble Mean', color = 'black', linestyle='--')
ax2.axhline(y=0, color = 'black', linestyle='--')
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels)
ax1.set_yticks(np.arange(-12,13,2))
ax2.set_yticks(np.arange(-12,13,2))
ax1.set_ylim([-12,12])
ax2.set_ylim([-12,12])
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels)
ax1.legend(bbox_to_anchor=(1, 1))
ax1.set_ylabel('$\Delta$SW CRE $W m^{-2}$')
ax2.set_ylabel('$\Delta$LW CRE $W m^{-2}$')
plt.savefig( location + 'Images/RRTGM/' + "SO_enemble_cre_bias_error_bars.pdf", format="pdf", bbox_inches='tight')
plt.show()