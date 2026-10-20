import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import matplotlib
import matplotlib.ticker as ticker
import sys
filename_out='spectra'
what_to_plot = ['FT_2.dat'      ,'FT_5.dat'    ,'FT_6.dat'     ,'FT_4.dat'   ,'FT_3.dat'         ,'FT_1_third.dat','FT_1_second.dat','FT_1.dat' ]
labels =       [ 'oct-H2O'      ,'tet1-H2O'    ,'tet2-H2O'     ,'OH$^{-}$'   , 'oct-OH'          ,'third'         ,'second'         ,'Total'    ]
linestyle =    ['-'             , '-'          ,'-'            ,'-'          ,'-'                ,'-.'            ,'--'             ,'-'        ]
colours=       ['tab:green'         ,'deepskyblue' ,'mediumorchid' ,'tab:blue'  ,'tab:orange'    ,'black'         ,'black'          ,'black'    ]

#write(54,*) '#    ,  1   ,     2      ,    3       ,    4      ,   5       ,  6       ,    7         ,   8         , 9           ,   10        ,   11       '
#write(54,*) '#time, full , oct-H2O-top, oct-H2O-bot, oct-OH-top, oct-OH-bot, OH- (free), H2O-tet1-top, H2O-tet2-top, H2O-tet1-bot, H2O-tet2-bot, rest'
#                1          2              3               4              5                 6              7               8               9               10            11
factors = np.ones(11)
factors[10]=0.5
y_scale=10**-5

x_min=2000
x_max=4000
y_min=-1.4
y_max=1

# plot settings
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
smallsize = 20
largesize = 24
plt.rcParams.update({'font.size': 32})
plt.rc('xtick', labelsize = smallsize, direction='in')
plt.rc('ytick', labelsize= smallsize, direction='in')
plt.rc('axes', labelsize = largesize)
plt.rc('axes', titlesize = largesize, linewidth=2)
#plt.rc('legend', fontsize=largesize)
plt.rc('lines', markersize=8, linewidth=2)
#plt.rc('legend', frameon=True,framealpha=1,)
plt.rcParams['figure.figsize'] = [6,4]
plt.rc('text', usetex=False)
#plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['mathtext.default'] = 'regular'

#rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
ticksize=15
labelsize=15
legendsize=15

all_data = {}
for i,what in enumerate(what_to_plot): 
    data=np.loadtxt('{}'.format(what))
    data[:,1]=-data[:,1]/150
    #Normalize step HARD CODED

    if what == "FT_4.dat" or what =="FT_1.dat" or what == "FT_1_third.dat" or what == "FT_1_second.dat":
        data[:,1]=data[:,1]/2
    if what == "FT_3.dat":
        data[:,1]=data[:,1]*2

    all_data[labels[i]]=data

#fmt = ['o','x','v','^','s']
for i, label in enumerate(labels):
  data=all_data[label]
  plt.plot(data[:,0], data[:,1]*factors[i]*y_scale,  color=colours[i], linestyle=linestyle[i],label='{}'.format(label))
  
# Customize the plot
#plt.title(r'H$_3$O$^+$(aq) ',fontsize=labelsize)
plt.xlabel(r'Wavenumber (cm$^{-1}$)',fontsize=labelsize)
plt.ylabel(r' Im$\chi^{(2)}$ (a.u.)',fontsize=labelsize)
#plt.ylabel(r"D (10$^{-9}$ (m$^2$/s)",fontsize=labelsize)
#plt.legend()
plt.grid(True)

plt.xlim([x_min,x_max]) 
plt.ylim([y_min,y_max]) 
plt.xticks(fontsize=ticksize)
#plt.yticks([],fontsize=ticksize)
plt.yticks(fontsize=ticksize)
#plt.legend(loc='upper center',ncol=3, fancybox=True, shadow=True,fontsize=legendsize, bbox_to_anchor=(0.5,1.40))
#plt.legend(loc='best',ncol=2, fancybox=True, shadow=True,fontsize=legendsize)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=6)
#plt.savefig('{}.pdf'.format(filename_out),bbox_inches='tight')
# Show the plot
plt.show()
print('{}.pdf'.format(filename_out))

