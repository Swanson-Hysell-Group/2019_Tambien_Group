# This file contains the functions used in the Seawater Model notebooks within this repository.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches





def dict_to_sorted_df(dictionary):
    """Convert a dictionary to a sorted pandas df.

    inputs:
    - dictionary = dictionary

    outputs:
    - df = sorted pandas dataframe
    """

    df = pd.DataFrame(list(dictionary.items()), columns=['var','val'])
    df.sort_values('var',inplace=True)
    df.reset_index(inplace=True,drop=True)

    return df





def plot_model_and_data(data, t, Sr87, Sr86, xy, t_start, axvlines, min8786, max8786, maxt, mint, save=None):
    """Plot model output and measured 87Sr/86Sr.

    inputs:
    - data = pandas dataframe containing measured 87Sr/86Sr
    - t = numpy array containing the time
    - Sr87 = numpy array containing Sr87 abundance in mols
    - Sr86 = numpy array containing Sr86 abundance in mols
    - xy = LOWESS fit
    - t_start = year that model starts (in years before 0Ma)
    - axvlines = list of axvlines
    - min8786 = minimum y-value to be used in the 87Sr/86Sr plot
    - max8786 = maximum y-value to be used in the 87Sr/86Sr plot
    - maxt = maximum/oldest time to plot (in Ma)
    - mint = minimum/youngest time to plot (in Ma)
    - save = optional save string
    """

    # shift the time vector, and convert it into Ma
    t = t * -1
    t = t + t_start
    t = t / 1e6

    # pull out the labels
    labels = np.array([])
    for i in range(len(data.index)):
        if data.loc[i,'label'] not in labels:
            labels = np.append(labels, data.loc[i,'label'])

    # generate colors
    color_idx = np.linspace(0, 1, len(labels))

    # plot
    fig, ax = plt.subplots()
    # the measured data
    for i in range(len(labels)):
        ax.scatter(data[data['label']==labels[i]]['age'],data[data['label']==labels[i]]['87Sr/86Sr_primary'],\
                   label=labels[i],color=plt.cm.cool(color_idx[i]))
    # the axvlines
    for i in range(len(axvlines)):
        if i==0:
            ax.axvline(axvlines[i], c='red', ls='--', lw=0.5, label='flux change')
        else:
            ax.axvline(axvlines[i], c='red', ls='--', lw=0.5)
    # the LOWESS fit
    ax.plot(xy[:,0], xy[:,1], label='LOWESS', color='mediumseagreen', lw=2)
    # the model output
    ax.plot(t,Sr87/Sr86,label='model',color='firebrick',lw=3)
    # the glacial
    ax.add_patch(patches.Rectangle((mint,min8786),\
                                   717-mint,\
                                   max8786-min8786,\
                                   facecolor='dodgerblue'))
    ax.set_xlim([mint,maxt])
    ax.set_ylim([min8786,max8786])
    ax.set_xlabel('age [Ma]')
    ax.set_ylabel('$^{87}$Sr/$^{86}$Sr')
    ax.invert_xaxis()
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, fontsize=8)
    if save!=None:
        plt.savefig(save)
    plt.show(fig)





def plot_model_and_data_2(data, t, Sr87, Sr86, min8786, max8786, maxt, mint):
    """Plot model output and measured 87Sr/86Sr.

    inputs:
    - data = pandas dataframe containing measured 87Sr/86Sr
    - t = numpy array containing the time
    - Sr87 = numpy array containing Sr87 abundance in mols
    - Sr86 = numpy array containing Sr86 abundance in mols
    - min8786 = minimum y-value to be used in the 87Sr/86Sr plot
    - max8786 = maximum y-value to be used in the 87Sr/86Sr plot
    - maxt = maximum/oldest time to plot (in Ma)
    - mint = minimum/youngest time to plot (in Ma)
    """
    # pull out the labels
    labels = np.array([])
    for i in range(len(data.index)):
        if data.loc[i,'label'] not in labels:
            labels = np.append(labels, data.loc[i,'label'])

    # generate colors
    color_idx = np.linspace(0, 1, len(labels))

    # plot
    fig, ax = plt.subplots()
    # the measured data
    for i in range(len(labels)):
        ax.scatter(data[data['label']==labels[i]]['age'],data[data['label']==labels[i]]['87Sr/86Sr_primary'],\
                   label=labels[i],color=plt.cm.cool(color_idx[i]))
    # the model output
    ax.plot(t,Sr87/Sr86,label='model',color='firebrick',lw=3)
    ax.set_xlim([mint,maxt])
    ax.set_ylim([min8786,max8786])
    ax.set_xlabel('age [Ma]')
    ax.set_ylabel('$^{87}$Sr/$^{86}$Sr')
    ax.invert_xaxis()
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show(fig)






def plot_seawater_model(t, Mg, Ca, Sr88, Sr87, Sr86, min8786, max8786, minMgCa, maxMgCa, mint, maxt, save=None):
    """Plot the results of the seawater model.

    inputs:
    - t = numpy array containing the time
    - Mg = numpy array containing Mg abundance in mols
    - Ca = numpy array containing Ca abundance in mols
    - Sr88 = numpy array containing Sr88 abundance in mols
    - Sr87 = numpy array containing Sr87 abundance in mols
    - Sr86 = numpy array containing Sr86 abundance in mols
    - min8786 = minimum y-value to be used in the 87Sr/86Sr plot
    - max8786 = maximum y-value to be used in the 87Sr/86Sr plot
    - minMgCa = minimum y-value to be used in the Mg/Ca plot
    - maxMgCa = maximum y-value to be used in the Mg/Ca plot
    - mint = minimum time to plot
    - maxt = maximum time to plot
    - save = optional save string

    ouputs:
    - fig = figure handle of output plot
    - ax = axes handles of output plot
    """

    fig, ax = plt.subplots(4,2,sharex=True,figsize=(10,10))
    ax[0,0].plot(t,Mg)
    ax[0,0].set_title('Mg')
    ax[0,0].set_ylabel('[mol]')
    ax[0,1].plot(t,Ca)
    ax[0,1].set_title('Ca')
    ax[0,1].set_ylabel('[mol]')
    ax[1,0].plot(t,Sr86)
    ax[1,0].set_title('$^{86}$Sr')
    ax[1,0].set_ylabel('[mol]')
    ax[1,1].plot(t,Sr87)
    ax[1,1].set_title('$^{87}$Sr')
    ax[1,1].set_ylabel('[mol]')
    ax[2,0].plot(t,Sr88)
    ax[2,0].set_title('$^{88}$Sr')
    ax[2,0].set_ylabel('[mol]')
    ax[2,1].plot(t,Sr88+Sr87+Sr86)
    ax[2,1].set_title('Sr')
    ax[2,1].set_ylabel('[mol]')
    ax[3,0].plot(t,Sr87/Sr86,color='firebrick',lw=3)
    ax[3,0].set_title('$^{87}$Sr/$^{86}$Sr')
    ax[3,0].set_xlabel('model duration [yrs]')
    ax[3,0].set_ylim(min8786,max8786)
    ax[3,1].plot(t,Mg/Ca,color='firebrick',lw=3)
    ax[3,1].set_title('Mg/Ca')
    ax[3,1].set_xlabel('model duration [yrs]')
    ax[3,1].set_ylim(minMgCa,maxMgCa)
    ax[3,1].set_xlim(mint,maxt)
    if save!=None:
        plt.savefig(save)
    plt.show(fig)
    return fig, ax





def plot_seawater_model_2(t, Mg, Ca, Sr88, Sr87, Sr86,\
                          bulk_W_carbs, bulk_W_ccs, bulk_W_LIPs,\
                          carb_DR_means, cc_DR_means, LIP_DR_means,\
                          min8786, max8786, minMgCa, maxMgCa, mint, maxt):
    """Plot the results of the seawater model.

    inputs:
    - t = numpy array containing the time
    - Mg = numpy array containing Mg abundance in mols
    - Ca = numpy array containing Ca abundance in mols
    - Sr88 = numpy array containing Sr88 abundance in mols
    - Sr87 = numpy array containing Sr87 abundance in mols
    - Sr86 = numpy array containing Sr86 abundance in mols
    - bulk_W_carbs = numpy array containing W_X_carb in kg/yr
    - bulk_W_ccs = numpy array containing W_X_cc in kg/yr
    - bulk_W_LIPs = numpy array containing W_X_LIP in kg/yr
    - carb_DR_means = numpy array containing carb_DR_mean in kg/km2/yr
    - cc_DR_means = numpy array containing cc_DR_mean in kg/km2/yr
    - LIP_DR_means = numpy array containing LIP_DR_mean in kg/km2/yr
    - min8786 = minimum y-value to be used in the 87Sr/86Sr plot
    - max8786 = maximum y-value to be used in the 87Sr/86Sr plot
    - minMgCa = minimum y-value to be used in the Mg/Ca plot
    - maxMgCa = maximum y-value to be used in the Mg/Ca plot
    - mint = oldest age to plot (in Ma)
    - maxt = youngest age to plot (in Ma)

    ouputs:
    - fig = figure handle of output plot
    - ax = axes handles of output plot
    """

    fig, ax = plt.subplots(5,3,sharex=True,figsize=(12.5,10))
    ax[0,0].plot(t,Mg)
    ax[0,0].set_title('Mg')
    ax[0,0].set_ylabel('[mol]')

    ax[0,1].plot(t,Ca)
    ax[0,1].set_title('Ca')
    ax[0,1].set_ylabel('[mol]')

    ax[0,2].plot(t,Sr88+Sr87+Sr86)
    ax[0,2].set_title('Sr')
    ax[0,2].set_ylabel('[mol]')

    ax[1,0].plot(t,Sr86)
    ax[1,0].set_title('$^{86}$Sr')
    ax[1,0].set_ylabel('[mol]')

    ax[1,1].plot(t,Sr87)
    ax[1,1].set_title('$^{87}$Sr')
    ax[1,1].set_ylabel('[mol]')

    ax[1,2].plot(t,Sr88)
    ax[1,2].set_title('$^{88}$Sr')
    ax[1,2].set_ylabel('[mol]')

    ax[2,0].plot(t,bulk_W_carbs)
    ax[2,0].set_title('W_X_carb')
    ax[2,0].set_ylabel('kg/yr')

    ax[2,1].plot(t,bulk_W_ccs)
    ax[2,1].set_title('W_X_cc')
    ax[2,1].set_ylabel('kg/yr')

    ax[2,2].plot(t,bulk_W_LIPs)
    ax[2,2].set_title('W_X_LIP')
    ax[2,2].set_ylabel('kg/yr')

    ax[3,0].plot(t,carb_DR_means)
    ax[3,0].set_title('carb_DR_mean')
    ax[3,0].set_ylabel('kg/km$^{2}$/yr')

    ax[3,1].plot(t,cc_DR_means)
    ax[3,1].set_title('cc_DR_mean')
    ax[3,1].set_ylabel('kg/km$^{2}$/yr')

    ax[3,2].plot(t,LIP_DR_means)
    ax[3,2].set_title('LIP_DR_mean')
    ax[3,2].set_ylabel('kg/km$^{2}$/yr')

    ax[4,0].plot(t,Sr87/Sr86,color='firebrick',lw=3)
    ax[4,0].set_title('$^{87}$Sr/$^{86}$Sr')
    ax[4,0].set_xlabel('time [Ma]')
    ax[4,0].set_ylim(min8786,max8786)

    ax[4,1].set_yticks([])
    ax[4,1].set_xlabel('time [Ma]')

    ax[4,2].plot(t,Mg/Ca,color='firebrick',lw=3)
    ax[4,2].set_title('Mg/Ca')
    ax[4,2].set_xlabel('time [Ma]')
    ax[4,2].set_ylim(minMgCa,maxMgCa)
    ax[4,2].set_xlim(mint,maxt)

    for i in range(ax.shape[0]):
        for j in range(ax.shape[1]):
            ax[i,j].yaxis.label.set_size(8)
            ax[i,j].tick_params(axis='both', which='major', labelsize=7)
            ax[i,j].get_yaxis().get_major_formatter().set_useOffset(False)

    plt.show(fig)
    return fig, ax
