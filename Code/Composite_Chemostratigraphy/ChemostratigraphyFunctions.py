# This file contains the functions used in the Chemostratigraphy notebooks within this repository.

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
import statistics
import jinja2
import json
import pandas as pd
#import mpld3





def age_model(section,constraints,indices):
    """Generate an age model for any given section.

    inputs:
    - section = pandas dataframe for a section
    - constraints = an ordered list of age constraints
    - indices = an ordered list of indices which correspond to the age constraints
    """

    # the number of input age constraints
    n = len(indices)

    # create an age column, and input the age constraints where appropriate
    for i in range(n):
        section.loc[indices[i],'age'] = constraints[i]

    # create a list containing the individual sedimentation rates
    sedrates = []
    for i in range(n-1):
        if indices[i+1]-indices[i] != 1:
            sedrate = -(section['strat_m'][indices[i+1]]-section['strat_m'][indices[i]])/\
            (section['age'][indices[i+1]]-section['age'][indices[i]])
            sedrates.append(sedrate)

    # fill the age column using the sedimentation rates
    j = 0
    for i in range(len(sedrates)):
        if indices[j+1]-indices[j] == 1:
            j = j + 1
        for k in range(indices[j]+1, indices[j+1]):
            meter_diff = section['strat_m'][k]-section['strat_m'][k-1]
            age_diff = (1/sedrates[i])*meter_diff
            section.loc[k,'age'] = section['age'][k-1]-age_diff
        j = j + 1





def c_sr_comp_plot(dataframes, colors, Sr_threshold, MnSr_threshold, save='No', ymin=-500, ymax=5000):
    """Plot d13C and 87Sr/86Sr for an arbitrary number of sections.

    inputs:
    note that all input lists must be ordered identically
    - dataframes = list of pandas dataframes
    - colors = list of strings of colours
    - Sr_threshold = list of [Sr] thresholds (make 0 if you do not wish to apply)
    - MnSr_threshold = list of [Mn]/[Sr] thresholds (make inf if you do not wish to apply)
    - save = string of save file name (if left as 'No', will not save)
    - ymin = optional - minimum ylim (default -500)
    - ymax = optional - maximum ylim (default 4400)
    """

    # the number of input dataframes
    n = len(dataframes)

    plt.figure(figsize=(6,15))

    # the C data
    plt.subplot(1,2,1)
    for i in range(n):
        # plot MTS data on top as triangles, if it exists
        if 'mts_mtx' in dataframes[i].columns:
            mtx_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mtx']
            mts_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mts']
            # plot
            plt.scatter(mtx_df['d13C'],mtx_df['cum_strat_height'],\
                        marker='o',c=colors[i],edgecolors='none')
            plt.scatter(mts_df['d13C'],mts_df['cum_strat_height'],\
                        marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
        else:
            plt.scatter(dataframes[i]['d13C'],dataframes[i]['cum_strat_height'],\
                        marker='o',c=colors[i],edgecolors='none')
    plt.xlim(-12, 10)
    plt.ylabel('stratigraphic height')
    plt.xlabel('$\delta^{13}$C')
    plt.gca().xaxis.grid(True)
    plt.ylim(ymin,ymax)

    # the Sr data
    plt.subplot(1,2,2)
    for i in range(n):
        if '87Sr/86Sr' in dataframes[i].columns:
            print(str(dataframes[i].section[0]) + ' has Sr data.')
            # create a Mn/Sr column
            for j in range(len(dataframes[i].index)):
                dataframes[i].loc[j,'Mn/Sr'] = dataframes[i]['Mn_ppm'][j] / dataframes[i]['Sr_ppm'][j]
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # also plot the data with no trace element (NTE) data
                mtx_df_NTE = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                           (np.isnan(dataframes[i]['Sr_ppm']))]
                mts_df_NTE = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                           (np.isnan(dataframes[i]['Sr_ppm']))]
                # plot
                plt.scatter(mtx_df_pass['87Sr/86Sr'],mtx_df_pass['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mtx_df_NTE['87Sr/86Sr'],mtx_df_NTE['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df_pass['87Sr/86Sr'],mts_df_pass['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
                plt.scatter(mts_df_NTE['87Sr/86Sr'],mts_df_NTE['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # also plot the data with no trace element (NTE) data
                df_NTE = dataframes[i][np.isnan(dataframes[i]['Sr_ppm'])]
                # plot
                plt.scatter(df_pass['87Sr/86Sr'],df_pass['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(df_NTE['87Sr/86Sr'],df_NTE['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.xlim(0.705, 0.7075)
    plt.ylim(ymin,ymax)
    axes = plt.gca()
    axes.get_yaxis().set_ticklabels([])
    axes.xaxis.grid(True)
    plt.xlabel('$^{87}$Sr/$^{86}$Sr')

    if save != 'No':
        plt.savefig('Output/' + save)

    plt.show()





def c_sr_o_ca_comp_plot(dataframes, colors, Sr_threshold, MnSr_threshold, ymin=-500, ymax=5000):
    """Plot d13C, 87Sr/86Sr, d18O, d44/40Ca for an arbitrary number of sections.

    inputs:
    note that all input lists must be ordered identically
    - dataframes = list of pandas dataframes
    - colors = list of strings of colours
    - Sr_threshold = list of [Sr] thresholds (make 0 if you do not wish to apply)
    - MnSr_threshold = list of [Mn]/[Sr] thresholds (make inf if you do not wish to apply)
    - ymin = optional - minimum ylim (default -500)
    - ymax = optional - maximum ylim (default 4400)
    """

    # number of input dataframes
    n = len(dataframes)

    plt.figure(figsize=(16,15))

    # the C data
    plt.subplot(1,4,1)
    for i in range(n):
        # plot MTS data on top as triangles, if it exists
        if 'mts_mtx' in dataframes[i].columns:
            mtx_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mtx']
            mts_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mts']
            # plot
            plt.scatter(mtx_df['d13C'],mtx_df['cum_strat_height'],\
                        marker='o',c=colors[i],edgecolors='none')
            plt.scatter(mts_df['d13C'],mts_df['cum_strat_height'],\
                        marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
        else:
            plt.scatter(dataframes[i]['d13C'],dataframes[i]['cum_strat_height'],\
                        marker='o',c=colors[i],edgecolors='none')
    plt.xlim(-12, 10)
    plt.ylabel('stratigraphic height')
    plt.xlabel('$\delta^{13}$C')
    plt.ylim(ymin,ymax)

    # the Sr data
    plt.subplot(1,4,2)
    for i in range(n):
        if '87Sr/86Sr' in dataframes[i].columns:
            print(str(dataframes[i].section[0]) + ' has Sr data.')
            # create a Mn/Sr column
            for j in range(len(dataframes[i].index)):
                dataframes[i].loc[j,'Mn/Sr'] = dataframes[i]['Mn_ppm'][j] / dataframes[i]['Sr_ppm'][j]
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # also plot the data with no trace element (NTE) data
                mtx_df_NTE = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                           (np.isnan(dataframes[i]['Sr_ppm']))]
                mts_df_NTE = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                           (np.isnan(dataframes[i]['Sr_ppm']))]
                # plot
                plt.scatter(mtx_df_pass['87Sr/86Sr'],mtx_df_pass['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mtx_df_NTE['87Sr/86Sr'],mtx_df_NTE['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df_pass['87Sr/86Sr'],mts_df_pass['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
                plt.scatter(mts_df_NTE['87Sr/86Sr'],mts_df_NTE['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # also plot the data with no trace element (NTE) data
                df_NTE = dataframes[i][np.isnan(dataframes[i]['Sr_ppm'])]
                # plot
                plt.scatter(df_pass['87Sr/86Sr'],df_pass['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(df_NTE['87Sr/86Sr'],df_NTE['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.xlim(0.705, 0.7075)
    plt.ylim(ymin,ymax)
    axes = plt.gca()
    axes.get_yaxis().set_ticklabels([])
    plt.xlabel('$^{87}$Sr/$^{86}$Sr')

    # the O data
    plt.subplot(1,4,3)
    for i in range(n):
        if 'd18O' in dataframes[i].columns:
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                mtx_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mtx']
                mts_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mts']
                # plot
                plt.scatter(mtx_df['d18O'],mtx_df['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df['d18O'],mts_df['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                plt.scatter(dataframes[i]['d18O'],dataframes[i]['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.ylim(ymin,ymax)
    axes = plt.gca()
    axes.get_yaxis().set_ticklabels([])
    plt.xlabel('$\delta^{18}$O')

    # the Ca data
    plt.subplot(1,4,4)
    for i in range(n):
        if 'd44_40Ca' in dataframes[i].columns:
            print(str(dataframes[i].section[0]) + ' has Ca data.')
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                mtx_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mtx']
                mts_df = dataframes[i][dataframes[i]['mts_mtx'] == 'mts']
                # plot
                plt.scatter(mtx_df['d44_40Ca'],mtx_df['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df['d44_40Ca'],mts_df['cum_strat_height'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                plt.scatter(dataframes[i]['d44_40Ca'],dataframes[i]['cum_strat_height'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.ylim(ymin,ymax)
    plt.xlim(-2,0)
    axes = plt.gca()
    axes.get_yaxis().set_ticklabels([])
    plt.xlabel('$\delta^{44/40}$Ca')
    plt.show()





def c_sr_plot(section_data, Sr_threshold=0, MnSr_threshold=float('inf'), save='No', xlim=None):
    """Plot d13C and 87Sr/86Sr for a single section.

    inputs:
    - section_data = pandas dataframe
    - Sr_threshold = [Sr] threshold (0 by default)
    - MnSr_threshold = [Mn]/[Sr] threshold (inf by default)
    - save = string of save file name (if left as 'No', will not save)
    - xlim = xlim of plot (default = None)
    """

    # if Sr data exists, plot both the C and Sr data
    if '87Sr/86Sr' in section_data.columns:
        plt.figure(figsize=(10, 5))

        # C data
        plt.subplot(2,1,1)
        # plot MTS data on top as triangles, if it exists
        if 'mts_mtx' in section_data.columns:
            mtx_df = section_data[section_data['mts_mtx'] == 'mtx']
            mts_df = section_data[section_data['mts_mtx'] == 'mts']
            # plot
            plt.scatter(mtx_df['strat_height'],mtx_df['d13C'],\
                        marker='o',c='dodgerblue',edgecolors='none')
            plt.scatter(mts_df['strat_height'],mts_df['d13C'],\
                        marker='^',facecolors='w',edgecolors='dodgerblue',linewidths=1)
        else:
            plt.scatter(section_data['strat_height'],section_data['d13C'],\
                        marker='o',c='dodgerblue',edgecolors='none')
        plt.ylabel('$\delta^{13}$C')
        axes = plt.gca()
        if xlim == None:
            xlim = axes.get_xlim()
        else:
            plt.xlim(xlim)
        axes.get_xaxis().set_ticklabels([])
        axes.yaxis.grid(True)
        plt.ylim(-12,10)

        # Sr data
        plt.subplot(2,1,2)
        # create a Mn/Sr column
        for i in range(len(section_data.index)):
            section_data.loc[i,'Mn/Sr'] = section_data['Mn_ppm'][i] / section_data['Sr_ppm'][i]
        # plot MTS data on top as triangles, if it exists
        if 'mts_mtx' in section_data.columns:
            # only select data that passes the thresholds
            mtx_df_pass = section_data[(section_data['mts_mtx']=='mtx') &\
                                        (section_data['Sr_ppm']>=Sr_threshold) &\
                                        (section_data['Mn/Sr']<=MnSr_threshold)]
            mts_df_pass = section_data[(section_data['mts_mtx']=='mts') &\
                                        (section_data['Sr_ppm']>=Sr_threshold) &\
                                        (section_data['Mn/Sr']<=MnSr_threshold)]
            # also plot the data with no trace element (NTE) data
            mtx_df_NTE = section_data[(section_data['mts_mtx']=='mtx') &\
                                       (np.isnan(section_data['Sr_ppm']))]
            mts_df_NTE = section_data[(section_data['mts_mtx']=='mts') &\
                                       (np.isnan(section_data['Sr_ppm']))]
            # plot
            plt.scatter(mtx_df_pass['strat_height'],mtx_df_pass['87Sr/86Sr'],\
                        marker='o',c='firebrick',edgecolors='none')
            plt.scatter(mtx_df_NTE['strat_height'],mtx_df_NTE['87Sr/86Sr'],\
                        marker='o',c='firebrick',edgecolors='none')
            plt.scatter(mts_df_pass['strat_height'],mts_df_pass['87Sr/86Sr'],\
                        marker='^',facecolors='w',edgecolors='firebrick',linewidths=1)
            plt.scatter(mts_df_NTE['strat_height'],mts_df_NTE['87Sr/86Sr'],\
                        marker='^',facecolors='w',edgecolors='firebrick',linewidths=1)
        else:
            # only select data that passes the thresholds
            df_pass = section_data[(section_data['Sr_ppm']>=Sr_threshold) &\
                                    (section_data['Mn/Sr']<=MnSr_threshold)]
            # also plot the data with no trace element (NTE) data
            df_NTE = section_data[np.isnan(section_data['Sr_ppm'])]
            # plot
            plt.scatter(df_pass['strat_height'],df_pass['87Sr/86Sr'],\
                        marker='o',c='firebrick',edgecolors='none')
            plt.scatter(df_NTE['strat_height'],df_NTE['87Sr/86Sr'],\
                        marker='o',c='firebrick',edgecolors='none')
        plt.ylabel('$^{87}$Sr/$^{86}$Sr')
        plt.xlabel('stratigraphic height [m]')
        plt.xlim(xlim)
        axes = plt.gca()
        axes.yaxis.grid(True)
        plt.ylim(0.706, 0.7075)

        plt.suptitle(str(section_data.section[0]), fontsize=14, fontweight='bold')

        if save != 'No':
            plt.savefig('Output/' + save)
        plt.show()

    # if Sr data doesn't exist, plot only the C data
    else:
        print('No Sr data.')

        plt.figure(figsize=(10, 2.5))
        if 'mts_mtx' in section_data.columns:
            mtx_df = section_data[section_data['mts_mtx'] == 'mtx']
            mts_df = section_data[section_data['mts_mtx'] == 'mts']
            # plot
            plt.scatter(mtx_df['strat_height'],mtx_df['d13C'],\
                        marker='o',c='dodgerblue',edgecolors='none')
            plt.scatter(mts_df['strat_height'],mts_df['d13C'],\
                        marker='^',facecolors='w',edgecolors='dodgerblue',linewidths=1)
        else:
            plt.scatter(section_data['strat_height'],section_data['d13C'],\
                        marker='o',c='dodgerblue',edgecolors='none')
        plt.xlabel('stratigraphic height [m]')
        plt.ylabel('$\delta^{13}$C')
        plt.title(str(section_data.section[0]), fontsize=14, fontweight='bold')
        plt.gca().yaxis.grid(True)
        plt.ylim(-12,10)
        if xlim != None:
            plt.xlim(xlim)

        if save != 'No':
            plt.savefig('Output/' + save)
        plt.show()





def clast_plot(dataframe, bins, slice_min=0, slice_max=0, save='No'):
    """Create a histogram for diamictite clast data.

    inputs:
    - dataframe = input dataframe
    - bins = number of bins to plot
    - slice_min = index to start data slice, inclusive (optional)
    - slice_max = index to end data slice, inclusive (optional)
    - save = string of save file name (if left as 'No', will not save)
    """

    # if no data slice is provided
    if slice_min==0 and slice_max==0:
        # get rid of the nan data
        d13C_data = []
        for i in range(len(dataframe.index)):
            if np.isfinite(dataframe['d13C'][i]):
                d13C_data.append(dataframe['d13C'][i])
        # plot
        plt.figure()
        plt.hist(d13C_data,bins=bins,color='goldenrod')
        plt.ylabel('n')
        plt.xlabel('$\delta^{13}C$')
        plt.xlim(-12,10)

    # if data slice is provided
    else:
        # get rid of the nan data
        d13C_data = []
        for i in range(slice_min,slice_max+1):
            if np.isfinite(dataframe['d13C'][i]):
                d13C_data.append(dataframe['d13C'][i])
        # plot
        plt.figure()
        plt.hist(d13C_data,bins=bins,color='goldenrod')
        plt.ylabel('n')
        plt.xlabel('$\delta^{13}C$')
        plt.xlim(-12,10)

    # save if desired
    if save != 'No':
        plt.savefig(save)
    plt.show()





def create_tambien_csv(sections, Sr_threshold, MnSr_threshold, save, out, plot, rename):
    """Take in selected sections, output collated data into a single sorted .csv, and plot the result.

    inputs:
    - sections = list of ordered pandas dataframes
    - Sr_threshold = list of ordered [Sr] thresholds - used for plotting only
    - MnSr_threshold = list of ordered [Mn]/[Sr] thresholds - used for plotting only
    - save = string of save path and file name
    - out = if True, return output as pandas dataframe
    - plot = if True, create plot
    - rename = if True, rename columns for input into Tonian_Composite.ipynb
    """

    # concatenate
    output = pd.concat(sections, ignore_index=True)

    # re-order columns
    if '87Sr/86Sr' in output.columns.values:
        cols = ['section','strat_height','d13C','d18O','cum_strat_height','mts_mtx','d44_40Ca',\
                '87Sr/86Sr','Ca_ppm','Mg_ppm','Al_ppm','Fe_ppm','K_ppm','Mn_ppm','Na_ppm','Sr_ppm']
    elif 'Sr_ppm' in output.columns.values:
        cols = ['section','strat_height','d13C','d18O','cum_strat_height','mts_mtx','Ca_ppm','Mg_ppm',\
                'Al_ppm','Fe_ppm','K_ppm','Mn_ppm','Na_ppm','Sr_ppm']
    else:
        cols = ['section','strat_height','d13C','d18O','cum_strat_height','mts_mtx']
    output = output[cols]

    # sort
    output.sort_values('cum_strat_height', inplace=True)
    output = output.reset_index(drop=True)

    # rename some columns
    if rename==True:
        output.rename(columns={'cum_strat_height':'strat_m', 'mts_mtx':'note'}, inplace=True)

    # round numbers:
    output['d13C'] = np.round(output['d13C'],2)
    output['d18O'] = np.round(output['d18O'],2)
    output['87Sr/86Sr'] = np.round(output['87Sr/86Sr'],6)
    output['Ca_ppm'] = np.round(output['Ca_ppm'])
    output['Mg_ppm'] = np.round(output['Mg_ppm'])
    output['Al_ppm'] = np.round(output['Al_ppm'])
    output['Fe_ppm'] = np.round(output['Fe_ppm'])
    output['K_ppm'] = np.round(output['K_ppm'])
    output['Mn_ppm'] = np.round(output['Mn_ppm'])
    output['Na_ppm'] = np.round(output['Na_ppm'])
    output['Sr_ppm'] = np.round(output['Sr_ppm'])

    # .csv output
    output.to_csv(save, index=False)

    # return output
    if out==True:
        return output

    # plot
    if plot==True:
        colors = []
        for i in range(len(sections)):
            colors.append('GoldenRod')
        c_sr_comp_plot(sections, colors, Sr_threshold, MnSr_threshold)





def create_tonian_csv(sections, labels, save, out):
    """Take in selected sections, output collated data into a single sorted .csv. Only critical columns are saved.

    inputs:
    - sections = ordered list of pandas dataframes
    - labels = ordered list of section labels
    - save = string of save path and file name
    - out = if True, return output as pandas dataframe
    """

    # append the labels column
    for i in range(len(sections)):
        for j in range(len(sections[i].index)):
            sections[i].loc[j,'label'] = labels[i]

    # concatenate
    output = pd.concat(sections, ignore_index=True)
    # sort
    output.sort_values('age', inplace=True)
    output = output.reset_index(drop=True)

    # re-order columns
    cols = ['label','age', 'strat_m', 'd13C', 'd18O','87Sr/86Sr_primary']
    output = output[cols]

    # round numbers:
    output['age'] = np.round(output['age'],3)
    output['strat_m'] = np.round(output['strat_m'],1)
    output['d13C'] = np.round(output['d13C'],2)
    output['d18O'] = np.round(output['d18O'],2)
    output['87Sr/86Sr_primary'] = np.round(output['87Sr/86Sr_primary'],6)

    # .csv output
    output.to_csv(save, index=False)

    # return output
    if out==True:
        return output





def o_c_comp_xplot(dataframes, colors, linregress=True):
    """Plot d18O vs d13C and perform a linear regression for an arbitrary number of sections.

    inputs:
    note that all input lists must be ordered identically
    - dataframes = list of pandas dataframes
    - colors = list of strings of colours
    - linregress = if True, perform linear regression on the data (default True)
    """

    # number of input dataframes
    n = len(dataframes)

    plt.figure()

    # lists to store all C and O data
    cumO = []
    cumC = []

    for i in range(n):
        # if we have both C and O data, plot it
        if 'd18O' in dataframes[i].columns and 'd13C' in dataframes[i].columns:
            plt.scatter(dataframes[i]['d18O'], dataframes[i]['d13C'], color=colors[i])
            # build cumulative C and O arrays
            for j in range(len(dataframes[i]['d13C'])):
                if math.isnan(dataframes[i]['d18O'][j]) and math.isnan(dataframes[i]['d18O'][j]):
                    pass
                else:
                    cumO.append(dataframes[i]['d18O'][j])
                    cumC.append(dataframes[i]['d13C'][j])

    if linregress == True:
        # perform the linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(cumO,cumC)

        # plot the linear regression
        x = np.linspace(min(cumO),max(cumO),num=3)
        y = slope*x + intercept
        plt.plot(x,y)

        # print statistics
        print('y = ' + str(slope) + 'x + ' + str(intercept))
        print('r = ' + str(r_value))
        print('r^2 = ' + str(r_value**2))
        print('p = ' + str(p_value))
        print('std_err = ' + str(std_err))

    plt.ylabel('$\delta^{13}$C')
    plt.xlabel('$\delta^{18}$O')
    plt.xlim(-20,5)
    plt.ylim(-12,12)
    plt.show()





def sr_comp_xplot(dataframes, colors, Sr_threshold, MnSr_threshold):
    """Plot [Sr] vs 87Sr/86Sr for an arbitrary number of sections.

    inputs:
    note that all input lists must be ordered identically
    - dataframes = list of pandas dataframes
    - colors = list of strings of colours
    - Sr_threshold = list of [Sr] thresholds (make 0 if you do not wish to apply)
    - MnSr_threshold = list of [Mn]/[Sr] thresholds (make inf if you do not wish to apply)
    """

    # the number of input dataframes
    n = len(dataframes)

    plt.figure()

    for i in range(n):
        if '87Sr/86Sr' in dataframes[i].columns:
            # create a Mn/Sr column
            for j in range(len(dataframes[i].index)):
                dataframes[i].loc[j,'Mn/Sr'] = dataframes[i]['Mn_ppm'][j] / dataframes[i]['Sr_ppm'][j]
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # plot
                plt.scatter(mtx_df_pass['Sr_ppm'],mtx_df_pass['87Sr/86Sr'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df_pass['Sr_ppm'],mts_df_pass['87Sr/86Sr'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # plot
                plt.scatter(df_pass['Sr_ppm'],df_pass['87Sr/86Sr'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.ylim(0.706, 0.7075)
    plt.xlim(0, 3500)
    plt.ylabel('$^{87}$Sr/$^{86}$Sr')
    plt.xlabel('[Sr] (ppm)')
    plt.show()





def sr_mn_comp_xplot(dataframes, colors, Sr_threshold, MnSr_threshold):
    """Plot [Mn]/[Sr] vs 87Sr/86Sr for an arbitrary number of sections.

    inputs:
    note that all input lists must be ordered identically
    - dataframes = list of pandas dataframes
    - colors = list of strings of colours
    - Sr_threshold = list of [Sr] thresholds (make 0 if you do not wish to apply)
    - MnSr_threshold = list of [Mn]/[Sr] thresholds (make inf if you do not wish to apply)
    """

    # the number of input dataframes
    n = len(dataframes)

    plt.figure()

    for i in range(n):
        if 'Sr_ppm' in dataframes[i].columns:
            # create a Mn/Sr column
            for j in range(len(dataframes[i].index)):
                dataframes[i].loc[j,'Mn/Sr'] = dataframes[i]['Mn_ppm'][j] / dataframes[i]['Sr_ppm'][j]
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # plot
                plt.scatter(mtx_df_pass['Mn/Sr'],mtx_df_pass['87Sr/86Sr'],\
                            marker='o',c=colors[i],edgecolors='none')
                plt.scatter(mts_df_pass['Mn/Sr'],mts_df_pass['87Sr/86Sr'],\
                            marker='^',facecolors='w',edgecolors=colors[i],linewidths=1)
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold[i]) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold[i])]
                # plot
                plt.scatter(df_pass['Mn/Sr'],df_pass['87Sr/86Sr'],\
                            marker='o',c=colors[i],edgecolors='none')
    plt.ylim(0.706, 0.7075)
    plt.xlim(0,5)
    plt.ylabel('$^{87}$Sr/$^{86}$Sr')
    plt.xlabel('[Mn]/[Sr]')
    plt.show()





def sr_variance(dataframes, interval, ymin, ymax):
    """For a given interval of [Sr], determine the variance in 87Sr/86Sr.

    inputs:
    - dataframes = pandas dataframes to be included
    - interval = interval over which variance will be calculated
    - ymin = minimum ylim
    - ymax = maximum ylim
    """

    print('Interval set at ' + str(interval) + 'ppm.')
    n = len(dataframes)

    # build up arrays with all Sr data
    total = np.array([])
    ratio = np.array([])
    for i in range(n):
        if '87Sr/86Sr' in dataframes[i].columns:
            for j in range(len(dataframes[i].index)):
                if np.isfinite(dataframes[i]['Sr_ppm'][j]) and np.isfinite(dataframes[i]['87Sr/86Sr'][j]):
                    total = np.append(total,dataframes[i]['Sr_ppm'][j])
                    ratio = np.append(ratio,dataframes[i]['87Sr/86Sr'][j])

    # sort the lists based on [Sr]
    p = total.argsort()
    total = total[p]
    ratio = ratio[p]

    # work through each interval and determine the variance
    var = np.array([])
    num = np.array([])
    mid = np.array([])
    bot = 0
    while bot <= max(total):
        top = bot + interval
        idx = np.where((total >= bot) * (total < top))
        num = np.append(num, len(ratio[idx]))
        mid = np.append(mid, ((top + bot) / 2))
        if len(ratio[idx]) >= 2:
            var = np.append(var, statistics.variance(ratio[idx]))
        else:
            var = np.append(var, float('nan'))
        bot = top

    # plot [Sr] vs 87Sr/86Sr
    fig, ax1 = plt.subplots(figsize=(4,4))
    ax1.scatter(total,ratio,c='GoldenRod')
    ax1.set_ylabel('$^{87}$Sr/$^{86}$Sr', color='GoldenRod')
    ax1.set_xlabel('[Sr] (ppm)')
    plt.ylim(ymin, ymax)
    for tl in ax1.get_yticklabels():
        tl.set_color('GoldenRod')

    # plot [Sr] vs variance
    ax2 = ax1.twinx()
    ax2.plot(mid,var,'r+-')
    ax2.set_ylabel('variance', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    plt.show()





def sr_xplot(section_data, Sr_threshold=0, MnSr_threshold=float('inf')):
    """Plot [Sr] vs 87Sr/86Sr for a single section.

    inputs:
    - section_data = pandas dataframe
    - Sr_threshold = [Sr] threshold (0 by default)
    - MnSr_threshold = [Mn]/[Sr] threshold (inf by default)
    """

    # create a Mn/Sr column
    for i in range(len(section_data.index)):
        section_data.loc[i,'Mn/Sr'] = section_data['Mn_ppm'][i] / section_data['Sr_ppm'][i]
    # plot MTS data on top as triangles, if it exists
    if 'mts_mtx' in section_data.columns:
        # only select data that passes the thresholds
        mtx_df_pass = section_data[(section_data['mts_mtx']=='mtx') &\
                                    (section_data['Sr_ppm']>=Sr_threshold) &\
                                    (section_data['Mn/Sr']<=MnSr_threshold)]
        mts_df_pass = section_data[(section_data['mts_mtx']=='mts') &\
                                    (section_data['Sr_ppm']>=Sr_threshold) &\
                                    (section_data['Mn/Sr']<=MnSr_threshold)]
        # plot
        plt.scatter(mtx_df_pass['Sr_ppm'],mtx_df_pass['87Sr/86Sr'],\
                    marker='o',c='firebrick',edgecolors='none')
        plt.scatter(mts_df_pass['Sr_ppm'],mts_df_pass['87Sr/86Sr'],\
                    marker='^',facecolors='w',edgecolors='firebrick',linewidths=1)
    else:
        # only select data that passes the thresholds
        df_pass = section_data[(section_data['Sr_ppm']>=Sr_threshold) &\
                                (section_data['Mn/Sr']<=MnSr_threshold)]
        # plot
        plt.scatter(df_pass['Sr_ppm'],df_pass['87Sr/86Sr'],\
                    marker='o',c='firebrick',edgecolors='none')
    plt.ylabel('$^{87}$Sr/$^{86}$Sr')
    plt.xlabel('[Sr] (ppm)')
    plt.ylim(0.706, 0.7075)
    plt.xlim(0, 3500)
    plt.title(str(section_data.section[0]), fontsize=14, fontweight='bold')
    plt.show()





# def tambien_sed_rates(dataframe, kind, ymin=0, ymax=0.1):
#     """Calculate sedimentation rates for each section in the input dataframe, and output as interactive html.
#
#     inputs:
#     - dataframe = pandas dataframe that contains the section data. Must have 'age' column
#     - kind = 'height' or 'age'
#     - ymin = minimum sedimentation rate (mm/yr) to plot (default = 0)
#     - ymax = maximum sedimentation rate(mm/yr) to plot (default = 0.1)
#
#     outputs:
#     - fig = mpld3 generated html figure handle for plotting within Jupyter notebook
#           = e.g. mpld3.display(fig)
#     """
#
#     # create mpld3 plugin to highlight lines - used only for interactive plotting
#     # adapted from: github.com/jakevdp/mpld3/issues/203
#     class HighlightLines(mpld3.plugins.PluginBase):
#         """A plugin to highlight lines on hover"""
#
#         JAVASCRIPT = """
#         mpld3.register_plugin("linehighlight", LineHighlightPlugin);
#         LineHighlightPlugin.prototype = Object.create(mpld3.Plugin.prototype);
#         LineHighlightPlugin.prototype.constructor = LineHighlightPlugin;
#         LineHighlightPlugin.prototype.requiredProps = ["line_ids"];
#         LineHighlightPlugin.prototype.defaultProps = {alpha_bg:0.3, alpha_fg:1.0}
#         function LineHighlightPlugin(fig, props){
#             mpld3.Plugin.call(this, fig, props);
#         };
#
#         LineHighlightPlugin.prototype.draw = function(){
#           for(var i=0; i<this.props.line_ids.length; i++){
#              var obj = mpld3.get_element(this.props.line_ids[i], this.fig),
#                  alpha_fg = this.props.alpha_fg;
#                  alpha_bg = this.props.alpha_bg;
#              obj.elements()
#                  .on("mouseover.highlight", function(d, i){
#                                 d3.select(this).transition().duration(50)
#                                   .style("stroke-opacity", alpha_fg); })
#                  .on("mouseout.highlight", function(d, i){
#                                 d3.select(this).transition().duration(200)
#                                   .style("stroke-opacity", alpha_bg); });
#           }
#         };
#         """
#
#         def __init__(self, line):
#             self.line = line
#             self.dict_ = {"type": "linehighlight",
#                           "line_ids": [mpld3.utils.get_id(line)],
#                           "alpha_bg": line.get_alpha(),
#                           "alpha_fg": 1.0}
#
#     # get the names of sections in the dataframe
#     sections = []
#     for i in range(len(dataframe.index)):
#         if dataframe['section'][i] not in sections:
#             sections.append(dataframe['section'][i])
#     N = len(sections)
#
#     # initialize list which will store lists for each section
#     sed_rates = []
#     mid_heights = []
#     mid_cum_heights = []
#     mid_ages = []
#     for i in range(N):
#         # initialize list which will store the values for each section
#         sed_rate = []
#         mid_height = []
#         mid_cum_height = []
#         mid_age = []
#         # make a temporary dataframe for the section in this iteration
#         section_data = dataframe[dataframe['section']==sections[i]]
#         section_data = section_data.reset_index(drop=True)
#         # remove duplicate data (molar tooth samples, duplicate samples)
#         to_drop = []
#         for j in range(len(section_data.index)-1):
#             if section_data['strat_height'][j+1] == section_data['strat_height'][j]:
#                 to_drop.append(j)
#         section_data = section_data.drop(to_drop)
#         section_data = section_data.reset_index(drop=True)
#         # fill lists
#         for j in range(len(section_data.index)-1):
#             sed_rate.append((section_data['strat_height'][j+1] - section_data['strat_height'][j]) /\
#                             (section_data['age'][j] - section_data['age'][j+1]))
#             sed_rate[j] = sed_rate[j] * 1e-3 #convert units to mm/yr
#             mid_height.append((section_data['strat_height'][j+1] + section_data['strat_height'][j]) / 2)
#             mid_cum_height.append((section_data['strat_m'][j+1] + section_data['strat_m'][j]) / 2)
#             mid_age.append((section_data['age'][j] + section_data['age'][j+1]) / 2)
#         # fill lists that store the lists
#         sed_rates.append(sed_rate)
#         mid_heights.append(mid_height)
#         mid_cum_heights.append(mid_cum_height)
#         mid_ages.append(mid_age)
#
#     # plot the data using mpld3
#     fig, ax = plt.subplots(figsize=(12,5))
#     if kind == 'height':
#         for i in range(N):
#             l, = ax.plot(mid_cum_heights[i],sed_rates[i],color='blueviolet',linewidth=4,alpha=0.1)
#             tooltip = mpld3.plugins.LineLabelTooltip(l, sections[i])
#             mpld3.plugins.connect(fig, HighlightLines(l))
#             mpld3.plugins.connect(fig, tooltip)
#         ax.set_ylabel('sedimentation rate [mm/yr]')
#         ax.set_xlabel('cumulative stratigraphic height [m]')
#         ax.set_ylim(ymin, ymax)
#     elif kind == 'age':
#         for i in range(N):
#             l, = ax.plot(mid_ages[i],sed_rates[i],color='blueviolet',linewidth=4,alpha=0.1)
#             tooltip = mpld3.plugins.LineLabelTooltip(l, sections[i])
#             mpld3.plugins.connect(fig, HighlightLines(l))
#             mpld3.plugins.connect(fig, tooltip)
#         ax.set_ylabel('sedimentation rate [mm/yr]')
#         ax.set_xlabel('age [Ma]')
#         ax.set_ylim(ymin, ymax)
#
#     return fig





def threshold_visualize(dataframes, Sr_threshold, MnSr_threshold, ymin, ymax):
    """Plot [Sr] vs 87Sr/86Sr and [Mn]/[Sr] vs 87Sr/86Sr, with the thresholds marked.

    inputs:
    - dataframes = pandas dataframes to be included
    - Sr_threshold = threshold [Sr] to be visualized
    - MnSr_threshold = threshold [Mn]/[Sr] to be visualized
    - ymin = minimum ylim
    - ymax = maximum ylim
    """

    print('[Sr] threshold set at ' + str(Sr_threshold) + 'ppm.')
    print('[Mn]/[Sr] threshold set at ' + str(MnSr_threshold))
    n = len(dataframes)
    fig, axs = plt.subplots(1, 2, figsize=(9,4), sharey=True)

    # [Sr] vs 87Sr/86Sr
    for i in range(n):
        if '87Sr/86Sr' in dataframes[i].columns:
            # create a Mn/Sr column
            if 'Mn/Sr' not in dataframes[i].columns:
                for j in range(len(dataframes[i].index)):
                    dataframes[i].loc[j,'Mn/Sr'] = dataframes[i]['Mn_ppm'][j] / dataframes[i]['Sr_ppm'][j]
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                mtx_df_fail = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            ((dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                            (dataframes[i]['Mn/Sr']>MnSr_threshold))]
                mts_df_fail = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            ((dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                            (dataframes[i]['Mn/Sr']>MnSr_threshold))]
                # plot
                axs[0].scatter(mtx_df_pass['Sr_ppm'],mtx_df_pass['87Sr/86Sr'],\
                               marker='o',c='g',edgecolors='none',label='primary')
                axs[0].scatter(mtx_df_fail['Sr_ppm'],mtx_df_fail['87Sr/86Sr'],\
                               marker='o',c='r',edgecolors='none',label='altered')
                axs[0].scatter(mts_df_pass['Sr_ppm'],mts_df_pass['87Sr/86Sr'],\
                               marker='^',facecolors='w',edgecolors='g',linewidths=1,label='primary (MTS)')
                axs[0].scatter(mts_df_fail['Sr_ppm'],mts_df_fail['87Sr/86Sr'],\
                               marker='^',facecolors='w',edgecolors='r',linewidths=1,label='altered (MTS)')
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                df_fail = dataframes[i][(dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                        (dataframes[i]['Mn/Sr']>MnSr_threshold)]
                # plot
                axs[0].scatter(df_pass['Sr_ppm'],df_pass['87Sr/86Sr'],\
                               marker='o',c='g',edgecolors='none',label='primary')
                axs[0].scatter(df_fail['Sr_ppm'],df_fail['87Sr/86Sr'],\
                               marker='o',c='r',edgecolors='none',label='altered')
    axs[0].axvline(Sr_threshold,color='b')
    axs[0].set_xlim(0, 3500)
    axs[0].set_ylabel("$^{87}$Sr/$^{86}$Sr")
    axs[0].set_xlabel("[Sr] (ppm)")

    # [Mn]/[Sr] vs 87Sr/86Sr

    # legend checks
    legend_pri = False
    legend_alt = False
    legend_pri_mts = False
    legend_alt_mts = False

    for i in range(n):
        if 'Sr_ppm' in dataframes[i].columns:
            # plot MTS data on top as triangles, if it exists
            if 'mts_mtx' in dataframes[i].columns:
                # only select data that passes the thresholds
                mtx_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                mts_df_pass = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            (dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                            (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                mtx_df_fail = dataframes[i][(dataframes[i]['mts_mtx']=='mtx') &\
                                            ((dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                            (dataframes[i]['Mn/Sr']>MnSr_threshold))]
                mts_df_fail = dataframes[i][(dataframes[i]['mts_mtx']=='mts') &\
                                            ((dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                            (dataframes[i]['Mn/Sr']>MnSr_threshold))]
                # plot
                if legend_pri == False:
                    axs[1].scatter(mtx_df_pass['Mn/Sr'],mtx_df_pass['87Sr/86Sr'],\
                                   marker='o',c='g',edgecolors='none',label='primary')
                    legend_pri = True
                else:
                    axs[1].scatter(mtx_df_pass['Mn/Sr'],mtx_df_pass['87Sr/86Sr'],\
                                   marker='o',c='g',edgecolors='none',label='_nolegend_')
                if legend_alt == False:
                    axs[1].scatter(mtx_df_fail['Mn/Sr'],mtx_df_fail['87Sr/86Sr'],\
                                   marker='o',c='r',edgecolors='none',label='altered')
                    legend_alt = True
                else:
                    axs[1].scatter(mtx_df_fail['Mn/Sr'],mtx_df_fail['87Sr/86Sr'],\
                                   marker='o',c='r',edgecolors='none',label='_nolegend_')
                if legend_pri_mts == False:
                    axs[1].scatter(mts_df_pass['Mn/Sr'],mts_df_pass['87Sr/86Sr'],\
                                   marker='^',facecolors='w',edgecolors='g',linewidths=1,label='primary (MTS)')
                    legend_pri_mts = True
                else:
                    axs[1].scatter(mts_df_pass['Mn/Sr'],mts_df_pass['87Sr/86Sr'],\
                                   marker='^',facecolors='w',edgecolors='g',linewidths=1,label='_nolegend_')
                if legend_alt_mts == False:
                    axs[1].scatter(mts_df_fail['Mn/Sr'],mts_df_fail['87Sr/86Sr'],\
                                   marker='^',facecolors='w',edgecolors='r',linewidths=1,label='altered (MTS)')
                    legend_alt_mts = True
                else:
                    axs[1].scatter(mts_df_fail['Mn/Sr'],mts_df_fail['87Sr/86Sr'],\
                                   marker='^',facecolors='w',edgecolors='r',linewidths=1,label='_nolegend_')
            else:
                # only select data that passes the thresholds
                df_pass = dataframes[i][(dataframes[i]['Sr_ppm']>=Sr_threshold) &\
                                        (dataframes[i]['Mn/Sr']<=MnSr_threshold)]
                df_fail = dataframes[i][(dataframes[i]['Sr_ppm']<Sr_threshold) |\
                                        (dataframes[i]['Mn/Sr']>MnSr_threshold)]
                # plot
                if legend_pri == False:
                    axs[1].scatter(df_pass['Mn/Sr'],df_pass['87Sr/86Sr'],\
                                   marker='o',c='g',edgecolors='none',label='primary')
                    legend_pri = True
                else:
                    axs[1].scatter(df_pass['Mn/Sr'],df_pass['87Sr/86Sr'],\
                                   marker='o',c='g',edgecolors='none',label='_nolegend_')
                if legend_alt == False:
                    axs[1].scatter(df_fail['Mn/Sr'],df_fail['87Sr/86Sr'],\
                                   marker='o',c='r',edgecolors='none',label='altered')
                    legend_alt = True
                else:
                    axs[1].scatter(df_fail['Mn/Sr'],df_fail['87Sr/86Sr'],\
                                   marker='o',c='r',edgecolors='none',label='_nolegend_')
    axs[1].axvline(MnSr_threshold,color='b')
    axs[1].set_xlim(0, 5)
    axs[1].set_ylim(ymin, ymax)
    axs[1].set_xlabel('[Mn]/[Sr]')
    axs[1].legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0)

    plt.show()





def tonian_comp_plot(kind, sections, colors, labels, xlim=None):
    """Plot d13C and 87Sr/86Sr for datasets in the Tonian Composite, for either height or age.

    inputs:
    - kind = 'height' or 'age'
    - sections = ordered list of pandas dataframes
    - colors = ordered list of colors
    - labels = ordered list of section labels
    """

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,6), sharex=True)

    # the d13C data
    for i in range(len(sections)):
        if kind == 'height':
            ax[0].scatter(sections[i]['strat_m'],sections[i]['d13C'],\
                          marker='o',c=colors[i],edgecolors='none',label=labels[i])
        elif kind == 'age':
            ax[0].scatter(sections[i]['age'],sections[i]['d13C'],\
                          marker='o',c=colors[i],edgecolors='none',label=labels[i])
    ax[0].set_ylim(-15,15)
    ax[0].set_ylabel('$\delta^{13}\mathrm{C}$')
    ax[0].legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0)

    # the 87Sr/86Sr data
    for i in range(len(sections)):
        if '87Sr/86Sr' in sections[i].columns or '87Sr/86Sr_primary' in sections[i].columns:
            if kind == 'height':
                ax[1].scatter(sections[i]['strat_m'],sections[i]['87Sr/86Sr_primary'],\
                              marker='o',c=colors[i],edgecolors='none',label=labels[i])
                if xlim != None:
                    ax[1].set_xlim(xlim)
            elif kind == 'age':
                ax[1].scatter(sections[i]['age'],sections[i]['87Sr/86Sr_primary'],\
                              marker='o',c=colors[i],edgecolors='none',label=labels[i])
                if xlim != None:
                    ax[1].set_xlim(xlim)
                else:
                    ax[1].invert_xaxis()
        else:
            print('Sr data not available for ' + labels[i])
    ax[1].set_ylim(0.705,0.709)
    ax[1].set_ylabel('$^{87}\mathrm{Sr}/^{86}\mathrm{Sr}$')
    if kind == 'height':
        ax[1].set_xlabel('composite stratigraphic height (m)')
    elif kind == 'age':
        ax[1].set_xlabel('age (Ma)')

    plt.show()





def section_sed_rates(dataframe, ylims=None):
    """Calculate sedimentation rates for a single section.

    inputs:
    - dataframe = pandas dataframe that contains the section data. Must have 'age' column.
    - ylims = y limits, as [ymin, ymax]
    """
    # initialize list which will store the values
    sed_rate = []
    mid_height = []
    mid_age = []
    # deep copy the dataframe
    dataframe = dataframe.copy()
    # remove duplicate data (molar tooth samples, duplicate samples)
    to_drop = []
    for j in range(len(dataframe.index)-1):
        if dataframe['strat_m'][j+1] == dataframe['strat_m'][j]:
            to_drop.append(j)
    dataframe = dataframe.drop(to_drop)
    dataframe = dataframe.reset_index(drop=True)
    # fill lists
    for j in range(len(dataframe.index)-1):
        sed_rate.append((dataframe['strat_m'][j+1] - dataframe['strat_m'][j]) /\
                        (dataframe['age'][j] - dataframe['age'][j+1]))
        sed_rate[j] = sed_rate[j] * 1e-3 #convert units to mm/yr
        mid_height.append((dataframe['strat_m'][j+1] + dataframe['strat_m'][j]) / 2)
        mid_age.append((dataframe['age'][j] + dataframe['age'][j+1]) / 2)

    # plot
    fig, ax = plt.subplots(figsize=(9,3))
    ax.plot(mid_height, sed_rate)
    ax.set_xlabel('stratigraphic height [m]')
    ax.set_ylabel('sedimentation rate [mm/yr]')
    if ylims != None:
        ax.set_ylim(ylims)
    plt.show(fig)





# def tonian_sed_rates(dataframes, labels, kind):
#     """Calculate sedimentation rates for each dataframe, and output as interactive html.
#
#     inputs:
#     - dataframes = list of pandas dataframes that contains the section data. Must have 'age' column
#     - labels = list of labels for each input dataframe
#     - kind = 'height' or 'age'
#
#     outputs:
#     - fig = mpld3 generated html figure handle for plotting within Jupyter notebook
#           = e.g. mpld3.display(fig)
#     """
#
#     # create mpld3 plugin to highlight lines - used only for interactive plotting
#     # adapted from: github.com/jakevdp/mpld3/issues/203
#     class HighlightLines(mpld3.plugins.PluginBase):
#         """A plugin to highlight lines on hover"""
#
#         JAVASCRIPT = """
#         mpld3.register_plugin("linehighlight", LineHighlightPlugin);
#         LineHighlightPlugin.prototype = Object.create(mpld3.Plugin.prototype);
#         LineHighlightPlugin.prototype.constructor = LineHighlightPlugin;
#         LineHighlightPlugin.prototype.requiredProps = ["line_ids"];
#         LineHighlightPlugin.prototype.defaultProps = {alpha_bg:0.3, alpha_fg:1.0}
#         function LineHighlightPlugin(fig, props){
#             mpld3.Plugin.call(this, fig, props);
#         };
#
#         LineHighlightPlugin.prototype.draw = function(){
#           for(var i=0; i<this.props.line_ids.length; i++){
#              var obj = mpld3.get_element(this.props.line_ids[i], this.fig),
#                  alpha_fg = this.props.alpha_fg;
#                  alpha_bg = this.props.alpha_bg;
#              obj.elements()
#                  .on("mouseover.highlight", function(d, i){
#                                 d3.select(this).transition().duration(50)
#                                   .style("stroke-opacity", alpha_fg); })
#                  .on("mouseout.highlight", function(d, i){
#                                 d3.select(this).transition().duration(200)
#                                   .style("stroke-opacity", alpha_bg); });
#           }
#         };
#         """
#
#         def __init__(self, line):
#             self.line = line
#             self.dict_ = {"type": "linehighlight",
#                           "line_ids": [mpld3.utils.get_id(line)],
#                           "alpha_bg": line.get_alpha(),
#                           "alpha_fg": 1.0}
#
#     # number of sections
#     N = len(dataframes)
#
#     # initialize list which will store lists for each section
#     sed_rates = []
#     mid_heights = []
#     mid_ages = []
#     for i in range(N):
#         # initialize list which will store the values for each section
#         sed_rate = []
#         mid_height = []
#         mid_age = []
#         # get dataframe for the section in this iteration
#         section_data = dataframes[i]
#         # remove duplicate data (molar tooth samples, duplicate samples)
#         to_drop = []
#         for j in range(len(section_data.index)-1):
#             if section_data['strat_m'][j+1] == section_data['strat_m'][j]:
#                 to_drop.append(j)
#         section_data = section_data.drop(to_drop)
#         section_data = section_data.reset_index(drop=True)
#         # fill lists
#         for j in range(len(section_data.index)-1):
#             sed_rate.append((section_data['strat_m'][j+1] - section_data['strat_m'][j]) /\
#                             (section_data['age'][j] - section_data['age'][j+1]))
#             sed_rate[j] = sed_rate[j] * 1e-3 #convert units to mm/yr
#             mid_height.append((section_data['strat_m'][j+1] + section_data['strat_m'][j]) / 2)
#             mid_age.append((section_data['age'][j] + section_data['age'][j+1]) / 2)
#         # fill lists that store the lists
#         sed_rates.append(sed_rate)
#         mid_heights.append(mid_height)
#         mid_ages.append(mid_age)
#
#     # plot the data using mpld3
#     fig, ax = plt.subplots(figsize=(12,5))
#     if kind == 'height':
#         for i in range(N):
#             l, = ax.plot(mid_heights[i],sed_rates[i],color='blueviolet',linewidth=4,alpha=0.1)
#             tooltip = mpld3.plugins.LineLabelTooltip(l, labels[i])
#             mpld3.plugins.connect(fig, HighlightLines(l))
#             mpld3.plugins.connect(fig, tooltip)
#         ax.set_ylabel('sedimentation rate [mm/yr]')
#         ax.set_xlabel('cumulative stratigraphic height [m]')
#         ax.invert_xaxis()
#     elif kind == 'age':
#         for i in range(N):
#             l, = ax.plot(mid_ages[i],sed_rates[i],color='blueviolet',linewidth=4,alpha=0.1)
#             tooltip = mpld3.plugins.LineLabelTooltip(l, labels[i])
#             mpld3.plugins.connect(fig, HighlightLines(l))
#             mpld3.plugins.connect(fig, tooltip)
#         ax.set_ylabel('sedimentation rate [mm/yr]')
#         ax.set_xlabel('age [Ma]')
#         ax.invert_xaxis()
#
#     return fig





def shift_chemostratigraphy(df, x_name, ref_ind, anchor1_ind, anchor2_ind, new_ref):
    """
    Move a data point to a new location, linearly shifting all points between two anchors.

    inputs:
    - df = dataframe which contains the data (edited in place)
    - x_name = string of the column name that has the values to be shifted
    - ref_ind = the index (0-indexed) of the reference point to be shifted (must be in between the anchors)
    - anchor1_ind = the index (0-indexed) of the first anchor point (must be smaller than the second anchor)
    - anchor2_ind = the index (0-indexed) of the second anchor point (must be larger than the first anchor)
    - new_ref = the new x value of the reference point
    """
    # get the original values
    old_ref = df.loc[ref_ind, x_name]
    old_anchor1 = df.loc[anchor1_ind, x_name]
    old_anchor2 = df.loc[anchor2_ind, x_name]

    # difference between old and new reference point
    ref_diff = new_ref - old_ref

    # iterate between our anchors (endpoints inclusive)
    for i in range(anchor1_ind, anchor2_ind+1):

        # denominators
        low_den = old_ref - old_anchor1
        high_den = old_anchor2 - old_ref

        # case where we are on the smaller side of (or on) the reference point
        if i <= ref_ind:
            to_add = ((df.loc[i, x_name] - old_anchor1) / low_den) * ref_diff

        # case where we are on the larger side of the reference point
        else:
            to_add = ((old_anchor2 - df.loc[i, x_name]) / high_den) * ref_diff

        # replace the value
        new_x = df.loc[i, x_name] + to_add
        df.loc[i, x_name] = new_x
