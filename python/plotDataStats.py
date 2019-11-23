#!/usr/bin/python2

import getopt
from itertools import izip
from math import ceil
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from time import clock

if os.path.isdir(os.path.dirname(os.path.realpath(__file__)) + "/../build/pyMods/"):
    sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../build/pyMods/")
if os.path.isdir(os.path.dirname(os.path.realpath(__file__)) + "/../lib/"):
    sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../lib/")
import DataStats

text_size = 20
colour_scheme = ["#D92120","#488BC2","#7FB972","#E6642C","#781C81","#D9AD3C","#BBBBBB","#4065B1"]
#colour_scheme = ['#D92120', '#4065B1', '#7FB972']
#colour_scheme = ['#D92120', '#7FB972', '#4065B1']
fill_colour_scheme = ["#C8100F","#377AB1","#6EA861","#D5531B","#670B7)","#C89C2B","#AAAAAA","#2F54A0"]
#fill_colour_scheme = ['#E6642C', '#488BC2', '#63AD99']
#fill_colour_scheme = ['#E6642C', '#63AD99', '#488BC2']

nucleotide_colour_scheme = ['#E6642C','#488BC2','#781C81','#B5BD4C','#BEBEBE']
baseline_colour = '#781C81';

plot_legend = True;
#line_width = 2
line_width = 4

def getMinMaxX(hists,zero_padding=True, cov_plot=False):
    min_x = sys.maxint
    max_x = 0
    
    if cov_plot:
        # Use only the first data (real data not simulated) for x values
        min_x = hists[0][0]
        max_x = hists[0][0]+len(hists[0][1])
        
        # Cut at 5 times the mean
        val = 0
        median_count = sum(hists[0][1])/2
        count = 0
        median = -1
        while count < median_count:
            median += 1
            count += hists[0][1][median]
        max_x = min(max_x, min_x+5*int(median) )

    else:
        for hist in hists:
            min_x = min( min_x, hist[0] )
            max_x = max( max_x, hist[0]+len(hist[1]) )
            pass
    
    # Pad the hists with a zero at the beginning and the end
    if zero_padding:
        min_x -= 1
        max_x += 1
        pass
    
    return (min_x,max_x)

def getMinMaxY(hists, cov_plot=False):
    min_y = sys.maxint
    max_y = 0
    
    if cov_plot:
        # Do not count x=0
        for hist in hists:
            min_y = min( min_y, min(hist[1][1:]) )
            max_y = max( max_y, max(hist[1][1:]) )
            pass
    else:
        for hist in hists:
            min_y = min( min_y, min(hist[1]) )
            max_y = max( max_y, max(hist[1]) )
            pass
    
    return (min_y * 0.9,max_y * 1.1)

def padList(data, min_x, max_x, normalize=False, normVector=False):
    if(0 == len(data[1])):
        return [0]*(max_x-min_x)
    else:
        if normalize:
            normSum = float(sum(data[1]))/100
            center_list = [x/normSum for x in data[1]]
            pass
        elif normVector:
            center_list = list(data[1])
            for x in range(len(data[1])):
                if normVector[x]:
                    center_list[x] = float(center_list[x])*100/normVector[x]
                    pass
                else:
                    center_list[x] = 0.0
                    pass
                pass
            pass
        else:
            center_list = list(data[1])
            pass
        
        if min_x <= data[0] and max_x >= data[0]+len(data[1]):
            # Pad
            return [0]*(data[0]-min_x) + center_list + [0]*(max_x-data[0]-len(data[1]))
        else:
            # Cut + Pad
            return [0]*max(0,data[0]-min_x) + center_list[max(0,min_x-data[0]):min(max_x-data[0], len(center_list))] + [0]*max(0,max_x-data[0]-len(data[1]))
            

def num_subplots( num_plots ):
    num_x_subplots = int( ceil( np.sqrt( num_plots ) ) )
    num_y_subplots = int( ceil( float(num_plots)/num_x_subplots) )
    return (num_x_subplots, num_y_subplots)

def finalizePlot(xtitle, ytitle):
    plt.xlabel(xtitle, fontsize=text_size)
    plt.ylabel(ytitle, fontsize=text_size)
    plt.tight_layout()
    return

def finalizeMultiPlot(fig, axs, xtitle, ytitle):
    fig.tight_layout(rect=[0.04, 0.04, 1, 1])

    fig.text(0.54, 0.04, xtitle, ha='center', va='center', fontsize=text_size)
    fig.text(0.04, 0.54, ytitle, ha='center', va='center', rotation='vertical', fontsize=text_size)
    return

def finalizeMultiPlotWithColorbar(fig, xtitle, ytitle, cmap, norm, horizontal):
    # Add colorbar and axis titles
    if horizontal:
        fig.subplots_adjust(bottom=0.2)
        cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.05])
        cb_orientation = 'horizontal'
        fig.text(0.5, 0.14, xtitle, ha='center', va='center', fontsize=text_size)
        pass
    else:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb_orientation = 'vertical'
        fig.text(0.5, 0.04, xtitle, ha='center', va='center', fontsize=text_size)
        pass
    ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation=cb_orientation)
    fig.text(0.06, 0.5, ytitle, ha='center', va='center', rotation='vertical', fontsize=text_size)
    return

def setAxis(ax, log, ylimits=None):
        # In case no axis is defined get the standard one
    if plt == ax:
        ax = plt.gca()
        pass

    if ylimits:
        ax.set_ylim(ylimits)

    if log:            
        ax.set_yscale('log', nonposy='clip')
        pass
    else:
        # Set at what range of exponents they are not plotted in exponential format for example (-3,4): [0.001-1000[
        ax.get_yaxis().get_major_formatter().set_powerlimits((-3,4))
        pass
    
    return

def plotMinimal(ax, names, hists, min_x, max_x, x_vals, log=False, marker=None, linewidth=line_width, linestyle='-', normalize=False, normVectors=False, cov_plot=False):
    if not normVectors:
        normVectors = [False] * len(names)
        pass
    
    for (name, hist, col, norm) in izip(names, hists, colour_scheme, normVectors):
        if(0 == len(hist[1]) or 0 == max(hist[1])):
            log = False
            pass
        ax.plot(x_vals, padList(hist, min_x, max_x, normalize=normalize, normVector=norm), linewidth=linewidth, marker=marker, linestyle=linestyle, label=name, color=col)
        pass

    if cov_plot:
        ylimits = getMinMaxY(hists, cov_plot=True)
    else:
        ylimits = None

    setAxis(ax, log, ylimits)

    return

def plotBackend(xtitle, ytitle, names, hists, ax=plt, zero_padding=True, log=False, shift_x=0, baseline=None, normalize=False, normVectors=False, cov_plot=False):
    (min_x, max_x) = getMinMaxX(hists, zero_padding, cov_plot=cov_plot)

    if 1 < max_x-min_x: # Otherwise the plot has only one data point and no lines can be plotted
        x_vals = range(min_x+shift_x, max_x+shift_x)

        if baseline:
            ax.fill_between(x_vals, padList(baseline[1], min_x, max_x), label=baseline[0], color=baseline_colour)
            pass
        
        plotMinimal(ax, names, hists, min_x, max_x, x_vals, log=log, normalize=normalize, normVectors=normVectors, cov_plot=cov_plot)
    
        finalizePlot(xtitle, ytitle)
    
        return True
    else:
        return False

def plot(pdf, xtitle, ytitle, names, hists, zero_padding=True, log=False, shift_x=0, baseline=None, legend='upper right', normalize=False, normVectors=False, cov_plot=False):
    plt.close()
    
    if plotBackend(xtitle, ytitle, names, hists, zero_padding=zero_padding, log=log, shift_x=shift_x, baseline=baseline, normalize=normalize, normVectors=normVectors, cov_plot=cov_plot):
        if plot_legend:
            plt.legend(loc=legend, shadow=True)
            pass

        pdf.savefig()
        pass
    return

def plotBaseQualities(pdf, xtitle, ytitle, names, means, minima, first_quartiles, medians, third_quartiles, maxima):
    boxplot_width=0.4 # Full width of boxplot is 2*boxplot_width
    
    plt.close()
    fig, ax = plt.subplots()
     
    (min_x,max_x) = getMinMaxX(means, zero_padding=False)
    x_range = xrange(min_x+1,max_x+1)
    
    for n in xrange( len(means) ):
        # Shift x based on n, so boxes only half overlap
        x_left = (2*float(n)/(len(means)+1)-1)*boxplot_width
        x_width = float(2)/(len(means)+1)*2*boxplot_width

        for (x_i, ymin_i, ymax_i) in izip(x_range, padList(first_quartiles[n], min_x, max_x), padList(third_quartiles[n], min_x, max_x)):
            ax.add_patch(
                patches.Rectangle(
                    (x_i+x_left, ymin_i),   # (x,y)
                    x_width,                # width
                    ymax_i - ymin_i,        # height
                    facecolor=fill_colour_scheme[n],
                    linewidth=0
                )
            )
            pass
        pass
     
    for n in xrange( len(means) ):
        # Horizontal lines for minima and maxima reach out to the same size as the boxes
        x_left = (2*float(n)/(len(means)+1)-1)*boxplot_width
        x_right = (2*float(n+2)/(len(means)+1)-1)*boxplot_width

        for x,y in izip(x_range, padList(minima[n], min_x, max_x)):
            ax.hlines(y=y, xmin=x+x_left, xmax=x+x_right, linewidth=1)
            pass
        
        for x,y in izip(x_range, padList(maxima[n], min_x, max_x)):
            ax.hlines(y=y, xmin=x+x_left, xmax=x+x_right, linewidth=1)
            pass
        pass
  
    for n in xrange( len(means) ):
        # Horizontal lines for medians reach out to the same size as the boxes
        x_left = (2*float(n)/(len(means)+1)-1)*boxplot_width
        x_right = (2*float(n+2)/(len(means)+1)-1)*boxplot_width

        for x,y in izip(x_range, padList(medians[n], min_x, max_x)):
            ax.hlines(y=y, xmin=x+x_left, xmax=x+x_right, linewidth=2, color=colour_scheme[n])
            pass
        pass
 
    # matplotlib patches don't automatically impact the scale of the ax, so we manually autoscale the x and y axes
    ax.autoscale_view()

    plotBackend(xtitle, ytitle, names, means, ax=ax, zero_padding=False, shift_x=1)
    if plot_legend:
        plt.legend(loc='center left', shadow=True)
        pass

    pdf.savefig()
    return

def plotTileAbundance(pdf, xtitle, ytitle, names, abundances, tiles):
    plt.close()

    common_tiles = {}
    for n, single_tiles in enumerate(tiles):
        for tile, ab_value in izip(single_tiles, abundances[n]):
            if tile not in common_tiles.keys():
                common_tiles[tile] = [0]*len(names);
                pass
            common_tiles[tile][n] = ab_value;
            pass
        pass

    if( 1 < len( common_tiles ) ): # Do not plot if it is just one tile
        tile_list = [];
        sorted_abundances = [[] for x in xrange(0,len(names))]
        for tile, abundance in common_tiles.items():
            tile_list.append(tile);
            for n, ab_value in enumerate(abundance):
                sorted_abundances[n].append(ab_value)
                pass
            pass
        
        x_vals = range(0,len(tile_list))
        
        for (name, abundance, col) in izip(names, sorted_abundances, colour_scheme):
            plt.plot(x_vals, abundance, linewidth=2, marker=None, linestyle='-', label=name, color=col)
            pass
        
        plt.xticks( x_vals, tile_list, rotation='vertical' )
        
        plt.autoscale()
        
        finalizePlot(xtitle, ytitle)
        if plot_legend:
            plt.legend(loc='lower center', shadow=True)
            pass
        
        pdf.savefig()
        pass

    return

def plotTileQuality(pdf, xtitle, ytitle, ztitle, names, tile_dicts):
    plt.close()
    min_x = sys.maxint
    max_x = 0
    max_z = 0

    tile_names = set()
    for ndict, dict in enumerate(tile_dicts):
        tile_names = tile_names.union(dict.keys())
 
        for key, hist in dict.items():
            min_x = min( min_x, hist[0] )
            max_x = max( max_x, hist[0]+len(hist[1]) )

            hist = (hist[0], [z for z in hist[1]])
            max_z = max([abs(z) for z in hist[1]+[max_z]])
            
            tile_dicts[ndict][key] = hist # Store convertion for later
            pass
        pass
    tile_names = sorted(tile_names)
    
    if 1 < len(tile_names):
        # Equalize min_z and max_z (except of the sign), so that 0 is in the middle and white
        min_z = -max_z
        z_value_range = 2*max_z
    
        x_width = float(1) # Read position is spaced in integeres
        y_width = float(1) # The spacing between tiles is the spacing between there ids which is 1
    
        # Define colour map
        RdWeGn = LinearSegmentedColormap.from_list('RdWeGn', ['#D92120', '#FFFFFF','#7FB972'])
        plt.cm.register_cmap(cmap=RdWeGn)
        z_axis_norm = Normalize(vmin=min_z, vmax=max_z)
    
        (num_x_subplots, num_y_subplots) = num_subplots( len(tile_dicts) )
        fig, axs = plt.subplots(num_y_subplots, num_x_subplots, sharex='col', sharey='row')
    
        subplot_x = 0
        subplot_y = 0
        for dict, name in izip(tile_dicts, names):
            if 1 == len(tile_dicts):
                ax = axs
                pass
            elif 2 == len(tile_dicts):
                ax = axs[subplot_x]
                pass
            else:
                ax = axs[subplot_y][subplot_x]
                pass
            ax.set_title(name)
            
            z_vals = np.zeros( len(tile_names)*(max_x-min_x) ).reshape( (len(tile_names), max_x-min_x) )
            patch_y_coords = []
            for tile_id, tile in enumerate(reversed(tile_names)):
                if tile in dict:
                    z_vals[tile_id] = np.array(padList(dict[tile], min_x, max_x))
                    pass
                else:
                    if not patch_y_coords or patch_y_coords[-1][0] != tile_id - patch_y_coords[-1][1]:
                        patch_y_coords.append( [tile_id,1] )
                        pass
                    else:
                        patch_y_coords[-1][1] += 1
                        pass
                    pass
                pass
            
            for (y, n_tiles) in patch_y_coords:
                ax.add_patch(
                    patches.Rectangle(
                        (min_x-x_width/2, y-y_width/2),   # (x,y)
                        max_x-min_x + x_width,            # width
                        y_width*n_tiles,                  # height
                        facecolor="grey",
                        linewidth=0
                    )
                )
                pass
    
            # Plot 2d histogram
            ax.imshow(z_vals, interpolation='nearest', aspect='auto', norm=z_axis_norm, cmap=RdWeGn) #
    
            ax.set_yticks( xrange( 0, len(tile_names) ) )
            if 0 == subplot_x:
                ax.set_yticklabels(reversed(tile_names))
                pass
            
            # matplotlib patches don't automatically impact the scale of the ax, so we manually autoscale the x and y axes
            #ax.autoscale_view()
        
            subplot_x += 1
            subplot_y += subplot_x/num_x_subplots
            subplot_x %= num_x_subplots
            pass
        
        # Add colorbar and axis titles
        finalizeMultiPlotWithColorbar(fig, xtitle, ytitle, RdWeGn, z_axis_norm, (num_x_subplots > num_y_subplots))

        pdf.savefig()
    return

def plotSequenceQualities(pdf, xtitle, ytitle, names, means, minima, first_quartiles, medians, third_quartiles, maxima):
    plt.close()
    
    min_x = getMinMaxX(minima,zero_padding=False)[0]
    max_x = getMinMaxX(maxima,zero_padding=False)[1]
    x_range = xrange(min_x,max_x)
    
    fig, ((ax_mean, ax_first_quartile, ax_minimum), (ax_median, ax_third_quartile, ax_maximum)) = plt.subplots(2,3, sharex='col', sharey='row')
    
    plotMinimal(ax_mean, names, means, min_x, max_x, x_range, normalize=True)
    if plot_legend:
        ax_mean.legend(loc='upper center',frameon=False)
        pass
    plotMinimal(ax_first_quartile, names, first_quartiles, min_x, max_x, x_range, normalize=True)
    plotMinimal(ax_minimum, names, minima, min_x, max_x, x_range, normalize=True)
    plotMinimal(ax_median, names, medians, min_x, max_x, x_range, normalize=True)
    plotMinimal(ax_third_quartile, names, third_quartiles, min_x, max_x, x_range, normalize=True)
    plotMinimal(ax_maximum, names, maxima, min_x, max_x, x_range, normalize=True)
    
    ax_mean.set_title("Mean", fontsize=text_size)
    ax_first_quartile.set_title("First quartile", fontsize=text_size)
    ax_minimum.set_title("Minimum", fontsize=text_size)
    ax_median.set_title("Median", fontsize=text_size)
    ax_third_quartile.set_title("Third quartile", fontsize=text_size)
    ax_maximum.set_title("Maximum", fontsize=text_size)
    
    finalizeMultiPlot(fig, [ax_mean, ax_first_quartile, ax_minimum, ax_median, ax_third_quartile, ax_maximum], xtitle, ytitle)
    
    pdf.savefig()
    return

def plot2dQuality(pdf, xtitle, ytitle, names, quals, axis_lable_dist = 5):
    zero_log_val = 0.01
    plt.close()
    min_x = sys.maxint
    max_x = 0
    min_y = sys.maxint
    max_y = 0
    for qual2d in quals:
        min_x = min( min_x, qual2d[0] )
        max_x = max( max_x, qual2d[0]+len(qual2d[1]) )
        
        for qual in qual2d[1]:
            if( len(qual[1]) ): # If there are no entries, don't set min_y to 0 because offset default is 0
                min_y = min( min_y, qual[0] )
                max_y = max( max_y, qual[0]+len(qual[1]) )
                pass
            pass
        pass

    if sys.maxint != min_y: # Otherwise the whole plot is empty        
        qual_maps = []
        
        min_z = sys.maxint
        max_z = 0
        for qual2d in quals:
            z_vals = np.zeros( (max_y-min_y)*(max_x-min_x) ).reshape( (max_y-min_y, max_x-min_x) ) + zero_log_val
    
            x_offset = qual2d[0]-min_x
            for x, qual in enumerate(qual2d[1]):
    
                #y_offset = qual[0]-min_y
                y_offset = max_y-qual[0]-len(qual[1]) # Offset is reversed
                for y, q in enumerate(reversed(qual[1])):
                    if 0 == q:
                        q = zero_log_val # Overwrite all zeros with a super small value, that cannot be a normal count
                        pass
                    z_vals[y+y_offset][x+x_offset] = q
    
                    min_z = min( min_z, q )
                    max_z = max( max_z, q )
                    pass
                pass
            
            qual_maps.append(z_vals)
            pass
    
        (num_x_subplots, num_y_subplots) = num_subplots( len(quals) )
        fig, axs = plt.subplots(num_y_subplots, num_x_subplots, sharex='col', sharey='row')
    
        subplot_x = 0
        subplot_y = 0
        for map, name in izip(qual_maps, names):
            if 1 == len(qual_maps):
                ax = axs
                pass
            elif 2 == len(qual_maps):
                ax = axs[subplot_x]
                pass
            else:
                ax = axs[subplot_y, subplot_x]
                pass
            ax.set_title(name)
    
            # Plot 2d histogram
            ax.imshow(map, interpolation='nearest', aspect='auto', vmin=min_z, vmax=max_z, cmap=plt.cm.magma, norm=LogNorm()) #
    
            # y Ticks have to be reversed so they increase from bottom to top and shifted since they always start at 0 and should start at the first multiple of axis_lable_dist above min_y
            ax.set_yticks( xrange( (max_y-1)%axis_lable_dist, max_y-min_y, axis_lable_dist ) )
            if 0 == subplot_x:
                ax.set_yticklabels( reversed( xrange( ((min_y-1)/axis_lable_dist+1)*axis_lable_dist, max_y, axis_lable_dist ) ) ) 
                pass
    
            # x Ticks only have to be shifted by min_x, since they already go from left to right
            ax.set_xticks( xrange( axis_lable_dist-((min_x-1)%axis_lable_dist+1), max_x-min_x, axis_lable_dist ) )
            if num_y_subplots-1 == subplot_y:
                ax.set_xticklabels( xrange( ((min_x-1)/axis_lable_dist+1)*axis_lable_dist, max_x, axis_lable_dist ) ) 
                pass
        
            subplot_x += 1
            subplot_y += subplot_x/num_x_subplots
            subplot_x %= num_x_subplots
            pass
    
        finalizeMultiPlot(fig, axs, xtitle, ytitle)

        pdf.savefig()
        pass
    return

def nucleotidePlot(pdf, xtitle, ytitle, names, hists, zero_padding=True, log=False, shift_x=0, normalize=False, pos_normalize=False, Nlegend=0, legend='upper center'):
    plt.close()

    min_x = sys.maxint
    max_x = 0
    
    for nuc_hists in hists:
        for hist in nuc_hists:
            min_x = min( min_x, hist[0] )
            max_x = max( max_x, hist[0]+len(hist[1]) )
            pass
        pass
    
    if normalize:
        normed_hists=[]
        for nuc_hists in hists:
            normed_hists.append([])
            for hist in nuc_hists:
                total = float(sum(hist[1]))
                normed_hists[-1].append( (hist[0], [x/total*100.0 for x in hist[1]]) )
                pass
            pass
        pass
    elif pos_normalize:
        # Update max_x to remove fluctuations at the end due to indels
        min_total = 0
        while max_x >= min_x and 50 > min_total:
            max_x -= 1
            min_total = sys.maxint
            for dataset in range(len(names)):
                total=0

                for hist in hists:
                    if max_x >= hist[dataset][0] and max_x < hist[dataset][0]+len(hist[dataset][1]):
                        total += hist[dataset][1][max_x-hist[dataset][0]]
                        pass
                    pass
                
                if total < min_total:
                    min_total = total
                    pass
                pass
            pass
        max_x += 1
                
        normed_hists=[]
        for nuc_hists in hists:
            normed_hists.append([])
            
            for hist in nuc_hists:
                normed_hists[-1].append( (min_x, [0.0] * (max_x-min_x)) )
                pass
            pass
        pass

        for dataset in range(len(names)):
            for x in xrange(min_x,max_x):
                total = 0;

                for hist in hists:
                    if x >= hist[dataset][0] and x < hist[dataset][0]+len(hist[dataset][1]):
                        total += hist[dataset][1][x-hist[dataset][0]]
                        pass
                    pass
                
                for hist, n_hist in izip(hists,normed_hists):
                    if x >= hist[dataset][0] and x < hist[dataset][0]+len(hist[dataset][1]) and 0 < total:
                        n_hist[dataset][1][x-min_x] = float(hist[dataset][1][x-hist[dataset][0]])*100/total
                        pass
                    pass
                pass
            pass
        pass
    else:
        normed_hists=hists
        pass

    
    if 1 < max_x-min_x: # Otherwise the plot is empty
        # Pad the hists with a zero at the beginning and the end
        if zero_padding:
            min_x -= 1
            max_x += 1
            pass
    
        x_vals = range(min_x+shift_x, max_x+shift_x)

        min_y = sys.maxint
        max_y = 0

        fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
        for N, nuc in enumerate(['A','C','G','T']):
            ax = axs[N/2, N%2]

            plotMinimal(ax, names, normed_hists[N], min_x, max_x, x_vals, log=log)
            
            tmp_min, tmp_max = ax.get_ylim()
            if tmp_min < min_y:
                min_y = tmp_min;
                pass
            if tmp_max > max_y:
                max_y = tmp_max;
                pass
        
            ax.set_title(nuc)
            pass
        
        for N in range(4):
            axs[N/2, N%2].set_ylim([min_y, max_y])
            pass
    
        finalizeMultiPlot(fig, axs.flatten, xtitle, ytitle)
        if plot_legend:
            if( "" == legend ):
                axs[Nlegend/2, Nlegend%2].legend(bbox_to_anchor=(1.1, 1.1))
                pass
            else:
                axs[Nlegend/2, Nlegend%2].legend(loc=legend)
                pass
            pass
    
        pdf.savefig()
        pass
    return

def plotCalledBases( pdf, xtitle, ytitle, names, called_bases, total_bases, include_unknown=True ):
    plt.close()
    
    nucleotides = ['A','C','G','T','N']

    min_x = sys.maxint
    max_x = 0
    for stats in total_bases:
        for call_nucs in stats:
            if(len(call_nucs[1])):
                min_x = min(min_x, call_nucs[0])
                max_x = max(max_x, call_nucs[0]+len(call_nucs[1]))
                pass
            pass
        pass

    normed_called_bases = [];
    ref_bases_unmapped = [];
    for stats, tot in izip(called_bases,total_bases):

        normed_stats = [ [None] * 5 for i in range(4) ]
        unmapped = [ np.zeros(max_x-min_x) for i in range(4) ]
        for call_base in xrange(5):
            total = np.zeros(max_x-min_x)
            for ref_base in xrange(4):
                normed_stats[ref_base][call_base] = np.array(padList(stats[ref_base][call_base], min_x, max_x), dtype=float)
                total += normed_stats[ref_base][call_base]
                pass
             
             
            left_over = np.array(padList(tot[call_base], min_x, max_x)) - total
 
            # Replace all 0 by 1 to avoid division by zero as 0 values should anyways just stay 0
            for i in range(len(total)):
                if 0 == total[i]:
                    total[i] = 1
                    pass
                pass
             
            for ref_base in xrange(4):
                # Split called bases from not mapped reads according to the proportions in the mapped reads along the ref_base
                unmapped[ref_base] += left_over * normed_stats[ref_base][call_base] / total
                pass
            pass
         
        normed_called_bases.append(normed_stats)
        ref_bases_unmapped.append(unmapped)
        pass


    for stats, bases_unmapped in izip(normed_called_bases, ref_bases_unmapped):
        for ref_nucs, unmapped in izip(stats, bases_unmapped):
            if( include_unknown ):
                total = unmapped #unmapped will be changed afterwards, but don't care
                pass
            else:
                total = np.zeros(max_x-min_x)
                pass
            
            for call_nucs in ref_nucs:
                total += call_nucs
                pass
            
            # Replace all 0 by 1 to avoid division by zero as 0 values should anyways just stay 0
            for i in range(len(total)):
                if 0 == total[i]:
                    total[i] = 1
                    pass
                pass

            for call_nucs in ref_nucs:
                call_nucs /= total / 100;
                pass

            pass
        pass
    
    min_y_value = 100; # Use the lowest value of correct calling as the minimum for the y axis, so erroneous callings can be better seen
    for stats in normed_called_bases:
        for base in range(4):
            for val in stats[base][base]:
                if val < min_y_value and 0 != val:
                    min_y_value = val
                    pass
                pass
            pass
        pass

    x_vals = np.insert(np.arange(min_x,max_x)+0.5, np.arange(max_x-min_x), np.arange(min_x,max_x)-0.5)
    x_vals[0] = min_x
    x_vals[-1] = max_x-1

    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    for ref_base in range(4):
        ax = axs[ref_base/2, ref_base%2]

        for n, nc_bases in enumerate(normed_called_bases):
            total = np.zeros( len(nc_bases[ref_base][0]) )
            for call_base in range(0,5):
                total += nc_bases[ref_base][call_base]
                pass
            for call_base in reversed([ref_base] + range(ref_base) + range(ref_base+1,5)): # Put current ref_base at last position on call_base
                if 0 == n: # First stats taken as reference and plotted as filled
                    ax.fill_between(x_vals, np.repeat(total, 2), color=nucleotide_colour_scheme[call_base])
                    pass
                else: # Other stats only plotted as lines
                    ax.plot(x_vals, np.repeat(total, 2), color=colour_scheme[n])
                    pass
                total -= nc_bases[ref_base][call_base]
                pass         
            pass
        pass
    
        ax.set_title(nucleotides[ref_base])
        ax.title.set_color(nucleotide_colour_scheme[ref_base])
        ax.set_ylim([min_y_value-1,100])
        pass

    finalizeMultiPlot(fig, axs, xtitle, ytitle)
    pdf.savefig()

    return

def plotMeanVsDispersion(pdf, xtitle, ytitle, names, mean_lists, dispersion_lists):
    plt.close()
        
    (num_x_subplots, num_y_subplots) = num_subplots( len(names) )
    fig, axs = plt.subplots(num_y_subplots, num_x_subplots, sharex='col', sharey='row')

    mean95 = 0;
    for mean in mean_lists:
        mean95 = max(mean95, sorted(mean)[-len(mean)/20])
    dispersion95 = 0;
    for dispersion in dispersion_lists:
        dispersion95 = max(dispersion95, sorted(dispersion)[-len(dispersion)/20])

    subplot_x = 0
    subplot_y = 0
    for mean, dispersion, name in izip(mean_lists, dispersion_lists, names):
        if 1 == len(names):
            ax = axs
            pass
        elif 2 == len(names):
            ax = axs[subplot_x]
            pass
        else:
            ax = axs[subplot_y, subplot_x]
            pass
        ax.set_title(name)

        # Plot 2d histogram
        ax.hist2d(mean, dispersion, bins=100, range=[[0.0, mean95], [0.0, dispersion95]], cmap=plt.cm.magma, norm=LogNorm())
    
        subplot_x += 1
        subplot_y += subplot_x/num_x_subplots
        subplot_x %= num_x_subplots
        pass

    finalizeMultiPlot(fig, axs, xtitle, ytitle)

    pdf.savefig()

def plotDataStats(statsFiles, oFile):
    names = []
    for sf in statsFiles:
        if 'gz' == sf.split('.')[-1]:
            names.append( sf.split('.',2)[0] )
            pass
        else:
            names.append( sf.rsplit('.',1)[0] )
            pass
        pass
    
    if not oFile:
        oFile = names[0] + ".pdf"
        pass

    for i, name in enumerate(names):
        if '/' in name:
            names[i] = name.rsplit('/',1)[1].split('.')[0][0:10]
            pass
        pass

    print "Start loading files: ", clock()
    stats = [];
    for sf in statsFiles:
        st = DataStats.DataStatsInterface(None)
        if st.Load(sf):
            stats.append(st)
            pass
        else:
            print "Error loading statistics from {0}".format(sf)
        pass

    if stats:
        with PdfPages(oFile) as pdf:
            plt.ioff()
            print "Start plotting files: ", clock()
            
            plot( pdf, "Coverage", "# bases of reference", names, [st.Coverage() for st in stats], zero_padding=False, cov_plot=True )
            #plot( pdf, "Coverage + strand", "# bases of reference", names, [st.CoverageStranded(False) for st in stats], zero_padding=False )
            #plot( pdf, "Coverage - strand", "# bases of reference", names, [st.CoverageStranded(True) for st in stats], zero_padding=False )
            plot( pdf, "Coverage + strand / total Coverage [%]", "# bases of reference", names, [st.CoverageStrandedPercent(False) for st in stats], zero_padding=False )
            plot( pdf, "Coverage - strand / total Coverage [%]", "# bases of reference", names, [st.CoverageStrandedPercent(True) for st in stats], zero_padding=False )
            #plot( pdf, "Coverage + strand / total Coverage [%]", "# bases of reference (min cov 10)", names, [st.CoverageStrandedPercentMinCov10(False) for st in stats], zero_padding=False )
            #plot( pdf, "Coverage - strand / total Coverage [%]", "# bases of reference (min cov 10)", names, [st.CoverageStrandedPercentMinCov10(True) for st in stats], zero_padding=False )
            #plot( pdf, "Coverage + strand / total Coverage [%]", "# bases of reference (min cov 20)", names, [st.CoverageStrandedPercentMinCov20(False) for st in stats], zero_padding=False )
            #plot( pdf, "Coverage - strand / total Coverage [%]", "# bases of reference (min cov 20)", names, [st.CoverageStrandedPercentMinCov20(True) for st in stats], zero_padding=False )
            
            plot( pdf, "Fragment duplication number", "% fragments", names, [st.FragmentDuplicationNumber() for st in stats], log=True, normalize=True )
            plotMeanVsDispersion( pdf, "Mean", "Dispersion", names, [st.MeanList() for st in stats], [st.DispersionList() for st in stats] )
            
            plot( pdf, "Reference sequence", "# pairs", names, [(0,st.Abundance()) for st in stats] )
            plot( pdf, "Insert length", "% read pairs", names, [st.InsertLengths() for st in stats], normalize=True )
            
            plot( pdf, "GC read content[%] (first)", "% reads", names, [st.GCReadContent(0) for st in stats], zero_padding=False, normalize=True )
            plot( pdf, "GC read content[%] (second)", "% reads", names, [st.GCReadContent(1) for st in stats], zero_padding=False, normalize=True )
            #plot( pdf, "GC read content on reference[%] (first)", "% reads", names, [st.GCReadContentReference(0) for st in stats], zero_padding=False, normalize=True )
            #plot( pdf, "GC read content on reference[%] (second)", "% reads", names, [st.GCReadContentReference(1) for st in stats], zero_padding=False, normalize=True )
            #plot( pdf, "GC read content mappped[%] (first)", "% reads", names, [st.GCReadContentMapped(0) for st in stats], zero_padding=False, normalize=True )
            #plot( pdf, "GC read content mapped[%] (second)", "% reads", names, [st.GCReadContentMapped(1) for st in stats], zero_padding=False, normalize=True )
            plot( pdf, "GC fragment content[%]", "% fragments", names, [st.GCFragmentContent() for st in stats], zero_padding=False, normalize=True )
            plot( pdf, "GC fragment content[%]", "% fragments / # sites", names, [st.GCFragmentContentBias() for st in stats], zero_padding=False, normalize=True )
            #plot( pdf, "N content (first)", "# reads", names, [st.NContent(0) for st in stats], zero_padding=False, log=True )
            #plot( pdf, "N content (second)", "# reads", names, [st.NContent(1) for st in stats], zero_padding=False, log=True )
            nucleotidePlot( pdf, "Read position (first)", "% A/C/G/T at position", names, [ [st.SequenceContent(0,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            nucleotidePlot( pdf, "Read position (second)", "% A/C/G/T at position", names, [ [st.SequenceContent(1,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            nucleotidePlot( pdf, "Read Position", "Preference for A/C/G/T", names, [ [(0, st.FragmentSurroundingBiasByBase(n)) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=-10)
            #nucleotidePlot( pdf, "Read position (first forward)", "% A/C/G/T reference at pos", names, [ [st.SequenceContentReference(0,0,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            #nucleotidePlot( pdf, "Read position (first reverse)", "% A/C/G/T reference at pos", names, [ [st.SequenceContentReference(0,1,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            #nucleotidePlot( pdf, "Read position (second forward)", "% A/C/G/T reference at pos", names, [ [st.SequenceContentReference(1,0,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            #nucleotidePlot( pdf, "Read position (second reverse)", "% A/C/G/T reference at pos", names, [ [st.SequenceContentReference(1,1,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=1, pos_normalize=True )
            nucleotidePlot( pdf, "position before(-)/after(+) fragment (forward)", "% A/C/G/T at position", names, [ [st.OutskirtContent(0,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=-20, pos_normalize=True )
            nucleotidePlot( pdf, "position before(-)/after(+) fragment (reverse)", "% A/C/G/T at position", names, [ [st.OutskirtContent(1,n) for st in stats] for n in xrange(4)], zero_padding=False, shift_x=-20, pos_normalize=True )
            plot( pdf, "Read position (first)", "# N", names, [st.SequenceContent(0,4) for st in stats], zero_padding=False, shift_x=1, legend='upper left' )
            plot( pdf, "Read position (second)", "# N", names, [st.SequenceContent(1,4) for st in stats], zero_padding=False, shift_x=1, legend='upper left' )
            
            print "Plotted coverage information: ", clock()
            
            
            plotBaseQualities( pdf, "Read position (first)", "Quality", names, [st.BaseQualityMean(0) for st in stats], [st.BaseQualityMinimum(0) for st in stats], [st.BaseQualityFirstQuartile(0) for st in stats], [st.BaseQualityMedian(0) for st in stats], [st.BaseQualityThirdQuartile(0) for st in stats], [st.BaseQualityMaximum(0) for st in stats] )
            plotBaseQualities( pdf, "Read position (second)", "Quality", names, [st.BaseQualityMean(1) for st in stats], [st.BaseQualityMinimum(1) for st in stats], [st.BaseQualityFirstQuartile(1) for st in stats], [st.BaseQualityMedian(1) for st in stats], [st.BaseQualityThirdQuartile(1) for st in stats], [st.BaseQualityMaximum(1) for st in stats] )
            #plotBaseQualities( pdf, "Read position (first)", "Quality (mapped reads)", names, [st.BaseQualityMeanReference(0) for st in stats], [st.BaseQualityMinimumReference(0) for st in stats], [st.BaseQualityFirstQuartileReference(0) for st in stats], [st.BaseQualityMedianReference(0) for st in stats], [st.BaseQualityThirdQuartileReference(0) for st in stats], [st.BaseQualityMaximumReference(0) for st in stats] )
            #plotBaseQualities( pdf, "Read position (second)", "Quality (mapped reads)", names, [st.BaseQualityMeanReference(1) for st in stats], [st.BaseQualityMinimumReference(1) for st in stats], [st.BaseQualityFirstQuartileReference(1) for st in stats], [st.BaseQualityMedianReference(1) for st in stats], [st.BaseQualityThirdQuartileReference(1) for st in stats], [st.BaseQualityMaximumReference(1) for st in stats] )
            #plot(pdf, "Read position (+ strand)", "Quality", names, [st.BaseQualityMeanPerStrand(0) for st in stats], zero_padding=False, shift_x=1)
            #plot(pdf, "Read position (- strand)", "Quality", names, [st.BaseQualityMeanPerStrand(1) for st in stats], zero_padding=False, shift_x=1)
            
            plot( pdf, "Mean sequence quality (first)", "% reads", names, [st.SequenceQualityMean(0) for st in stats], zero_padding=False, normalize=True, legend='upper left' )
            plot( pdf, "Mean sequence quality (second)", "% reads", names, [st.SequenceQualityMean(1) for st in stats], zero_padding=False, normalize=True, legend='upper left' )
            plotSequenceQualities( pdf, "Sequence quality (first)", "% reads", names, [st.SequenceQualityMean(0) for st in stats], [st.SequenceQualityMinimum(0) for st in stats], [st.SequenceQualityFirstQuartile(0) for st in stats], [st.SequenceQualityMedian(0) for st in stats], [st.SequenceQualityThirdQuartile(0) for st in stats], [st.SequenceQualityMaximum(0) for st in stats] )
            plotSequenceQualities( pdf, "Sequence quality (second)", "% reads", names, [st.SequenceQualityMean(1) for st in stats], [st.SequenceQualityMinimum(1) for st in stats], [st.SequenceQualityFirstQuartile(1) for st in stats], [st.SequenceQualityMedian(1) for st in stats], [st.SequenceQualityThirdQuartile(1) for st in stats], [st.SequenceQualityMaximum(1) for st in stats] )
            plot( pdf, "Sequence error probability mean (first)", "% reads", names, [st.SequenceQualityProbabilityMean(0) for st in stats], zero_padding=False, normalize=True )
            plot( pdf, "Sequence error probability mean (second)", "% reads", names, [st.SequenceQualityProbabilityMean(1) for st in stats], zero_padding=False, normalize=True )           
 
            plot2dQuality( pdf, "Sequence quality (first)", "Base quality", names, [ ( st.BaseQualityForSequenceStart(0), [st.BaseQualityForSequence(0, sq) for sq in xrange( st.BaseQualityForSequenceStart(0), st.BaseQualityForSequenceEnd(0) )] ) for st in stats] )
            plot2dQuality( pdf, "Sequence quality (second)", "Base quality", names, [ ( st.BaseQualityForSequenceStart(1), [st.BaseQualityForSequence(1, sq) for sq in xrange( st.BaseQualityForSequenceStart(1), st.BaseQualityForSequenceEnd(1) )] ) for st in stats] )
            plot2dQuality( pdf, "Preceding quality (first)", "Following quality", names, [ ( st.BaseQualityForPrecedingQualityStart(0), [st.BaseQualityForPrecedingQuality(0, sq) for sq in xrange( st.BaseQualityForPrecedingQualityStart(0), st.BaseQualityForPrecedingQualityEnd(0) )] ) for st in stats] )
            plot2dQuality( pdf, "Preceding quality (second)", "Following quality", names, [ ( st.BaseQualityForPrecedingQualityStart(1), [st.BaseQualityForPrecedingQuality(1, sq) for sq in xrange( st.BaseQualityForPrecedingQualityStart(1), st.BaseQualityForPrecedingQualityEnd(1) )] ) for st in stats] )
            plot2dQuality( pdf, "Sequence quality (first)", "Sequence quality (second)", names, [ ( st.SequenceQualityPairsStart(), [st.SequenceQualityPairs(sq) for sq in xrange( st.SequenceQualityPairsStart(), st.SequenceQualityPairsEnd() )] ) for st in stats] )
            plot(pdf, "Fragment length", "Mean Sequence Quality (first)", names, [st.MeanSequenceQualityMeanByFragmentLength(0) for st in stats], legend='lower right' )
            plot(pdf, "Fragment length", "Mean Sequence Quality (second)", names, [st.MeanSequenceQualityMeanByFragmentLength(1) for st in stats], legend='lower right' )
 
            plot2dQuality( pdf, "Distance to start of error region", "Error rate", names, [ ( st.ErrorRatesByDistanceStart(), [st.ErrorRatesByDistance(dist) for dist in xrange( st.ErrorRatesByDistanceStart(), st.ErrorRatesByDistanceEnd() )] ) for st in stats] )
            plot2dQuality( pdf, "GC", "Error rate", names, [ ( st.ErrorRatesByGCStart(), [st.ErrorRatesByGC(gc) for gc in xrange( st.ErrorRatesByGCStart(), st.ErrorRatesByGCEnd() )] ) for st in stats] )
 
            nucleotidePlot( pdf, "Nucleotide content (first) [%]", "Average sequence quality", names, [ [st.AverageSequenceQualityForBase(0, n) for st in stats] for n in xrange(4)], Nlegend=1, legend='' )
            nucleotidePlot( pdf, "Nucleotide content (second) [%]", "Average sequence quality", names, [ [st.AverageSequenceQualityForBase(1, n) for st in stats] for n in xrange(4)], Nlegend=1, legend='' )
            plot( pdf, "GC content (first)", "Average sequence quality", names, [st.AverageSequenceQualityForGC(0) for st in stats], zero_padding=False)
            plot( pdf, "GC content (second)", "Average sequence quality", names, [st.AverageSequenceQualityForGC(1) for st in stats], zero_padding=False)
            
            nucleotidePlot( pdf, "Base quality (first)", "% of A/C/G/T bases", names, [ [st.NucleotideQuality(0,n) for st in stats] for n in xrange(4)], normalize=True, zero_padding=False )
            nucleotidePlot( pdf, "Base quality (second)", "% of A/C/G/T bases", names, [ [st.NucleotideQuality(1,n) for st in stats] for n in xrange(4)], normalize=True, zero_padding=False )
            
            #plot( pdf, "Base quality (first)", "# N", names, [st.NucleotideQuality(0,4) for st in stats] )
            #plot( pdf, "Base quality (second)", "# N", names, [st.NucleotideQuality(1,4) for st in stats] )
            
            #for qual in [2,4,20,39,40]:
            #    plot( pdf, "Quality {} in read (first)".format(qual), "% reads", names, [st.SequenceQualityContent(0,qual) for st in stats], log=True, normalize=True )
            #    plot( pdf, "Quality {} in read (second)".format(qual), "% reads", names, [st.SequenceQualityContent(1,qual) for st in stats], log=True, normalize=True )
            #    pass
            
            plotTileQuality( pdf, "Read position (first)", "Tile", "Quality mean difference", names, [{ tile : st.TileQualityMeanDifference(0, tile_id) for tile_id, tile in enumerate( st.TileNames() ) } for st in stats ] )
            plotTileQuality( pdf, "Read position (second)", "Tile", "Quality mean difference", names, [{ tile : st.TileQualityMeanDifference(1, tile_id) for tile_id, tile in enumerate( st.TileNames() ) } for st in stats ] )
 
            print "Plotted quality information: ", clock()
            
            plot( pdf, "Error Coverage", "# bases of reference", names, [st.ErrorCoverage() for st in stats], zero_padding=False, log=True )
            plot( pdf, "Error Coverage / Base Coverage [%]", "# bases of reference", names, [st.ErrorCoveragePercent() for st in stats], zero_padding=False, log=True )
            #plot( pdf, "Error Coverage / Base Coverage [%]", "# bases of reference (min cov 10)", names, [st.ErrorCoveragePercentMinCov10() for st in stats], zero_padding=False, log=True )
            #plot( pdf, "Error Coverage / Base Coverage [%]", "# bases of reference (min cov 20)", names, [st.ErrorCoveragePercentMinCov20() for st in stats], zero_padding=False, log=True )
            plot2dQuality( pdf, "Error Coverage + strand", "Error Coverage - strand", names, [ ( 0, [( 0, [st.ErrorCoveragePercentStranded(i, j) for j in xrange( 101 )]) for i in xrange( 101 )] ) for st in stats], axis_lable_dist = 10 )
            #plot2dQuality( pdf, "Error Coverage + strand (min cov 10)", "Error Coverage - strand (min cov 10)", names, [ ( 0, [( 0, [st.ErrorCoveragePercentStrandedMinCov10(i, j) for j in xrange( 101 )]) for i in xrange( 101 )] ) for st in stats], axis_lable_dist = 10 )
            #plot2dQuality( pdf, "Error Coverage + strand (min cov 20)", "Error Coverage - strand (min cov 20)", names, [ ( 0, [( 0, [st.ErrorCoveragePercentStrandedMinCov20(i, j) for j in xrange( 101 )]) for i in xrange( 101 )] ) for st in stats], axis_lable_dist = 10 )
            plot( pdf, "Errors in read (first)", "% reads", names, [st.ErrorsPerRead(0) for st in stats], log=True, normalize=True )
            plot( pdf, "Errors in read (second)", "% reads", names, [st.ErrorsPerRead(1) for st in stats], log=True, normalize=True )
            
            plotCalledBases( pdf, "Base Quality (first)", "% Called", names, [ [ [st.CalledBasesByBaseQuality(0, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.NucleotideQuality(0,n) for n in xrange(5)] for st in stats], False )
            #plotCalledBases( pdf, "Base Quality (first)", "% Called", names, [ [ [st.CalledBasesByBaseQuality(0, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.NucleotideQuality(0,n) for n in xrange(5)] for st in stats] )
            plotCalledBases( pdf, "Base Quality (second)", "% Called", names, [ [ [st.CalledBasesByBaseQuality(1, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.NucleotideQuality(1,n) for n in xrange(5)] for st in stats], False )
            #plotCalledBases( pdf, "Base Quality (second)", "% Called", names, [ [ [st.CalledBasesByBaseQuality(1, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.NucleotideQuality(1,n) for n in xrange(5)] for st in stats] )
            plotCalledBases( pdf, "Position (first)", "% Called", names, [ [ [st.CalledBasesByPosition(0, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.SequenceContent(0,n) for n in xrange(5)] for st in stats], False )
            #plotCalledBases( pdf, "Position (first)", "% Called", names, [ [ [st.CalledBasesByPosition(0, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.SequenceContent(0,n) for n in xrange(5)] for st in stats] )
            plotCalledBases( pdf, "Position (second)", "% Called", names, [ [ [st.CalledBasesByPosition(1, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.SequenceContent(1,n) for n in xrange(5)] for st in stats], False )
            #plotCalledBases( pdf, "Position (second)", "% Called", names, [ [ [st.CalledBasesByPosition(1, i, j) for j in range(5)] for i in range(4)] for st in stats], [ [st.SequenceContent(1,n) for n in xrange(5)] for st in stats] )
            
            #for prev_nuc in range(4):
            #    plotCalledBases( pdf, "Base Quality (first) call {}".format(prev_nuc), "% Called", names, [ [ [st.CalledBasesByBaseQualityPerPreviousCalledBase(0, i, j, prev_nuc) for j in range(5)] for i in range(4)] for st in stats], [ [st.NucleotideQuality(0,n) for n in xrange(5)] for st in stats], False )
            #    pass
            
            nucleotidePlot( pdf, "Homopolymer length", "# A-/C-/G-/T-mers", names, [ [st.HomopolymerDistribution(n) for st in stats] for n in xrange(4)], log=True )
            #plot( pdf, "Homopolymer length", "# N-mers", names, [st.HomopolymerDistribution(4) for st in stats], log=True )
            #for qual in [2,4,20,39,40]:
            #    plot( pdf, "Homoqualimer length (quality={})".format(qual), "# qualimers", names, [st.HomoqualityDistribution(qual) for st in stats], log=True )
            #    pass
            
            plot( pdf, "Length of Insertion", "# insertions", names, [st.InDelErrorByLength(0) for st in stats] )
            plot( pdf, "Read position", "# inserted bases", names, [st.InDelErrorByPosition(0) for st in stats] )
            plot( pdf, "Read GC", "# inserted bases", names, [st.InDelErrorByGC(0) for st in stats] )
            plot( pdf, "Length of Deletion", "# deletions", names, [st.InDelErrorByLength(1) for st in stats] )
            plot( pdf, "Read position", "# deleted bases", names, [st.InDelErrorByPosition(1) for st in stats] )
            plot( pdf, "Read GC", "# deleted bases", names, [st.InDelErrorByGC(1) for st in stats] )
            
            print "Plotted base calling information: ", clock()
            
            #plot( pdf, "Read length (first)", "# reads", names, [st.ReadLengths(0) for st in stats] )
            #plot( pdf, "Read length (second)", "# reads", names, [st.ReadLengths(1) for st in stats] )
            
            #plot( pdf, "Proper pair mapping quality", "# reads", names, [st.ProperPairMappingQuality() for st in stats], legend='upper left' )
            #plot( pdf, "Improper pair mapping quality", "# reads", names, [st.ImproperPairMappingQuality() for st in stats], legend='upper left' )
            #plot( pdf, "Single read mapping quality", "# reads", names, [st.SingleReadMappingQuality() for st in stats], legend='upper left' )
            
            plotTileAbundance( pdf, "Tile", "# read pairs", names, [st.TileAbundance() for st in stats], [st.TileNames() for st in stats] )
            
            plotTileAbundance( pdf, "Adapters (first)", "# adapters", names, [[counts for counts in st.AdapterCount(0) if 0 < counts] for st in stats], [[st.AdapterName(0, i) for i,counts in enumerate(st.AdapterCount(0)) if 0 < counts] for st in stats] )
            plotTileAbundance( pdf, "Adapters (second)", "# adapters", names, [[counts for counts in st.AdapterCount(1) if 0 < counts] for st in stats], [[st.AdapterName(1, i) for i,counts in enumerate(st.AdapterCount(1)) if 0 < counts] for st in stats] )

            print "Plotted other information: ", clock()
            pass
        pass

    pass

def usage():
    print "Usage: python plotDataStats.py [OPTIONS] File [File2 File3 ...]"
    print "Plots the DataStats from an boost archive from readar."
    print "  -h, --help            display this help and exit"
    print "  -o, --output          define plotting output file [File with ending pdf]"
    pass

def main(argv):
    try:
        optlist, args = getopt.getopt(argv, 'ho:', ['help','output='])
        pass
    except getopt.GetoptError:
        print "Unknown option\n"
        usage()
        sys.exit(2)
        pass

    oFile = ''
    for opt, par in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
            pass
        elif opt in ("-o", "--output"):
            oFile = par
            pass
        pass

    if 1 > len(args) or len(args) > 8:
        print "Wrong number of arguments. Only one to eight files are supported.\n"
        usage()
        sys.exit(2)
        pass

    plotDataStats(args,oFile)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])
    pass
