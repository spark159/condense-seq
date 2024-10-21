import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import seaborn as sns
import copy
import math
from scipy import signal
from scipy import stats
from scipy.interpolate import interpn
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import statis_edit as statis

def tuple_cmp(a, b):
    if a[0] <= b[0]:
        return -1
    else:
        return 1

def norm(L):
    total = sum(L)
    return [L[i]/float(total) for i in range(len(L))]


# plot density scatter plot
def density_scatter (X,
                     Y,
                     xlim=[None,None],
                     ylim=[None,None],
                     bins=20,
                     density=False,
                     sort=True,
                     cbar=True,
                     cmap=None,
                     s=3,
                     xlabel=None,
                     ylabel=None,
                     title=None,
                     fig_width=4,
                     fig_height=3,
                     ax=None,
                     save=False,
                     note='',
                     **kwargs):
    
    if ax == None:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    # filter out data nan or outside of limits
    newX, newY = [], []
    for x, y in zip(X, Y):
        if np.isnan(x) or np.isnan(y):
            continue
        if xlim[0] != None and x < xlim[0]:
            continue
        if xlim[1] != None and x > xlim[1]:
            continue
        if ylim[0] != None and y < ylim[0]:
            continue
        if ylim[1] != None and y > ylim[1]:
            continue
        newX.append(x)
        newY.append(y)

    X, Y = newX, newY
    del newX
    del newY

    if cmap == None:
        # "jet-like" colormap with white background
        pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                                       [(0, '#ffffff'),
                                                        (0.03, 'tab:cyan'),
                                                        (0.1, 'tab:blue'),
                                                        (0.3, 'tab:green'),
                                                        (0.5, 'yellow'),
                                                        (0.7, 'tab:orange'),
                                                        (0.9, 'tab:red'),
                                                        (1, 'darkred')],
                                                       N=256)
        cmap = pastel_jet


    # make 2d histogram
    data, X_e, Y_e = np.histogram2d(X,
                                    Y,
                                    bins = bins,
                                    density=density)

    # interpolate the 2d histogram to a continuous density function
    Z = interpn((0.5*(X_e[1:]+X_e[:-1]),
                 0.5*(Y_e[1:]+Y_e[:-1])),
                data,
                np.vstack([X, Y]).T,
                method = "splinef2d",
                bounds_error = False)

    # convert nan to zero
    Z[np.where(np.isnan(Z))] = 0.0

    # sort the points by density, so that the densest points are plotted last
    if sort :
        idx = Z.argsort()
        X, Y, Z = np.asarray(X)[idx], np.asarray(Y)[idx], Z[idx]
        
    img = ax.scatter(X,
                     Y,
                     c=Z,
                     s=s,
                     cmap=cmap,
                     **kwargs)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    
    if cbar:
        cbar = plt.colorbar(img)
        cbar.ax.tick_params(labelsize=5)
    #cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=img)
    #cbar.ax.set_ylabel('Density')

    if make_fig:
        if save:
            plt.savefig("_".join(['DensityScatter', note]) + ".png",
                        format='png',
                        dpi=300,
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax

# pair-wise data correlation plot
def plot_corr_matrix (id_data,
                      id_label=None,
                      ids=None,
                      xlim=[None,None],
                      ylim=[None,None],
                      pair_corr=None,
                      corr='Spearman',
                      fig_scale=1,
                      cell_size=1,
                      label_color='black',
                      text_color='black',
                      scatter_style='dot',
                      ms=1,
                      mfc='k',
                      mec='k',
                      alpha=0.5,
                      bins=20,
                      xscale='linear',
                      yscale='linear',
                      basex=None,
                      basey=None,
                      cbar=True,
                      cmap='Reds',
                      vmin=0.1,
                      vmax=0.9,
                      save=False,
                      title=None,
                      cbar_label=None,
                      note=''):

    if ids == None:
        ids = sorted(id_data.keys())
    if id_label == None:
        id_label = {id:str(id) for id in ids}

    data_num = len(ids)

    # set figure and subplot sizes
    cell_width = cell_size * fig_scale  # subplot axes width in inch 
    cell_height = cell_size * fig_scale # subplot axes height in inch

    left  = 0.1 * fig_scale    # left space of the figure in inch
    right = 0.1 * fig_scale    # right space of the figure in inch
    bottom = 0.1 * fig_scale   # bottom space of the figure in inch
    top = 0.1 * fig_scale      # top space of the figure in inch
    wspace = 0.2 * fig_scale   # width space between subplots in inch
    hspace = 0.2 * fig_scale   # height space between subplots in inch

    if cbar:
        right += wspace + 0.3 * cell_size * fig_scale # give more space for adding color bar

    if title != None:
        top += hspace + 0.3 * cell_size * fig_scale # give more space for adding title

    nrows = data_num # row number
    ncols = data_num # column number

    fig_width = cell_width*ncols + wspace*(ncols-1) + left + right
    fig_height = cell_height*nrows + hspace*(nrows-1) + top + bottom
    
    fig, axes = plt.subplots(figsize=(fig_width, fig_height),
                             nrows=nrows,
                             ncols=ncols)

    # adjust margin spaces
    fig.subplots_adjust(left=left/fig_width,
                        bottom=bottom/fig_height,
                        right=1.0-right/fig_width,
                        top=1.0-top/fig_height,
                        wspace=wspace/cell_width,
                        hspace=hspace/cell_height)

    for i in range(data_num):
        for j in range(data_num):
            id1, id2 = ids[i], ids[j]
            data1, data2 = id_data[id1], id_data[id2]
            label1, label2 = id_label[id1], id_label[id2]

            if i > j:
                if scatter_style == 'dot':
                    axes[i,j].plot(data1,
                                   data2,
                                   'k.',
                                   ms=ms,
                                   mfc=mfc,
                                   mec=mec,
                                   alpha=alpha)

                elif scatter_style == 'histogram':
                    axes[i,j].hist2d(data1,
                                     data2,
                                     range=[xlim,ylim],
                                     bins=bins)

                elif scatter_style == 'density':
                    density_scatter(data1,
                                    data2,
                                    cbar=False,
                                    xlim=xlim,
                                    ylim=ylim,
                                    ax=axes[i,j])

                #wspace = 0.1*(max(X) - min(X))
                #hspace = 0.1*(max(Y) - min(Y))
                #axes[i,j].set_xticks([min(data1), max(data1)])
                #axes[i,j].set_xticklabels([str(round(min(data1),1)),
                #                           str(round(max(data1),1))], rotation=45)
                #axes[i,j].set_yticks([min(data2), max(data2)])
                #axes[i,j].set_yticklabels([str(round(min(data2),1)),
                #                           str(round(max(data2),1))])

                axes[i,j].set_xlim(xlim)
                axes[i,j].set_ylim(ylim)

                axes[i,j].set_xscale(xscale, basex=basex)
                axes[i,j].set_yscale(yscale, basey=basey)

                if j > 0 and i < data_num -1:
                    axes[i,j].tick_params(axis='both',
                                          which='both',
                                          labelbottom=False,
                                          labelleft=False)
                if j == 0 and i < data_num -1:
                    axes[i,j].tick_params(axis='x',
                                          which='both',
                                          labelbottom=False)
                if j > 0 and i == data_num - 1:
                    axes[i,j].tick_params(axis='y',
                                          which='both',
                                          labelleft=False)

            elif i == j:
                assert label1 == label2
                axes[i,j].text(4.5,
                               4.5,
                               label1,
                               color=label_color,
                               ha="center",
                               va="center",
                               fontsize=10,
                               weight='bold')
                axes[i,j].set_xlim([0, 9])
                axes[i,j].set_ylim([0, 9])
                axes[i,j].set_axis_off()
            else:
                assert i < j

                if pair_corr == None:
                    if corr == 'Spearman':
                        value = stats.spearmanr(data1,
                                                data2,
                                                nan_policy='omit')[0]
                    elif corr == 'Pearson':
                        value = stats.pearsonr(data1,
                                               data2,
                                               nan_policy='omit')[0]
                else:
                    value = pair_corr[(id1, id2)]

                matrix = np.zeros((10, 10))
                matrix[:] = value
                img = axes[i,j].imshow(matrix,
                                       cmap=cmap,
                                       vmin=vmin,
                                       vmax=vmax,
                                       origin='lower')
                axes[i,j].text(4.5,
                               4.5,
                               str(round(value,2)),
                               ha="center",
                               va="center",
                               fontsize=10,
                               color=text_color,
                               weight='bold')
                axes[i,j].set_xlim([0, 9])
                axes[i,j].set_ylim([0, 9])
                axes[i,j].tick_params(axis='both',
                                      which='both',
                                      bottom=False,
                                      top=False,
                                      left=False,
                                      right=False,
                                      labelbottom=False,
                                      labeltop=False,
                                      labelleft=False,
                                      labelright=False)
    
    if cbar:
        
        gs = gridspec.GridSpec(nrows=3,
                               ncols=2,
                               width_ratios=[1.0-0.5*right/fig_width,
                                             0.5*right/fig_width],
                               height_ratios=[0.1, 2, 1],
                               left=left/fig_width,
                               bottom=bottom/fig_height,
                               right=1.0-left/fig_width,
                               top=1.0-top/fig_height,
                               wspace=0,
                               hspace=0)

        cax = fig.add_subplot(gs[1, 1])
        cbar=fig.colorbar(img, cax=cax)

        if cbar_label == None:
            cbar_label = corr + ' correlation'
        
        cax.set_ylabel(cbar_label,
                       rotation=-90,
                       va="bottom",
                       fontsize=10)

    if title !=None:

        gs = gridspec.GridSpec(nrows=2,
                               ncols=1,
                               height_ratios=[top/fig_height, 1.0-top/fig_height],
                               left=0,
                               bottom=0,
                               right=1,
                               top=1,
                               wspace=0,
                               hspace=0)

        tax = fig.add_subplot(gs[0, 0])

        tax.text(0,
                 0,
                 title,
                 ha="center",
                 va="center",
                 fontsize=20)
        
        tax.set_xlim([-10, 10])
        tax.set_ylim([-5, 5])
        tax.set_axis_off()
    
    if save:
        plt.savefig("Corr_matrix_%s.png" % (note), dpi=1000, bbox_inches='tight')
    plt.show()
    plt.close()
    return

# plot chromosome ideogram
def plot_ideogram (Gtype_ideogram,
                   Gtick_locs=[],
                   Gtick_labels=[],
                   fig_width=15,
                   fig_height=2,
                   aspect='auto',
                   save=False,
                   ax=None,
                   note=''):

    ideograms = Gtype_ideogram.values()

    if len(ideograms) <=0:
        return
    
    xaxis_len = None
    for i in range(len(ideograms)):
        ideogram = ideograms[i]
        if i == 0:
            xaxis_len = len(ideogram)
            continue
        assert len(ideogram) == xaxis_len

    if ax == None:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    ax.imshow(np.transpose(Gtype_ideogram['num']),
              cmap='Greys',
              aspect=aspect)

    ax.imshow(np.transpose(Gtype_ideogram['acen']),
              cmap ='Reds',
              vmin=0,
              vmax=20,
              aspect=aspect)

    ax.imshow(np.transpose(Gtype_ideogram['var']),
              cmap ='Purples',
              vmin=0,
              vmax=20,
              aspect=aspect)

    ax.set_yticks([])
    ax.set_xticks(Gtick_locs)
    ax.set_xticklabels(Gtick_labels)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.set_xlim([0, xaxis_len+1])

    if make_fig:
        if save:
            plt.savefig("_".join(['Ideogram', note]) + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot data in genome-wide with chromosome ideogram
def plot_genome_wide (side_names,
                      name_sig,
                      name_color=None,
                      name_alpha=None,
                      name_lw=None,
                      name_linestyle=None,
                      name_label=None,
                      side_ylabel=None,
                      side_ycolor=None,
                      side_ylim=None,
                      side_yscale=None,
                      side_basey=None,
                      xtick_locs=[],
                      xtick_labels=[],
                      Gtype_ideogram=None,
                      Gtick_locs=[],
                      Gtick_labels=[],
                      fig_width=15,
                      fig_height=4,
                      spine_option=None,
                      height_ratios=[8,1],
                      hspace=0.65,
                      axes=None,
                      xlabel="Position (Mb)",
                      legend_loc='best',
                      shade_wins=None,
                      save=False,
                      note=''):

    names = []
    for side in side_names:
        names += side_names[side]

    if len(names) <=0:
        return

    xaxis_len = None
    for i in range(len(names)):
        name = names[i]
        sig = name_sig[name]
        if i == 0:
            xaxis_len = len(sig)
            continue
        assert len(sig) == xaxis_len

    if Gtype_ideogram:
        for ideogram in Gtype_ideogram.values():
            assert len(ideogram) == xaxis_len
        
    if axes != None:
        if Gtype_ideogram:
            assert len(axes) == 2
        else:
            assert len(axes) == 1

        make_fig = False
        
    else:
        if Gtype_ideogram:
            fig, axes = plt.subplots(nrows=2,
                                     ncols=1,
                                     figsize=(fig_width, fig_height),
                                     gridspec_kw = {'hspace':hspace,
                                                    'height_ratios':height_ratios})
        else:
            fig, ax = plt.subplots(nrows=1,
                                   ncols=1,
                                   figsize=(fig_width, fig_height))
            axes = [ax]
            
        make_fig = True

    side_ax = {'left':axes[0], 'right':axes[0].twinx()}

    for side in side_names:
        ax = side_ax[side]
        names = side_names[side]

        for name in names:
            try:
                color = name_color[name]
            except:
                color = None
            try:
                alpha = name_alpha[name]
            except:
                alpha = None
            try:
                lw = name_lw[name]
            except:
                lw = None
            try:
                linestyle = name_linestyle[name]
            except:
                linestyle = None
            try:
                label = name_label[name]
            except:
                label = None

            ax.plot(name_sig[name],
                    color=color,
                    alpha=alpha,
                    lw=lw,
                    linestyle=linestyle,
                    label=label)

        try:
            ylabel = side_ylabel[side]
        except:
            ylabel = None
        try:
            ycolor = side_ycolor[side]
        except:
            ycolor = 'k'
        try:
            ylim = side_ylim[side]
        except:
            ylim = [None, None]
        try:
            yscale = side_yscale[side]
        except:
            yscale = 'linear'
        try:
            basey = side_basey[side]
        except:
            basey = None

        ax.set_ylabel(ylabel, color=ycolor, fontsize=12)
        ax.tick_params('y', colors=ycolor)
        ax.set_ylim(ylim)
        ax.set_yscale(yscale, basey=basey)

        for spine in ['top', 'bottom', 'left', 'right']:
            try:
                option = spine_option[spine]
            except:
                continue
            ax.spines[spine].set_visible(option)

    if shade_wins:
        for win in shade_wins:
            st, ed = win
            axes[0].axvspan(st, ed-1, alpha=0.15, facecolor='grey')

    axes[0].set_xlim([0, xaxis_len+1])
    axes[0].set_xticks(xtick_locs)
    axes[0].set_xticklabels(xtick_labels, fontsize=10)
    axes[0].set_xlabel(xlabel, fontsize=15)

    if name_label:
        leg = axes[0].legend(framealpha=1, loc=legend_loc)
        for legobj in leg.legendHandles:
            legobj.set_alpha(1.0)
            legobj.set_linewidth(2.0)

    axes[0].set_zorder(1)  # default zorder is 0 for ax1 and ax2
    axes[0].patch.set_visible(False)  # prevents ax1 from hiding ax2

    if Gtype_ideogram:        
        plot_ideogram (Gtype_ideogram,
                       Gtick_locs,
                       Gtick_labels,
                       ax=axes[1])

    if make_fig:
        if save:
            plt.savefig("_".join(['Gwide', note]) + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    if Gtype_ideogram:
        return side_ax, axes[1]
    
    return side_ax

# plot data in genome-wide with chromosome ideogram (mulitple)
def plot_genome_wide_multiple (side_names_list,
                               name_sig,
                               name_color=None,
                               name_alpha=None,
                               name_lw=None,
                               name_linestyle=None,
                               name_label=None,
                               side_ylabel_list=None,
                               side_ycolor_list=None,
                               side_ylim_list=None,
                               side_yscale_list=None,
                               side_basey_list=None,
                               xtick_locs=[],
                               xtick_labels=[],
                               Gtype_ideogram=None,
                               Gtick_locs=None,
                               Gtick_labels=None,
                               fig_width=None,
                               fig_height=None,
                               spine_option={},
                               xlabel="Position (Mb)",
                               legend_loc='best',
                               shade_wins=None,
                               save=False,
                               note=''):

    names = []
    for side_names in side_names_list:
        for side in side_names:
            names += side_names[side]

    if len(names) <=0:
        return

    xaxis_len = None
    for i in range(len(names)):
        name = names[i]
        sig = name_sig[name]
        if i == 0:
            xaxis_len = len(sig)
            continue
        assert len(sig) == xaxis_len

    if Gtype_ideogram:
        for ideogram in Gtype_ideogram.values():
            assert len(ideogram) == xaxis_len

    if fig_width == None:
        fig_width = round(6.4*float(xaxis_len)/10000, 1)

    if fig_height == None:
        fig_height = 1.4*len(side_names_list)
    
    nrows = len(side_names_list)
    ncols = 1
    
    height_ratios = [7]*len(side_names_list)

    if Gtype_ideogram:
        fig_height +=1
        nrows +=1
        height_ratios += [1]

    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             figsize=(fig_width, fig_height),
                             sharex=True,
                             gridspec_kw = {'hspace':0.12,
                                            'wspace':0.12,
                                            'height_ratios':height_ratios})

    for k in range(len(side_names_list)):
        side_names = side_names_list[k]

        try:
            side_ylabel = side_ylabel_list[k]
        except:
            side_ylabel = None
        try:
            side_ycolor = side_ycolor_list[k]
        except:
            side_ycolor = None
        try:
            side_ylim = side_ylim_list[k]
        except:
            side_ylim = None
        try:
            side_yscale = side_yscale_list[k]
        except:
            side_yscale = None
        try:
            side_basey = side_basey_list[k]
        except:
            side_basey = None

        side_ax = plot_genome_wide (side_names,
                                    name_sig,
                                    name_color=name_color,
                                    name_alpha=name_alpha,
                                    name_lw=name_lw,
                                    name_linestyle=name_linestyle,
                                    name_label=name_label,
                                    side_ylabel=side_ylabel,
                                    side_ycolor=side_ycolor,
                                    side_ylim=side_ylim,
                                    side_yscale=side_yscale,
                                    side_basey=side_basey,
                                    xtick_locs=xtick_locs,
                                    xtick_labels=xtick_labels,
                                    xlabel=xlabel,
                                    legend_loc=legend_loc,
                                    shade_wins=shade_wins,
                                    spine_option=spine_option,
                                    axes=[axes[k]])

        axes[k].tick_params(top='off',
                            bottom='off',
                            left='on',
                            right='on',
                            labelbottom='off')

    if Gtype_ideogram:
        plot_ideogram (Gtype_ideogram,
                       Gtick_locs,
                       Gtick_labels,
                       ax=axes[-1])

        axes[-1].set_xticks([])
        axes[-1].set_yticks([])

    if save:
        plt.savefig("_".join(['Gwide_multi', note]) + ".svg",
                    format='svg',
                    bbox_inches='tight')
    else:
        plt.tight_layout()
        plt.show()    
    plt.close()

    return axes

# plot boxplot
def plot_boxplot (key_values,
                  key_name=None,
                  keys=None,
                  ylabel='',
                  title=None,
                  rotation=None,
                  ylim=[None, None],
                  color='white',
                  ycolor='black',
                  save=False,
                  fig_width=8,
                  fig_height=6,
                  ax=None,
                  note=""):

    if keys:
        new_keys = []
        for key in keys:
            if key in key_values:
                new_keys.append(key)
        keys = new_keys
    else:
        keys = key_values.keys()

    if not ax:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    pos_list = [i for i in range(len(keys))]

    bp = ax.boxplot([key_values[key] for key in keys],
                    0, "",
                    positions=pos_list,
                    widths=0.5,
                    patch_artist=True,
                    boxprops=dict(facecolor=color))

    ax.set_ylabel(ylabel, color=ycolor)
    ax.tick_params('y', colors=ycolor)

    if title:
        ax.set_title(title)

    ax.set_xticks(range(len(keys)))

    if key_name:
        xticklabels = [key_name[key] for key in keys]
    else:
        xticklabels = keys
        
    ax.set_xticklabels(xticklabels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=rotation)

    ax.set_xlim([-0.5, len(keys)-0.5])
    ax.set_ylim(ylim)
    
    if make_fig:
        if save:
            plt.savefig("boxplot_" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax

# plot boxplot for two groups
def plot_boxplot_pair (key_values1,
                       key_values2,
                       key_name=None,
                       keys=None,
                       color1='pink',
                       color2='lightblue',
                       ylim1=[None, None],
                       ylim2=[None, None],
                       ycolor1='red',
                       ycolor2='blue',
                       ylabel1='',
                       ylabel2='',
                       title=None,
                       label1='',
                       label2='',
                       legend_loc='upper right',
                       rotation=None,
                       sharey = False,
                       fig_width=8,
                       fig_height=6,
                       save=False,
                       axes=None,
                       note=""):

    if keys:
        new_keys = []
        for key in keys:
            if key in key_values1 and key in key_values2:
                new_keys.append(key)
        keys = new_keys
    else:
        keys = set(key_values1.keys()) & set(key_values2.keys())


    if not axes:
        fig, ax1 = plt.subplots(figsize=(fig_width, fig_height))

        if sharey:
            ax2 = ax1
            assert ylabel1 == ylabel2
            assert ycolor1 == ycolor2
            assert ylim1 == ylim2
        else:
            ax2 = ax1.twinx()

        make_fig = True

    else:
        assert len(axes) == 2
        ax1, ax2 = axes
        make_fig = False

    pos_list1 = [i - 0.2 for i in range(len(keys))]
    bp1 = ax1.boxplot([key_values1[key] for key in keys],
                      0, "",
                      positions=pos_list1,
                      widths=0.3,
                      patch_artist=True,
                      boxprops=dict(facecolor=color1))

    pos_list2 = [i + 0.2 for i in range(len(keys))]
    bp2 = ax2.boxplot([key_values2[key] for key in keys],
                      0, "",
                      positions=pos_list2,
                      widths=0.3,
                      patch_artist=True,
                      boxprops=dict(facecolor=color2))

    
    ax1.set_xticks(range(len(keys)))

    if key_name:
        xticklabels = [key_name[key] for key in keys]
    else:
        xticklabels = keys

    ax1.set_xticklabels(xticklabels,
                        ha="right",
                        rotation_mode="anchor",
                        rotation=rotation)
    
    ax1.set_xlim([-0.5, len(keys)-0.5])

    ax1.set_ylabel(ylabel1, color=ycolor1)
    ax1.tick_params('y', colors=ycolor1)
    ax1.set_ylim(ylim1)

    if not sharey:
        ax2.set_ylabel(ylabel2, color=ycolor2)
        ax2.tick_params('y', colors=ycolor2)
        ax2.set_ylim(ylim2)

    if title:
        plt.title(title)

    ax1.legend([bp1["boxes"][0],
                bp2["boxes"][0]],
               [label1, label2],
               loc='upper right')

    if make_fig:
        if save:
            plt.savefig("boxplot_pair" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    if sharey:
        return [ax1]
    else:
        return [ax1, ax2]

# plot boxplot of multiple groups
def plot_boxplot_multiple (key_values_list,
                           key_name=None,
                           axhline=None,
                           keys=None,
                           ylabel='',
                           title=None,
                           rotation=None,
                           ylim=[None, None],
                           colors=[],
                           ycolor='black',
                           labels=[],
                           legend_loc='best',
                           ax=None,
                           save=False,
                           note=""):

    common_keys = set([])
    for i in range(len(key_values_list)):
        key_values = key_values_list[i]
        if i == 0:
            common_keys |= set(key_values.keys())
            continue
        common_keys &= set(key_values.keys())
    common_keys = list(common_keys)

    if keys:
        new_keys = []
        for key in keys:
            if key in common_keys:
                new_keys.append(key)
        keys = new_keys
    else:
        keys = common_keys

    if not labels:
        labels = [None]*len(key_values_list)

    if not colors:
        colors = ['white']*len(key_values_list)

    offset = 0.6/len(key_values_list)

    if not ax:
        fig, ax = plt.subplots()
        make_fig = True
    else:
        make_fig = False

    if axhline != None:
        ax.axhline(y=axhline, linestyle='--', color='k', alpha=0.5)

    bp_list = []
    for i in range(len(key_values_list)):
        key_values = key_values_list[i]
        pos_list = [k + offset*i for k in range(len(keys))]
        bp = ax.boxplot([key_values[key] for key in keys],
                         positions=pos_list,
                         showfliers=False,
                         notch=True,
                         widths=offset*0.7,
                         patch_artist=True,
                         boxprops=dict(facecolor=colors[i]))
        bp_list.append(bp)

    for bp in bp_list:
        for median in bp['medians']:
            median.set_color('red')


    xtick_locs = [k + 0.5*offset*(len(key_values_list)-1) for k in range(len(keys))]
    ax.set_xticks(xtick_locs)

    if key_name:
        xticklabels = [key_name[key] for key in keys]
    else:
        xticklabels = keys

    ax.set_xticklabels(xticklabels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=rotation)
        
    ax.set_xticklabels(xticklabels, rotation=rotation)

    ax.set_xlim([xtick_locs[0]-offset*3, xtick_locs[-1]+offset*3])
    ax.set_ylim(ylim)

    ax.legend([bp["boxes"][0] for bp in bp_list],
              labels,
              loc=legend_loc)

    ax.set_ylabel(ylabel, color=ycolor)
    ax.tick_params('y', colors=ycolor)

    if title:
        ax.set_title(title)

    if make_fig:
        if save:
            plt.savefig("boxplot_multi" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax

# plot metagene profile
def plot_profile (profile,
                  offset=0,
                  pad_len=0,
                  color='black',
                  lw=3,
                  alpha=1,
                  xtick_locs=[],
                  xtick_labels=[],
                  rotation=45,
                  xlabel='',
                  ylabel='',
                  ylim=[None,None],
                  title='',
                  ax=None,
                  fig_width=5,
                  fig_height=3.5,
                  save=False,
                  note=''):

    if not ax:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    profile[:pad_len] = [np.NaN]*pad_len
    profile[len(profile)-pad_len:] = [np.NaN]*pad_len

    X = [k + offset for k in range(len(profile))]
    ax.plot(X, profile, lw=lw, alpha=alpha, color=color)

    ax.set_xticks(xtick_locs)
    ax.set_xticklabels(xtick_labels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=rotation)

    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)

    if title:
        ax.set_title(title, fontsize=20)

    ax.set_ylim(ylim)

    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)

    if make_fig:
        if save:
            plt.savefig("profile" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot metagene profile for pair
def plot_profile_pair (profile1,
                       profile2,
                       offset=0,
                       pad_len=0,
                       color1='black',
                       color2='black',
                       lw1=3,
                       lw2=3,
                       alpha1=1,
                       alpha2=1,
                       xtick_locs=[],
                       xtick_labels=[],
                       rotation=45,
                       xlabel='',
                       sharey=False,
                       ylabel1='',
                       ylabel2='',
                       ycolor1='black',
                       ycolor2='black',
                       ylim1=[None,None],
                       ylim2=[None,None],
                       title='',
                       axes=None,
                       fig_width=5,
                       fig_height=3.5,
                       save=False,
                       note=''):

    assert len(profile1) == len(profile2)

    if not axes:
        fig, ax1 = plt.subplots(figsize=(fig_width, fig_height))

        if sharey:
            ax2 = ax1
            assert ylabel1 == ylabel2
            assert ycolor1 == ycolor2
            assert ylim1 == ylim2
        else:
            ax2 = ax1.twinx()

        make_fig = True

    else:
        assert len(axes) == 2
        ax1, ax2 = axes
        make_fig = False

    profile1[:pad_len] = [np.NaN]*pad_len
    profile1[len(profile1)-pad_len:] = [np.NaN]*pad_len

    profile2[:pad_len] = [np.NaN]*pad_len
    profile2[len(profile2)-pad_len:] = [np.NaN]*pad_len

    X = [k + offset for k in range(len(profile1))]
    ax1.plot(X, profile1, lw=lw1, alpha=alpha1, color=color1)
    ax2.plot(X, profile2, lw=lw2, alpha=alpha2, color=color2)

    ax1.set_xticks(xtick_locs)
    ax1.set_xticklabels(xtick_labels,
                        ha="right",
                        rotation_mode="anchor",
                        rotation=rotation)
    ax1.set_xlabel(xlabel, fontsize=15)

    ax1.set_ylabel(ylabel1, color=ycolor1, fontsize=15)
    ax1.tick_params('y', colors=ycolor1)
    ax1.set_ylim(ylim1)

    if not sharey:
        ax2.set_ylabel(ylabel2, color=ycolor2, fontsize=15)
        ax2.tick_params('y', colors=ycolor2)
        ax2.set_ylim(ylim2)

    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.tick_params(axis='both', which='minor', labelsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='minor', labelsize=15)

    if title:
        ax.set_title(title, fontsize=20)

    if make_fig:
        if save:
            plt.savefig("profile_pair" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    if sharey:
        return [ax1]
    else:
        return [ax1, ax2]


# plot metagene profiles for multiple input
def plot_profile_multiple (profiles,
                           offset=0,
                           pad_len=0,
                           colors=[],
                           lws=[],
                           alphas=[],
                           xtick_locs=[],
                           xtick_labels=[],
                           labels=[],
                           legend_loc='best',
                           rotation=45,
                           xlabel='',
                           ylabel='',
                           ylim=[None,None],
                           title='',
                           ax=None,
                           fig_width=5,
                           fig_height=3.5,
                           save=False,
                           note=''):

    profile_len = len(profiles[0])
    for profile in profiles:
        assert len(profile) == profile_len

    if not ax:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    X = [k + offset for k in range(profile_len)]

    for i in range(len(profiles)):
        profile = profiles[i]
        profile[:pad_len] = [np.NaN]*pad_len
        profile[profile_len - pad_len:] = [np.NaN]*pad_len

        color = colors[i]
        lw = lws[i]
        alpha = alphas[i]
        label = labels[i]

        ax.plot(X,
                profile,
                lw=lw,
                alpha=alpha,
                color=color,
                label=label)

    ax.set_ylim(ylim)

    ax.set_xticks(xtick_locs)
    ax.set_xticklabels(xtick_labels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=rotation)

    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)

    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)

    ax.legend(loc=legend_loc, fontsize=14, frameon=False)

    if title:
        ax.set_title(title, fontsize=20)

    if make_fig:
        if save:
            plt.savefig("profile_multiple" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot heatmap for metagene profiles
def plot_profile_heatmap (profiles,
                          offset=0,
                          pad_len=0,
                          cmap='jet',
                          vmin=None,
                          vmax=None,
                          log_scale=None,
                          xtick_locs=[],
                          xtick_labels=[],
                          rotation=45,
                          xlabel='',
                          ax=None,
                          fig_width=3,
                          fig_height=5,
                          aspect='auto',
                          save=False,
                          note=''):

    profile_len = len(profiles[0])
    for profile in profiles:
        assert len(profile) == profile_len

    if not ax:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    X = [k + offset for k in range(profile_len)]

    img = []
    for i in range(len(profiles)):
        profile = profiles[i]
        left_pad = profile[pad_len]
        profile[:pad_len] = [left_pad]*pad_len
        right_pad = profile[profile_len - pad_len-1]
        profile[profile_len - pad_len:] = [right_pad]*pad_len
        img.append(profile)
    img = np.asarray(img)

    if log_scale:
        minimum = np.min(img)
        img = img - minimum + 1
        if log_scale == 'log2':
            img = np.log2(img)
        elif log_scale == 'log10':
            img = np.log10(img)
    
    ax.imshow(img,
              aspect=aspect,
              cmap=cmap,
              vmin=vmin,
              vmax=vmax)
    
    ax.set_xticks([loc - offset for loc in xtick_locs])
    ax.set_xticklabels(xtick_labels,
                       fontsize=10,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=45)
    ax.set_yticks([])
    ax.set_yticklabels([])

    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=10)

    ax.set_xlabel(xlabel, fontsize=15)

    if make_fig:
        if save:
            plt.savefig("profile_heatmap" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot dendrogram
def plot_dendrogram(Z,
                    idx_name,
                    node_color=None,
                    name_color=None,
                    fig_width=8,
                    fig_height=8,
                    ax=None,
                    save=False,
                    note=''):
    
    if not ax:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False
        
    if node_color != None:
        dendrogram(Z,
                   link_color_func=lambda k: node_color[k],
                   orientation='right',
                   ax=ax)
    else:
        dendrogram(Z,
                   orientation='right',
                   ax=ax)

    new_labels = []
    old_labels = ax.get_yticklabels()
    for label in old_labels:
        idx = int(label.get_text())
        name = idx_name[idx]
        label.set_text(name)
        if name_color:
            label.set_color(name_color[name])
        new_labels.append(label)
    ax.set_yticklabels(new_labels, weight='bold')

    if make_fig:
        if save:
            plt.savefig("dendrogram_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot GSEA output
def plot_GSEA (gene_value,
               gs_list,
               gsname_label=None,
               gsname_color=None,
               fig_width=3,
               fig_height=5,
               fontsize=10,
               save=False,
               note=''):

    min_value = min(gene_value.values())
    max_value = max(gene_value.values())
    lim_value = min(abs(min_value), abs(max_value))
    nes_list = [gs['nes'] for gs in gs_list]
    min_nes = min(nes_list)
    max_nes = max(nes_list)

    cmap = mpl.cm.get_cmap("binary")
    bg_cmap = mpl.cm.get_cmap("Spectral_r")

    nrows = len(gs_list)+1
    ncols = 2
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             figsize=(fig_width, fig_height),
                             sharey=True,
                             gridspec_kw={'width_ratios': [0.6, 0.4],
                                          'wspace':0.03})

    for i in range(nrows):
        for j in range(ncols):

            if i == 0 and j == 0:
                # plot ranked values
                axes[i][j].bar(range(len(gene_value)),
                               sorted(gene_value.values(), reverse=True),
                               color='k')
                axes[i][j].set_xlim([0, len(gene_value)-1])

            if i > 0:
                gs = gs_list[i-1]
                nes = gs['nes']
                gsname = gs['name']

                try:
                    label = gsname_label[gsname]
                except:
                    label = gsname
                
                try:
                    color = gsname_color[gsname]
                except:
                    color = 'k'

                if j == 1:
                    # put gene-set name
                    axes[i][j].annotate(label,
                                        (0,0),
                                        ha='left',
                                        va='center',
                                        weight='bold',
                                        color=color,
                                        fontsize=fontsize)
                    axes[i][j].set_xlim([0, max_nes])
                    axes[i][j].set_ylim([-1,1])

                elif j == 0:
                    # plot background color gradient
                    binnum = 10
                    binsize = len(gene_value)/10
                    for k in range(binnum):
                        color = float(k+0.5)/binnum
                        axes[i][j].axvspan(k*binsize,
                                           (k+1)*binsize-1,
                                           color=bg_cmap(color),
                                           alpha=0.3,
                                           lw=0,
                                           zorder=1)
                    axes[i][j].axhspan(0, 1, color='white', alpha=1, lw=0, zorder=1)

                    # plot gene set ranks
                    gene_info = gs['genes']
                    max_es = max([abs(gene_info[gene]['es']) for gene in gene_info])
                    for gene in gene_info:
                        rank = gene_info[gene]['rank']
                        pos = rank - 1
                        #value = gene_value[gene]
                        es = gene_info[gene]['es']

                        color = statis.rescale(rank,
                                               old_st=1,
                                               old_ed=len(gene_value),
                                               new_st=0,
                                               new_ed=1) #color by rank
                    
                        alpha = statis.rescale(abs(es),
                                               old_st=0,
                                               old_ed=max_es,
                                               new_st=0.05,
                                               new_ed=1) # emphasis by running enrichment

                        axes[i][j].axvline(x=pos, color=cmap(alpha), linewidth=0.25, alpha=1)
                        axes[i][j].set_xlim([0, len(gene_value)-1])
                        axes[i][j].set_ylim([-1, 1])

            # remove ticks
            axes[i][j].tick_params(top='off',
                                   bottom='off',
                                   left='off',
                                   right='off',
                                   labelleft='off',
                                   labelright='off',
                                   labeltop='off',
                                   labelbottom='off')

            if i <= 0 or j != 0:
                # remove frames
                axes[i][j].spines['top'].set_visible(False)
                axes[i][j].spines['bottom'].set_visible(False)
                axes[i][j].spines['left'].set_visible(False)
                axes[i][j].spines['right'].set_visible(False)

    if save:
        plt.savefig("GSEA_%s.png" % (note),
                    dpi=500,
                    bbox_inches='tight')
    else:
        plt.tight_layout()
        plt.show()    
    plt.close()

    return axes


def plot_partition (ID_value,
                    p_wins,
                    hist_bins=1000,
                    label_fontsize=12,
                    xlabel='',
                    ylabel='',
                    title='',
                    xlim=None,
                    ylim=None,
                    fig_width=2.4,
                    fig_height=2,
                    ax=None,
                    save=False,
                    note=''):

    if not ax:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    ax.hist(ID_value.values(), bins=hist_bins)

    if xlim == None:
        xlim = ax.get_xlim()
    if ylim == None:
        ylim = ax.get_ylim()

    p_num = len(p_wins)

    p_lines = set([])
    for p_win in p_wins:
        p_lines |= set(p_win)
    
    for p_line in p_lines:
        plt.axvline(x=p_line, color='k', linestyle='--')

    for i in range(p_num):
        st, ed = p_wins[i]
        if i == 0:
            x = np.mean([max([xlim[0], st]), ed])
        elif i == p_num - 1:
            x = np.mean([st, min([xlim[1], ed])])
        else:
            x = np.mean([st, ed])
        ax.text(x,
                np.mean(ylim),
                statis.print_Roman(i+1),
                fontsize=label_fontsize,
                va='center',
                ha='center')
        
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.tick_params(axis='both', which='minor', labelsize=5)

    if make_fig:
        if save:
            plt.savefig("partition_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


def plot_partition_heatmap (img,
                            xtick_labels=None,
                            ytick_labels=None,
                            fig_width=2.5,
                            fig_height=8,
                            vmin=None,
                            vmax=None,
                            cmap='coolwarm',
                            ax=None,
                            save=False,
                            note=''):

    if not ax:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    im = ax.imshow(img,
                   vmin=vmin,
                   vmax=vmax,
                   cmap=cmap,
                   aspect='auto')

    if xtick_labels == None:
        xtick_labels = range(1, len(img[0])+1)
    if ytick_labels == None:
        ytick_labels = range(1, len(img)+1)
    
    ax.set_xticks(range(len(xtick_labels)))
    ax.set_xticklabels(xtick_labels,
                       fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(ytick_labels)))    
    ax.set_yticklabels(ytick_labels,
                       horizontalalignment='right',
                       fontname='monospace',
                       fontsize=9)

    if make_fig:
        if save:
            plt.savefig("partition_heatmap_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


def plot_partition_heatmap_multiple (img_list,
                                     xtick_labels_list=None,
                                     ytick_labels_list=None,
                                     fig_width=2.5,
                                     fig_height=8,
                                     vmin=None,
                                     vmax=None,
                                     cmap='coolwarm',
                                     axes=None,
                                     save=False,
                                     note=''):

    if not axes:
        fig, axes = plt.subplots(nrows=1,
                                 ncols=len(img_list),
                                 figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    for i in range(len(img_list)):
        img = img_list[i]
        
        try:
            xtick_labels = xtick_labels_list[i]
        except:
            xtick_labels = None
        try:
            ytick_labels = ytick_labels_list[i]
        except:
            ytick_labels = None
        
        plot_partition_heatmap (img,
                                xtick_labels=xtick_labels,
                                ytick_labels=ytick_labels,
                                vmin=vmin,
                                vmax=vmax,
                                cmap=cmap,
                                ax=axes[i])

    return axes

def plot_ATGC_periodicity (AT_sig,
                           GC_sig,
                           AT_color='r',
                           GC_color='b',
                           AT_ylabel='AA/AT/TA/TT freqeuncy',
                           GC_ylabel='CC/CG/GC/GG freqeuncy',
                           xlabel="Super Helical Location",
                           NCPlen=147,
                           SHL_lines = range(-6, 7, 2),
                           SHL_xticks = range(-7, 8, 2),
                           fig_width=2.6,
                           fig_height=2,
                           axes=None,
                           save=False,
                           note=''):

    if not axes:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        ax1 = ax
        ax2 = ax1.twinx()
        make_fig = True
    else:
        assert len(axes) == 2
        ax1, ax2 = axes
        make_fig = False

    ax1.plot(AT_sig, color=AT_color, lw=2)
    ax2.plot(GC_sig, color=GC_color, lw=2)

    lines = [NCPlen/2 + i*10 for i in SHL_lines]
    for line in lines:
        ax1.axvline(x=line, color='k', linestyle='--', alpha=0.25)
    
    ax1.set_ylabel(AT_ylabel, color='r', fontsize=8)
    ax2.set_ylabel(GC_ylabel, color='b', fontsize=8)
    ax1.tick_params('y', colors='r', labelsize=8)
    ax2.tick_params('y', colors='b', labelsize=8)
    ax1.set_xticks([NCPlen/2 + 10*i for i in SHL_xticks])
    ax1.set_xticklabels([str(10*i) for i in SHL_xticks], fontsize=5)
    ax1.set_xlabel(xlabel, fontsize=8)

    if make_fig:
        if save:
            plt.savefig("ATGCperiod_" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return [ax1, ax2]


def plot_ATGC_periodicity_multiple (AT_sigs,
                                    GC_sigs,
                                    AT_colors=[],
                                    GC_colors=[],
                                    AT_ylabel='AA/AT/TA/TT freqeuncy',
                                    GC_ylabel='CC/CG/GC/GG freqeuncy',
                                    xlabel="Super Helical Location",
                                    NCPlen=147,
                                    SHL_lines = range(-6, 7, 2),
                                    SHL_xticks = range(-7, 8, 2),
                                    fig_width=2.6,
                                    fig_height=2,
                                    AT_labels=[],
                                    GC_labels=[],
                                    AT_legend_loc='upper left',
                                    GC_legend_loc='lower right',
                                    axes=None,
                                    save=False,
                                    note=''):

    if not axes:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        ax1 = ax
        ax2 = ax1.twinx()
        make_fig = True
    else:
        assert len(axes) == 2
        ax1, ax2 = axes
        make_fig = False

    assert len(AT_sigs) == len(GC_sigs)
        
    if not AT_colors:
        AT_colors = ['r'] * len(AT_sigs)
    if not GC_colors:
        GC_colors = ['b'] * len(GC_sigs)

    if not AT_labels:
        AT_labels = [None] * len(AT_sigs)
        AT_legend = False
    else:
        AT_legend = True
    if not GC_labels:
        GC_labels = [None] * len(GC_sigs)
        GC_legend = False
    else:
        GC_legend = True
        
    for AT_sig, AT_color, AT_label in zip(AT_sigs, AT_colors, AT_labels):
        ax1.plot(AT_sig, color=AT_color, lw=2, label=AT_label)

    for GC_sig, GC_color, GC_label in zip(GC_sigs, GC_colors, GC_labels):
        ax2.plot(GC_sig, color=GC_color, lw=2, label=GC_label)

    lines = [NCPlen/2 + i*10 for i in SHL_lines]
    for line in lines:
        ax1.axvline(x=line, color='k', linestyle='--', alpha=0.25)
    
    ax1.set_ylabel(AT_ylabel, color='r', fontsize=8)
    ax2.set_ylabel(GC_ylabel, color='b', fontsize=8)
    ax1.tick_params('y', colors='r', labelsize=8)
    ax2.tick_params('y', colors='b', labelsize=8)
    ax1.set_xticks([NCPlen/2 + 10*i for i in SHL_xticks])
    ax1.set_xticklabels([str(10*i) for i in SHL_xticks], fontsize=5)
    ax1.set_xlabel(xlabel, fontsize=8)

    if AT_legend:
        ax1.legend(loc=AT_legend_loc, fontsize=14)
    if GC_legend:
        ax2.legend(loc=GC_legend_loc, fontsize=14)

    if make_fig:
        if save:
            plt.savefig("ATGCperiod_multi_" + note + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return [ax1, ax2]


def plot_polar (phases,
                amplts,
                phase_offset=0,
                colors=None,
                alphas=None,
                marker_sizes=16,
                texts=None,
                text_colors=None,
                text_sizes=7,
                text_ha='center',
                text_va='center',
                rticks=None,
                rtick_labels=[],
                rlabel_pos=135,
                fig_width=2.2,
                fig_height=2.2,
                labels=None,
                legend_loc='best',
                title=None,
                ax=None,
                note='',
                save=False):

    if not ax:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height),
                               subplot_kw={'projection': 'polar'})
        make_fig = True
    else:
        make_fig = False

    assert len(phases) == len(amplts)

    if colors == None:
        colors = [None] * len(phases)
    else:
        if type(colors) != list:
            colors = [colors] * len(phases)
            
    if alphas == None:
        alphas = [1] * len(phases)
    else:
        if type(alphas) != list:
            alphas = [alphas] * len(phases)

    if marker_sizes == None:
        marker_sizes = [16] * len(phases)
    else:
        if type(marker_sizes) != list:
            marker_sizes = [marker_sizes] * len(phases)
            
    if texts == None:
        texts = [None] * len(phases)
    else:
        if type(texts) != list:
            texts = [texts] * len(phases)
            
    if text_colors == None:
        text_colors = [None] * len(phases)
    else:
        if type(text_colors) != list:
            text_colors = [text_colors] * len(phases)

    if text_sizes == None:
        text_sizes = [16] * len(phases)
    else:
        if type(text_sizes) != list:
            text_sizes = [text_sizes] * len(phases)


    if labels == None:
        labels = [None] * len(phases)
        legend = False
    else:
        legend = True

    for i in range(len(phases)):
        phase = phases[i] + phase_offset
        amplt = amplts[i]
        color = colors[i]
        alpha = alphas[i]
        marker_size = marker_sizes[i]
        label = labels[i]
        text = texts[i]
        text_color = text_colors[i]
        text_size = text_sizes[i]
        
        ax.plot(phase,
                amplt,
                '.',
                markersize=marker_size,
                color=color,
                alpha=alpha,
                label=label)

        if text != None:
            ax.text(phase,
                    amplt,
                    text,
                    ha=text_ha,
                    va=text_va,
                    color=text_color,
                    size=text_size)

    ax.set_rlabel_position(rlabel_pos)
    ax.set_rticks(rticks) 
    ax.set_yticklabels(rtick_labels, fontsize=5)
    ax.tick_params(axis='both', which='major', labelsize=5, pad=-5)
    ax.tick_params(axis='both', which='minor', labelsize=5, pad=-5)

    if legend:
        ax.legend(loc=legend_loc,
                  fontsize=14)

    if title:
        ax.set_title(title)

    if make_fig:
        if save:
            plt.savefig("polar_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax

def plot_rlen_dist (rlen_counts,
                    colors=None,
                    labels=None,
                    alphas=None,
                    fig_width=3,
                    fig_height=2,
                    xlim=[None, None],
                    xticks=range(0, 500, 50),
                    xlabel="Read length (bp)",
                    ylabel="Read counts",
                    title="Read length distribution",
                    legend_loc='best',
                    save=False,
                    note=''):

    if ax == None:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    for i in range(len(rlen_counts)):
        rlen_count = rlen_counts[i]
        X = sorted(rlen_count.keys())
        Y = [rlen_count[x] for x in X]

        try:
            color = colors[i]
        except:
            color = 'black'

        try:
            label = labels[i]
        except:
            label = None

        try:
            alpha = alphas[i]
        except:
            alpha = 1
        
        ax.plot(X,
                Y,
                color=color,
                label=label,
                alpha=alpha)

    ax.set_xlim(xlim)
    ax.set_grid(True)
    
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(xtick) for xtick in xticks])
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if labels:
        ax.legend(loc=legend_loc)

    if make_fig:
        if save:
            plt.savefig("_".join(['rlen', note]) + ".svg",
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot NMF basis matrix
def plot_NMF_basis_matrix (basis_matrix,
                           basis_idxs,
                           features,
                           feature_cmaps=None,
                           xlabel='Property class',
                           title=None,
                           fig_width=5,
                           fig_height=6,
                           ax=None,
                           save=False,
                           note=''):

    if not ax:
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    assert len(basis_idxs) == len(basis_matrix)
    assert len(features) == len(basis_matrix[0])

    if not feature_cmaps:
        feature_cmaps = ['Greys'] * len(features)
        
    elif type(feature_cmaps) == list:
        if len(feature_cmaps) < len(features):
            repeat = float(len(features)) / len(feature_cmaps)
            repeat = int(math.ceil(repeat))
            feature_cmaps = feature_cmaps * repeat
            feature_cmaps = feature_cmaps[:len(features)]

    elif type(feature_cmaps) == str:
        feature_cmaps = [feature_cmaps] * len(features)
        
    imgs = []
    for i in range(len(features)):
        img = np.zeros((len(features), len(basis_idxs)))
        img[:] = np.nan
        for j in range(len(basis_idxs)):
            idx = basis_idxs[j]
            img[i][j] = basis_matrix[idx][i]
        imgs.append(img)

    for img, cmap in zip(imgs, feature_cmaps):
        plt.imshow(img, cmap=cmap, aspect='auto')

    for i in range(len(basis_idxs)):
        idx = basis_idxs[i]
        for j in range(len(features)):
            mean = np.mean(basis_matrix[:,j])
            std = np.std(basis_matrix[:,j])
            value = basis_matrix[idx][j]
            if value > mean + std:
                color = 'white'
            else:
                color = 'black'
            ax.text(i,
                    j,
                    str(round(value, 2)),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=color)

    ax.set_xticks(range(len(basis_idxs)))
    ax.set_xticklabels(range(1, len(basis_idxs)+1),
                       fontsize=8)

    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features,
                       fontsize=8)

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=10)

    if title:
        ax.set_title(title, fontsize=12)

    if make_fig:
        if save:
            plt.savefig("NMF_basis_matrix_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax


# plot p-value matrix for pair-wise pvalue data
def plot_pvalue_matrix (pair_pvalue,
                        keys=None,
                        key_label=None,
                        take_neglog10=False,
                        dummy=0.0,
                        rotation=45,
                        cmap='viridis',
                        vmin=None,
                        vmax=None,
                        cbar=True,
                        cbar_label=None,
                        fig_width=5,
                        fig_height=5,
                        save=False,
                        note='',
                        ax=None):

    if not ax:        
        fig, ax = plt.subplots(nrows=1,
                               ncols=1,
                               figsize=(fig_width, fig_height))
        make_fig = True
    else:
        make_fig = False

    if not keys:
        keys = sorted(pair_pvalue.keys())

    if not key_label:
        labels = keys
    else:
        labels = [key_label[key] for key in keys]

    matrix = np.zeros((len(keys), len(keys)))
    matrix[:] = np.nan

    for i in range(len(keys)-1):
        for j in range(i+1, len(keys)):
            pvalue = pair_pvalue[keys[i]][keys[j]]
            if take_neglog10:
                value = -np.log10(pvalue + dummy)
            else:
                value = pvalue

            matrix[i][j] = value
            matrix[j][i] = value

    img = ax.imshow(matrix,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels,
                       ha="right",
                       rotation_mode="anchor",
                       rotation=rotation)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)

    if cbar:
        if cbar_label == None:
            if take_neglog10:
                cbar_label = '-log10 p-value'
            else:
                cbar_label = 'p-value'
                
        cbar = plt.colorbar(img,
                            shrink=0.6)

        cbar.ax.set_ylabel(cbar_label,
                           rotation=-90,
                           va="bottom")

    if make_fig:
        if save:
            plt.savefig("pvalue_matrix_%s.svg" % (note),
                        format='svg',
                        bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()    
        plt.close()

    return ax
