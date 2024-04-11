import matplotlib.pyplot as plt
from matplotlib.patches import Patch, RegularPolygon
from matplotlib import cm, colors
import matplotlib

import numpy as np
import pandas as pd

import os

from cytospace.common import read_file



def format_label(label,max_length=14,max_lines=3):
    label_words = label.split(' ')
    label_chunks = []
    chunk = ''
    for word in label_words:
        if len(chunk)==0:
            if len(word)>max_length:
                word = word[:max_length]+'...'
                label_chunks.append(word)
            else:
                chunk += word
        elif len(chunk)+len(word)<=(max_length-1):
            chunk = chunk+' '+word
        else:
            label_chunks.append(chunk)
            chunk = word
    if len(chunk)>0:
        label_chunks.append(chunk)
    
    if len(label_chunks)>max_lines:
        label_chunks = label_chunks[:max_lines]
        last_line = label_chunks[-1]
        if len(last_line)>(max_length-3):
            label_chunks[-1] = last_line[:max_length]+'...'
        else: 
            label_chunks[-1] = last_line+'...'

    return '\n'.join(label_chunks)

def rand_jitter(arr,interval):
    return arr + np.random.uniform(-interval/2,interval/2,len(arr))

def plot_results_bulk_ST_by_spot(coordinates, metadata, dir_out, output_prefix, geometry='honeycomb', num_cols=3):
 
    # Define output files
    fout_png_all = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot.png")
    fout_pdf_all = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot.pdf")
    
    coordinates = coordinates.loc[metadata.index,:]
    cell_types = list(metadata.columns)[:-1]
    cell_types = list(np.sort(cell_types))
    cell_types.insert(0,'Total cells')

    X = coordinates.iloc[:,0].astype(float)
    Y = coordinates.iloc[:,1].astype(float)
    
    # distinguish between row/col indices and coordinates
    # based on range (500) and type (int vs. float) of values
    scale = Y.max() < 500 and ((Y - Y.round()).abs() < 1e-5).all()

    # define representative interval between each adjacent spot/point
    y_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(Y)))))
    x_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(X)))))
    
    if geometry == 'honeycomb' and scale:
        print('Detecting row and column indexing of Visium data; rescaling for coordinates')
        
        # Rotate
        X_prev = X
        Y_prev = Y
        X = Y_prev
        Y = 1-X_prev
        
        # Rescale
        Y = 1.75*Y

    elif geometry == 'square' and scale:
        print('Detecting row and column indexing of legacy ST data; rotating for coordinates')
        
        # Rotate
        X_prev = X
        Y_prev = Y
        X = Y_prev
        Y = 1-X_prev

    else:        
        # Rotate 
        Y = 1-Y
        
    if geometry == 'honeycomb':
        hex_vert = 6

        if scale:
            hex_rot = 0
            hex_rad = y_int
            hex_rad = hex_rad + 0.2*hex_rad

        else:
            hex_rot = 0
            hex_rad = x_int
            hex_rad = hex_rad + 0.28*hex_rad

    elif geometry == 'square':
        hex_vert = 4
        hex_rot = 45
        interval = y_int
        hex_rad = 0.5*np.sqrt(2*interval**2)

    else:
        print("Unknown geometry specified.")
        exit()

    num_rows = int(len(cell_types)/num_cols)
    if num_rows*num_cols < len(cell_types):
        num_rows = num_rows+1
    width = max(X)-min(X)
    height = max(Y)-min(Y)

    plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'], 'size':'12'})
    plt.rcParams['figure.dpi'] = 450
    plt.rcParams['savefig.dpi'] = 450

    fig, axes = plt.subplots(num_rows,num_cols,figsize=(width/height*3*num_cols,3*num_rows))
    full_frac = 0.047 * (3*num_rows / (width/height*3*num_cols))
    k = 0
    for i in range(num_rows):
        for j in range(num_cols):

            ax = axes[i,j]

            if k >= len(cell_types):
                ax.axis('off')
            else:
                ct = cell_types[k]
                ax.set_aspect('equal')

                node_assignment = metadata.loc[:,ct]

                viridis = cm.get_cmap('viridis')
                norm = matplotlib.colors.Normalize(min(node_assignment), max(node_assignment))

                node_assignment = (1/max(node_assignment))*node_assignment
                colors = viridis(node_assignment)

                for x, y, c in zip(X, Y, colors):
                    hex = RegularPolygon((x, y), numVertices=hex_vert, radius=hex_rad, 
                                         orientation=np.radians(hex_rot), 
                                         facecolor=c, edgecolor=None)
                    ax.add_patch(hex)

                # Also add scatter points in hexagon centres - not sure why this line has to be here
                ax.scatter(X, Y, c=[c[0] for c in colors],alpha=0)
                ct_label = ct
                # Reformat labels that are too long
                if len(ct_label)>14:
                    ct_label = format_label(ct_label)
                ax.set_title(ct_label)
                cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=viridis),ax=ax,label='',fraction=0.036, pad = 0.04)
                ax.axis('off')

            k += 1

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.suptitle('Number of cells per spot mapped by CytoSPACE')        
    fig.savefig(fout_png_all, facecolor="w", bbox_inches='tight')
    fig.savefig(fout_pdf_all, facecolor="w", bbox_inches='tight')


def plot_results_bulk_ST_jitter(assigned_locations, dir_out, output_prefix, geometry="honeycomb",max_num_cells=50000):
    
    # Define output files
    fout_png_jitter = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot_jitter.png")
    fout_pdf_jitter = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot_jitter.pdf")

    if assigned_locations.shape[0] > max_num_cells:
        assigned_locations = assigned_locations.sample(max_num_cells)

    X = assigned_locations.iloc[:,-2]
    Y = assigned_locations.iloc[:,-1]
    cell_types = assigned_locations['CellType'].values.astype(str)

    # distinguish between row/col indices and coordinates
    # based on range (500) and type (int vs. float) of values
    scale = Y.max() < 500 and ((Y - Y.round()).abs() < 1e-5).all()

    # define representative interval between each adjacent spot/point
    y_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(Y)))))
    x_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(X)))))


    if geometry == 'honeycomb' and scale:
        print('Detecting row and column indexing of Visium data; rescaling for coordinates')
        
        # Rotate
        X_prev = X
        Y_prev = Y
        X = Y_prev
        Y = 1-X_prev
        
        # Rescale
        Y = 1.75*Y
        y_interval = 1.75*x_int
        x_interval = 1.75*y_int

    elif geometry == 'square' and scale:
        print('Detecting row and column indexing of legacy ST data; rotating for coordinates')
        
        # Rotate
        X_prev = X
        Y_prev = Y
        X = Y_prev
        Y = 1-X_prev
        y_interval = 1.75*x_int
        x_interval = y_int

    else:
        # Rotate 
        Y = 1-Y
        y_interval = y_int
        x_interval = x_int

    X = rand_jitter(X.values,x_interval)
    Y = rand_jitter(Y.values,y_interval)


    plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'], 'size':'12'})
    plt.rcParams['figure.dpi'] = 450
    plt.rcParams['savefig.dpi'] = 450
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_aspect('equal')
    colors = ["#222222", "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
                  "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", 
                  "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", "#5A5156", "#E4E1E3", 
                  "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
                  "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
                  "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6",
                  "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79",
                  "#66B0FF", "#3B00FB"]
    
    unique_cell_types = np.unique(cell_types)
    color_per_ct = dict(zip(unique_cell_types,colors[:len(unique_cell_types)]))
    cell_type_colors = [color_per_ct[ct] for ct in cell_types]
    ax.scatter(X,Y,c=cell_type_colors,s=2)
    ax.axis('off')
    plt.title('Single cells mapped to tissue by CytoSPACE')
    fig.tight_layout()

    legend_elements = []
    for ct in unique_cell_types:
        legend_elements.append(Patch(facecolor=color_per_ct[ct],label=format_label(ct)))
        
    plt.legend(bbox_to_anchor=(1.05,1.0),handles=legend_elements)
    fig.savefig(fout_png_jitter, facecolor="w", bbox_inches='tight')
    fig.savefig(fout_pdf_jitter, facecolor="w", bbox_inches='tight')


def plot_results_single_cell_ST(assigned_locations, dir_out, output_prefix, max_num_cells=50000):
    
    # Define output files
    fout_png = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot_single_cell.png")
    fout_pdf = os.path.join(dir_out, f"{output_prefix}cell_type_assignments_by_spot_single_cell.pdf")

    if assigned_locations.shape[0] > max_num_cells:
        assigned_locations = assigned_locations.sample(max_num_cells)

    X = assigned_locations.iloc[:,-2].values
    Y = assigned_locations.iloc[:,-1].values
    cell_types = assigned_locations['CellType'].values.astype(str)

    Y = 1-Y # Convert to top left origin - this is standard

    plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial'], 'size':'12'})
    plt.rcParams['figure.dpi'] = 450
    plt.rcParams['savefig.dpi'] = 450
    
    fig = plt.figure()
    ax = fig.add_subplot(111)


    ax.set_aspect('equal')

    colors = ["#222222", "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
                  "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", 
                  "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", "#5A5156", "#E4E1E3", 
                  "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
                  "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
                  "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6",
                  "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79",
                  "#66B0FF", "#3B00FB"]
        
    unique_cell_types = np.unique(cell_types)
    color_per_ct = dict(zip(unique_cell_types,colors[:len(unique_cell_types)]))
    cell_type_colors = [color_per_ct[ct] for ct in cell_types]
    ax.scatter(X,Y,c=cell_type_colors,s=2)
    ax.axis('off')
    plt.title('Single cells mapped to tissue by CytoSPACE')
    fig.tight_layout()

    legend_elements = []
    for ct in unique_cell_types:
        legend_elements.append(Patch(facecolor=color_per_ct[ct],label=format_label(ct)))
        
    plt.legend(bbox_to_anchor=(1.05,1.0),handles=legend_elements)
    fig.savefig(fout_png, facecolor="w", bbox_inches='tight')
    fig.savefig(fout_pdf, facecolor="w", bbox_inches='tight')

def plot_results(output_directory, output_prefix, assigned_locations_path=None, metadata_path=None, coordinates_data=None, 
                        geometry='honeycomb',num_cols=3, single_cell_ST_mode=False,max_num_cells=50000):
    # check input files
    if assigned_locations_path is None:
        assigned_locations_path = os.path.join(output_directory, f"{output_prefix}assigned_locations.csv")
        if not os.path.exists(assigned_locations_path):
            raise ValueError("Cannot locate assigned_locations.csv file required for plotting")
    
    if (not single_cell_ST_mode) and (metadata_path is None):
        metadata_path = os.path.join(output_directory, f"{output_prefix}cell_type_assignments_by_spot.csv")
        if not os.path.exists(metadata_path):
            raise ValueError("Cannot locate metadata_path required for plotting bulk ST data")
    
    if (not single_cell_ST_mode) and (coordinates_data is None):
        raise ValueError("coordinates_data is required unless in single_cell_ST_mode")

    assigned_locations = read_file(assigned_locations_path)

    if single_cell_ST_mode:
        plot_results_single_cell_ST(assigned_locations, output_directory, output_prefix, max_num_cells=max_num_cells)
    else:
        metadata = read_file(metadata_path)
        metadata.index = ['SPOT_'+str(idx) for idx in metadata.index]
        plot_results_bulk_ST_by_spot(coordinates_data, metadata, output_directory, output_prefix, geometry=geometry, num_cols=num_cols)
        plot_results_bulk_ST_jitter(assigned_locations, output_directory, output_prefix, geometry=geometry, max_num_cells=max_num_cells)



