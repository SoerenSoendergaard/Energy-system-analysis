# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 10:24:16 2024

@author: soere
"""


import PyMGA
from PyMGA.utilities.plot import near_optimal_space_2D, near_optimal_space_matrix,set_options
import numpy as np
import yaml
import pandas as pd
import statistics 
import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import warnings
import dill
import os
from scipy.spatial import ConvexHull
from matplotlib.lines import Line2D

import PyMAA
from PyMAA.utilities.plot import near_optimal_space_matrix
import numpy as np
import yaml
import pandas as pd

from matplotlib.path import Path
import matplotlib.cm as cm

#%% Plotting statisticks 
def get_subfolder_names(directory_path):
    # Get a list of all items in the directory
    items = os.listdir(directory_path)

    # Filter out subfolders
    subfolders = [item for item in items if os.path.isdir(os.path.join(directory_path, item))]

    # Return the list of subfolder names
    return subfolders

def PlotContoursAndPoints(PathToContours,PathToPoints,chosen_variables,MAA_variables):
    
    ## First read and plot the points
    MAA_Samples_df = pd.read_csv(PathToPoints)
    # Append 'Cost' as a separate element to the chosen_variables list
    chosen_variables_Cost = chosen_variables.copy()
    chosen_variables_Cost.append('Cost')
    variable_samples = MAA_Samples_df[chosen_variables_Cost]
    
    # plt.scatter(x=variable_samples[f'{chosen_variables[0]}'], 
    #             y=variable_samples[f'{chosen_variables[1]}'], 
    #             c=variable_samples['Cost'], cmap='viridis',alpha=0.6)  # Adjust the colormap as needed
    
    # # Add a color bar
    # cbar = plt.colorbar()
    # cbar.set_label('Cost')  # Set the label for the color bar
    
    ## Read and plot the contours
    Contours = get_subfolder_names(PathToContours)
    NrContours = len(Contours)
    
    for i in range(NrContours):
    # for i in range(0):
        slack = Contours[i]
        file_name = f'Verticies'
        file_path_v = os.path.join(PathToContours, slack, file_name)
        verticies = np.loadtxt(file_path_v, delimiter=',', dtype=int)
        
        All_verticies = pd.DataFrame(verticies,
                                    columns = MAA_variables)
        
        variable_verticies = All_verticies[chosen_variables]
        

        hull = ConvexHull(variable_verticies.values)
        
        # plot simplexes
        for simplex in hull.simplices:
            l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                    color = 'k',
                    linewidth = 2, zorder = 1)
        
    
    plt.xlabel(f'{chosen_variables[0]} [MW]')
    plt.ylabel(f'{chosen_variables[1]} [MW]')
    plt.title('Near optimal space - Greenfield')

def PlotContours(PathToContours,chosen_variables,MAA_variables):
    
    ## First read and plot the points
  
    ## Read and plot the contours
    
    for i in range(len(PathToContours)):
    # for i in range(0):
       
        verticies = np.loadtxt(PathToContours[i], delimiter=',', dtype=int)
        
        All_verticies = pd.DataFrame(verticies,
                                    columns = MAA_variables)
        
        variable_verticies = All_verticies[chosen_variables]
        

        hull = ConvexHull(variable_verticies.values)
        
        # plot simplexes
        for simplex in hull.simplices:
            l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                    color = 'k',
                    linewidth = 2, zorder = 1)
        
    
    plt.xlabel(f'{chosen_variables[0]} [MW]')
    plt.ylabel(f'{chosen_variables[1]} [MW]')
    plt.title('Near optimal space - Greenfield')
    
    
def PlotContoursAndPoints_2(PathToContours,PathToPoints,chosen_variables,MAA_variables):
    
    ## First read and plot the points
    MAA_Samples_df = pd.read_csv(PathToPoints)
    # Append 'Cost' as a separate element to the chosen_variables list
    chosen_variables_Cost = chosen_variables.copy()
    chosen_variables_Cost.append('Cost')
    variable_samples = MAA_Samples_df[chosen_variables_Cost]
    
    # plt.scatter(x=variable_samples[f'{chosen_variables[0]}'], 
    #             y=variable_samples[f'{chosen_variables[1]}'], 
    #             c=variable_samples['Cost'], cmap='viridis')  # Adjust the colormap as needed
    
    ## Read and plot the contours
    Contours = get_subfolder_names(PathToContours)
    NrContours = len(Contours)
    assigned_points = set()
    for i in range(NrContours):
        slack = Contours[i]
        file_name = f'Verticies'
        file_path_v = os.path.join(PathToContours, slack, file_name)
        verticies = np.loadtxt(file_path_v, delimiter=',', dtype=int)

        All_verticies = pd.DataFrame(verticies, columns=MAA_variables)
        variable_verticies = All_verticies[chosen_variables]

        hull = ConvexHull(variable_verticies.values)

        # Create a path object for the convex hull
        path = Path(variable_verticies.values[hull.vertices])

        # Check which points are inside the convex hull
        points_inside = path.contains_points(variable_samples[[chosen_variables[0], chosen_variables[1]]].values)
        points_inside = np.logical_and(points_inside, ~np.isin(range(len(points_inside)), list(assigned_points)))

        # Assign a unique color to each contour
        color = cm.viridis(i / NrContours)  # Adjust the colormap as needed

        # Scatter plot for points inside the convex hull with unique color
        plt.scatter(x=variable_samples.loc[points_inside, chosen_variables[0]], 
                    y=variable_samples.loc[points_inside, chosen_variables[1]], 
                    c=color,  # Use a different color for points inside each contour
                    marker='x',  # Use a different marker for points inside the contour
                    s=50,  # Adjust marker size as needed
                    label=f'Points inside Contour {slack}')

        # Update the set of assigned points
        assigned_points.update(np.where(points_inside)[0])
        
        # Plot simplexes for convex hull
        for simplex in hull.simplices:
            plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                    color='k', linewidth=2, zorder=2)
        
    
    plt.xlabel(f'{chosen_variables[0]} [MW]')
    plt.ylabel(f'{chosen_variables[1]} [MW]')
    plt.title('Near optimal spaces - Greenfield')


    
def PlotOuterCountoursAndPoints(chosen_variables,OuterContours ,Points,MAA_variables):
   # combined_samples_list = []  # Initialize an empty list to store DataFrames
    ColoScheme = ['spring','winter']
    for i in range(len(Points)):
        MAA_Samples_df_1 = pd.read_csv(Points[i])
        # Append 'Cost' as a separate element to the chosen_variables list
        chosen_variables_Cost = chosen_variables.copy()
        chosen_variables_Cost.append('Cost')
        variable_samples = MAA_Samples_df_1[chosen_variables_Cost]
        # Append the DataFrame to the list
      #  combined_samples_list.append(variable_samples)
    
        # Concatenate the i DataFrames
       # combined_samples = pd.concat(combined_samples_list, ignore_index=True)
        
        plt.scatter(x=variable_samples[f'{chosen_variables[0]}'], 
                    y=variable_samples[f'{chosen_variables[1]}'], 
                    c=variable_samples['Cost'], cmap=ColoScheme[i], alpha=0.1)  # Adjust the colormap as needed
        
        # Add a color bar
        cbar = plt.colorbar()
        #cbar.set_label('Cost 60')  # Set the label for the color bar
        
    LineStyle = ['k-','k--']
    for i in range(len(OuterContours)):
        verticies = np.loadtxt(OuterContours[i], delimiter=',', dtype=int)
        
        All_verticies = pd.DataFrame(verticies,
                                    columns = MAA_variables)
        
        variable_verticies = All_verticies[chosen_variables]
        
    
        hull = ConvexHull(variable_verticies.values)
        
        # plot simplexes
        for simplex in hull.simplices:
            l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 
                           LineStyle[i], color = 'k',
                    linewidth = 2, zorder = 1)
        
def Trans_Histogram(Samples):
    Samples_df = pd.read_csv(Samples[0])
    # Create a single plot with multiple bar plots
    

    fig, axs = plt.subplots(len(Samples_df.columns), 1, figsize=(5, 5 * len(Samples_df.columns)))

    # Loop through each column
    for j, col in enumerate(Samples_df.columns):
        axs[j].hist(Samples_df[col], bins=20, orientation='horizontal', color='skyblue', edgecolor='black')
        axs[j].set_title(col)
        axs[j].set_ylabel(col)
        axs[j].set_xlabel('Frequency')

    plt.tight_layout()
    plt.show()
    
    # for i in range(len(Samples)):
    #     Samples_df = pd.read_csv(Samples[i])
        
        
def Trans_Histogram_2(Sample, MAA_variables ,Type=None):

    # Filter columns based on the specified component type
    if Type:
        columns_to_plot = [col for col, values in MAA_variables.items() if values[0] == Type]
    else:
        columns_to_plot = Sample.columns
    
    
    # Create a single plot with overlaid histograms
    fig, ax = plt.subplots(figsize=(8, 6))

    # Loop through each column
    for j, col in enumerate(columns_to_plot):
        ax.hist(Sample[col], bins=20, alpha=0.5, edgecolor='black', label=col, color=f'C{j}')

    ax.set_title('Overlayed Histograms of Sample Columns')
    ax.set_xlabel('Value')
    ax.set_ylabel('Frequency')
    ax.legend()
    #plt.show()
        
        
def plot_histograms_shared_x(MAA_variables, Type=None, bins=20, alpha=0.5, 
                              edgecolor='black', linewidth=1.2, *samples, labels=None):
    
    # Filter columns based on the specified component type
    if Type:
        columns_to_plot = [col for col, values in MAA_variables.items() if values[0] == Type]
    else:
        columns_to_plot = samples[0].columns
    
    # Create a single plot with overlaid histograms
    fig, axs = plt.subplots(len(samples), 1, figsize=(8, 6*len(samples)), sharex=True)

    # Loop through each sample
    for i, sample in enumerate(samples):
        # Determine label for the legend
        label = labels[i] if labels else f'Sample {i+1}'

        # Loop through each column
        for j, col in enumerate(columns_to_plot):
            axs[i].hist(sample[col], bins=bins, alpha=alpha, edgecolor=edgecolor,
                        linewidth=linewidth, label=f'{col}')
            
            # Convert y-axis to percentage occurrence
            axs[i].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{(x / len(sample[col])) * 100:.0f}%'))


        axs[i].set_title(f'Histogram - {label} reduction scenario')
        axs[i].set_ylabel('Percentage of samples %')
        axs[i].legend()

    axs[-1].set_xlabel('Installed Capacity [MW]')
    plt.show()   
        
        
def boxplots_shared_x(sample, MAA_variables, Type=None):
    # Create a single plot with overlaid boxplots
    fig, ax = plt.subplots(figsize=(12, 6))

    # Filter columns based on the specified component type
    if Type:
        columns_to_plot = [col for col, values in MAA_variables.items() if values[0] == Type]
    else:
        columns_to_plot = sample.columns

    # Select columns individually using iloc
    boxplots = sample.iloc[:, [sample.columns.get_loc(col) for col in columns_to_plot]].boxplot(grid=False, ax=ax, vert=False, patch_artist=True)

    # Add labels and title
    ax.set_xlabel('Values')
    ax.set_title(f'Boxplots of {", ".join(columns_to_plot)} Columns')

# Add information or annotations
    info_text = f"Type: {Type}\n"
    info_text += "Boxplot Explanations:\n"
    info_text += " - Box: Middle 50% of values (Interquartile Range)\n"
    info_text += f" - Lower Whisker: {round(0.25 * 100)}th percentile\n"
    info_text += " - Line inside box: Median (50th percentile)\n"
    info_text += f" - Upper Whisker: {round(0.75 * 100)}th percentile\n"
    info_text += " - Dots: Outliers\n"
    ax.text(0.7, 0.9, info_text, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))



    plt.show()    
    
    
def boxplots_shared_x_multiple(dataframes_list,Reductions, MAA_variables, Type=None):
    # Create subplots with shared x-axis
    
    
    fig, axs = plt.subplots(len(dataframes_list), 1, figsize=(12, 8), sharex=True, sharey=True)

    for i, samples_df in enumerate(dataframes_list):
        # Filter columns based on the specified component type
        if Type:
            columns_to_plot = [col for col, values in MAA_variables.items() if values[0] == Type]
        else:
            columns_to_plot = samples_df.columns

        # Select columns individually using iloc
        # boxplots = samples_df.iloc[:, [samples_df.columns.get_loc(col) for col in columns_to_plot]].boxplot(grid=False, ax=axs[i], vert=False, patch_artist=True)

# Select columns individually using iloc
        boxplots = samples_df.iloc[:, [samples_df.columns.get_loc(col) for col in columns_to_plot]].boxplot(
            grid=False, 
            ax=axs[i], 
            vert=False, 
            patch_artist=False,
            showfliers=False,  # Exclude outliers
            whis=[5, 95]  # Set whiskers to the 5th and 95th percentiles
        )
        
    

        # Add labels and title
        # axs[i].set_title(f'Boxplots of {", ".join(columns_to_plot)} Columns - scenario')
        # axs[i].set_title(rf'{Reductions[i]}% - Reduction')
        # axs[i].legend(rf'{Reductions[i]}% - Reduction')

    fig.suptitle('Boxplots of generator capacity', fontsize=16)
    plt.xlabel('Capacity [MW]')
    plt.tight_layout()
    plt.savefig("Boxplts.pdf", format="pdf", bbox_inches="tight")
    plt.show()


def CapacityModes(dataframes_list, Reductions, MAA_variables, Type=None):
    # Create subplots with shared x-axis
    plt.figure(figsize=(10, 6))

    # if len(dataframes_list) == 1:
    #     axs = [axs]  # Ensure axs is a list even when there's only one subplot

    for i, samples_df in enumerate(dataframes_list):
        # Filter columns based on the specified component type
        if Type:
            columns_to_plot = [col for col, values in MAA_variables.items() if values[0] == Type]
        else:
            columns_to_plot = samples_df.columns
            
            
        # Make line plot of each mode
    Colors = ["lightblue", "darkblue", "yellow", "green"]
    j=-1
    for column in columns_to_plot:
        j=j+1
        Mode = np.zeros(len(Reductions))
        LowerBound = np.zeros(len(Reductions))
        UpperBound = np.zeros(len(Reductions))
        for i, samples_df in enumerate(dataframes_list):
            counts, bin_edges = np.histogram(samples_df[column], bins=50)
            max_bin_index = np.argmax(counts)
            Mode[i] = (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]) / 2
            
            
            # Calculate the top 5% range
            total_counts = np.sum(counts)
            top_5_percent_count = 0.25 * total_counts
            accumulated_counts = 0
            lower_index, upper_index = max_bin_index, max_bin_index

            while accumulated_counts < top_5_percent_count:
                if lower_index > 0:
                    lower_index -= 1
                    accumulated_counts += counts[lower_index]
                if accumulated_counts < top_5_percent_count and upper_index < len(counts) - 1:
                    upper_index += 1
                    accumulated_counts += counts[upper_index]
            
            LowerBound[i] = bin_edges[lower_index]
            UpperBound[i] = bin_edges[upper_index + 1]
            
            
        plt.plot(Reductions,Mode,Colors[j],linewidth = 3, label=column)
        
        # Plot the top 5% range
        plt.fill_between(Reductions, LowerBound, UpperBound, color=Colors[j % len(Colors)], alpha=0.3)

        
    plt.legend()
    plt.xlabel("CO2 reduction [%]")
    plt.ylabel("Installed capacity [MW]")
    plt.title("Most frequent optimal capacity\n Shaded area = bounds of the top 10% most frequent capacities")
    plt.savefig("Frequencies.pdf", format="pdf", bbox_inches="tight")
    plt.show()

        


#%% Computing pathways
# def GeneratePaths(dataframes_list,NrPaths,Decay,StartingPoint):
    
#     LengthOfPath = len(dataframes_list)
    
#     # Create the desired number of paths
#     TransitionPath = []
#     TransitionPath.append(StartingPoint)
    
#     AllPaths = []
#     while len(AllPaths)<NrPaths:
#         print(f'{len(AllPaths)}')
#         Point = StartingPoint
#         TransitionPath = []
#         TransitionPath.append(StartingPoint)

#         for j in range(LengthOfPath):
#         #while k < (LengthOfPath):
#             # If starting point is defined (it needs to be a dataframe)
#             if Decay == "no":
                
#                 # Find legal paths
#                 ValidPaths = dataframes_list[j][(Point.squeeze().values <= dataframes_list[j].values).all(axis=1)]
                
#                 Point = ValidPaths.sample(n=1)
                
                
#                 # Only accept point if it has valid path options
#                 HasValidConnections = "no"
#                 while HasValidConnections == "no" and j<LengthOfPath-1:
                    
#                     # Check if this point has valid paths in all spaces to come
#                     k = 0
#                     for i in range(LengthOfPath-1):
#                         ValidPaths_temp = dataframes_list[j+i][(Point.squeeze().values <= dataframes_list[j+i].values).all(axis=1)]
                   
#                         if len(ValidPaths_temp) == 0:
#                             # No valid paths, pick a new random point
#                             Point = ValidPaths.sample(n=1)
                            
#                         elif j<LengthOfPath-2: # Maybe it increases robustness to check two steps into the future
#                             for z in range(len(ValidPaths_temp)):
#                                 # Check if one of the points in the valid path has a valid path
#                                 tempPoint = ValidPaths_temp.iloc[z]
#                                 ValidPaths_temp_2 = dataframes_list[j+i+1][(tempPoint.squeeze().values <= dataframes_list[j+i+1].values).all(axis=1)]
                                
#                                 if len(ValidPaths_temp_2) != 0:
#                                     k = 1
#                                 else:
#                                     # Ensure that this point is not used
#                                     2+2
                                    
#                         else:
#                             k = 1

#                         if k == 1:                               
#                             HasValidConnections = "yes"
                    
                    
#                 TransitionPath.append(Point)
                
#         result_df = pd.concat(TransitionPath, ignore_index=True)
#         TransitionPath = []
#         AllPaths.append(result_df)   
        
#     return AllPaths


def GenerateRobustPaths(dataframes_list, NrPaths, Decay, StartingPoint,MAA_variables):
    if Decay == "no":

        LengthOfPath = len(dataframes_list)

        # Create the desired number of paths
        TransitionPath = []
        TransitionPath.append(StartingPoint)
        AllPaths = []
        
        while len(AllPaths) < NrPaths:
            print(f"{len(AllPaths)}")
            Point = StartingPoint
            TransitionPath = []
            TransitionPath.append(StartingPoint)
            
            ## First, try without "Pick a random valid endpoint"
            # Create a new 6-dimensional space representing the intersection
            intersected_space = np.zeros((LengthOfPath, 6))
            
            for j in range(LengthOfPath):
                
                # Find valid solution space from the verticies and the relevant point
                tmp = dataframes_list[j].copy()
                verticies = tmp
                
                for z in range(len(verticies)):
                    for k in range(verticies.shape[1]):
                        # Check each dimension 'k' of each point 'z'
                        
                        if Point.iloc[0,k] > verticies.iloc[z,k]:
                            # Move the dimension of the vertice, to the point
                            verticies.iloc[z,k] = Point.iloc[0,k]
                            
                            
                ### Plotting for debugging
                chosen_variables = ['OnWind', 'OffWind']
                variable_verticies = verticies[chosen_variables]
                

                hull = ConvexHull(variable_verticies.values)
                
                # plot simplexes
                for simplex in hull.simplices:
                    l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                            color = 'k',
                            linewidth = 2, zorder = 1)
                
            
                # Also plot original shape
                verticies_orig = dataframes_list[j]
                chosen_variables = ['OnWind', 'OffWind']
                variable_verticies = verticies_orig[chosen_variables]
                

                hull = ConvexHull(variable_verticies.values)
                
                # plot simplexes
                for simplex in hull.simplices:
                    l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                            color = 'r', alpha = 0.6,
                            linewidth = 2, zorder = 1)
                
                
                plt.scatter(x=Point[f'{chosen_variables[0]}'], 
                            y=Point[f'{chosen_variables[1]}'])  # Adjust the colormap as needed
                
                plt.xlabel(f'{chosen_variables[0]} [MW]')
                plt.ylabel(f'{chosen_variables[1]} [MW]')
                plt.title('Near optimal spaces - Greenfield')
                plt.show()
                    
                        
                
            

    

def GeneratePaths(dataframes_list, NrPaths, Decay, StartingPoint):

    if Decay == "no":

        LengthOfPath = len(dataframes_list)

        # Create the desired number of paths
        TransitionPath = []
        TransitionPath.append(StartingPoint)

        AllPaths = []
        while len(AllPaths) < NrPaths:
            print(f"{len(AllPaths)}")
            Point = StartingPoint
            TransitionPath = []
            TransitionPath.append(StartingPoint)
            
            ## Pick a random valid endpoint
            ValidEndPoints = dataframes_list[-1].loc[
                (Point.squeeze().values <= dataframes_list[-1].values).all(axis=1)]
            EndPoint = ValidEndPoints.sample(n=1)

            # Find upper bounds for each variable in the range of all dataframes
            max_values = [df.max().values for df in dataframes_list]

            # Choose the minimum value for each variable among the found maxes
            max_bound_values = np.min(max_values, axis=0)
            min_indices = np.argmin(max_values, axis=0)

            # Create a dataframe with one row and the same column names as the dataframes
            max_bound = pd.DataFrame([max_bound_values], columns=dataframes_list[0].columns)
         
            #NrMiddlePoint = LengthOfPath - 1
            for j in range(LengthOfPath):
                # If the boundary is defined from the j'th dataframe, move the boundary
                
                if any(min_indices==j) and j < LengthOfPath-1:
                    # Just move all the boundaries (boundaries not at j'th, will be recalc not moved)
                    dataframes_list_moved = dataframes_list[j+1:]
                    max_values = [df.max().values for df in dataframes_list_moved]

                    max_bound_values = np.min(max_values, axis=0)
                    min_indices = np.argmin(max_values, axis=0)+j+1

                    max_bound = pd.DataFrame([max_bound_values], columns=dataframes_list[0].columns)

                

                # Find legal paths within the defined range, both using the end point
                # and the moving boundaries
                
                # First find all points below the boundaries
                ValidPaths = dataframes_list[j][
                    (Point.squeeze().values <= dataframes_list[j].values).all(axis=1)
                    & (dataframes_list[j].values <= max_bound.values).all(axis=1)
                ]
                
                # Then exclude those points which have at least one dimension above the endpoint
                ValidPaths_2 = ValidPaths[(ValidPaths.values <= EndPoint.squeeze().values).all(axis=1)]
                if len(ValidPaths_2) == 0:
                    # Pick a new random endpoint
                    NrValidEndPoints = len(ValidEndPoints)
                    print(f'Starting to try and move endpoint')
                    # for i in range(round(NrValidEndPoints * 0.0001)):
                    for i in range(4):
                        # Print a message when a quarter of the dataframes have been processed
                        if (i + 1) % (4 // 4) == 0:
                            print(f"{i + 1}")
                        ## Pick a random valid endpoint
                        ValidEndPoints = dataframes_list[-1].loc[
                            (Point.squeeze().values <= dataframes_list[-1].values).all(axis=1)]
                        EndPoint = ValidEndPoints.sample(n=1)
                        ValidPaths_2 = ValidPaths[(ValidPaths.values <= EndPoint.squeeze().values).all(axis=1)]
                        
                        if len(ValidPaths_2) > LengthOfPath-j:
                            print(f'End point was moved in step {j+1}')
                            break
                
                # # Well the random approach might have missed some
                # # So if the random approach fails maybe systematically go through the
                # # possible endpoint movements
                # if len(ValidPaths_2) == 0:
                #     # Pick a new random endpoint
                #     print(f'Starting systematic check')
                #     ValidEndPoints = dataframes_list[-1].loc[
                #         (Point.squeeze().values <= dataframes_list[-1].values).all(axis=1)]
                #     NrValidEndPoints = len(ValidEndPoints)
                #     for i in range(NrValidEndPoints):
                        
                #         # Print a message when a quarter of the dataframes have been processed
                #         if (i + 1) % (len(NrValidEndPoints) // 4) == 0:
                #             print(f"{i + 1} dataframes processed out of {len(dataframes_list)}")
                        
                #         ## Pick a random valid endpoint
                #         EndPoint = ValidEndPoints.iloc[i,:]
                #         ValidPaths_2 = ValidPaths[(ValidPaths.values <= EndPoint.squeeze().values).all(axis=1)]
                        
                #         if len(ValidPaths_2) != 0:
                #             print(f'End point was moved in step {j+1} (systematic check)')
                #             break        
                
                if len(ValidPaths_2) == 0:
                     print(f'failed after 1000 endpoint moves, returning any completed paths')
                     

                     return AllPaths
 
                         
                        
                
                ## Plot for debugging:
                    # I have 6 columns of data and so i want to make a 3 times 3 subplot
                    # The title of each subplut can be the column name. In each subplot all
                    # the datapoints of that row must be plotted for the starting point and the 3 dataframes
                    # So that in subplot 1,1 at x = 0, the y axis must have the value of the first column of the
                    # starting point, and at x = 1 all the datapoints from the first column of dataframe 1 must
                    # be plotted, and so on
                    # Then i also want to plot the upperboundary from "max_bound" of that plot as a vertical line.

                # Sample a valid point
                Point = ValidPaths_2.sample(n=1)

                # # Plot for debugging
                # fig, axes = plt.subplots(2, 3, figsize=(15, 10))
                # fig.suptitle(f"Column Values - Path: {len(AllPaths) + 1}, Step: {j + 1}")
                
                # for col in range(3):
                #     ax = axes[0, col]
                #     col_name = StartingPoint.columns[col]
                #     df_values_start = StartingPoint.iloc[:, col].values
                #     df_values_point = Point.iloc[:, col].values
                #     df_values_endpoint = EndPoint.iloc[:, col].values
                #     df_values_data = [df[col_name].values for df in dataframes_list]
                #     ax.scatter([0] * len(df_values_start), df_values_start, color='red', label='Starting Point')
                #     ax.scatter([j+1] * len(df_values_point), df_values_point, color='black', label='Chosen Point',zorder=10)
                #     ax.scatter([LengthOfPath] * len(df_values_endpoint), df_values_endpoint, color='blue', label='Chosen end point',zorder=10)
                    
                #     for i, values in enumerate(df_values_data, start=1):
                #         ax.scatter([i] * len(values), values, label=f'Dataframe {i}')
                
                    
                #     ax.scatter([j+1] * len(ValidPaths), ValidPaths[col_name], color='pink', marker='o', label='Valid Points (moving)')
                #     #ax.scatter([j+1] * len(ValidPaths), ValidPaths_2[col_name], color='red', marker='o', label='Valid Points (endpoint)')

                #     ax.axhline(max_bound.iloc[0, col], color='green', linestyle='--', label='Upper Boundary')
                #     ax.set_title(f'Column {col_name}')
                #     ax.legend()
                
                # for col in range(3):
                #     ax = axes[1, col]
                #     col_name = StartingPoint.columns[col + 3]  # Use columns 4, 5, and 6
                #     df_values_start = StartingPoint.iloc[:, col + 3].values
                #     df_values_point = Point.iloc[:, col+3].values
                #     df_values_endpoint = EndPoint.iloc[:, col+3].values
                #     df_values_data = [df[col_name].values for df in dataframes_list]
                #     ax.scatter([0] * len(df_values_start), df_values_start, color='red', label='Starting Point')
                #     ax.scatter([j+1] * len(df_values_point), df_values_point, color='black', label='Chosen Point',zorder=10)
                #     ax.scatter([LengthOfPath] * len(df_values_endpoint), df_values_endpoint, color='blue', label='Chosen end point',zorder=10)
                #     for i, values in enumerate(df_values_data, start=1):
                #         ax.scatter([i] * len(values), values, label=f'Dataframe {i}')
                
                #     ax.scatter([j+1] * len(ValidPaths), ValidPaths[col_name], color='pink', marker='o', label='Valid Points')
                #     #ax.scatter([j+1] * len(ValidPaths), ValidPaths_2[col_name], color='red', marker='o', label='Valid Points (endpoint)')

                #     ax.axhline(max_bound.iloc[0, col + 3], color='green', linestyle='--', label='Upper Boundary')
                #     ax.set_title(f'Column {col_name}')
                #     ax.legend()
                
                # plt.show()

            
                
                ## Then add the new point, so that it is distinsguishable

                TransitionPath.append(Point)

            # Store the generated path
            #TransitionPath.append(StartingPoint)
            result_df = pd.concat(TransitionPath, ignore_index=True)
            TransitionPath = []
            AllPaths.append(result_df)

    return AllPaths


def GenerateBestPaths(dataframes_list, NrPaths, Decay, StartingPoint):

    if Decay == "no":

        LengthOfPath = len(dataframes_list)

        # Create the desired number of paths
        TransitionPath = []
        TransitionPath.append(StartingPoint)

        AllPaths = []
        while len(AllPaths) < NrPaths:
            print(f"{len(AllPaths)}")
            Point = StartingPoint
            TransitionPath = []
            TransitionPath.append(StartingPoint)
            
            ## Pick a random valid endpoint
            ValidEndPoints = dataframes_list[-1].loc[
                (Point.squeeze().values <= dataframes_list[-1].values).all(axis=1)]
            EndPoint = ValidEndPoints.sample(n=1)

            # Find upper bounds for each variable in the range of all dataframes
            max_values = [df.max().values for df in dataframes_list]

            # Choose the minimum value for each variable among the found maxes
            max_bound_values = np.min(max_values, axis=0)
            min_indices = np.argmin(max_values, axis=0)

            # Create a dataframe with one row and the same column names as the dataframes
            max_bound = pd.DataFrame([max_bound_values], columns=dataframes_list[0].columns)
         
            #NrMiddlePoint = LengthOfPath - 1
            for j in range(LengthOfPath):
                # If the boundary is defined from the j'th dataframe, move the boundary
                
                if any(min_indices==j) and j < LengthOfPath-1:
                    # Just move all the boundaries (boundaries not at j'th, will be recalc not moved)
                    dataframes_list_moved = dataframes_list[j+1:]
                    max_values = [df.max().values for df in dataframes_list_moved]

                    max_bound_values = np.min(max_values, axis=0)
                    min_indices = np.argmin(max_values, axis=0)+j+1

                    max_bound = pd.DataFrame([max_bound_values], columns=dataframes_list[0].columns)

                

                # Find legal paths within the defined range, both using the end point
                # and the moving boundaries
                
                # First find all points below the boundaries
                ValidPaths = dataframes_list[j][
                    (Point.squeeze().values <= dataframes_list[j].values).all(axis=1)
                    & (dataframes_list[j].values <= max_bound.values).all(axis=1)
                ]
                
                # Then exclude those points which have at least one dimension above the endpoint
                ValidPaths_2 = ValidPaths[(ValidPaths.values <= EndPoint.squeeze().values).all(axis=1)]
                if len(ValidPaths_2) == 0:
                    # Pick a new random endpoint
                    NrValidEndPoints = len(ValidEndPoints)
                    print(f'Starting to try and move endpoint')
                    # for i in range(round(NrValidEndPoints * 0.0001)):
                    for i in range(4):
                        # Print a message when a quarter of the dataframes have been processed
                        if (i + 1) % (4 // 4) == 0:
                            print(f"{i + 1}")
                        ## Pick a random valid endpoint
                        ValidEndPoints = dataframes_list[-1].loc[
                            (Point.squeeze().values <= dataframes_list[-1].values).all(axis=1)]
                        EndPoint = ValidEndPoints.sample(n=1)
                        ValidPaths_2 = ValidPaths[(ValidPaths.values <= EndPoint.squeeze().values).all(axis=1)]
                        
                        if len(ValidPaths_2) > LengthOfPath-j:
                            print(f'End point was moved in step {j+1}')
                            break
                       
                
                if len(ValidPaths_2) == 0:
                     print(f'failed after 1000 endpoint moves, returning any completed paths')
                     

                     return AllPaths

                Point = ValidPaths_2.sample(n=1)

                TransitionPath.append(Point)

            # Store the generated path
            #TransitionPath.append(StartingPoint)
            result_df = pd.concat(TransitionPath, ignore_index=True)
            TransitionPath = []
            AllPaths.append(result_df)

    return AllPaths

import pandas as pd
import random

def GenerateRandomPaths_2(dataframes_list, NrPaths):
    """
    Generate random paths by selecting a random row from each dataframe.

    Parameters:
    - dataframes_list: List of dataframes.
    - NrPaths: Number of paths to generate.

    Returns:
    - List of dataframes, each with a random row from each dataframe.
    """
    paths_list = []

    for i in range(NrPaths):
        Path = []
        # Make i'th path
        for j in range(len(dataframes_list)):
            Path.append(dataframes_list[j].sample(n=1))
        
        # random_paths = [df.sample(n=1, random_state=random.seed()).iloc[0] for df in dataframes_list]
        # paths_list.append(pd.DataFrame([random_paths], columns=sum([df.columns.tolist() for df in dataframes_list], [])))
        paths_list.append( pd.concat(Path, ignore_index=True))
        
        # Print a message when a quarter of the paths have been generated
        # if (i + 1) % (NrPaths // 4) == 0:
        #     print(f"{i + 1} paths generated out of {NrPaths}")
    return paths_list
def GenerateRandomPaths(dataframes_list, NrPaths):
    """
    Generate random paths by selecting a random row from each dataframe.

    Parameters:
    - dataframes_list: List of dataframes.
    - NrPaths: Number of paths to generate.

    Returns:
    - List of dataframes, each with a random row from each dataframe.
    """
    paths_list = []

    # Make i'th path
    for i in range(NrPaths):
        Path = [df.sample(n=1).iloc[0] for df in dataframes_list]
        paths_list.append(Path)

        # Print a message when a quarter of the paths have been generated
        if (i + 1) % (NrPaths // 4) == 0:
            print(f"{i + 1} paths generated out of {NrPaths}")

    # Concatenate paths outside the loop
    concatenated_paths = [pd.concat(path, ignore_index=True) for path in paths_list]
    return concatenated_paths

def filter_non_decreasing(dataframes_list):
    """
    Filter dataframes where no columns have decreasing values.

    Parameters:
    - dataframes_list: List of dataframes.

    Returns:
    - List of dataframes where no columns have decreasing values.
    """
    #non_decreasing_paths = [df for df in dataframes_list if all(df[col].is_monotonic_increasing for col in df.columns)]
    
    non_decreasing_paths = []
    indices = []

    for i, df in enumerate(dataframes_list):
        # Check if all columns are non-decreasing
        if all(df[col].is_monotonic_increasing or col in ["GasToEl", "GasToSector","El To Sectors"] for col in df.columns):
        # if all(df[col].is_monotonic_increasing for col in df.columns):
            non_decreasing_paths.append(df)
            indices.append(i)

        # Print a message when a quarter of the dataframes have been processed
        if (i + 1) % (len(dataframes_list) // 4) == 0:
            print(f"{i + 1} dataframes processed out of {len(dataframes_list)}")

    
    return non_decreasing_paths,indices


def PlottingPaths(dataframes_list, Paths, labels=None):
    
    NrPaths = len(Paths)
    
    if labels is None:
       # Default labels if not provided
       labels = [f"Dataframe {i + 1}" for i in range(len(dataframes_list))]
   
        
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"{NrPaths} near optimal paths")
    
    for col in range(3):
        ax = axes[0, col]
        col_name = dataframes_list[0].columns[col]  # Get the column name
        df_values_data = [df[col_name].values for df in dataframes_list]  # Get the solution points of column dimension
        df_min = [min(values) for values in df_values_data]
        df_max = [max(values) for values in df_values_data]
    
        ax.plot(labels[1:], df_min, marker='o', color='black', label='Min')
        ax.plot(labels[1:], df_max, marker='o', color='black', label='Max')
        ax.set_title(f'Column {col_name}')
        ax.set_xlabel('Dataframe Index')
        #ax.set_ylabel('Column Values')
        ax.legend()
       

    
    
    for col in range(3):
        ax = axes[1, col]
        col_name = dataframes_list[0].columns[col + 3]  # Get the column name
        df_values_data = [df[col_name].values for df in dataframes_list]  # Get the solution points of column dimension
        df_min = [min(values) for values in df_values_data]
        df_max = [max(values) for values in df_values_data]
    
        # ax.plot(range(1, len(dataframes_list) + 1), df_min, marker='o', 
        #         color = 'black' ,label='Min')
        # ax.plot(range(1, len(dataframes_list) + 1), df_max, marker='o', 
        #         color = 'black', label='Max')
        
        ax.plot(labels[1:], df_min, marker='o', color='black', label='Min')
        ax.plot(labels[1:], df_max, marker='o', color='black', label='Max')
        
        ax.set_title(f'{col_name}')
        ax.set_xlabel('Dataframe Index')
        #ax.set_ylabel('Column Values')
        ax.legend()
    
  # Adding paths to the plot
    for i, path in enumerate(Paths):
        for col in range(6):
            ax = axes[col // 3, col % 3]  # Use the subplot index for each dimension
            col_name = path.columns[col]  # Get the column name
            values = path[col_name].values
            ax.plot(labels, values, marker='o', label=f'Path {i + 1}')
            ax.set_title(f'{col_name}')
            ax.set_xlabel('% CO2 reduction')
            # if col == 0:  # Set legend only for the first dimension
            #     ax.legend()
           
    
    # # Plot for debugging
    # fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    # fig.suptitle(f"{NrPaths} near optimal paths")
    
    # for col in range(3):
    #     ax = axes[0, col]
    #     col_name = Path[0].columns[col] # Get the column name
    #     df_values_data = [df[col_name].values for df in dataframes_list] # Get the solution points of column dimension
    #     df_min = min(df_values_data)
    #     df_max = max(df_values_data)
        
    # for col in range(3):
    #     ax = axes[1, col]
    #     col_name = Path[0].columns[col+3] # Get the column name
    #     df_values_data = [df[col_name].values for df in dataframes_list] # Get the solution points of column dimension
    #     df_min = min(df_values_data)
    #     df_max = max(df_values_data)
        
    
    
# # Plot for debugging
# fig, axes = plt.subplots(2, 3, figsize=(15, 10))
# fig.suptitle(f"Column Values - Path: {len(AllPaths) + 1}, Step: {j + 1}")

# for col in range(3):
#     ax = axes[0, col]
#     col_name = StartingPoint.columns[col]
#     df_values_start = StartingPoint.iloc[:, col].values
#     df_values_point = Point.iloc[:, col].values
#     df_values_endpoint = EndPoint.iloc[:, col].values
#     df_values_data = [df[col_name].values for df in dataframes_list]
#     ax.scatter([0] * len(df_values_start), df_values_start, color='red', label='Starting Point')
#     ax.scatter([j+1] * len(df_values_point), df_values_point, color='black', label='Chosen Point',zorder=10)
#     ax.scatter([LengthOfPath] * len(df_values_endpoint), df_values_endpoint, color='blue', label='Chosen end point',zorder=10)
    
#     for i, values in enumerate(df_values_data, start=1):
#         ax.scatter([i] * len(values), values, label=f'Dataframe {i}')

    
#     ax.scatter([j+1] * len(ValidPaths), ValidPaths[col_name], color='pink', marker='o', label='Valid Points (moving)')
#     #ax.scatter([j+1] * len(ValidPaths), ValidPaths_2[col_name], color='red', marker='o', label='Valid Points (endpoint)')

#     ax.axhline(max_bound.iloc[0, col], color='green', linestyle='--', label='Upper Boundary')
#     ax.set_title(f'Column {col_name}')
#     ax.legend()

# for col in range(3):
#     ax = axes[1, col]
#     col_name = StartingPoint.columns[col + 3]  # Use columns 4, 5, and 6
#     df_values_start = StartingPoint.iloc[:, col + 3].values
#     df_values_point = Point.iloc[:, col+3].values
#     df_values_endpoint = EndPoint.iloc[:, col+3].values
#     df_values_data = [df[col_name].values for df in dataframes_list]
#     ax.scatter([0] * len(df_values_start), df_values_start, color='red', label='Starting Point')
#     ax.scatter([j+1] * len(df_values_point), df_values_point, color='black', label='Chosen Point',zorder=10)
#     ax.scatter([LengthOfPath] * len(df_values_endpoint), df_values_endpoint, color='blue', label='Chosen end point',zorder=10)
#     for i, values in enumerate(df_values_data, start=1):
#         ax.scatter([i] * len(values), values, label=f'Dataframe {i}')

#     ax.scatter([j+1] * len(ValidPaths), ValidPaths[col_name], color='pink', marker='o', label='Valid Points')
#     #ax.scatter([j+1] * len(ValidPaths), ValidPaths_2[col_name], color='red', marker='o', label='Valid Points (endpoint)')

#     ax.axhline(max_bound.iloc[0, col + 3], color='green', linestyle='--', label='Upper Boundary')
#     ax.set_title(f'Column {col_name}')
#     ax.legend()

# plt.show()
    
def remove_duplicate_dataframes(dfs):
    # Function to convert a dataframe to a hashable representation
    def dataframe_to_tuple(df):
        return tuple(tuple(row) for _, row in df.iterrows())

    # Create a set to store hashable representations
    unique_representations = set()

    # List to store unique dataframes
    unique_dfs = []

    # Iterate over dataframes
    for df in dfs:
        # Convert the dataframe to a hashable representation
        df_tuple = dataframe_to_tuple(df)
        
        # Check if the representation is unique
        if df_tuple not in unique_representations:
            unique_representations.add(df_tuple)
            unique_dfs.append(df)

    return unique_dfs    
