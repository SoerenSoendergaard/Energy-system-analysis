# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 20:54:12 2024

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
import random
import PathToPyPSA

from sklearn.metrics import mean_squared_error

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def find_maximally_different_dataframes(dataframes_list, n=1, exclude_last_n_columns=0):
    """
    Find the n maximally different dataframes from a list of dataframes.

    Parameters:
    - dataframes_list (list): List of pandas DataFrames.
    - n (int): Number of maximally different dataframes to return (default is 1).
    - exclude_last_n_columns (int): Number of last columns to exclude during dissimilarity computation (default is 0).

    Returns:
    - List of tuples containing (index1, dataframe1, index2, dataframe2, dissimilarity_score)
      for the n maximally different dataframes.
    """

    def compute_dissimilarity(df1, df2):
        """
        Compute dissimilarity between two dataframes using mean squared error.

        Parameters:
        - df1, df2 (DataFrame): Two pandas DataFrames to compare.

        Returns:
        - Dissimilarity score between the two dataframes.
        """
        return mean_squared_error(df1.values.flatten(), df2.values.flatten())

    num_dataframes = len(dataframes_list)
    max_dissimilarities = []

    # Compare each pair of dataframes
    for i in range(num_dataframes):
        for j in range(i + 1, num_dataframes):
            df1 = dataframes_list[i].iloc[:, :-exclude_last_n_columns]
            df2 = dataframes_list[j].iloc[:, :-exclude_last_n_columns]
            dissimilarity = compute_dissimilarity(df1, df2)
            max_dissimilarities.append((i, dataframes_list[i], j, dataframes_list[j], dissimilarity))

    # Sort the dissimilarities list in descending order of dissimilarity scores
    max_dissimilarities.sort(key=lambda x: x[4], reverse=True)

    # Return the n maximally different dataframes
    return max_dissimilarities[:n]



def AssignYearlyCosts(SourceFolder,Reductions,Paths, MAA_variables):
    
    Slack = [0.02, 0.04, 0.06 ,0.08]
    
    # Check that there are as many columns in paths, as there are reductions
    if len(Reductions) != Paths[0].shape[0]-1:
        print(f'Number of reductions must match number of dataframe rows')
    
    
    for df in Paths:
        df['Cost'] = None
    
    # Assign costs to the i'th period
    
    
    # Get all combinations
    [combination_keys,num_combinations,combinations_dictionary] = GetCombinations(MAA_variables)
    
    for i in range(len(Reductions)):
        print(f'Assigning cost to {Reductions[i]} - case')
        # Get optimal cost
        folder_name = f'{SourceFolder}_{Reductions[i]}'
        file_name = f'OptimalCost'
        file_path = rf'C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\{folder_name}\{file_name}'
        
        obj = np.loadtxt(file_path, delimiter=',', dtype=float)
        
        # If nothing else, then a point will have the maximum value
        for df in Paths:
            df.loc[i+1, 'Cost'] = obj*1.1
            
        
        # Make a temporary dataframe with all the point of this reduction
        rows_df = pd.concat([df.iloc[[i+1]] for df in Paths], ignore_index=True)

        
        # Assign costs be going through all 2d slice combinations of the contour plot
        for j in range(num_combinations):
            Combination = combinations_dictionary[combination_keys[j]]
            chosen_variables = [Combination['Names'][0],Combination['Names'][1]]
            # Create simplexes and check 
            
            
            variable_samples = rows_df[chosen_variables]
            
            assigned_points = set()
            
            for k in range(len(Slack)):
               
                # Load verticies
                slack = Slack[k]
                folder_name = f'{SourceFolder}_MultiSlack_{Reductions[i]}\{slack}'
                file_name = f'Verticies'
                file_path = rf'C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\{folder_name}\{file_name}'
                verticies = np.loadtxt(file_path, delimiter=',', dtype=int)
                
                All_verticies = pd.DataFrame(verticies,
                                            columns = MAA_variables)
                
                
                variable_verticies = pd.DataFrame(All_verticies,
                                            columns = chosen_variables)
                
                Points = variable_verticies.values
                hull = ConvexHull(Points)
                path = Path(variable_verticies.values[hull.vertices])

                # Check which points are inside the convex hull
                points_inside = path.contains_points(variable_samples[[chosen_variables[0], chosen_variables[1]]].values)
                points_inside = np.logical_and(points_inside, ~np.isin(range(len(points_inside)), list(assigned_points)))
                # Update the set of assigned points
                assigned_points.update(np.where(points_inside)[0])
                
                # Assign a unique cost to each contour
                cost = obj*(1+slack)  # You can adjust the cost assignment as needed
                
                if j==0:
                    # Update the 'Cost' column for points inside the convex hull
                    rows_df.loc[points_inside, 'Cost'] = cost
                else: 
                    
                    tmpCost = rows_df.copy()
                    tmpCost.loc[points_inside, 'Cost'] = cost
                    # Update 'Cost' column only if the new cost is higher than the existing cost
                    
                #     # Plotting for debugging
                
                # for simplex in hull.simplices:
                #     l0, = plt.plot(variable_verticies.values[simplex, 0], variable_verticies.values[simplex, 1], 'k-', 
                #             color = 'k',
                #             linewidth = 2, zorder = 1)
                    
                    
           
            if j > 0:
            # Update 'Cost' column in Sample_With_Cost only for the points where tmpCost has higher costs
                higher_cost_points = tmpCost['Cost'] > rows_df['Cost']
                rows_df.loc[higher_cost_points, 'Cost'] = tmpCost.loc[higher_cost_points, 'Cost']
                
        for z, df in enumerate(Paths):
            cost = rows_df.loc[z, 'Cost']  # Get the cost for the current dataframe
            df.loc[i+1, 'Cost'] = cost  # Assign the cost to the 'Cost' column of the first row of the current dataframe

        
    # plt.scatter(x=variable_samples[f'{chosen_variables[0]}'], 
    #             y=variable_samples[f'{chosen_variables[1]}'],
    #             c = rows_df['Cost'], cmap='viridis')
        
    
    # plt.xlabel(f'{chosen_variables[0]} [MW]')
    # plt.ylabel(f'{chosen_variables[1]} [MW]')
    # plt.show()
    
    return Paths
        
        
def GetCombinations(MAA_variables):

    technologies = list(MAA_variables.keys())
    
    # Create a new dictionary with 21 sub-dictionaries
    combinations_dictionary = {}
    
    ###########
    combinations_dictionary = {}
    
    # Iterate through all combinations
    for i in range(len(technologies)):
        for j in range(i + 1, len(technologies)):
            combination_key = f"{technologies[i]}_{technologies[j]}"
            combined_names = [technologies[i], technologies[j]]
    
            combination_values = {
                'Names': combined_names,
                'types': [MAA_variables[tech][0] for tech in combined_names],
                'carrier': [component for tech in combined_names for component in MAA_variables[tech][1]],
                'exploring': [MAA_variables[tech][2] for tech in combined_names]
            }
    
            combinations_dictionary[combination_key] = combination_values
    #########
    
    combination_keys = list(combinations_dictionary.keys())
    
    num_combinations = len(combinations_dictionary)   
    
    return combination_keys,num_combinations,combinations_dictionary
        
        
def TotalCost(PathsWithCost,Reductions_moderate,Years_moderate,Reductions_RapidStart,Years_RapidStart,
              Reductions_RapidEnd,Years_RapidEnd):

    # Go through the 3 scenarios: 1 = RapidStart, 2 = moderate, 3 = RapidEnd
    
    # Normalize according to greenfield
    TransitionPath = [0.5,0.6,0.7,0.8,0.9,0.98]
    Years = [2027,2028.5,2030,2035,2040,2045]

    OptimumCostArray = []
    for i in range(len(TransitionPath)):
        file_path = rf"C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\FlexDNK_Ons_Offs_Sol_Biom_Bat_X_{TransitionPath[i]}\OptimalCost"
        obj = np.loadtxt(file_path, delimiter=',', dtype=float)
        OptimumCostArray.append(obj)

    Greenfield_Trans_Cost_moderate = np.trapz(OptimumCostArray, x=Years)
    Greenfield_Trans_Cost_RapidStart = np.trapz(OptimumCostArray, x=Years_RapidStart)
    Greenfield_Trans_Cost_RapidEnd = np.trapz(OptimumCostArray, x=Years_RapidEnd)
    
    for z in range(3):
        
        if z == 0:
            Years = Years_RapidStart
            Reductions = Reductions_RapidStart
            
            
        if z == 1:
            Years = Years_moderate
            Reductions = Reductions_moderate
            
        if z == 2:
            Years = Years_RapidEnd
            Reductions = Reductions_RapidEnd
    
    # First plot the transition function
    
          
        TotalCostList = []
        # Calculate the cost of each path
        for i in range(len(PathsWithCost)):
            
            YearlyCosts = PathsWithCost[i]['Cost'][1:] # The first entry is the starting year for some reason
            
            TotalCost = np.trapz(YearlyCosts, x=Years)
            TotalCostList.append(TotalCost)
            
            if z==0:
                TotalCostList_RapidStart = TotalCostList/Greenfield_Trans_Cost_RapidStart
                
            
            if z == 1:
                TotalCostList_moderate = TotalCostList/Greenfield_Trans_Cost_moderate
                
                
            if z == 2:
                TotalCostList_RapidEnd = TotalCostList/Greenfield_Trans_Cost_RapidEnd
                
    
            
    # Make subplot with transition strategy and histogram
    
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    
    axes[0, 0].plot(Years_RapidStart,Reductions_RapidStart,linewidth = 3)
    axes[0, 0].set_title("Rapid start",fontsize = 16)
    axes[0, 0].set_ylabel("CO2 - reduction",fontsize = 14)
    axes[0, 1].plot(Years_moderate,Reductions_moderate,linewidth = 3)
    axes[0, 1].set_title("Moderate timing",fontsize = 16)
    axes[0, 2].plot(Years_RapidEnd,Reductions_RapidEnd,linewidth = 3)
    axes[0, 2].set_title("Rapid end",fontsize = 16)
    
    axes[1, 0].hist(TotalCostList_RapidStart, bins=30, color='skyblue', edgecolor='black')
    axes[1, 0].set_ylabel("Occurences",fontsize = 14)
    axes[1, 1].hist(TotalCostList_moderate, bins=30, color='skyblue', edgecolor='black')
    axes[1, 1].set_xlabel("Transition cost normalized to greenfield transition",fontsize = 14)
    axes[1, 2].hist(TotalCostList_RapidEnd, bins=30, color='skyblue', edgecolor='black')
    
    for i, total_cost_list in enumerate([TotalCostList_RapidStart, TotalCostList_moderate, TotalCostList_RapidEnd]):
        axes[1, i].hist(total_cost_list, bins=30, color='skyblue', edgecolor='black')
        
        max_value = np.max(total_cost_list)
        min_value = np.min(total_cost_list)
        difference = (max_value - min_value) / min_value * 100
        difference = round(difference, 2)
        
        textbox_text = f'max/min difference:\n {difference}%'
        axes[1, i].text(0.65, 0.95, textbox_text, transform=axes[1, i].transAxes, fontsize=14, verticalalignment='top', horizontalalignment='right')
    
    
    plt.tight_layout()

    # Show the plot
    plt.show()
    
    return TotalCostList_RapidStart, TotalCostList_moderate, TotalCostList_RapidEnd  


def PlotYearlyCosts(PathsWithCost, Years):
    """Plot the yearly costs for each path in PathsWithCost."""
    
    # Create a new figure
    # plt.figure(figsize=(10, 6))
    
    # Iterate over each dataframe in PathsWithCost
    for i, df in enumerate(PathsWithCost):
        # Extract the data from the 7th column (assuming it starts from index 0)
        path_data = df.iloc[:, 6]  # 7th column
        
        # Plot the data as a blue line plot with dots
        plt.plot(Years, path_data, marker='o', color='blue' )
    
    # Add labels and title
    plt.xlabel('Year')
    plt.ylabel('Yearly Cost')
    plt.title('Yearly Costs for Each Path')
    
    # Add legend
    plt.legend()
    
    # Show plot
    # plt.show()
    

def Discount(n,r):
    """Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6"""

    if r > 0:
        return r/(1. - 1./(1.+r)**n)
    else:
        return 1/n      
        
def SDR(PathsWithCost,Years,DiscountRate):
    
    DiscountedYearlyCosts = []
    
    indexes = [index for index, year in enumerate(Years) if year <= 2050]
    vector_discount_rates = [DiscountRate[0] if index in indexes else DiscountRate[1] for index in range(len(Years))]
    
    # Define the starting year
    starting_year = 2020
    
    # Calculate the number of years since the starting year for each entry
    n = [year - starting_year for year in Years]
    
    # Discount each path
    for i in range(len(PathsWithCost)):
       YearlyCosts = PathsWithCost[i]['Cost'][1:]
       
       # Get the index
       # DiscountedYearlyCosts_tmp = [Discount(n_val, discount_rate) * cost for n_val, 
       #                              discount_rate, cost in zip(n, vector_discount_rates, YearlyCosts)]
       
       DiscountedYearlyCosts_tmp = [cost/((1+discount_rate)**n_val) for n_val, 
                                    discount_rate, cost in zip(n, vector_discount_rates, YearlyCosts)]
        
       DiscountedYearlyCosts.append(DiscountedYearlyCosts_tmp)
       
    return DiscountedYearlyCosts
 

       
def GetTopAndBottomPaths(NrTopAndBottomPaths,TotalCost_moderate,PathsWithCost):
    
    # Get indices of sorted costs
    sorted_indices = np.argsort(TotalCost_moderate)
    
    # Get indices of top and bottom paths
    top_indices = sorted_indices[-NrTopAndBottomPaths:]
    bottom_indices = sorted_indices[:NrTopAndBottomPaths]
    
    # Get corresponding dataframes
    top_paths = [PathsWithCost[i] for i in top_indices]
    bottom_paths = [PathsWithCost[i] for i in bottom_indices]
    
    # Get corresponding costs
    top_costs = [TotalCost_moderate[i] for i in top_indices]
    bottom_costs = [TotalCost_moderate[i] for i in bottom_indices]
    
    # Get corresponding YearlyCosts
    
    return top_paths, bottom_paths, top_costs, bottom_costs, top_indices, bottom_indices

def Plotting_NOS_Bounadaries(dataframes_list_list, labels=None,transparency=1):
    
    
    
    # if labels is None:
    #    # Default labels if not provided
    #    labels = [f"Dataframe {i + 1}" for i in range(len(dataframes_list))]
       

   
        
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f" Near optimal solution spaces")
    
    # Create a ScalarMappable object for the colorbar
    colors = ['black', 'red', 'blue', 'green']
    slacks = ["2% slack","4% slack"]
    for z, dataframes_list in enumerate(dataframes_list_list):
        color = colors[z % len(colors)]
        dataframes_list = dataframes_list_list[z]
        
    
        for col in range(3):
            ax = axes[0, col]
            col_name = dataframes_list[0].columns[col]  # Get the column name
            df_values_data = [df[col_name].values for df in dataframes_list]  # Get the solution points of column dimension
            df_min = [min(values) for values in df_values_data]
            df_max = [max(values) for values in df_values_data]
        
            ax.plot(labels[1:], df_min, marker='o', color=color)
            ax.plot(labels[1:], df_max, marker='o', color=color)
            ax.set_title(f'{col_name}')
            # ax.set_xlabel('CO2 reduction [%]')
            ax.set_ylabel('Capacity [MW]')
            ax.tick_params(axis='x', rotation=45)
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
            
            if col == 2:
                ax.plot(labels[1:], df_min, marker='o', color=color, label=rf'Min and max of {slacks[z]}')
                ax.plot(labels[1:], df_max, marker='o', color=color)
            else:
                ax.plot(labels[1:], df_min, marker='o', color=color)
                ax.plot(labels[1:], df_max, marker='o', color=color)
                
            ax.set_title(f'{col_name}')
            ax.set_xlabel('Years')
            ax.set_ylabel('Capacity [MW]')
            ax.tick_params(axis='x', rotation=45)
            ax.legend()    
    plt.tight_layout()
    plt.show()
        
        
        
        
        
        

def PlottingPaths(dataframes_list, Paths, labels=None, total_costs=None,transparency=1):
    
    NrPaths = len(Paths)
    
    if labels is None:
       # Default labels if not provided
       labels = [f"Dataframe {i + 1}" for i in range(len(dataframes_list))]
       
    if total_costs is None:
     total_costs = [0] * len(Paths)  # Default total costs if not provided

   
        
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"{NrPaths} near optimal paths")
    
    # Create a ScalarMappable object for the colorbar
    sm = plt.cm.ScalarMappable(cmap='winter', norm=plt.Normalize(vmin=min(total_costs), vmax=max(total_costs)))
    sm.set_array([])  # Dummy array required for the ScalarMappable
    
    for col in range(3):
        ax = axes[0, col]
        col_name = dataframes_list[0].columns[col]  # Get the column name
        df_values_data = [df[col_name].values for df in dataframes_list]  # Get the solution points of column dimension
        df_min = [min(values) for values in df_values_data]
        df_max = [max(values) for values in df_values_data]
    
        ax.plot(labels[1:], df_min, marker='o', color='black', label='Min')
        ax.plot(labels[1:], df_max, marker='o', color='black', label='Max')
        ax.set_title(f'{col_name}')
        ax.set_xlabel('CO2 reduction [%]')
        ax.set_ylabel('Capacity [MW]')
        ax.tick_params(axis='x', rotation=45)
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
    for i, (path, cost) in enumerate(zip(Paths, total_costs)):
        for col in range(6):
            ax = axes[col // 3, col % 3]  # Use the subplot index for each dimension
            col_name = path.columns[col]  # Get the column name
            values = path[col_name].values
            ax.plot(labels, values, marker='o', label=f'Path {i + 1}', alpha=transparency,color=plt.cm.winter(np.interp(cost, [min(total_costs), max(total_costs)], [0, 1])))
            ax.set_title(f'{col_name}')
            ax.set_xlabel('CO2 reduction [%]')
            ax.set_ylabel('Capacity [MW]')
            ax.tick_params(axis='x', rotation=45)
            # if col == 0:  # Set legend only for the first dimension
            #     ax.legend()
    # # Add colorbar

    plt.tight_layout()
    cbar_ax = fig.add_axes([0.2, -0.02, 0.6, 0.02])  # Position and size of the colorbar
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Total Cost (Normalized relative to greenfield optimum transition)')
          

def RunningTotalCost(YearlyCosts, Years, StartYear, Plotting):
    AllRunningCosts = []
    Years.insert(0, StartYear)
    for j in range(len(YearlyCosts)):
        running_total_cost = []
        costs = YearlyCosts[j].fillna(0)
        costs_list = costs.tolist()
        for i in range(len(YearlyCosts[j])):
            area_under_curve = np.trapz(costs_list[:i+1], x= Years[:i+1])
            running_total_cost.append(area_under_curve)
        
        AllRunningCosts.append(running_total_cost)
        
    if Plotting == 'yes':
        # Do the plot
        2+2
    
    return AllRunningCosts



def Max_Tech_Paths(PathsWithCost,Technologies):
    
    # Read all paths and for each desired technology:
    # - Find all paths where this technology has the majority of installed capacity
    # - If multiple paths exists pick the cheapest
    # - If none exists just pick the path with the most end capacity of the given tech
    
    InterestingPaths = []
    for tech in Technologies:
        
        # Get all end configurations of this tech
        # 
    # Initialize an empty list to store the values
        Tech_values = []
        
        # Iterate over each DataFrame in the list
        for df in PathsWithCost:
            # Extract the last value of the 'OnWind' column and append it to the list
            Tech_values.append(df[tech].iloc[-1])
        
        # Create a DataFrame from the list of values
        result_df = pd.DataFrame(Tech_values, columns=['Tech_Last_Value'])
        
        MaxPathId = result_df.idxmax()
        
        InterestingPaths.append(PathsWithCost[MaxPathId[0]])
    return InterestingPaths
  

def SavePathData(PathCapacities,PathTotalCost,PathTotalCO2,PathPartOfEURPerWatt,
                                OffWind_Curtail,OnWind_Curtail,Solar_Curtail,FilePath):
    
    if not os.path.exists(FilePath):
        os.makedirs(FilePath)
    
    # Capacities
    FileName = 'PathCapacities'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}.h5"

    # Save the list of dataframes to an HDF5 file
    with pd.HDFStore(FullFilePath) as store:
        for i, df in enumerate(PathCapacities):
            store[f'df_{i}'] = df

    # Total cost
    # FilePath = rf"C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\PathResults"
    FileName = 'PathTotalCost.npy'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}"

    np.save(FullFilePath, PathTotalCost)


            
    # PathTotalCO2
    # FilePath = rf"C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\PathResults"
    FileName = 'PathTotalCO2.npy'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}"
    
    np.save(FullFilePath, PathTotalCO2)
    
    # Curtailment data
    
    # OnWind
    FileName = 'OnWind_Curtail.npy'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}"
    
    np.save(FullFilePath, OnWind_Curtail)
    
    
    # OffWind
    FileName = 'OffWind_Curtail.npy'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}"
    
    np.save(FullFilePath, OffWind_Curtail)
    
    # Solar
    FileName = 'Solar_Curtail.npy'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}"
    
    np.save(FullFilePath, Solar_Curtail)
    
    

            
    # PathPartOfEURPerWatt
    # FilePath = rf"C:\Users\soere\OneDrive\UNI\Kandidat\Semester 4\Results\Green\PathResults"
    FileName = 'PathPartOfEURPerWatt'
    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}.h5"

    # Save the list of dataframes to an HDF5 file
    with pd.HDFStore(FullFilePath) as store:
        for i, df in enumerate(PathPartOfEURPerWatt):
            store[f'df_{i}'] = df


def LoadPathData(PathResultsFolder):
    
    ## Load the databases
    FilePath = PathResultsFolder
    FileName = "PathCapacities"

    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}.h5"

    # Load the list of dataframes from the HDF5 file
    PathCapacities = []
    with pd.HDFStore(FullFilePath) as store:
        for key in store.keys():
            PathCapacities.append(store[key])
    
    
    FilePath = PathResultsFolder
    FileName = "PathPartOfEURPerWatt"

    # Combine FilePath and FileName to create the full file path
    FullFilePath = rf"{FilePath}\{FileName}.h5"

    # Load the list of dataframes from the HDF5 file
    PathPartOfEURPerWatt = []
    with pd.HDFStore(FullFilePath) as store:
        for key in store.keys():
            PathPartOfEURPerWatt.append(store[key])
            
    
    ## Load the numpy arrays
    FileName = "PathTotalCost"
    FullFilePath = rf"{FilePath}\{FileName}.npy"
    PathTotalCost =  np.load(FullFilePath, allow_pickle=True)
    # Reshaping
    PathTotalCost = PathTotalCost.reshape((PathTotalCost.shape[0], PathTotalCost.shape[2]))


    FileName = "PathTotalCO2"
    FullFilePath = rf"{FilePath}\{FileName}.npy"
    PathTotalCO2 =  np.load(FullFilePath, allow_pickle=True)
    # Reshaping
    PathTotalCO2 = PathTotalCO2.reshape((PathTotalCO2.shape[0], PathTotalCO2.shape[2]))
    
    FileName = "OnWind_Curtail"
    FullFilePath = rf"{FilePath}\{FileName}.npy"
    PathOnWind_Curtail =  np.load(FullFilePath, allow_pickle=True)
    # Reshaping
    PathOnWind_Curtail = PathOnWind_Curtail.reshape((PathOnWind_Curtail.shape[0], PathOnWind_Curtail.shape[2]))
    
    FileName = "OffWind_Curtail"
    FullFilePath = rf"{FilePath}\{FileName}.npy"
    PathOffWind_Curtail =  np.load(FullFilePath, allow_pickle=True)
    # Reshaping
    PathOffWind_Curtail = PathOffWind_Curtail.reshape((PathOffWind_Curtail.shape[0], PathOffWind_Curtail.shape[2]))
    
    
    FileName = "Solar_Curtail"
    FullFilePath = rf"{FilePath}\{FileName}.npy"
    PathSolar_Curtail =  np.load(FullFilePath, allow_pickle=True)
    # Reshaping
    PathSolar_Curtail = PathOffWind_Curtail.reshape((PathSolar_Curtail.shape[0], PathSolar_Curtail.shape[2]))
    
    
    
    return PathCapacities, PathTotalCost, PathTotalCO2, PathPartOfEURPerWatt, PathOnWind_Curtail, PathOffWind_Curtail,PathSolar_Curtail
    
    
  
def StackedPathPlot(Reduction, Capacities, TotalCost, TotalCO2, PartOfEURPerWatt):
    
    fig,ax = plt.subplots()
    plt.stackplot(Reduction,Capacities.iloc[:,0],
                  Capacities.iloc[:,1],
                  Capacities.iloc[:,2],
                  Capacities.iloc[:,3],
                  Capacities.iloc[:,4],
                  Capacities.iloc[:,5],
                  # Capacities.iloc[:,6],
                  #Capacities.iloc[:,7],
                  Capacities.iloc[:,8],

                  labels = ["Onshore","Offshore","Solar","Biomass (El)","Gas (El)","Gas (sectors)",
                            #"battery storage",
                            # "X storage",
                            "El-sector-link"],
                  colors = ["lightblue","darkblue","yellow","green","black","brown","purple"])
    plt.legend(fancybox=True, loc='upper center', bbox_to_anchor=(0.5, -0.2),shadow=True, ncol=4)
    plt.xlabel("CO2 reduction",fontsize=14)
    plt.ylabel("Installed capacity [MW]",fontsize=14)
    plt.title("Optimal capcity as a function of CO2 reduction \n Greenfield")
    
    # Create a twin Axes for the secondary y-axis
    ax2 = ax.twinx()
    ax2.scatter(Reduction, TotalCost, color='red', label='Total Cost')
    ax2.set_ylabel("Annual system cost (eur)", fontsize=14, color='red')
    plt.show()
    
    
    fig,ax = plt.subplots()
    plt.stackplot(Reduction,PartOfEURPerWatt.iloc[:,0],
                  PartOfEURPerWatt.iloc[:,1],
                  PartOfEURPerWatt.iloc[:,2],
                  PartOfEURPerWatt.iloc[:,3],
                  PartOfEURPerWatt.iloc[:,4],
                  PartOfEURPerWatt.iloc[:,5],
                  PartOfEURPerWatt.iloc[:,6],
                  PartOfEURPerWatt.iloc[:,7],
                  # PartOfEURPerWatt.iloc[:,8],
                  
                  labels = ["Onshore","Offshore","Solar","Biomass (El)","Gas (El)","Gas (sectors)",
                            "battery storage",
                            "X storage",
                            # "El-sector-link"])
                            ],
                  colors = ["lightblue","darkblue","yellow","green","black","brown","pink","grey"])
    plt.legend(fancybox=True, loc='upper center', bbox_to_anchor=(0.5, -0.2),shadow=True, ncol=4)
    plt.xlabel("CO2 reduction",fontsize=14)
    plt.ylabel("Part of energy price [EUR/MWh]",fontsize=14)
    plt.title("System cost as function of CO2 reduction \n Greenfield")
    
    
    Co2_1990 = 51.1*1e6
    CO2_Targets = [Co2_1990*(1-0.40),(1-0.50)*Co2_1990,(1-0.60)*Co2_1990,
                   (1-0.70)*Co2_1990,(1-0.80)*Co2_1990,(1-0.90)*Co2_1990,(1-0.98)*Co2_1990]
    
    # Create a twin Axes for the secondary y-axis
    ax2 = ax.twinx()
    ax2.scatter(Reduction, TotalCO2, color='black', edgecolor='red',s=70,
                label='Yearly CO2 emmissions',marker = 's')
    ax2.scatter(Reduction, CO2_Targets, color='red', edgecolor='black',s=70,
                label='Yearly CO2 targets',marker = '^')
    ax2.set_ylabel("Annual CO2 emissions (Tonnes)", fontsize=14, color='red')
    
    ax2.legend(loc='upper center')
    
    # ax2.legend('CO2 emmisions','CO2 Target')
    
    plt.show()
    
    
def FullySolvePath(Paths,Order,Costs,NrPaths,SavePath,AppendFolder = None):
    
    if Order == 'Top':
        
        ## Actually i broke this after making crucial plots, so i have not had 
        # enough time and incentive ti fix it.
        # Sort paths
        sorted_indices = np.argsort(Costs)
        
        # Get idx of top paths
        # top_indices = sorted_indices[-NrPaths:]
        # PathsToSolve = [Paths[i] for i in top_indices]
        # Actually the above is not the right way. I need to just pick them out
        # one by one so that i can get NrPaths even though some fail.
        
    if Order == 'Bottom':
        
        sorted_indices = np.argsort(Costs)
        
        # bottom_indices = sorted_indices[:NrPaths]
        # PathsToSolve = [Paths[i] for i in bottom_indices]
        
        
    if Order == 'Random':

        # random_paths = random.sample(Paths, NrPaths)
        indices = np.arange(len(Costs))  # Generate indices
        np.random.shuffle(indices)  # Shuffle the indices randomly
        sorted_indices = indices
        
    if AppendFolder is None:
        
        PathTransitionCost, PathCapacities, PathTotalCost, PathTotalCO2, PathPartOfEURPerWatt, PathOffWind_Curtail, PathOnWind_Curtail, PathSolar_Curtail = [],[], [], [], [], [], [], []
        
    # Solve the paths
    
    i = -1
    while len(PathTransitionCost)<NrPaths:
        i = i+1
        
        Path = Paths[sorted_indices[i]]
        
        # if Order == 'Top':
        #     Path = Paths[sorted_indices[i]]
            
        # if Order == 'Bottom':
        #     Path = Paths[sorted_indices[i]]
            
        
        # if Order == 'Random':
        #     Path = Paths[sorted_indices[i]]
        
        Capacities, TotalCost, TotalCO2, PartOfEURPerWatt,OffWind_Curtail,OnWind_Curtail,Solar_Curtail = PathToPyPSA.Run_PyPSA_For_Paths(Path)
        
        
        if AppendFolder is None:
            
            if len(Capacities) != 0:
        
                # Save the results
                PathCapacities.append(Capacities)
                TotalCost_reshaped = TotalCost.reshape(1, -1)
                PathTotalCost.append(TotalCost_reshaped)
    
                TotalCO2_reshaped = TotalCO2.reshape(1, -1)
                PathTotalCO2.append(TotalCO2_reshaped)
                PathPartOfEURPerWatt.append(PartOfEURPerWatt)
                
                TransitionCost = np.trapz(TotalCost, x=[2025,2027,2028.5,2030,2035,2040,2045])
                PathTransitionCost = np.append(PathTransitionCost,TransitionCost)
                
                OffWind_Curtail_reshaped = OffWind_Curtail.reshape(1, -1)
                PathOffWind_Curtail.append(OffWind_Curtail_reshaped)
                
                OnWind_Curtail_reshaped = OnWind_Curtail.reshape(1, -1)
                PathOnWind_Curtail.append(OnWind_Curtail_reshaped)
                
                Solar_Curtail_reshaped = Solar_Curtail.reshape(1, -1)
                PathSolar_Curtail.append(Solar_Curtail_reshaped)
                
                ## Really save the results

                SavePathData(PathCapacities,PathTotalCost,PathTotalCO2,PathPartOfEURPerWatt,
                                                PathOffWind_Curtail,PathOnWind_Curtail,PathSolar_Curtail,SavePath)
                
        
            
    
    
def Exogen(Capacities,TotalCost,PathPartOfEURPerWatt,Years):
    
    ## Correct prices by exogen changes
    # Years = [2025, 2027, 2028.5, 2030, 2035, 2040, 2045]
    ## New prices
    Discount_rate = 0.07
    # Onshore
    Lifetime = [28.5, 29, 29.5, 30, 30, 30,30]
    investment = [1077.17*1000,1060526, 1048043,1035.56*1000,1006.56*1000, 977.57*1000,970.32*1000 ]
    FOM = [1.23,1,23,1.22,1.22, 1.2, 1.19, 1.18]
    VOM = [1.42,1.4 ,1.4 ,1.35, 1.3, 1.24, 1.23]
    
    DemandHours = 214813939.1300003
    
    
    
    
    
    Exo_DiscountedYearlyCost_arr = np.zeros([700,7])
    for i in range(len(Capacities)):
        
        # Back track the Network.generators_t.p
        
    
        # Capex = np.zeros(len(Years))
        for j in range(len(Years)):
        
            Capex = annuity(Lifetime[j],Discount_rate)*investment[j]*(1+(FOM[j]/100))
            Marg = VOM[j]
            
            PriceOfTech = Capex*Capacities[i]["OnWind"][j] + sum(Network.generators.marginal_cost[j]*Network.generators_t.p.iloc[:,j])
            
            PriceOfTech
            
            PartOfEURPerWatt.iloc[i,j] = PriceOfTech/DemandHours
        
        
    a=b/c
    
    
    
    # interpolation template
    Lifetime_25 = 0.25
    Lifetime_30 = 0.34
    
    # Years for interpolation
    years = [2025, 2030]
    lifetimes = [Lifetime_25, Lifetime_30]
    
    # Years for which to interpolate
    years_to_interpolate = [2027, 2028.5]
    
    # Perform linear interpolation
    interpolated_lifetimes = np.interp(years_to_interpolate, years, lifetimes)
    
    
    

def annuity(n,r):
    """Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6"""

    if r > 0:
        return r/(1. - 1./(1.+r)**n)
    else:
        return 1/n
    
    
def StackedExoCompare(Reduction, Capacities, PathPartOfEURPerWatt, Pareto):

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharey=True)  # Create subplot with 2 rows and 3 columns
    
    # Iterate over each subplot
    for i, ax in enumerate(axes.flatten()):
        # Plot stackplot
        # if i<5:
        ax.stackplot(Reduction, Capacities[Pareto[i]].iloc[:, 0], Capacities[Pareto[i]].iloc[:, 1], Capacities[Pareto[i]].iloc[:, 2],
                     Capacities[Pareto[i]].iloc[:, 3], Capacities[Pareto[i]].iloc[:, 4], Capacities[Pareto[i]].iloc[:, 5],
                     Capacities[Pareto[i]].iloc[:, 6],Capacities[Pareto[i]].iloc[:, 7], Capacities[Pareto[i]].iloc[:, 8],  # Include only required columns
                     labels=["Onshore", "Offshore", "Solar", "Biomass (El)", "Gas (El)", "Gas (sectors)",
                         "Battery storage","X-storage","El-sector-link"],
                     colors=["lightblue", "darkblue", "yellow", "green", "black", "brown","pink","grey", "purple"])
        ax.set_xlabel("Years", fontsize=14)
        ax.set_ylabel("Installed capacity [MW]", fontsize=14)
        ax.tick_params(axis='x', rotation=45)
        # ax.set_title(f'p{i + 1}', fontsize=14)
        ax.set_title(f'random {i + 1}', fontsize=14)# Set subplot title
        if i==5:
            ax.set_title(f'Myopic greenfield', fontsize=14)# Set subplot title

    
    
    
    # Get handles and labels for all subplots
    handles, labels = axes[0, 0].get_legend_handles_labels()
    
    # Create a single legend for all subplots and place it outside the plot
    fig.legend(handles, labels, fancybox=True, loc='upper center', 
               bbox_to_anchor=(0.5, -0.05), fontsize='large',shadow=True, ncol=4)
    
    plt.tight_layout()  # Adjust layout to prevent overlap
    plt.savefig("StackedCap.pdf", format="pdf", bbox_inches="tight")
    plt.show()
    
    
    # ## Do it again but different
    # fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharey=True)  # Create subplot with 2 rows and 3 columns
    
    # # Iterate over each subplot
    # for i, ax in enumerate(axes.flatten()):
    #     # Plot stackplot
    #     # if i<5:
    #     ax.stackplot(Reduction, PathPartOfEURPerWatt[Pareto[i]].iloc[:, 0], PathPartOfEURPerWatt[Pareto[i]].iloc[:, 1], PathPartOfEURPerWatt[Pareto[i]].iloc[:, 2],
    #                  PathPartOfEURPerWatt[Pareto[i]].iloc[:, 3], PathPartOfEURPerWatt[Pareto[i]].iloc[:, 4], PathPartOfEURPerWatt[Pareto[i]].iloc[:, 5], 
    #                  PathPartOfEURPerWatt[Pareto[i]].iloc[:, 6],PathPartOfEURPerWatt[Pareto[i]].iloc[:, 7],  # Include only required columns
    #                  labels = ["Onshore","Offshore","Solar","Biomass (El)","Gas (El)","Gas (sectors)",
    #                            "battery storage",
    #                            "X storage",
    #                            # "El-sector-link"])
    #                            ],
    #                  colors = ["lightblue","darkblue","yellow","green","black","brown","pink","grey"])
    #     ax.set_xlabel("Years", fontsize=14)
    #     ax.set_ylabel("Part of energy price [EUR/MWh]", fontsize=14)
    #     ax.tick_params(axis='x', rotation=45)
    #     ax.set_title(f'p{i + 1}', fontsize=14)  # Set subplot title

    
    
    
    # # Get handles and labels for all subplots
    # handles, labels = axes[0, 0].get_legend_handles_labels()
    
    # # Create a single legend for all subplots and place it outside the plot
    # fig.legend(handles, labels, fancybox=True, loc='upper center', 
    #            bbox_to_anchor=(0.5, -0.05), fontsize='large',shadow=True, ncol=4)
    
    # plt.tight_layout()  # Adjust layout to prevent overlap
    # plt.savefig("PartOfPrice.pdf", format="pdf", bbox_inches="tight")
    # plt.show()
    
    


    
    

    
