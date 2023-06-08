# -*- coding: UTF-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import timedelta
import sys

csv_filename = sys.argv[1]
flux_tower_csv = sys.argv[2] if len(sys.argv) > 2 else 'NA'

#csv_file = 'SWB2_variable_values__col_236__row_373__x_386474__y_1289444.csv'
bits = csv_filename.replace('___','_').replace('__','_').replace('.','_').split('_')
output_file_prefix = csv_filename.split('.')[0]
locator_text = 'x: ' + bits[8] + ', y: ' + bits[10]
df = pd.read_csv(csv_filename)

flux_df = None

if flux_tower_csv != 'NA':
    flux_df_in = pd.read_csv(flux_tower_csv,header=0,comment='#')
    flux_df = munge_flux_tower_data(flux_df_in)

# eliminate unwanted spaces before and after items in the column name list
df.columns = map(str.strip, df.columns)
df['date'] = pd.to_datetime(df['date'])
df.set_index('date', inplace=True)

years = np.unique(df.index.year)

output_filename = 'daily_mass_balance_plot__' + output_file_prefix + '.pdf'

with PdfPages(output_filename) as pdf:

    for n in range(0,len(years)):
        fig, ax = plt.subplots(10,1,figsize=(15, 25),sharex=False)
        df_ss = df.loc[df.index.year==years[n]]
        y0 = np.zeros(len(df_ss))

        p_temp           = 0
        p_inp            = 1
        p_snow_stor      = 2
        p_int_stor       = 3
        p_soil_stor      = 4
        p_fao56          = 5
        p_fao56_evap     = 6
        p_et             = 7
        p_out            = 8
        p_mb             = 9

    
        lu = str(int(df_ss['landuse_code'][0]))
        soil_group = str(int(df_ss['soil_group'][0]))
    
        ax.flat[p_temp].plot(df_ss.index, y0+32, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_inp].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_snow_stor].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_int_stor].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_soil_stor].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_et].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_out].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
        ax.flat[p_mb].plot(df_ss.index, y0, marker=' ', linestyle='-', color='black', linewidth=0.15)
    
        # Temperature plot
        ax.flat[p_temp].fill_between(df_ss.index, df_ss['tmin'],df_ss['tmax'], color='grey', alpha=0.4)
        ax.flat[p_temp].fill_between(df_ss.index, df_ss['tmin'],df_ss['tmax'], color='cyan', alpha=0.5,
            where=df_ss['tmax'] <= 32)
        ax.flat[p_temp].set(ylabel='Air temp (°F)')
    
        # Inputs plot
        ax.flat[p_inp].bar(df_ss.index,df_ss['rainfall'],color='blue',alpha=0.7, label='rainfall')
        ax.flat[p_inp].bar(df_ss.index, df_ss['snowmelt'], bottom=df_ss['rainfall'],color='green',alpha=0.6,
            label='snowmelt')
        ax.flat[p_inp].bar(df_ss.index, df_ss['irrigation'], bottom=df_ss['snowmelt']+df_ss['rainfall'],color='rebeccapurple',                alpha=0.6,label='irrigation')
    
        ax2 = ax.flat[p_inp].twinx()
    
        gross_precip = df_ss['rainfall'] + df_ss['snowfall']
    
        ax2.plot(df_ss.index, np.cumsum(gross_precip), color="darkblue", label='cumulative gross precip')   
        ax2.plot(df_ss.index, np.cumsum(df_ss['irrigation']), color="darkviolet", label='cumulative irrigation')   
    
        ax.flat[p_inp].legend()
        ax2.legend()
        ax2.set(ylabel='Cumulative input (inches)')
        ax.flat[p_inp].set(ylabel='Input (inches)')
    
        # Snow storage plot
        ax.flat[p_snow_stor].fill_between(df_ss.index, y0,df_ss['snow_storage'],color='blue',alpha=0.4,
           where=df_ss['snow_storage'] > 0, interpolate=True)
        ax.flat[p_snow_stor].set(ylabel='Snow storage (inches)')   
    
        # Interception storage plot
        ax.flat[p_int_stor].fill_between(df_ss.index, y0,df_ss['interception_storage'],color='green',
           alpha=0.8, where=df_ss['interception_storage'] > 0, interpolate=True)
        ax.flat[p_int_stor].set(ylabel='Interception storage (inches)')   
    
        # Soil storage plot
        ax.flat[p_soil_stor].fill_between(df_ss.index, y0,df_ss['soil_storage'],color='tan',alpha=0.4)
        ax.flat[p_soil_stor].plot(df_ss.index, df_ss['soil_storage_max'],color='brown')
        ax.flat[p_soil_stor].set(ylabel='Soil storage\n (inches)')
    
        # FAO-56 params plot
        ax.flat[p_fao56].plot(df_ss.index, df_ss['plant_stress_coef_ks'],color='green',
            marker=' ', linestyle='--', alpha=0.4, label='plant stress coefficient (Ks)')
        ax.flat[p_fao56].plot(df_ss.index, df_ss['crop_coefficient_kcb'],color='orange',
            marker='.', linestyle="-", alpha=0.5, label='crop coefficient (Kcb)')

        ax2 = ax.flat[p_fao56].twinx()
        ax2.plot(df_ss.index, df_ss['current_rooting_depth'],color='brown',
            marker=' ', alpha=0.7, label='current_rooting_depth')

        ax.flat[p_fao56].legend()    
        ax.flat[p_fao56].set(ylabel='Coefficient value\n(dimensionless)')
        ax2.legend()    
        ax2.set(ylabel='Rooting depth (feet)')

        # FAO-56 evap params plot
        ax.flat[p_fao56_evap].plot(df_ss.index, df_ss['evaporable_water_storage'],color='green',
            marker=' ', linestyle='--', alpha=0.4, label='evaporable water storage')
        ax.flat[p_fao56_evap].plot(df_ss.index, df_ss['evaporable_water_deficit'],color='red',
            marker=' ', alpha=0.7, label='evaporable_water_deficit')

        ax2 = ax.flat[p_fao56_evap].twinx()
        ax2.plot(df_ss.index, df_ss['evap_reduction_coef_kr'],color='purple',
            marker=' ', linestyle="-.", alpha=0.5, label='evap reduction coef (Kr)')
        ax2.plot(df_ss.index, df_ss['surf_evap_coef_ke'],color='blue',
            marker=' ', linestyle="-", alpha=0.5, label='surface evap coef (Ke)')
        ax2.plot(df_ss.index, df_ss['crop_coefficient_kcb'],color='orange',
            marker='.', linestyle="-", alpha=0.5, label='crop coefficient (Kcb)')
        ax.flat[p_fao56_evap].legend()    
        ax2.legend()
        ax.flat[p_fao56_evap].set(ylabel='Evaporation storage')
        ax2.set(ylabel='Coefficient value\n(dimensionless)')

        # Actual ET plot
        ax.flat[p_et].bar(df_ss.index, -df_ss['bare_soil_evap'], 
            color='purple', label='bare_soil_evap')
        ax.flat[p_et].bar(df_ss.index, -df_ss['crop_etc'], color='gold', 
            label='crop_etc', bottom= -df_ss['bare_soil_evap'])    
        ax.flat[p_et].bar(df_ss.index, -df_ss['actual_et_interception'], 
            color='green', bottom = -df_ss['actual_et_soil'], label='actual_et_interception')
        ax.flat[p_et].plot(df_ss.index, -df_ss['actual_ET'],marker=' ',linestyle='dashdot',color='black', 
            label='total actual ET')    
        ax.flat[p_et].plot(df_ss.index, -df_ss['reference_ET0'],marker=' ',linestyle='--',color='orange', 
            label='reference ET0')    
        ax.flat[p_et].set(ylabel='Evapotranspiration\n from soil (inches)')

        if flux_df is not None:
            flux_df_ss = flux_df.loc[flux_df.index.year==years[n]]
            ax.flat[p_et].plot(flux_df_ss.index, -flux_df_ss['AET_in'], marker='x',linestyle='none',
                color='red',label='flux tower')

        ax2 = ax.flat[p_et].twinx()
        ax2.plot(df_ss.index, np.cumsum(-df_ss['actual_ET']), color="darkgreen", label='cumulative actual ET')   
        ax2.legend()
        ax.flat[p_et].legend()
    
        # Ouflows
        ax.flat[p_out].bar(df_ss.index,-df_ss['runoff'],color='red',alpha=0.5, label='runoff')
        ax.flat[p_out].bar(df_ss.index, -df_ss['net_infiltration'],color='blue',alpha=0.7,
            label='net infiltration', bottom = -df_ss['runoff'])
        ax.flat[p_out].bar(df_ss.index, -df_ss['rejected_net_infiltration'],color='orange',alpha=0.7,
            label='rejected net infiltration', bottom = -df_ss['net_infiltration'])
        ax2 = ax.flat[p_out].twinx()
        ax2.plot(df_ss.index, np.cumsum(-df_ss['runoff']), color="darkred", label='cumulative runoff')   
        ax2.plot(df_ss.index, np.cumsum(-df_ss['net_infiltration']), color="darkblue", 
            label='cumulative net infiltration')   
        ax2.legend()
        ax.flat[p_out].legend()
        
        mb = df_ss['rainfall'] + df_ss['snowmelt'] + df_ss['runon']  + df_ss['irrigation']  \
             - df_ss['actual_ET']                                                           \
             - df_ss['runoff'] - df_ss['delta_soil_storage'] - df_ss['net_infiltration']    \
             - df_ss['rejected_net_infiltration']
    
        ax.flat[p_mb].plot(df_ss.index,mb)
        
        ax.flat[p_mb].set(ylabel="Mass balance (inches)")
        ax.flat[0].set(title='water balance for ' + str(years[n]) + ';  land-use: ' + lu + '  soil group: ' +     soil_group + '  ' + locator_text)
    
        fig.tight_layout()
        #fig.show()   
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()