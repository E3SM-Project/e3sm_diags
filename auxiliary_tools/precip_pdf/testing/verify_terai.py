import sys                                                                                            
sys.path.insert(0, '/global/u2/c/chengzhu/e3sm_diags')                                                
from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter                              
from e3sm_diags.driver.utils.dataset_xr import Dataset                                                
import numpy as np                                                                                    
                                                                                                      
param = PrecipPDFParameter()                                                                          
param.reference_data_path ='/global/cfs/cdirs/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_Daily_1DD/time_series/'                 
param.ref_name = 'GPCP_Daily_1DD'                                                                     
param.ref_start_yr = '2001'                                                                           
param.ref_end_yr = '2010'                                                                             
param.ref_timeseries_input = True                                                                     
param.run_type = 'model_vs_obs'                                                                       
                                                                                                      
ref_data = Dataset(param, data_type='ref')                                                            
ref_ds = ref_data.get_time_series_dataset('PRECT', single_point=True)                                 
var = ref_ds['PRECT']
units = var.attrs.get('units', 'unknown')
print('=== Raw GPCP Data Check ===')
print(f'Units: {units}')
print(f'Shape: {var.shape}')
print(f'Time range: {var.time.dt.year.values[0]}-{var.time.dt.year.values[-1]}')
print(f'Min value: {np.nanmin(var.values):.6f}')
print(f'Max value: {np.nanmax(var.values):.6f}')
print(f'Mean value: {np.nanmean(var.values):.6f}')
print(f'Median value: {np.nanmedian(var.values):.6f}')
print(f'Percent zeros: {100 * np.sum(var.values == 0) / var.values.size:.2f}%')
print(f'Percent < 0.1: {100 * np.sum(var.values < 0.1) / var.values.size:.2f}%')
