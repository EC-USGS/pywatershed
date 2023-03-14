import pandas as pd
import xarray as xr

df = pd.read_csv('SWB2_variable_values__col_43__row_162__x_150393__y_1179331.csv',index_col='date')
df.columns = df.columns.str.strip()

df.index = pd.to_datetime(df.index)

ds = df.to_xarray()