# %%
from dask.distributed import Client
client = Client('tcp://134.121.21.204:8786')
client

# %%
import xarray as xr
import hvplot.xarray
import zarr

# %%
mpl = xr.open_mfdataset('/mnt/disk2/data/mpl/smt*2010*.cdf', concat_dim='time', combine='nested', parallel=True, chunks={'time': 5000})

# %%
mpl['time'] = mpl.base_time + mpl.time_offset

# %%
#mpl.energy.hvplot()

# %%
mpl.to_zarr('/mnt/disk2/data/mpl/zarr/')

# %%
client.close()
