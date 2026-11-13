"""Portable posterior IO without installing h5py.
New files are standard NetCDF3 64-bit-offset, with group-prefixed variable names.
Archived NetCDF4/HDF5 files are read through the existing Homebrew h5dump.
load() reconstructs the same xarray DataTree used by ArviZ and analysis scripts.
"""
from pathlib import Path
import json, re, subprocess, tempfile, os
import numpy as np
import pandas as pd
import xarray as xr
ROOT=Path(__file__).resolve().parents[1]

def save(tree, path):
    path=Path(path)
    groups={}
    for group in ['posterior','sample_stats']:
        if group not in tree: continue
        ds=tree[group].to_dataset().copy()
        ds.attrs={}
        for name in ds.variables: ds[name].attrs={};ds[name].encoding={}
        groups[group]=ds.rename({v:group+'__'+v for v in ds.data_vars})
    flat=xr.merge(list(groups.values()),compat='override')
    flat.attrs={'format_note':'NetCDF3: variable prefix encodes DataTree group; load with src/nc_io.py', 'groups':','.join(groups)}
    tmp=path.with_suffix('.partial.nc')
    flat.to_netcdf(tmp,engine='scipy',format='NETCDF3_64BIT')
    os.replace(tmp,path)

def hdf_array(path,var):
    cmd='/opt/homebrew/bin/h5dump'
    h=subprocess.check_output([cmd,'-H','-d',var,str(path)],text=True,stderr=subprocess.DEVNULL)
    shape=tuple(map(int,re.search(r'SIMPLE \{ \( ([\d, ]+) \)',h)[1].split(',')))
    dtype='<f4' if 'H5T_IEEE_F32LE' in h else '<f8'
    with tempfile.NamedTemporaryFile(dir=ROOT/'tmp',suffix='.bin') as f:
        subprocess.run([cmd,'-d',var,'-b','LE','-o',f.name,str(path)],stdout=subprocess.DEVNULL,check=True)
        return np.fromfile(f.name,dtype=dtype).reshape(shape)

def load(path):
    path=Path(path)
    with path.open('rb') as f: magic=f.read(8)
    if magic[:3]==b'CDF':
        flat=xr.load_dataset(path,engine='scipy')
        d={}
        for group in flat.attrs['groups'].split(','):
            names=[v for v in flat.data_vars if v.startswith(group+'__')]
            ds=flat[names].rename({v:v.split('__',1)[1] for v in names})
            # Drop coordinates that belong only to other groups.
            ds=ds.drop_vars([c for c in ds.coords if c not in ds.dims],errors='ignore')
            ds.attrs={};d[group]=ds
        return xr.DataTree.from_dict(d)
    if magic!=b'\x89HDF\r\n\x1a\n': raise ValueError(f'Unrecognised file {path}')
    name=path.stem.replace('idata_','')
    info=path.parent/f'fit_{name}.json'
    cov=json.loads(info.read_text())['spec']['X'] if info.exists() else ['idhm_renda','idhm_educ','logpop','sh_age_lt20','sh_age_ge35','sh_esc_le7','sh_nonwhite']
    header=subprocess.check_output(['/opt/homebrew/bin/h5dump','-H','-g','/posterior',str(path)],text=True)
    present=re.findall(r'^   DATASET "([^"]+)"',header,re.M)
    ds={}
    for v in ['b0','beta','sigma_u','rho','sigma_hr','gamma_hr','u']:
        if v not in present: continue
        a=hdf_array(path,'/posterior/'+v)
        dims=['chain','draw']+(['cov'] if v=='beta' else ['muni'] if v=='u' else ['hr'] if v=='gamma_hr' else [])
        ds[v]=(dims,a)
    base=pd.read_csv(ROOT/'data/analysis_muni.csv',dtype={'CODMUNRES':str,'regsaud_id':str})
    coords={'cov':cov,'muni':base.CODMUNRES.to_numpy()}
    if 'gamma_hr' in ds: coords['hr']=np.unique(base.regsaud_id)
    post=xr.Dataset(ds,coords=coords)
    return xr.DataTree.from_dict({'posterior':post})
