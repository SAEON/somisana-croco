from datetime import datetime
import crocotools_py.preprocess as pre
    
wasa_grid = '/media/external_1/geo_em.d03.nc'
wasa_zarr_dir = '/media/external_1/uv_ds_10m_v2'
croco_grd = '/home/gfearon/code/somisana-croco/configs/algoa_01/croco_v1.3.1/GRID/croco_grd.nc'
croco_blk_dir = '/home/gfearon/code/somisana-croco/configs/algoa_01/croco_v1.3.1/ERA5'
out_wasa_dir = '/home/gfearon/code/somisana-croco/configs/algoa_01/croco_v1.3.1/WASA3'
ref_date = datetime(1993,1,1)
start_date = datetime(2011,1,1)
end_date = datetime(2019,12,1)

pre.make_WASA3_from_blk(wasa_grid, 
                  wasa_zarr_dir, 
                  croco_grd,
                  croco_blk_dir,
                  out_wasa_dir,
                  ref_date,
                  start_date,
                  end_date,
                  croco_atmos='ERA5',
                  interp_method='linear'
                  )    

# after running make_WASA3_from_blk we found some missing values in the WASA data for some months
# see months_with_missing in this directory
# so for now we're just interpolating over some of the gaps (if it's not too much missing data) 
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2012M07_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2012M07.nc')
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2012M08_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2012M08.nc')
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2017M08_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2017M08.nc')
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2017M12_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2017M12.nc')
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2000M11_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2000M11.nc')
# pre.fill_blk(croco_grd, out_wasa_dir + '/croco_blk_WASA3_Y2000M12_corrupted.nc', out_wasa_dir + '/croco_blk_WASA3_Y2000M12.nc')