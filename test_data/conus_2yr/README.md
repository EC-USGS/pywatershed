## conus_2yr

This dir contains the version controlled files for this test domain.

The actual run directory is configured as shown:

```
(pynhm_tg) tallgrass ~/pynhm/test_data > pwd
/home/jmccreight/pynhm/test_data

(pynhm_tg) tallgrass ~/pynhm/test_data > ls -lah conus_2yr_run_domain
lrwxrwxrwx 1 jmccreight users 66 Mar 29 13:37 conus_2yr_run_domain -> /caldera/projects/usgs/water/impd//jmccreight/conus_2yr_run_domain

(pynhm_tg) tallgrass ~/pynhm/test_data > ls -lah conus_2yr_run_domain/
total 2.7G
drwxr-sr-x 3 jmccreight impd 4.0K Apr  1 12:32 .
drwxr-sr-x 4 jmccreight impd 4.0K Mar 29 13:30 ..
lrwxrwxrwx 1 jmccreight impd   55 Apr  1 12:32 control.test -> /home/jmccreight/pynhm/test_data/conus_2yr/control.test
lrwxrwxrwx 1 jmccreight impd   57 Apr  1 12:32 conus_2yr.yaml -> /home/jmccreight/pynhm/test_data/conus_2yr/conus_2yr.yaml
-rw-r--r-- 1 jmccreight impd 3.6K Mar 29 21:18 model.out
-r-xr-xr-x 1 jmccreight impd 211M Mar 29 13:33 myparam.param
drwxr-sr-x 2 jmccreight impd 4.0K Mar 29 20:58 output
-r-xr-xr-x 1 jmccreight impd 401M Mar 29 13:35 prcp.cbh
-r-xr-xr-x 1 jmccreight impd 481M Mar 29 13:35 rhavg.cbh
-r-xr-xr-x 1 jmccreight impd  56M Mar 29 13:34 sf_data
-rw-r--r-- 1 jmccreight impd 655M Mar 29 20:58 soltab_debug
-r-xr-xr-x 1 jmccreight impd 482M Mar 29 13:35 tmax.cbh
-r-xr-xr-x 1 jmccreight impd 479M Mar 29 13:35 tmin.cbh
```