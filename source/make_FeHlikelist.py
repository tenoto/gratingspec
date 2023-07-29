#!/usr/bin/env python

import os
cmd = 'fselect /Users/enoto/work/tmp/atomdb_v3.0.9/apec_v3.0.9_linelist.fits+2  apec_v3.0.9_linelist_fe_Hlike.fits "Element==26 .and. ion==26"'
print(cmd)
os.system(cmd)