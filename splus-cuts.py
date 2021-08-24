import splusdata
import pandas as pd

conn = splusdata.connect('user', 'pass')

df = pd.read_csv('selected-gals.csv')

for key, value in df.iterrows():
    hdu = conn.get_cut(value.RA, value.DEC, 128, 'R')
    hdu.writeto('stamps-splus/%.6f_%.6f.fz' % (value.RA, value.DEC))