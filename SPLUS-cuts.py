import splusdata
import pandas as pd

conn = splusdata.connect('juliamoliveira', '10203040')

df = pd.read_csv('selected-gals.csv')

for key, value in df.iterrows():
    hdu = conn.get_cut(value.RA, value.DEC, 128, 'R')
    print(hdu)
    print(value.RA, value.DEC)
    hdu.writeto(f'{value.ID}.fits')