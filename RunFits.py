from config import * 
from FitFunction import FitPoint

nSlots = 4

### Create the inputs
inputs = []
for iy in gal:
  for ix in x:
    if not CheckFit(ix,iy):
      inputs.append([ix,iy,meth])

### Run!
from multiprocessing import Pool
pool = Pool(nSlots)
results = pool.map(FitPoint, inputs)
pool.close()
pool.join()
