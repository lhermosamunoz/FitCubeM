from configMix import * 
from FitFunctionMix import FitPoint

nSlots = 4

### Create the inputs
inputs = []
#gal = [5]
#x = [21]
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
