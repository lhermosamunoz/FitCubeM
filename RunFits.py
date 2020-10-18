from config import * 
from FitFunction import FitPoint

nSlots = 4

### Create the inputs
inputs = []
for iy in igal:
  for ix in x:
    if CheckFit(x,y):
      inputs.append([ix,iy])

### Run!
from multiprocessing import Pool
pool = Pool(nSlots)
results = pool.map(FitPoint, inputs)
pool.close()
pool.join()

