#!/usr/bin/env python
'''
Test script downloading and running cell models from the Allen Brain
Institute's Cell Type Database
'''
import os
import numpy as np
import neuron
import LFPy
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
try:
    from allensdk.model.biophysical_perisomatic import runner 
    from allensdk.api.queries.biophysical_perisomatic_api import BiophysicalPerisomaticApi
except ImportError as ie:
    raise ie, 'install AllenSDK from http://alleninstitute.github.io/AllenSDK/'

#fetch files from the Cell Type Database
working_directory = 'neuronal_model'
bp = BiophysicalPerisomaticApi('http://api.brain-map.org')
bp.cache_stimulus = False # set to True to download the large stimulus NWB file
neuronal_model_id = 472451419    # get this from the web site as above
bp.cache_data(neuronal_model_id, working_directory=working_directory)

#change directory
os.chdir(working_directory)

#compile and load NEURON NMODL files
os.system('nrnivmodl modfiles')    
neuron.load_mechanisms('.')


#load model in neuron using the Allen SDK
description = runner.load_description('manifest.json')

# configure NEURON
utils = runner.Utils(description)
h = utils.h

# configure model
manifest = description.manifest
morphology_path = description.manifest.get_path('MORPHOLOGY')
utils.generate_morphology(morphology_path.encode('ascii', 'ignore'))
utils.load_cell_parameters()

#go back to root folder
os.chdir('..')

#load cell into LFPy from NEURON namespace
cell = LFPy.Cell(morphology=None, delete_sections=False,
                 tstartms=0, tstopms=3000,
                 passive=False, nsegs_method=None,
                 extracellular=False, v_init=-100)
cell.set_rotation(x=np.pi)

#perform a single sweep
PointProcParams = {
    'idx' : cell.get_idx('soma'),
    'record_current' : True,
    'pptype' : 'IClamp',
    'amp' : 0.1,
    'dur' : 2000,
    'delay' : 500,
}

pointProcess = LFPy.StimIntElectrode(cell, **PointProcParams)

#run simulation, record input current and soma potential (default)
cell.simulate(rec_istim=True, rec_variables = [])


#plot
plt.rcParams.update({'font.size': 9.0})
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(3, 3)

ax = fig.add_subplot(gs[:, 0])
zips = []
for x, z in cell.get_idx_polygons(projection=('x', 'y')):
    zips.append(zip(x, z))

polycol = PolyCollection(zips,
                         edgecolors='none',
                         facecolors='k')
ax.add_collection(polycol)
ax.axis(ax.axis('equal'))
ax.set_title('morphology')
ax.set_ylabel('y (mum)')
ax.set_xlabel('x (mum)')

#plot response
ax = fig.add_subplot(gs[0:2, 1:])
ax.plot(cell.tvec, cell.somav)
ax.set_xticks([])
ax.set_ylabel('V (mV)')
ax.set_title('somatic response')

#plot stimulus current
ax = fig.add_subplot(gs[2, 1:])
ax.plot(cell.tvec, pointProcess.i)
ax.set_ylabel('I (nA)')
ax.set_xlabel('t (ms)')
ax.set_title('stimulus current')

plt.show()
