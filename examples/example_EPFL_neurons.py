#!/usr/bin/env python
'''
Test implementation using cell models of the Blue Brain Project with LFPy.
The example assumes that the complete set of cell models available from
https://bbpnmc.epfl.ch/nmc-portal/downloads is unzipped in this folder. 
'''

import os
import sys
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import PolyCollection
from glob import glob
import numpy as np
import neuron
import LFPy

#load some required neuron-interface files
neuron.h.load_file("stdrun.hoc")
neuron.h.load_file("import3d.hoc")

#load only some layer 5 pyramidal cell types
neurons = glob(os.path.join('hoc_combos_syn.1_0_10.allzips', 'L5_*PC*'))[:1]

#flag for cell template file to switch on (inactive) synapses
add_synapses = False

#attempt to set up a folder with all unique mechanism mod files, compile, and
#load them all
if not os.path.isdir('hoc_combos_syn.1_0_10.allmods'):
    os.mkdir('hoc_combos_syn.1_0_10.allmods')
for nrn in neurons:
    for nmodl in glob(os.path.join(nrn, 'mechanisms', '*.mod')):
        while not os.path.isfile(os.path.join('hoc_combos_syn.1_0_10.allmods',
                                              os.path.split(nmodl)[-1])):
            os.system('cp {} {}'.format(nmodl,
                                        os.path.join('hoc_combos_syn.1_0_10.allmods',
                                                     '.')))
os.chdir('hoc_combos_syn.1_0_10.allmods')
os.system('nrnivmodl')
neuron.load_mechanisms('.')
os.chdir('..')

#load the LFPy SinSyn mechanism for stimulus
neuron.load_mechanisms(LFPy.__path__[0])


def get_templatename(f):
    '''
    Assess from hoc file the templatename being specified within
    
    Arguments
    ---------
    f : file, mode 'r'
    
    Returns
    -------
    templatename : str
    
    '''    
    f = file("template.hoc", 'r')
    for line in f.readlines():
        if 'begintemplate' in line.split():
            templatename = line.split()[-1]
            print 'template {} found!'.format(templatename)
            continue
    
    return templatename
    

for nrn in neurons:
    os.chdir(nrn)

    #get the template name
    f = file("template.hoc", 'r')
    templatename = get_templatename(f)
    f.close()
    
    #get biophys template name
    f = file("biophysics.hoc", 'r')
    biophysics = get_templatename(f)
    f.close()
    
    #get morphology template name
    f = file("morphology.hoc", 'r')
    morphology = get_templatename(f)
    f.close()
    
    #get synapses template name
    f = file(os.path.join("synapses", "synapses.hoc"), 'r')
    synapses = get_templatename(f)
    f.close()
    

    print('Loading constants')
    neuron.h.load_file('constants.hoc')

    if not hasattr(neuron.h, morphology): 
        """Create the cell model"""
        # Load morphology
        neuron.h.load_file(1, "morphology.hoc")
    if not hasattr(neuron.h, biophysics): 
        # Load biophysics
        neuron.h.load_file(1, "biophysics.hoc")
    if not hasattr(neuron.h, synapses):
        # load synapses
        neuron.h.load_file(1, os.path.join('synapses', 'synapses.hoc'))
    if not hasattr(neuron.h, templatename): 
        # Load main cell template
        neuron.h.load_file(1, "template.hoc")


    for morphologyfile in glob('morphology/*'):
        # Instantiate the cell(s) using LFPy
        cell = LFPy.TemplateCell(morphology=morphologyfile,
                         templatefile=os.path.join(nrn, 'template.hoc'),
                         templatename=templatename,
                         templateargs=1 if add_synapses else 0,
                         tstopms=500.,)
    
        #set view as in most other examples
        cell.set_rotation(x=np.pi/2)
        
        ##some stimuli
        #PointProcParams = {
        #    'idx' : 0,
        #    'record_current' : False,
        #    'pptype' : 'IClamp',
        #    'amp' : 0.5,
        #    'dur' : 200,
        #    'delay' : 200,
        #}

        PointProcParams = {
            'idx' : 0,
            'pptype' : 'SinSyn',
            'delay' : 200.,
            'dur' : 200.,
            'pkamp' : 0.5,
            'freq' : 0.,
            'phase' : np.pi/2,
            'bias' : 0.,
            'record_current' : False
        }
        
        pointProcess = LFPy.StimIntElectrode(cell, **PointProcParams)
        
        #run simulation
        cell.simulate(rec_imem=True)

        electrode = LFPy.RecExtElectrode(cell, x = np.array([20.]))
        electrode.calc_lfp()
        
        #plot
        gs = GridSpec(2, 2)
        fig = plt.figure()
        fig.suptitle(nrn + '\n' + morphologyfile)
        zips = []
        for x, z in cell.get_idx_polygons(projection=('x', 'z')):
            zips.append(zip(x, z))    
        polycol = PolyCollection(zips,
                                 edgecolors='none',
                                 facecolors='k')
        ax = fig.add_subplot(gs[:, 0])
        ax.add_collection(polycol)
        plt.plot(electrode.x, electrode.z, 'ro')
        ax.axis(ax.axis('equal'))
    
        ax = fig.add_subplot(gs[0, 1])
        ax.plot(cell.tvec, cell.somav)

        ax = fig.add_subplot(gs[1, 1])
        ax.plot(cell.tvec, electrode.LFP[0, ])
    
    os.chdir(os.path.join('..', '..'))

plt.show()

    
