import time

import matplotlib
matplotlib.use('Qt4Agg')

import numpy

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units

from galaxies import fi_utils
#===========================================================================================================
home = '/home/mher'

params = {
         'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-test', # the path of the dir containing the simulation
         }

snap_index = 4
#===========================================================================================================

######################################## loading the original snapshot #####################################

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

# setting the filename
suffix = '%06d' % snap_index
snapName = 'fiout.%s' % suffix 
filename = params['rundir'] + '/firun/' + snapName 

#loading the sph simulation data
print 'loading snapshot %s : ' % filename
gas_fi, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

# getting the gas particles in cgs units
gas1 = fi_utils.convert_units_to_pdr_units(gas_fi, metallicity)

print 'done reading fi snapshot : %s' % filename
print 'number of sph particles in snapshot = %d' %  len(gas1)

################################ loading the fi info from the .npz.ext snapshot file ########################
################################ loading the fi info from the .npz.ext snapshot file ########################
################################ loading the fi info from the .npz.ext snapshot file ########################
################################ loading the fi info from the .npz.ext snapshot file ########################
################################ loading the fi info from the .npz.ext snapshot file ########################
################################ loading the fi info from the .npz.ext snapshot file ########################

################################ copying the attributes #####################################################
if True:

    # the name of the sampled snapshot
    sampled_snap_filename = filename + '.ext.npz'
    
    gas2, attr_read = fi_utils.load_gas_particle_info(sampled_snap_filename)

    ## making gas1 in the same order as gas2
    inds1 = numpy.where(numpy.log(gas1.n) < numpy.log(1e2))[0]
    inds2 = numpy.where(numpy.log(gas1.n) > numpy.log(1e2))[0]
    inds2s = inds2[numpy.argsort(gas1.n[inds2])]
    gas1_ordered = gas1[numpy.hstack((inds1, inds2s))]
    
    ## assigning the attributes as empty arrays to gas2
    attr2_id = numpy.zeros(len(gas2), 'i4')
    attr2_radius = numpy.zeros(len(gas2), 'f8')
    
    inds_has_children = gas2.get_inds_has_children()
    inds_children = gas2.get_inds_children()
    
    npp = numpy.int32(float(inds_children.size)/float(inds_has_children.size))
    no = len(gas1_ordered)
    print 'number of children per parent', npp
    print 'number of particles in the original set', no
    
    # the index of the child (from zero to number of sampled points - 1)
    child_indicies = gas2.children[inds_has_children]
    
    # the index of the parent
    pinds = gas2.parent 
    
    ## the data of the attributes to be assigned from the original set 
    attr1_id = gas1_ordered.id
    attr1_radius = gas1_ordered.radius
    
    ## copying the attributes to the parents (the are at the beginning of the block)
    attr2_id[0:no] = attr1_id 
    attr2_radius[0:no] = attr1_radius
    
    ids_parents = numpy.int32(gas1_ordered.id[inds_has_children])
    radii_parents = gas1_ordered.radius[inds_has_children]
    
    #for i, child in enumerate(child_indicies[::100]):
    for i, child in enumerate(child_indicies):
            
        cid = numpy.int32(child)
        
        print 'proccessing child index', cid

        # the range of indicies of the childrent
        indx_rng = [no + i*npp, no + (i+1)*npp]

        # unique indecies in suspect bunch of children (all should have the same index)
        u_pinds = numpy.unique(pinds[indx_rng[0]:indx_rng[1]])
        if u_pinds.size > 1:
            raise ValueError('more than one unqiue index found for this bunch')

        parent_index = inds_has_children[i]
        #print gas2.children[parent_index] 
        
        ## copying the attributes to the childrent
        attr2_id[indx_rng[0]:indx_rng[1]] = attr1_id[parent_index]
        attr2_radius[indx_rng[0]:indx_rng[1]] = attr1_radius[parent_index]

        print '\tid set = %d, radius set = %.2f' % (attr1_id[parent_index], attr1_radius[parent_index])
        
        '''        
        inds_children_in_gas2 = numpy.where(pinds == cid)

        if inds_children_in_gas2[0].size == 0:
            raise ValueError('children not found...smoething wrong')
        if inds_children_in_gas2[0].size != npp:
            raise ValueError('children not found...smoething wrong')
         
        attr2_id[inds_children_in_gas2] = ids_parents[i]
        attr2_radius[inds_children_in_gas2] = radii_parents[i]
        
        print 'id set = %d, radius set = %.2f' % (ids_parents[i], radii_parents[i])
        '''
        
        '''
        # the range of indicies of the childrent
        indx_rng = [no + i*npp, no + (i+1)*npp]
        
        # unique indecies in suspect bunch of children (all should have the same index)
        u_pinds = numpy.unique(pinds[indx_rng[0]:indx_rng[1]])
        if u_pinds.size > 1:
            raise ValueError('more than one unqiue index found for this bunch')
        
        ## copying the attributes to the childrent
        attr2_id[indx_rng[0]:indx_rng[1]] = attr1_id[u_pinds[0]]
        attr2_radius[indx_rng[0]:indx_rng[1]] = attr1_radius[u_pinds[0]]
    
        ############ misc stuff ###########
        print gas1.x[u_pinds[0]]
        print gas2.x[indx_rng[0]:indx_rng[1]]
        ###################################
        '''
    gas2.id = attr2_id
    gas2.radius = attr2_radius
    
    
    attr_list = ('Av', 'G0', 'Pe', 'T', 'gmech', 'mass', 'n', 'vdisp', 
                 'vx', 'vy', 'vz', 'x', 'y', 'z', 'id', 'radius',
                 'weights', 'children', 'parent')
    fi_utils.save_gas_particle_info(sampled_snap_filename, gas2, attr_list)
    print 'saved snapshot with sampled data to \n\t %s' % sampled_snap_filename 

    if True:
        '''doing some random checks to see if the ids match.  The positions of a parent and
        its childrent should be overlapping, so the difference should be zero, which is what
        is checked here, we could have checked other attributes which are the same as well'''
        
        q1 = gas1_ordered.y
        id1 = gas1_ordered.id
         
        q2 = gas2.y
        id2 = gas2.id
        
        gas2_ids_with_children = gas2.id[gas2.get_inds_has_children()]
        
        for i, this_id in enumerate(gas2_ids_with_children):
            
            q1_this = q1[ id1 == this_id]
            q2_this = q2[ id2 == this_id]

            print q1_this, numpy.unique(q2_this)            
            print '--------------', i

###################### loading the .npz.ext info and copying some attr info to the states.npz###############
###################### loading the .npz.ext info and copying some attr info to the states.npz###############
###################### loading the .npz.ext info and copying some attr info to the states.npz###############
###################### loading the .npz.ext info and copying some attr info to the states.npz###############
###################### loading the .npz.ext info and copying some attr info to the states.npz###############
###################### loading the .npz.ext info and copying some attr info to the states.npz###############
        
###########################################################################################
if True:

    # the name of the sampled snapshot
    sampled_snap_filename = filename + '.ext.npz'    
    gas2, attr_read = fi_utils.load_gas_particle_info(sampled_snap_filename)
    
    # the name of the states snapshot
    states_snap_filename = filename + '.states.npz'
    gas3, attr_read = fi_utils.load_gas_particle_info(states_snap_filename)

    # getting the original gas set
    gaso = gas2[gas2.get_inds_original_set()]
    
    ## assigning the attributes as empty arrays to gas2
    attr3_id = numpy.zeros(len(gas3), 'i4')
    attr3_radius = numpy.zeros(len(gas3), 'f8')

    ## determinig the particles based on their (x,v) location (must be unique for each particle
    ## in the original set

    # the state of the original particles
    qo = gaso.x**2 + gaso.y**2 + gaso.z**2 + gaso.vx**2 + gaso.vy**2 + gaso.vz**2
    if qo.size != len(len(gaso)):
        raise ValueError("not all particles can be identified with a unique state")
    
    # the state of the sampled particles of gas3 
    q3 = gas3.x**2 + gas3.y**2 + gas3.z**2 + gas3.vx**2 + gas3.vy**2 + gas3.vz**2
    


###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############
###################### checking the assigned information by selecting ids z###############

if True:
    
    ids_children = numpy.unique(gas3.id[ gas3.get_inds_children() ])
    
    q1 = gas1.x
    ids1 = gas1.id
    
    q3 = gas3.x
    ids3 = gas3.id
    
    for i, cid in enumerate(ids_children):
        
        q1_this = q1[ ids1 == cid]
        q3_this = q3[ ids3 == cid]
        
        if numpy.unique(q1_this) - numpy.unique(q3_this) != 0:
            raise ValueError('difference detected')
        
        if i % 100 == 0:
            print '------------------%d-------------------' % i
        