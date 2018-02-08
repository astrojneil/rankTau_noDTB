import yt
import numpy as np
import trident as tri
import matplotlib.pyplot as plt
import subprocess
from yt.data_objects.particle_filters import add_particle_filter
import h5py
from yt.units import centimeter, gram, second, Kelvin, erg

kpc = 3.086e+21*centimeter
c_speed = 3.0e10  #cm/s
mp = 1.6726e-24*gram #grams
kb = 1.3806e-16*erg/Kelvin   #egs/K


#add metallicity to dataset, constant Z = 1 Zsun
def _metallicity(field, data):
    v = data['ones']  #sets metallicity to 1 Zsun
    return data.apply_units(v, "Zsun")

#add absorption to dataset
def _absorption(field, data):
    absorb = np.exp(-sigma*data[('all', 'O_p5_number_density')])
    return absorb

def no_ColumnDepth(pfilter, data):
    colDepth = data['O_p5_number_density']
    filter = np.logical_or(colDepth != 0.0, colDepth > 0.0)
    return filter

def deleteProjections(listofFiles):
    for num in listofFiles:
        subprocess.call("rm KH_hdf5_chk_"+num+"_proj.h5")

######## CHANGE THESE #######
ProjectionFolder = 'Projections_Uneg2'
#ion = 'O VI'
#fieldname = 'O_p5_number_density'
#ionfolder = '/OVI/'
########


#function to read in spectrum, convert to velocity and return
#wavelength and flux numpy arrays
def convert_to_vel(filename, rest_wave):
    #load spectrum; will need to add a case when this is a txt file rather than h5
    wavelength = []
    flux = []
    f = open(filename, 'r')
    i = 1
    for line in f:
        if i > 3:
            splitline = line.split()
            wavelength.append(float(splitline[0]))
            flux.append(float(splitline[1]))
        i=i+1

    wavelength = np.array(wavelength)
    flux = np.array(flux)
    vel = c_speed*(wavelength/rest_wave -1)
    #returns with velocities in cm/s
    return vel, flux

#function to determine the velocity of the cloud for a particular run
#will return an array of velocities (km/s) for the appropriate velocity bins to use
#in the observational data.
def findCloudVel(directory, runName, f_list):
    velList = []
    velFrames = []
    for i in f_list:
        data = yt.load(directory+runName+'/KH_hdf5_chk_'+i)
        allDataRegion = data.all_data()
        cloudRegion = allDataRegion.cut_region(['obj["density"] >= 3.33e-25'])  #at or above 1/3 original density (1e-24)
        avg_vx_cloud = cloudRegion.quantities.weighted_average_quantity('velx', 'ones') #velx is the velocity in the radial direction (towards/away obs)

        #need to add the frame velocity!
        #get the frame velocity  - not needed in x direciton
        #f = h5py.File(directory+runName+'/KH_hdf5_chk_'+i, 'r')
        #velframe = f['real scalars'][7][1] #cm/s
        #velFrames.append(velframe/1.0e5) #append frame vel in km/s
        #f.close()

        vx_cloud = (avg_vx_cloud.value)/1.0e5  #convert to km/s
        velList.append(vx_cloud)

    #print('Frame vels:')
    #print(velFrames)
    print('Cloud vels:')
    print(velList)
    return velList #, velFrames


#given one data file sphere size, and an ion, return average absorption  (given velocity in km/s)
def calcRankTau(directory, runName, runNumber, ions, velocityBin):
    #load and add ions to the dataset
    velocityBin = velocityBin*1.0e5  #convert to cm/s
    #frameVel = frameVel*1.0e5  #convert to cm/s
    data = yt.load(directory+runName+'/KH_hdf5_chk_'+runNumber)
    data.add_field(('gas', 'metallicity'), function=_metallicity, display_name="Metallicity", units='Zsun')

    absorbFracs = []
    #add ion fields to the dataset
    for ion in ions:
        tri.add_ion_fields(data, ions=[ion['ion']])

        #define the function for b(x) for this particular ion
        ### function definition in a for loop... Ew... ###
        def _soundSpeed(field, data):
            topFrac = 2.0*kb*data['temperature']
            botFrac = ion['massNum']*mp
            b = np.sqrt(topFrac/botFrac)
            return b

        #add b_soundSpeed to the dataset
        data.add_field(('gas', 'b_soundSpeed'), function=_soundSpeed, display_name="B Sound Speed", units = 'cm/s', force_override=True)

        #define the function for dTau for this particular ion/velocity
        ### function definition in a for loop... Ew... ###
        def _dTau(field, data):
            term1 = data[ion['fieldname']]/data['b_soundSpeed']
            #assumes velocity is now cm/s
            expTerm = (((data['velocity_x']) - velocityBin*(centimeter/second))/data['b_soundSpeed'])**2.0
            dtau = term1*np.exp(-1.0*expTerm)
            return dtau

        #add dTau to the dataset
        data.add_field(('gas', 'tau'), function=_dTau, display_name="dTau", units = 's/cm**4', force_override=True)

        #make cylinder, project through the y axis
        #reg = data.disk([0.0, 0.2*kpc, 0.0], [0,1,0], (0.25, 'kpc'), (1.6, 'kpc'))

        #select entire region
        reg = data.all_data()

	#plot = yt.ProjectionPlot(data, 'x', 'density')
	#plot.set_zlim('density', 1e-5, 5e-5)
	#plot.save()

        #projection of the number density
        p = data.proj('tau', 0)
        #make fixed resolution buffer (same as when a figure is made)
        #to have an array with every value in each pixel
        frb = yt.FixedResolutionBuffer(p, (-2.464e21, 2.464e21, -2.464e21, 2.464e21), (800, 800))
        #flatten the frb to a 1d array
        flattened_tau = frb['tau'].flatten()


        tau_coef = c_speed*ion['sigma']/np.sqrt(np.pi)
        projectedTau = tau_coef*flattened_tau

        #rank sort the projected taus and then make plot
        sorted_tau = np.sort(projectedTau, kind='quicksort')
        totalpixel = len(sorted_tau)

        pixelNum = np.arange(0.0, totalpixel, 1.0)
        pixelFrac = pixelNum/totalpixel
        #print(pixelFrac)

        writeFile = open('../rankTau_noDTB'+ion['ionfolder']+'rankTau'+runName+'_v'+str(int(round(velocityBin, 0)))+'.txt', 'w')
        for i in range(len(pixelFrac)):
            writeString = str(sorted_tau[i].value)+'\n'
            writeFile.write(writeString)
        writeFile.close()

        #clip = 620800 #0.97
        #plt.plot(pixelFrac[clip:], sorted_tau[clip:])
        #fig = plt.gcf()
        #fig.savefig('rankTau/test_rankTau'+runName+'_v'+str(int(round(velocityBin/1e5, 0)))+'_'+ion['ion']+'.png')
        #fig.clear()


    return absorbFracs

    #test.add_field(('all', 'absorption'), function=_absorption, display_name="Absorption", particle_type=True)
    #add_particle_filter("no_colDepth", function=no_ColumnDepth, filtered_type='all', requires=[ion['fieldname']])
    #test.add_particle_filter('no_colDepth')

    #print("Filtered reloaded projection:")
    #print(testregion[('all', 'O_p5_number_density')][hasColumnDens].size) #this is the size of projection with the same filter prior to saving
    #print("Particle Filter")
    #print(testregion['no_colDepth', 'absorption'].size)

    #plot the filtered column density and absorptions
    #t = yt.ParticlePlot(test, ('no_colDepth', 'px'), ('no_colDepth', 'py'), ('no_colDepth', 'O_p5_number_density'))
    #t.set_ylim(-0.8*kpc, 0.8*kpc)
    #t.set_xlim(-0.8*kpc, 0.8*kpc)
    #t.save()

    #t2 = yt.ParticlePlot(test, ('no_colDepth', 'px'), ('no_colDepth', 'py'), ('no_colDepth', 'absorption'))
    #t2.set_ylim(-0.8*kpc, 0.8*kpc)
    #t2.set_xlim(-0.8*kpc, 0.8*kpc)
    #t2.save()

def runthroughRuns(ion, fieldname, ionfolder, runList):
    for run in runList:
        #savefolder = '../'+ProjectionFolder+'/'+run['Name']+ionfolder
        #calcAbsorb(run['Dir'], run['Name'], run['f_list'], savefolder, ion, fieldname)
        velocities = findCloudVel(run['Dir'], run['Name'], run['f_list'])
        print velocities

    print("Finished with Ion"+ion)

def main():
##### Runs to have spectra generated #######
    run1 = { 'Name':'T0.3_v1000_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.8,
        'tcc':1.7,
        'f_list':['0013', '0038', '0080', '0132']}
    run2 = { 'Name':'T3_v3000_chi3000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.6,
        'tcc':1.8,
        'f_list':['0001', '0004', '0007', '0010']}
    run3 = { 'Name':'T1_v1700_chi1000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.5,
        'tcc':1.8,
        'f_list':['0002', '0010', '0017', '0028']}
    run4 = { 'Name':'T0.3_v1000_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.8,
        'f_list':['0025', '0033', '0042', '0058']}
    run5 = { 'Name':'T3_v3000_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.6,
        'f_list':['0021', '0030', '0040', '0062']}
    run6 = { 'Name':'T1_v1700_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.5,
        'f_list':['0021', '0029', '0038', '0052']}
    run7 = { 'Name':'HC_v1000_chi300_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'f_list':['0054', '0060', '0080', '0107']}
    run8 = { 'Name':'HC_v1700_chi1000_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'f_list':['0024', '0050', '0082', '0083']}
    run9 = { 'Name':'HC_v3000_chi3000_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'f_list':['0007', '0015', '0026', '0049']}
    run10 = { 'Name':'LowCond_v1700_chi300_cond',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'f_list':['0016', '0075', '0115', '0184']}

    run11 = { 'Name':'T0.3_v1700_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':6.5,
        'tcc':1.0,
        'f_list':['0003', '0020', '0046', '0078']}
    run12 = { 'Name':'T0.3_v3000_chi300_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':11.4,
        'tcc':0.56,
        'f_list':['0001', '0004', '0014', '0035']}
    run13 = { 'Name':'T3_v860_chi3000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.2,
        'f_list':['0001', '0003', '0006', '0010']}
    run14 = { 'Name':'T10_v1500_chi10000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.5,
        'f_list':['0001', '0002', '0004', '0008']}
    run15 = { 'Name':'T1_v480_chi1000_cond',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.4,
        'f_list':['0002', '0009', '0017', '0031']}

    run16 = { 'Name':'T0.3_v1700_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.5,
        'f_list':['0022', '0032', '0053', '0085']}
    run17 = { 'Name':'T0.3_v3000_chi300',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':11.4,
        'f_list':['0028', '0044', '0065', '0110']}
    run18 = { 'Name':'T3_v430_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':0.5,
        'f_list':['0010', '0011', '0016', '0024']}
    run19 = { 'Name':'T3_v860_chi3000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0010', '0018', '0030', '0038']}
    run20 = { 'Name':'T1_v3000_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.2,
        'f_list':['0022', '0032', '0048', '0095']}
    run21 = { 'Name':'T10_v1500_chi10000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0014', '0021', '0029', '0044']}
    run22 = { 'Name':'T1_v480_chi1000',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0009', '0013', '0024', '0035']}


#add the runs to the list that will have spectra generated
    runList = []
    runList.append(run1)
    runList.append(run2)
    runList.append(run3)
    runList.append(run4)
    runList.append(run5)
    runList.append(run6)
    runList.append(run7)
    runList.append(run8)
    runList.append(run9)
    runList.append(run10)
    runList.append(run11)
    runList.append(run12)
    runList.append(run13)
    runList.append(run14)
    runList.append(run15)
    runList.append(run16)
    runList.append(run17)
    runList.append(run18)
    runList.append(run19)
    runList.append(run20)
    runList.append(run21)
    runList.append(run22)
    '''
    run6 = { 'Name':'T1_v1700_chi1000_lref4',
                'Dir':'../../Blob_paper1/Files/',
                'Mach':3.5,
                'f_list':['0020', '0030', '0040', '0050']}
    run7 = { 'Name':'T1_v1700_chi1000_lref6',
                'Dir':'../../Blob_paper1/Files/',
                'Mach':3.5,
                'f_list':['0021', '0029', '0038']}
    '''


### dictionaries of ion info
    ion1 = {'ion':'O VI',
        'fieldname':'O_p5_number_density',
        'ionfolder': '/OVI/',
        'rest_wave': 1031.91,
        'data_file': '../Files/S1226-o6-forJNeil',
        'sigma': 1.1776e-18,
        'massNum': 16.0}
    ion2 = {'ion':'C IV',
        'fieldname':'C_p3_number_density',
        'ionfolder': '/CIV/',
        'rest_wave': 1548.18,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 2.5347e-18,
        'massNum': 12.0}
    ion3 = {'ion':'N V',
        'fieldname':'N_p4_number_density',
        'ionfolder': '/NV/',
        'rest_wave': 1242.8,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 8.3181e-19,   #UM. is it actually e-19?   -12/7/17
        'massNum': 14.0}
    ion4 = {'ion':'C II',
        'fieldname':'C_p1_number_density',
        'ionfolder': '/CII/',
        'rest_wave': 1335.66,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 1.4555e-19,
        'massNum': 12.0}
    ion5 = {'ion': 'Ne VII',
            'fieldname': 'Ne_p6_number_density',
            'ionfolder':'/NeVII/',
            'rest_wave': 465.22,
            'data_file': '../Files/S1226-o6-forJNeil',
            'sigma': 2.1922e-19,
            'massNum': 20.0}
    ion6 = {'ion': 'C III',
            'fieldname': 'C_p2_number_density',
            'ionfolder':'/CIII/',
            'rest_wave': 977.02,
            'data_file': '../Files/S1226-o6-forJNeil',
            'sigma': 6.359e-18,
            'massNum': 12.0}
    ion7 = {'ion': 'Mg II',
            'fieldname': 'Mg_p1_number_density',
            'ionfolder':'/MgII/',
            'rest_wave': 1239.92,
            'data_file': '../Files/S1226-o6-forJNeil',
            'sigma': 6.60717e-21,
            'massNum': 24.0}
    ion8 = {'ion': 'Si III',
            'fieldname': 'Si_p2_number_density',
            'ionfolder':'/SiIII/',
            'rest_wave': 1206.5,
            'data_file': '../Files/S1226-o6-forJNeil',
            'sigma': 2.2258e-20,
            'massNum': 28.0}
    ion9 = {'ion': 'Si IV',
            'fieldname': 'Si_p3_number_density',
            'ionfolder':'/SiIV/',
            'rest_wave': 1402.77,
            'data_file': '../Files/S1226-redward-forJNeil',
            'sigma': 3.0694e-18,
            'massNum': 28.0}
    ion10 = {'ion': 'H I 1215.67',
            'fieldname': 'H_p0_number_density',
            'ionfolder':'/HI/',
            'rest_wave': 1215.67,
            'data_file': '../Files/S1226-o6-forJNeil',
            'sigma': 4.3394e-18,
            'massNum': 2.0}

    ionList = []
    ionList.append(ion1)
    ionList.append(ion2)
    ionList.append(ion3)
    ionList.append(ion4)
    ionList.append(ion5)
    ionList.append(ion6)
    ionList.append(ion7)
    ionList.append(ion8)
    ionList.append(ion9)
    ionList.append(ion10)



### run though functions!
    for run in runList:

        #find appropriate velocity bins
        velBins = findCloudVel(run['Dir'], run['Name'], run['f_list'])  #velocities are in km/s
        #velBins = [116.88]
        for v in range(len(velBins)): #velocity = velBins[v]
            #make the ranked tau plot for each ion/velocity bin
            tauList = calcRankTau(run['Dir'], run['Name'], run['f_list'][v], ionList, velBins[v])
            print('Done with velocity '+str(v+1)+'of '+str(len(velBins)))

        print('Finished with Run: '+run['Name'])

            #everything beyond this is chi square absorption fitting
        '''
            cList = np.logspace(-2, 4, 1000)
            #cList = [1.0]
            chi_squareList = []

            for c in cList: #covering fraction = c
                chi_square = 0.0

                for i in range(len(ionList)):
                    #find average Absorption for this ion at this velocity
                    averageAbs = np.exp(c*tauList[i])
                    #load obs data as velocities/normalized fluxes
                    velocities, flux = convert_to_vel(ionList[i]['data_file'], ionList[i]['rest_wave'])
                    #save flux and variance
                    #make velocity bin negative to try to fit to the majority of the absorption profile
                    vel_flux = np.interp(-1.0*velBins[v], velocities, flux)
                    vel_var = (0.05*vel_flux)**2.0

                    #calculate chi term
                    chi_term = (averageAbs-vel_flux)**2.0/vel_var
                    chi_square = chi_square+chi_term
                    #print(str(c)+', '+str(averageAbs)+', '+str(vel_flux))
                chi_squareList.append(chi_square)
            #find minmum of chiSquares
            minChi_index = np.argmin(chi_squareList)

            #print/save results
            bestFitAbsorb = []
            for i in range(len(ionList)):
                averageAbs = np.exp(cList[minChi_index]*tauList[i])
                bestFitAbsorb.append(averageAbs)

            absorbString =  ", ".join(map(str, bestFitAbsorb))
            saveData.write(str(-1.0*velBins[v])+', '+str(np.log10(cList[minChi_index]))+', '+absorbString+'\n')

        #delete projection files made - clean up!
        #deleteProjections(run['f_list'])
        '''




if __name__ =="__main__":
    main()
