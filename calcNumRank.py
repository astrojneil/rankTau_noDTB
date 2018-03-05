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
        avg_vy_cloud = cloudRegion.quantities.weighted_average_quantity('velx', 'ones') #vely is the velocity in the radial direction (towards/away obs)

        #need to add the frame velocity!
        #get the frame velocity
        f = h5py.File(directory+runName+'/KH_hdf5_chk_'+i, 'r')
        velframe = 0.0 #f['real scalars'][7][1] #cm/s
        velFrames.append(velframe/1.0e5) #append frame vel in km/s
        f.close()

        vy_cloud = (avg_vy_cloud.value+velframe)/1.0e5  #convert to km/s
        velList.append(vy_cloud)

    print('Frame vels:')
    print(velFrames)
    print('Cloud vels:')
    print(velList)
    return velList, velFrames


#given one data file sphere size, and an ion, return average absorption  (given velocity in km/s)
def calcRankTau(directory, runName, runNumber, ions, velocityBin, frameVel):
    #load and add ions to the dataset
    frameVel = frameVel*1.0e5*(centimeter/second)  #convert to cm/s
    data = yt.load(directory+runName+'/KH_hdf5_chk_'+runNumber)
    data.add_field(('gas', 'metallicity'), function=_metallicity, display_name="Metallicity", units='Zsun')

    bList = []
    vList = []
    #add ion fields to the dataset
    for ion in ions:
        tri.add_ion_fields(data, ions=[ion['ion']])

        #select entire region
        reg = data.all_data()

        #projection of the number density (= fieldname)
        p = data.proj(ion['fieldname'], 0)
        #make fixed resolution buffer (same as when a figure is made)
        #to have an array with every value in each pixel
        frb = yt.FixedResolutionBuffer(p, (-2.464e21, 2.464e21, -2.464e21, 2.464e21), (800, 800))
        #flatten the frb to a 1d array
        flattened_num = frb[ion['fieldname']].flatten()
        projectedNum = flattened_num

        #rank sort the projected taus and then make plot
        sorted_num = np.sort(projectedNum, kind='quicksort')

        #find the average velocity for the ion
        average_velx = reg.quantities.weighted_average_quantity('velx', ion['fieldname'])+frameVel


        #write the ranked column densities to a file
        writeFile = open('../rankNum_noDTB'+ion['ionfolder']+'rankNum'+runName+'_v'+str(velocityBin)+'.txt', 'w')
        for i in range(len(sorted_num)):
            writeString = str(sorted_num[i].value)+'\n'
            writeFile.write(writeString)
        writeFile.close()


        #add b^2 thermal as a field
        def _btherm(field, data):
            topFrac = 2.0*kb*data['temperature']
            botFrac = ion['massNum']*mp
            b = topFrac/botFrac
            return b

        data.add_field(('gas', 'b_therm'), function=_btherm, display_name="B thermal squared", units = 'cm**2/s**2', force_override=True)
        #project b thermal
        btherm = data.proj('b_therm', 0, weight_field=ion['fieldname'])
        frb_therm = yt.FixedResolutionBuffer(btherm, (-2.464e21, 2.464e21, -2.464e21, 2.464e21), (800, 800))
        flattened_btherm = frb_therm['b_therm'].flatten()

        #add b^2 doppler as a field
        def _bdop(field, data):
            b = (data['velx'] - average_velx)**2
            return b  #multiply by number density to weight the b value

        data.add_field(('gas', 'b_doppler'), function=_bdop, display_name="B doppler squared", units = 'cm**2/s**2', force_override=True)
        #project b doppler
        bdop = data.proj('b_doppler', 0, weight_field=ion['fieldname'])
        frb_dop = yt.FixedResolutionBuffer(bdop, (-2.464e21, 2.464e21, -2.464e21, 2.464e21), (800, 800))
        flattened_bdop = frb_dop['b_doppler'].flatten()

        good_btherm = flattened_btherm[~np.isnan(flattened_btherm)]
        good_bdop = flattened_bdop[~np.isnan(flattened_bdop)]

        #add the two b's and sqrt
        b_tot_Sq = good_btherm+good_bdop

        b_tot = np.sqrt(b_tot_Sq)

        #compute the average of the summation of b projections
        b_avg = np.average(b_tot)

        bList.append(b_avg)
        vList.append(average_velx)
        print(b_avg)
        print(average_velx)

    return bList, vList


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

    run23 = { 'Name':'T1_v1700_chi1000_lref6',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.5,
        'f_list':['0021', '0029', '0038']}


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
    writebv_file = open('../rankNum_noDTB/Totalrun_allIon_bv.txt', 'w')
    writebv_file.write('Run, frame, Ion, Average_vel(cm/s), b(cm/s)\n')
    for run in runList:

        #find appropriate velocity bins
        velBins, velFrames = findCloudVel(run['Dir'], run['Name'], run['f_list'])  #velocities are in km/s
        #velBins = [116.88]
        for v in range(len(velBins)): #velocity = velBins[v]
            #make the ranked tau plot for each ion/velocity bin
            bList, vList = calcRankTau(run['Dir'], run['Name'], run['f_list'][v], ionList, v, velFrames[v])

            for i in range(len(ionList)):
                writeString = run['Name']+', '+str(v)+', '+ionList[i]['ionfolder'][1:-1]+', '+str(vList[i].value)+', '+str(bList[i])+'\n'
                print(writeString)
                writebv_file.write(writeString)
    writebv_file.close()


if __name__ =="__main__":
    main()
