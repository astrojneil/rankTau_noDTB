import yt
import numpy as np
import trident as tri
import matplotlib.pyplot as plt
import subprocess
from yt.data_objects.particle_filters import add_particle_filter
import h5py
from yt.units import centimeter, gram, second, Kelvin, erg
import scipy

kpc = 3.086e+21*centimeter
c_speed = 3.0e10  #cm/s
mp = 1.6726e-24*gram #grams
kb = 1.3806e-16*erg/Kelvin   #egs/K

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

def openIonFile(ion, runName):
    openfile = open('../rankNum'+ion['ionfolder']+ion['ionfolder'][1:-1]+'_bestFitParameters.txt', 'r')
    singleIon = np.zeros((4,6))
    vel = 0
    for line in openfile:
        ls = line.split(', ')
        if ls[0] == runName:
            singleIon[vel][0] = float(ls[1]) #velocity
            singleIon[vel][1] = float(ls[2]) #b
            singleIon[vel][2] = float(ls[3]) #N
            singleIon[vel][3] = float(ls[6]) #q
            singleIon[vel][4] = float(ls[9])  #area
            singleIon[vel][5] = ion['sigma']  #sigma
            vel = vel+1

    return singleIon

def model(t, b, x):
    return t*(0.01)/(1.01-x**b)

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
    #runList.append(run1)
    #runList.append(run2)
    #runList.append(run3)
    runList.append(run4)
    runList.append(run5)
    runList.append(run6)
    #runList.append(run7)
    #runList.append(run8)
    #runList.append(run9)
    #runList.append(run10)
    #runList.append(run11)
    #runList.append(run12)
    #runList.append(run13)
    #runList.append(run14)
    #runList.append(run15)
    runList.append(run16)
    runList.append(run17)
    runList.append(run18)
    runList.append(run19)
    runList.append(run20)
    runList.append(run21)
    runList.append(run22)


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
    #ionList.append(ion1)
    ionList.append(ion2)
    #ionList.append(ion3)
    ionList.append(ion4)
    #ionList.append(ion5)
    #ionList.append(ion6)
    #ionList.append(ion7)
    #ionList.append(ion8)
    ionList.append(ion9)
    #ionList.append(ion10)

    #writefile = open('../rankTau/ranktau_outflowApp.txt', 'w')
    #writefile.write('RunName, averageChi2_alpha, averageChi2_beta, maximumChi2_alpha, maximumChi2_beta\n')

    otherwrite = open('../rankNum/rankNum_cgm_actualBestFit.txt', 'w')
    otherwrite.write('RunName, vel, beta, betaChi\n')

    alphas = np.logspace(-5, 2, 10000)
    betas = np.logspace(-5, 2, 10000)

    for run in runList:
        #save a lot of info
        runMultiCloud = np.zeros((len(ionList), 4))  #run[row][col] = run[ion][vel]
        runSingleCloud = np.zeros((len(ionList), 4))
        ionInfo = np.zeros((len(ionList), 4, 6))  #ion[ion][vel][param]
        obsFlux = np.zeros((len(ionList), 2)) #obs[ion][coldens/error]
        obsFlux[0] = [14.2, 0.04] #CIV   - these are log10(N) and error on log10(N)
        obsFlux[2] = [12.85, 0.12] #SiIII
        obsFlux[1] = [13.5, 0.5] #CII
        #obsFlux[0] = [14.5, 0.2] #OVI
        #obsFlux[1] = [14.34, 0.2] #CIV
        #obsFlux[2] = [13.17, 0.2] #SiIII

        for i in range(len(ionList)):
            #open ionfile, find the velocities=0 b=1 N=2 q=3 area=4 and sigma=5
            #returns a 4 by 4 array, one row for each velocity = singleIon[vel][param]
            singleIon = openIonFile(ionList[i], run['Name'])#open ion file and get info for this run
            ionInfo[i] = singleIon

            for v in range(4):
                velBin = singleIon[v][0]
                #calculate the observed flux at this velocity for this ion

                #load obs data as velocities/normalized fluxes
                velocities, flux = convert_to_vel(ionList[i]['data_file'], ionList[i]['rest_wave'])
                #save flux and variance
                #make velocity bin negative to try to fit to the majority of the absorption profile
                vel_flux = np.interp(-1.0*velBin, velocities, flux)
                vel_var = (0.05)

                #obsFlux[i][v][0] = vel_flux
                #obsFlux[i][v][1] = vel_var

        #now that information is saved, do the chisquare analysis
        alphaBest_in = []
        alphaBest_chi = []
        fluxList = []
        obsList = []
        for v in range(4):

            alphaChi = []
            for alpha in alphas:
                sumChia = 0.0
                for i in range(len(ionList)):
                    estN = np.log10(alpha*ionInfo[i][v][4]) #log10(area) log10(average N)
                    chi = (obsFlux[i][0] - estN)**2/obsFlux[i][1]
                    if i==2:
                        fluxList.append(estN)
                        #obsList.append(np.exp(-1*(10**(obsFlux[i][0]))*tauTerm*ionList[i]['sigma']))
                    #print(np.exp(-1*(10**(obsFlux[i][0]))*c_speed*ionList[i]['sigma']))
                    sumChia = sumChia + chi
                alphaChi.append(sumChia)

            #find minmum of chiSquares
            minChia_index = np.argmin(alphaChi)
            alphaBest_in.append(minChia_index)
            alphaBest_chi.append(alphaChi[minChia_index])

            plt.plot(alphas, alphaChi)
            plt.yscale('log')
            plt.xscale('log')
            fig = plt.gcf()
            fig.savefig('test.png')

            #print/save results
            otherwrite.write(run['Name']+', '+str(v)+', '+str(alphas[minChia_index])+', '+str(alphaChi[minChia_index])+'\n')

        #find the average and maximum chiSquares across velocities
        #save this information to a file
        #writefile.write(run['Name']+', '+str(np.average(alphaBest_chi))+', '+str(betaAvg)+', '+str(np.max(alphaBest_chi))+', '+str(betaMax)+'\n')




if __name__ =="__main__":
    main()
