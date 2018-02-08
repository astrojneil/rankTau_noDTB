import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
kpc = 3.086e21
c = 3.e5

#function to read in spectrum, convert to velocity and return
#wavelength and flux numpy arrays
def convert_to_vel(ion, rest_wave):
    #load spectrum; will need to add a case when this is a txt file rather than h5
    if ion=='OVI':
        filename = '../Files/S1226-o6-forJNeil'
    else:
        filename = '../Files/S1226-redward-forJNeil'
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
    vel = c*(wavelength/rest_wave -1)

    return vel, flux

def openRun(runName, filename):
    openfile = open(filename, 'r')
    singleIon = np.zeros((4,3))
    vel = 0
    for line in openfile:
        ls = line.split(', ')
        if ls[0] == runName:
            singleIon[vel][0] = float(ls[1]) #vel
            singleIon[vel][1] = float(ls[2]) #alpha
            if runName[-5:] != '_cond':
                singleIon[vel][2] = float(ls[3]) #b
            vel = vel+1

    return singleIon

def openIonFile(ion, runName):
    openfile = open('../rankTau'+ion['ionfolder']+ion['ionfolder'][1:-1]+'_bestFitParameters_3.txt', 'r')
    singleIon = np.zeros((4,4))
    vel = 0
    for line in openfile:
        ls = line.split(', ')
        if ls[0] == runName:
            singleIon[vel][0] = float(ls[1]) #velocity
            singleIon[vel][1] = float(ls[2]) #tau
            singleIon[vel][2] = float(ls[5]) #q
            singleIon[vel][3] = float(ls[8]) #area
            vel = vel+1

    return singleIon

def model(t, b, x):
    return t*(0.01)/(1.01-x**b)

### dictionaries of ion info
ion1 = {'ion':'OVI',
        'fieldname':'O_p5_number_density',
        'ionfolder': '/OVI/',
        'rest_wave': 1031.91,
        'data_file': '../Files/S1226-o6-forJNeil',
        'sigma': 1.1776e-18,
        'massNum': 16.0}
ion2 = {'ion':'CIV',
        'fieldname':'C_p3_number_density',
        'ionfolder': '/CIV/',
        'rest_wave': 1548.18,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 2.5347e-18,
        'massNum': 12.0}
ion3 = {'ion':'NV',
        'fieldname':'N_p4_number_density',
        'ionfolder': '/NV/',
        'rest_wave': 1242.8,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 8.3181e-19,
        'massNum': 14.0}

ion4 = {'ion':'CII',
        'fieldname':'C_p1_number_density',
        'ionfolder': '/CII/',
        'rest_wave': 1335.66,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 1.4555e-19,
        'massNum': 12.0}

ion9 = {'ion': 'Si IV',
        'fieldname': 'Si_p3_number_density',
        'ionfolder':'/SiIV/',
        'rest_wave': 1402.77,
        'data_file': '../Files/S1226-redward-forJNeil',
        'sigma': 3.0694e-18,
        'massNum': 28.0}

ionList = []
ionList.append(ion1)
ionList.append(ion2)
ionList.append(ion3)
ionList.append(ion4)
ionList.append(ion9)


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

bestFitFile = '../rankTau/ranktau_outflow_actualBestFit.txt'

for run in runList:
    runInfo = openRun(run['Name'], bestFitFile)

    f, ax = plt.subplots(5, 1, sharey=True, sharex=True)

    for i in range(len(ionList)):
        singleIon = openIonFile(ionList[i], run['Name'])
        estFlux_a = []
        estFlux_b = []
        fitvel = []
        for v in range(4):
            if run['Name'][-5:] !='_cond':
                def modelFun(x):
                    return np.exp(-runInfo[v][2]*model(singleIon[v][1], singleIon[v][2], x))
                Flux_b = scipy.integrate.quad(modelFun, 0.0, 1.0)[0]
                estFlux_b.append(Flux_b)

            Flux_a = np.exp(-runInfo[v][1]*singleIon[v][3])

            fitvel.append(-runInfo[v][0])
            estFlux_a.append(Flux_a)



        obsvel, obsflux = convert_to_vel(ionList[i]['ion'], ionList[i]['rest_wave'])
        ax[i].plot(obsvel, obsflux)
        ax[i].scatter(fitvel, estFlux_a, color = 'red')
        if run['Name'][-5:] !='_cond':
            ax[i].scatter(fitvel, estFlux_b, color = 'green')

        ax[i].set_xlim(-600, 200)
        ax[i].set_title(ionList[i]['ion'])

    ax[0].set_title('BestFit '+run['Name'])
    plt.ylim(0, 1.5)
    plt.xlabel('Velocity (km/s)')
    ax[1].set_ylabel('Normalized Flux')


    fig = plt.gcf()
    fig.set_size_inches(6.5, 10.5)
    fig.savefig('outflowApp/testOutflowApp_'+run['Name']+'.png')

    plt.close()
