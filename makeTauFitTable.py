import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import scipy
import yt
import h5py


##cloud pixels!!
cloudpix = 7880

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

        vx_cloud = (avg_vx_cloud.value)
        velList.append(vx_cloud)

    #print('Frame vels:')
    #print(velFrames)
    print('Cloud vels:')
    print(velList)
    return velList #, velFrames

#define functions for the emcee best fit
def model(t, b, x):
    return t*(0.01)/(1.01-x**b)

def lnprior(theta, maxtau):
    t, b = theta
    if 0.0 < t < maxtau and 0.0 < b < 1.e2:
        return 0.0
    return -np.inf

def lnlike(theta, taus, frac, err):
    t, b = theta
    model_test = model(t, b, frac)
    return -0.5*(np.sum((taus - model_test)**2/err**2))

def lnprob(theta, taus, frac, err, maxtau):
    lp = lnprior(theta, maxtau)
    if np.isfinite(lp):
        return lp + lnlike(theta, taus, frac, err)
    return -np.inf

#function to open the ranked tau file
def openRankedFile(filename):
    openfile = open(filename, 'r')
    rankedTaus = []
    for line in openfile:
        rankedTaus.append(float(line))
    openfile.close()
    return np.array(rankedTaus)

def runemcee(taus, pixelFrac):
    ndim, nwalkers = 2, 200
    burnin = 50
    r = np.zeros(ndim)+1.e-3
    pos = [r + 1.e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(taus, pixelFrac, 1.e-2, 2*max(taus)))

    print("Running MCMC...")
    sampler.run_mcmc(pos, 1000, rstate0=np.random.get_state())
    print("Done")

    samples = sampler.chain[:, burnin:,:].reshape((-1, ndim))
    tau_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis = 0)))

    def modelFun(x):
        return model(tau_mcmc[0], b_mcmc[0], x)

    model_area = scipy.integrate.quad(modelFun, 0.0, 1.0)

    return tau_mcmc, b_mcmc, model_area


def main():
##### Runs to work add to table #######
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

    ### dictionaries of ion info
    ion1 = {'ion':'O VI',
            'fieldname':'O_p5_number_density',
            'ionfolder': '/OVI/'}
    ion2 = {'ion':'Mg II',
            'fieldname':'Mg_p1_number_density',
            'ionfolder': '/MgII/'}
    ion3 = {'ion':'N V',
            'fieldname':'N_p4_number_density',
            'ionfolder': '/NV/'}
    ion4 = {'ion':'C IV',
            'fieldname':'C_p3_number_density',
            'ionfolder': '/CIV/'}
    ion5 = {'ion':'Si III',
            'fieldname':'Si_p2_number_density',
            'ionfolder': '/SiIII/'}
    ion6 = {'ion':'Si IV',
            'fieldname':'Si_p3_number_density',
            'ionfolder': '/SiIV/'}
    ion7 = {'ion':'Ne VII',
            'fieldname':'Ne_p6_number_density',
            'ionfolder': '/NeVII/'}
    ion8 = {'ion':'H I',
            'fieldname':'H_p0_number_density',
            'ionfolder': '/HI/'}
    ion9 = {'ion':'C II',
            'fieldname':'C_p1_number_density',
            'ionfolder': '/CII/'}
    ion10 = {'ion':'C III',
            'fieldname':'C_p2_number_density',
            'ionfolder': '/CIII/'}


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


    #now the real work
    for ion in ionList:
        print('Started '+ion['ionfolder'][1:-1])
        ionTable = open('../rankTau_noDTB'+ion['ionfolder']+ion['ionfolder'][1:-1]+'_bestFitParameters_3.txt', 'w')
        #write column headers
        ionTable.write('Run, velocity(cm/s), tau_fit, upper_tau_err, lower_tau_err, q_fit, upper_q_err, lower_q_err, area, area_err\n')

        #define priors here? Not great to have answers depend on the priors for each ion...

        for run in runList:
            print(run['Name'])
            velBins = findCloudVel(run['Dir'], run['Name'], run['f_list'])  #velocities are in cm/s
            print(velBins)

            for v in range(len(velBins)): #velocity = velBins[v]
                taus = openRankedFile('../rankTau_noDTB'+ion['ionfolder']+'rankTau'+run['Name']+'_v'+str(int(round(velBins[v], 0)))+'.txt')
                totalpixel = len(taus)
                pixelNum = np.arange(0.0, totalpixel, 1.0)
                pixelFrac = pixelNum/totalpixel
                clip = totalpixel-cloudpix
                print(ion['ion']+': '+run['Name']+': '+str(int(round(velBins[v], 0)))+':')

                pixelFrac_new = (pixelNum[clip:]-clip)/len(pixelNum[clip:])
                taus_new = taus[clip:]
                t_fit, b_fit, area_fit = runemcee(taus_new, pixelFrac_new)

                runinfo = run['Name']+', '+str(int(round(velBins[v], 0)))
                tauinfo = str(t_fit[0])+', '+str(t_fit[1])+', '+str(t_fit[2])
                binfo = str(b_fit[0])+', '+str(b_fit[1])+', '+str(b_fit[2])
                areainfo = str(area_fit[0])+', '+str(area_fit[1])
                writeString = runinfo+', '+tauinfo+', '+binfo+', '+areainfo+'\n'

                ionTable.write(writeString)
        print("Finished with ion: "+ion['ionfolder'][1:-1])
        ionTable.close()


if __name__ =="__main__":
    main()
