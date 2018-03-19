import numpy as np

blue1 = '#4c95cb'
blue2 = '#0068b6'
blue3 = '#00487f'
green1 = '#4ca85c'
green2 = '#008317'
green3 = '#005b10'
yellow1 = '#f4bd54'
yellow2 = '#f0a10b'
yellow3= '#a87007'
red1 = '#df4c4c'
red2 = '#d20000'
red3 = '#930000'
##### Runs to have spectra generated #######
run1 = { 'Name':'T0.3_v1000_chi300_cond',
        'Name_plot':'M3.6-v1000-T0.3-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.8,
        'tcc':1.7,
        'f_list':['0013', '0038', '0080', '0132'],
        'marker':'^',
        'color':green2}
run2 = { 'Name':'T0.3_v1700_chi300_cond',
        'Name_plot':'M6.5-v1700-T0.3-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':6.5,
        'tcc':1.0,
        'f_list':['0003', '0020', '0046', '0078'],
        'marker':'^',
        'color':yellow1}
run3 = { 'Name':'T0.3_v3000_chi300_cond',
        'Name_plot':'M11.4-v3000-T0.3-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':11.4,
        'tcc':0.56,
        'f_list':['0001', '0004', '0014', '0035'],
        'marker':'^',
        'color':red1}
run4 = { 'Name':'T3_v3000_chi3000_cond',
        'Name_plot':'M3.6-v3000-T3-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.6,
        'tcc':1.8,
        'f_list':['0001', '0004', '0007', '0010'],
        'marker':'^',
        'color':red2}
run5 = { 'Name':'T3_v860_chi3000_cond',
        'Name_plot':'M1.0-v860-T3-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.2,
        'f_list':['0001', '0003', '0006', '0010'],
        'marker':'^',
        'color':blue2}
run6 = { 'Name':'T1_v1700_chi1000_cond',
        'Name_plot':'M3.5-v1700-T1-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':3.5,
        'tcc':1.8,
        'f_list':['0002', '0010', '0017', '0028'],
        'marker':'^',
        'color':yellow2}
run7 = { 'Name':'T1_v480_chi1000_cond',
        'Name_plot':'M1.0-v480-T1-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.4,
        'f_list':['0002', '0009', '0017', '0031'],
        'marker':'^',
        'color':blue1}
run8 = { 'Name':'T10_v1500_chi10000_cond',
        'Name_plot':'M1.0-v1500-T10-c',
        'Dir':'../../Blob_paper2/Files/',
        'Mach':1.0,
        'tcc':6.5,
        'f_list':['0001', '0002', '0004', '0008'],
        'marker':'^',
        'color':yellow3}
run9 = { 'Name':'T0.3_v1000_chi300',
        'Name_plot':'M3.8-v1000-T0.3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.8,
        'f_list':['0025', '0033', '0042', '0058'],
        'marker':'o',
        'color':green2}
run10 = { 'Name':'T0.3_v1700_chi300',
        'Name_plot':'M6.5-v1700-T0.3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.5,
        'f_list':['0022', '0032', '0053', '0085'],
        'marker':'o',
        'color':yellow1}
run11 = { 'Name':'T0.3_v3000_chi300',
        'Name_plot':'M11.4-v3000-T0.3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':11.4,
        'f_list':['0028', '0044', '0065', '0110'],
        'marker':'o',
        'color':red1}
run12 = { 'Name':'T3_v3000_chi3000',
        'Name_plot':'M3.6-v3000-T3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.6,
        'f_list':['0021', '0030', '0040', '0062'],
        'marker':'o',
        'color':red3}
run13 = { 'Name':'T3_v430_chi3000',
        'Name_plot':'M0.5-v430-T3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':0.5,
        'f_list':['0010', '0011', '0016', '0024'],
        'marker':'o',
        'color':blue2}
run14 = { 'Name':'T3_v860_chi3000',
        'Name_plot':'M1.0-v860-T3',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0010', '0018', '0030', '0038'],
        'marker':'o',
        'color':blue3}
run15 = { 'Name':'T1_v1700_chi1000',
        'Name_plot':'M3.5-v1700-T1',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':3.5,
        'f_list':['0021', '0029', '0038', '0052'],
        'marker':'o',
        'color':yellow2}
run16 = { 'Name':'T1_v3000_chi1000',
        'Name_plot':'M6.2-v3000-T1',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':6.2,
        'f_list':['0022', '0032', '0048', '0095'],
        'marker':'o',
        'color':red2}
run17 = { 'Name':'T1_v480_chi1000',
        'Name_plot':'M1.0-v480-T1',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0009', '0013', '0024', '0035'],
        'marker':'o',
        'color':blue1}
run18 = { 'Name':'T10_v1500_chi10000',
        'Name_plot':'M1.0-v1500-T10',
        'Dir':'../../Blob_paper1/Files/',
        'Mach':1.0,
        'f_list':['0014', '0021', '0029', '0044'],
        'marker':'o',
        'color':yellow3}
run19 = { 'Name':'HC_v1000_chi300_cond',
        'Name_plot':'M3.8-v1000-T0.3-hc',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.8,
        'f_list':['0054', '0060', '0080', '0107'],
        'marker':'s',
        'color':green2}
run20 = { 'Name':'HC_v1700_chi1000_cond',
        'Name_plot':'M3.5-v1700-T1-hc',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.5,
        'f_list':['0024', '0050', '0082', '0083'],
        'marker':'s',
        'color':yellow2}
run21 = { 'Name':'HC_v3000_chi3000_cond',
        'Name_plot':'M3.6-v3000-T3-hc',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':3.6,
        'f_list':['0007', '0015', '0026', '0049'],
        'marker':'s',
        'color':red2}
run22 = { 'Name':'LowCond_v1700_chi300_cond',
        'Name_plot':'M6.5-v1700-T0.3-lc',
        'Dir':'../../Blob_paper3/Files/',
        'Mach':6.5,
        'f_list':['0016', '0075', '0115', '0184'],
        'marker':'.',
        'color':yellow2}

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

def writeFloat(number, prec):
    newnum = format(float(number), prec)
    newString = '$'+newnum+'$'
    return newString

tableFileLoc = '../rankNum/'
tableName = 'rankNum_outflowApp_guass.txt'
#tableName = 'rankNum_cgm_actualBestFit.txt'
tableFile = open(tableFileLoc+tableName, 'r')
Latex = open(tableFileLoc+tableName[:-4]+'Latex.txt', 'w')
p = '.3f'

i = 1
for line in tableFile:
    if i==1:
        ls = line.replace(', ', ' & ')
        writeString = ls
    else:
        ls = line.split(', ')
        if ls[0][-4:]=='cond':
            for run in runList:
                if run['Name']==ls[0]:
                    writeName = run['Name_plot']

            writeString = writeName+' & '+writeFloat(ls[1], p)+' & '+ls[2]+' & '+writeFloat(ls[3], p)+' & '+ls[4]+' & '+writeFloat(ls[5], p)+' & '+ls[6]
            #writeString = writeName+' & '+ls[1]+' & '+writeFloat(ls[2], p)+' & '+writeFloat(ls[3], p)+'\n'

        else:
            for run in runList:
                if run['Name']==ls[0]:
                    writeName = run['Name_plot']

            writeString = writeName+' & '+writeFloat(ls[1], p)+' & '+writeFloat(ls[2], p)+' & '+writeFloat(ls[3], p)+' & '+writeFloat(ls[4], p)+' & '+writeFloat(ls[5], p)+' & '+writeFloat(ls[6], p)+'\n'
           # writeString = writeName+' & '+ls[1]+' & '+writeFloat(ls[2], p)+' & '+writeFloat(ls[3], p)+'\n'

    Latex.write(writeString)
    i = i+1

tableFile.close()
Latex.close()
