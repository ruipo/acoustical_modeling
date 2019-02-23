
# Make and run OASN files for two varying bottoms. 
# To Do: define bottom paramters and intervals
#        define outFileName
#        define water sound speed layers
#        input settings for each block


import numpy as np
import types
import os
import subprocess
import time


# define options
options = 'N J 1 R'

# define frequencies
freqs = np.array([100,100,1,0]) # [min freq, max freq, # of freqs, integration contour offset]

# define sound speed layers
layers = np.array([[0., 0., 0., 0., 0., 0., 0.], # [depth, comp. speed, shear speed, alpha_p, alpha_s, density, RMS roughness]
                   [0., 1451., 0., 0., 0., 1., 0.]])


# define bottoms - each column is a different iteration of the subbottom; make an array for each subbottom layer
num_bottoms = 2 # number of bottoms
num_b1 = 3 # number of bottom1 types
num_b2 = 3 # number of bottom2 types

bottom1 = np.array([102.5*np.ones(num_b1), # depth
                    np.linspace(1500.,1600.,num_b1), # c_p
                    -1750.*np.ones(num_b1), # c_s
                    0.35*np.ones(num_b1), # alpha_p
                    np.zeros(num_b1), # alpha_s
                    1.75*np.ones(num_b1), # density
                    np.zeros(num_b1)]) # RMS of interface roughness

bottom2 = np.array([200*np.ones(num_b2), # depth
                    np.linspace(1700.,1800.,num_b2), # c_p
                    np.zeros(num_b2), # c_s
                    0.35*np.ones(num_b2), # alpha_p
                    np.zeros(num_b2), # alpha_s
                    1.75*np.ones(num_b2), # density
                    np.zeros(num_b2)]) # RMS of interface roughness

# define array
array = types.SimpleNamespace()

array.n = 15 # number of array elements
array.dist = 5 # array element spacing
array.L = array.dist * (array.n-1) # total array aperature

#vertical array
array.depth = np.round(np.linspace(10,80,array.n)) # depth of each array element
array.x = 0. # x-coord of each array element
array.y = 0. # y-coord of each array element

array.type = 1 # element type: OASES type 1 is hydrophone
array.gain = 0 # element gain (dB)


# define sources
src = types.SimpleNamespace()

src.src = np.array([0, 50, 0, 1]) # [surf. noise src strength, white noise level, deep src level, # of discr. srcs]

if src.src[0] != 0: #if there is sea surface noise
    src.seacs = np.array([1400,1e8]) # [cmins, cmaxs]
    src.seawav = np.array([400,400,100]) # [samples in cont., discr, evanes.]

if src.src[2] != 0: #if there is deep noise source
    src.dscs = np.array([1400,1e8]) # [cmins, cmaxs]
    src.dswav = np.array([400,400,100]) # [samples in cont., discr, evanes.]
    
if src.src[3] != 0: #if there is discrete source    
    src.ndsrc = np.array([50,1,0,120]) # [depth(m), x-coord(km), y-coord(km), source level]
    src.ndcs = np.array([1400,1e8]) # [cmins, cmaxs]
    src.ndwav = np.array([-1,0,0]) # [# of sampling points,first sampling point, last sampling point]; [-1,0,0] for auto
    
    #repeat for each discrete source
    
# define signal replicas  

if 'R' in options:
    reps = types.SimpleNamespace()

    reps.zs = np.array([50,50,1]) # depths of sampling of replicas
    reps.xs = np.array([1,1,1]) # x-coord of sampling 
    reps.ys = np.array([0,0,1]) # y-coord of sampling

    reps.cvals = np.array([1400,1e8]) # [cmins, cmaxs]
    reps.wavs = np.array([-1,0,0]) # [# of sampling points,first sampling point, last sampling point]; [-1,0,0] for auto

ii = 0

# define function to write a line
def line(var):
    s = ''
    for ii in range(len(var)):
        sii = '{'+str(ii)+'} '
        s = ''.join([s,sii])
    return s[0:-1].format(*var)


### start writing input files
allDirs = []

for ii in range(bottom1.shape[1]): # for each column in bottom1
    for jj in range(bottom2.shape[1]): # for each column in bottom2
        
        outFileName = 'b1' + str(ii+1) + '_' + 'b2' + str(jj+1) + '.dat'

        if not os.path.isdir(outFileName[0:-4]):
            os.mkdir(outFileName[0:-4])

        allDirs.append(outFileName[0:-4])

        outFile = open(outFileName[0:-4] + '/' + outFileName, 'w')

        # Begin writing file:

        # Block 1: Title
        outFile.write('workshop case ' + outFileName[0:-4] + '\n')

        # Block 2: Options
        outFile.write(options + '\n')

        # Block 3: Frequencies
        outFile.write(line(freqs) + '\n')

        outFile.write('\n')

        # Block 4: Environment

        #number of layers
        outFile.write(str(layers.shape[0] + num_bottoms) + '\n')

        #water layers
        layers[np.isnan(layers)] = 0.
        for kk in range(layers.shape[0]):
            outFile.write(line(layers[kk,:]) + '\n')

        #bottom layers
        outFile.write(str(line(bottom1[:,ii]) + '\n'))
        outFile.write(str(line(bottom2[:,jj]) + '\n'))

        outFile.write('\n')

        # Block 5: Array
        outFile.write(str(array.n) + '\n')
        for irec in range(len(array.depth)):
            outFile.write(line(np.array([array.depth[irec],array.x,array.y])))
            outFile.write(' ' + str(array.type))
            outFile.write(' ' + str(array.gain) + '\n')    

        outFile.write('\n')

        # Block 6: Sources
        outFile.write(line(src.src) + '\n')
        
        outFile.write('\n')

        # Block 7: Sea Surface
        if src.src[0] != 0: #if there is sea surface noise
            outFile.write(line(src.seacs) + '\n')
            outFile.write(line(src.seawav) + '\n')
            
            outFile.write('\n')

        # Block 8: Deep Noise
        if src.src[2] != 0: #if there is deep noise source
            outFile.write(line(src.dscs) + '\n')
            outFile.write(line(src.dswav) + '\n')
            
            outFile.write('\n')            
            
        # Block 9: Discrete Sources
        if src.src[3] != 0: #if there is discrete source   
            outFile.write(line(src.ndsrc) + '\n')
            outFile.write(line(src.ndcs) + '\n')
            outFile.write(line(src.ndwav) + '\n')

            outFile.write('\n')

        # Block 10: Signal replicas
        if 'R' in options:
            outFile.write(line(reps.zs[0:-1]) + ' ' + str(int(reps.zs[-1])) + '\n')
            outFile.write(line(reps.xs[0:-1]) + ' ' + str(int(reps.xs[-1])) + '\n')
            outFile.write(line(reps.ys[0:-1]) + ' ' + str(int(reps.ys[-1])) + '\n')

            outFile.write('\n')

            outFile.write(line(reps.cvals) + '\n')
            outFile.write(line(reps.wavs) + '\n')


# Run OASN

baseDir = os.getcwd()

os.chdir(baseDir)
if not os.path.isdir('chk_files'):
    os.mkdir('chk_files')

proc_list = []
logs = []

# print('Checking for OASN')
oasn_check = subprocess.Popen('which  oasn'.split())
oasn_check.wait()

# print(str(oasn_check.poll()==0))
# print(str(oasn_check.poll()))

if oasn_check.returncode == 0:
    print('OASN found!')
    for ii in range(len(allDirs)):
        os.chdir(baseDir + '/' + allDirs[ii])
        
        logs = open('oasn_log_'+str(ii+1)+'.txt','w')
        print(('oasn {'+str(ii)+'}').format(*allDirs))
        proc = subprocess.Popen(('oasn {'+str(ii)+'}').format(*allDirs).split(),stdout=logs,stderr=logs)
        proc.wait()
        print('  Return code: ' + str(proc.returncode))

        if not proc.returncode == 0:
            print('Trying once more: ')

            print(('oasn {'+str(ii)+'}').format(*allDirs))
            proc = subprocess.Popen(('oasn {'+str(ii)+'}').format(*allDirs).split(),stdout=logs,stderr=logs)
            proc.wait()
            print('  Return code: ' + str(proc.returncode))
	
        if proc.returncode == 0:
            os.rename(baseDir + '/' + allDirs[ii] + '/' + allDirs[ii] + '.chk', baseDir + '/chk_files/' + allDirs[ii] + '.chk')

        logs.close()
        time.sleep(0.05)

    
os.chdir(baseDir)
print('Done!')

