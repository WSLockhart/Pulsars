import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm  #color maps
import matplotlib.animation as animation
import sys #this is for accessing command line input

global Nt #number of time steps one rotation of the star is divided into
global Nen #number of energy bands the spectrum is divided into
global column1; global column2; global column3;

'''
Output.dat contains 3 columns. Depending on what's plotted, these are either:
[time, energy, flux] or [x-position , y-position, intensity].
'''

# # load global variables from definitions.h
# with open("../FullNS/definitions.h") as f:
#     content = f.readlines()
# # you may also want to remove whitespace characters like `\n` at the end of each line
# content = [x.strip() for x in content]
# print(content)


# load the data from the output.dat file into 3 array columns:
print('loading data...')
filename = "outputs/output_90-60_FLAT_PP.dat"
column1,column2,column3 = np.loadtxt(filename,skiprows=4,unpack=True)
print('data imported.')

#NOTE: can we pull these from 'definitions.h' instead?
Nt = 63    #NB: One less than the input in definitions (default 31)
Nen = 126  #NB: Two less than the input in definitions (default 126)

###############################################################################


def plotLightCurve(e):
    '''plot flux versus time at a given energy: I(t)|e '''

    # grab all the flux-time data for a given photon energy.
    # e.g. if Nen=128, we want data[e,128+e,256+e,384+e...]
    print('plotting light curve...')
    time_data,energy_data,flux_data = column1,column2,column3
    # we want to plot flux versus time _at a given energy_
    energy_value = energy_data[e]
    time_values = []
    flux_values = []
    for t in range(0,Nt-1):
        time_values.append(column1[Nen*t + e])
        flux_values.append(column3[Nen*t + e])
    plt.plot(time_values,flux_values)
    plt.title('energy = ' + str(energy_value) + ' KeV')
    plt.xlabel('Phase'); plt.ylabel('Flux (photons/cm^2/s/keV)')
    plt.show()


def plotPulseProfile():
    '''plot the total flux across all energy bands versus phase'''

    time_data,total_flux = column1,column3
    plt.plot(time_data,total_flux)
    plt.title('Pulse Profile')
    plt.xlabel('Phase'); plt.ylabel('Flux (photons/cm^2/s/keV)')
    plt.show()


def plotPulseProfiles():
    '''compare old and new profiles'''

    fileA = "outputs/old.dat"; fileB = "outputs/new.dat"

    c1a,c2a,c3a = np.loadtxt(fileA,skiprows=4,unpack=True)
    c1b,c2b,c3b = np.loadtxt(fileB,skiprows=4,unpack=True)
    time_data_A,total_flux_A = c1a,c3a
    time_data_B,total_flux_B = c1b,c3b

    plt.plot(time_data_A,total_flux_A,label="flat temp")
    plt.plot(time_data_B,total_flux_B,label="gralla dipole")
    plt.xlabel('Phase'); plt.ylabel('Flux (photons/cm^2/s/keV)')
    plt.title('Comparing dipole models')
    plt.legend()
    plt.show()


# def ppp(filename,spotsize):
#     c1,c2,c3 = np.loadtxt(filename,skiprows=4,unpack=True)
#     data = (c1,c3*(50/spotsize)**2)  # note, power is ** not ^
#     return data
#
# def michi():
#     '''flux scaled by (50º/pho)^2 to reproduce the plots in Michi's paper'''
#
#     fig = plt.figure(figsize=(12,5))
#     plt.subplot(1,2,1)
#     x,y = ppp("outputs/output_PP_m_50a.dat",50); plt.plot(x,y,label="rho=50")
#     x,y = ppp("outputs/output_PP_m_25a.dat",25); plt.plot(x,y,label="rho=25")
#     x,y = ppp("outputs/output_PP_m_10a.dat",10); plt.plot(x,y,label="rho=10")
#     # x,y = ppp("outputs/output_PP_m_2a.dat",2); plt.plot(x,y,label="rho=2")
#     plt.xlabel('Energy (keV)'); plt.ylabel('Scaled Flux (photons/cm^2/s/keV)')
#     plt.title('Pulse Profiles for i=90º, theta_s=80º')
#     plt.ylim(-10,400)
#     plt.legend()
#     plt.subplot(1,2,2)
#     x,y = ppp("outputs/output_PP_m_50b.dat",50); plt.plot(x,y,label="rho=50")
#     x,y = ppp("outputs/output_PP_m_25b.dat",25); plt.plot(x,y,label="rho=25")
#     x,y = ppp("outputs/output_PP_m_10b.dat",10); plt.plot(x,y,label="rho=10")
#     # x,y = ppp("outputs/output_PP_m_2b.dat",2); plt.plot(x,y,label="rho=2")
#     plt.xlabel('Phase'); plt.ylabel('Scaled Flux (photons/cm^2/s/keV)')
#     plt.title('Pulse Profiles for i=70º, theta_s=40º')
#     plt.ylim(-10,400)
#     plt.legend()
#     plt.show()

###############################################################################


def plotSpectrum(t):
    '''plot flux versus energy at a given time: I(e)|t '''

    #grab all the flux-energy data for a given phase (time)
    time_data,energy_data,flux_data = column1,column2,column3
    N = int(len(flux_data)/Nt) # N is the number of data points per time slice
    x = energy_data[N*t: N*(t+1)]
    y = flux_data[N*t: N*(t+1)]
    plt.plot(x,y)
    plt.title('phase = ' + str(t) + '/' + str(Nt))
    plt.xlabel('Energy (keV)'); plt.ylabel('Flux (photons/cm^2/s/keV)')
    # plt.legend()
    plt.show()


# def plotSpectra():
#     '''plot several spectra together on the same graph '''
#
#     times = (0,4,8,12,16)
#     for t in times:
#         x,y = getSpectrum(t)
#         plt.plot(x,y,label=t)
#         # plt.title('phase = ' + str(t))
#     plt.xlabel('Energy (keV)'); plt.ylabel('Flux (photons/cm^2/s/keV)')
#     plt.legend()
#     plt.show()


###############################################################################


def plotHotSpot():
    print('plotting hot spot...')
    myplot = plot(column1,column2,column3)
    plt.colorbar()
    plt.show()


def plot(x_pos,y_pos,intensity):
    x = np.array(x_pos)
    y = np.array(y_pos)
    I = np.array(intensity)
    # #figure out the shape of the data array:
    # data_size = len(intensity)
    # firstX = x_pos[0]
    # xnum = 0
    # for i in range(data_size):
    #     if x_pos[i] != firstX:
    #         break
    #     xnum += 1
    # ynum = int(data_size / xnum)
    # nrows, ncols = ynum, xnum
    # grid = I.reshape((nrows, ncols))
    # rangeI = np.max(I) - np.min(I)
    # frm = plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()),
    #        interpolation='nearest', cmap=cm.seismic, vmin=-rangeI,vmax=rangeI)
    max_value = np.max(I); min_value = np.min(I)
    # for val in I:  #BS for making cool plot with star in blue
    #     if val < 0:
    #         val = -max_value/10
    rangeI = max(max_value**2, min_value**2)**(1/2)
    frm = plt.scatter(x,y,c=I,cmap=cm.seismic,vmin=-0.035,vmax=0.035)
    #set vmin=0 for positive-only, vmin=-rangeI for total J^2 plots
    plt.axes().set_aspect('equal')
    return frm


def makeVideo():
    '''create a video of the hotspot as the star rotates'''

    num_periods = 3
    x_pos,y_pos,intensity = column1,column2,column3
    fig = plt.figure()
    movie_frames = []
    N = int(len(intensity)/Nt) # N is the number of data points per time slice
    print('creating video frames...')
    for cycle in range(num_periods):
        for t in range(Nt):
            print('t = ' + str(t+1) + '/' + str(Nt))
            # separate the data by time slice and plot
            newframe = plot(x_pos[N*t:N*(t+1)],y_pos[N*t:N*(t+1)],intensity[N*t:N*(t+1)])
            movie_frames.append([newframe])
    ani = animation.ArtistAnimation(fig, movie_frames, interval=150, blit=True, repeat_delay=1000)
    print('saving video...')
    videoname = filename + '.mp4'
    ani.save(videoname)
    print('complete.')
    # plt.show()


# ###############################################################################
# f = open('test_data.txt',"r")
# contents = f.read()
# print(contents)
#
# f = open("write.txt","w+")
# for i in range(10):
#      f.write("This is line %d\r\n" % (i+1))
# f.close()
