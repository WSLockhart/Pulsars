import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm  #color maps
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

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
#

# load the data from the output.dat file into 3 array columns:
print('loading data... ',end="\r")
sourcefile = "LC_1.50_11.0_ 12_090_075_0.35_0.0"
filename = "FINAL_PLOTS/FIG2-3_DipoleLightCurves/10Hz/" + sourcefile + ".dat"
column1,column2,column3 = np.loadtxt(filename,skiprows=0,unpack=True)
# filename2 = "FINAL_PLOTS/FIG8_QDLightCurve/80-60/HS_80-60-1.5_data/" + sourcefile + ".dat"
# mov1,mov2,mov3 = np.loadtxt(filename2,skiprows=0,unpack=True)
print('data imported.')

Nt = 64     #NB: One less than the input in definitions (default 64)
Nen = 126  #NB: Two less than the input in definitions (default 126)

###############################################################################


def temp_contours():  #FIGURE 7
    '''contour plot of the temperature ratio between the two hemispheres.
    Data from jtemp.c. For Figure 7 in the paper. '''

    Qlist = np.linspace(0,2,21)
    Zlist = np.linspace(0,90,46)  #NB: these are hard-coded
    Z,Q = np.meshgrid(Zlist,Qlist)
    T = column3[0:21*46]
    grid = T.reshape((21, 46))

    levels = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    cplot = plt.contour(Q,Z,grid,levels, colors='k')
    # plt.clabel(cplot, fmt = '%2.1f',inline=True, fontsize=11)
    cfplot = plt.contourf(Q,Z,grid,levels, cmap="Blues")
    cbar=plt.colorbar(cfplot)
    cbar.set_label(r'$T_1 / T_2$',rotation=270, fontsize=12, labelpad=16)
    plt.xlabel('q')
    plt.ylabel(r'$\zeta$ (degrees)')
    plt.show()


def plotPP():
    '''plot the total flux across all energy bands versus phase'''
    time_data,total_flux = column1,column3
    plt.plot(time_data,total_flux)
    plt.show()


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

    x1 = np.append(column1,1);
    y1 = np.append(column3,column3[0]) #repeat first value no there's no space
    # peaks = calcPeakData(y1)
    # time_data = x1
    # flux_data = y1

    time_data = np.append(x1,x1+1) #plot two periods
    flux_data = np.append(y1,y1)

    average_flux = np.sum(column3)/len(column3)
    flux_data *= 1/(average_flux) #normalize

    # amp = round(np.max(total_flux)/np.min(total_flux),2)
    plt.plot(time_data,flux_data,linewidth=2.2,color='#6B4C9A')#,label='Amplitude = '+str(amp))

    plt.xlim(0,2)
    # plt.ylim(0.5,1.5)
    plt.ylim(0,np.max(flux_data)*1.1)
    # plt.axes().set_aspect(1.1) #for FIG9
    # plt.set_xticks([0,0.5,1,1.5,2])
    # plt.set_yticks([0.5,0.75,1,1.25,1.5])
    plt.tick_params(direction='in')
    plt.axvline(x=1,color='gray',linestyle='dashed',linewidth=1)
    # plt.title(sourcefile)
    plt.xlabel('Phase')
    # plt.ylabel('Total Flux (photons/cm^2/s)')
    plt.ylabel('Normalized Flux')
    # for p in peaks:
    #     plt.axvline(x=p/Nt,color='green',linewidth=0.75)
    #     plt.axvline(x=1+(p/Nt),color='green',linewidth=0.75)
        # plt.plot(p/Nt,total_flux[p],color='red',marker='o',markersize=2.5)
    # if len(peaks) == 2:
    #     t1 = peaks[0]; t2 = peaks[1]
    #     A = max(column3[t1],column3[t2])/min(column3[t1],column3[t2]) #ratio of two peaks
    # plt.title(r'$\Delta \phi$ = ' + str(round(phase_shift,2)) + '; $A_1/A_2$ = ' + str(round(A,2)),color='green', fontsize=14)
    # plt.legend(loc=4)
    plt.show()



def calcPeakData(flux_data):
    '''returns the phase shift and amplitude difference of the two peaks'''

    print('spread = ' , str(round(np.max(flux_data)/np.min(flux_data),3)) )
    peaks = []
    for t in range(Nt):
        #it's a local peak if it's bigger than its neighbors. Note the mod (%)
        if (flux_data[t]>flux_data[(t-1)%Nt]) and (flux_data[t]>flux_data[(t+1)%Nt]):
            peaks.append(t)
    if len(peaks) == 2:
        t1 = peaks[0]; t2 = peaks[1]
        peaktopeak = flux_data[t1] / flux_data[t2]
        phase_shift = abs((t1-t2)/Nt)
        if phase_shift > 0.5:  #the real phase shift can't exceed 180º=0.5, if we wrap around
            phase_shift = 1-phase_shift
        print('phase shift = ' + str(round(phase_shift,3)) + ', peak-to-peak ratio = ' + str(round(peaktopeak,3)))
    elif len(peaks) == 1:
        print('only 1 peak detected: amplitude = ' + str(round(flux_data[peaks[0]],2)) +
            ', phase = ' + str(round(peaks[0]/Nt,2)))
    else: print('more than two peaks detected!')
    return peaks


# def ppp(filename,spotsize):
#     c1,c2,c3 = np.loadtxt(filename,unpack=True)
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
    myplot = plotHS(mov1,mov2,mov3)
    # max_value = np.max(mov3); min_value = np.min(mov3)
    # rangeI = max(max_value**2, min_value**2)**(1/2)
    # plt.scatter(mov1,mov2,c=mov3,cmap=cm.seismic,vmin=-rangeI,vmax=rangeI,s=0.3)
    # plt.colorbar()
    radius = 6.45 #NB:hard-coded, need to double-check each time! Here use xlim=ylim=7.5
    # radius = 7.25 #NB:hard-coded, need to double-check each time!   Here use xlim=ylim=8.43
    circle1=plt.Circle((0,0),radius,color='black',fill=False)
    plt.axes().add_artist(circle1)
    plt.show()


def plotHS(x_pos,y_pos,intensity):
    x = np.array(x_pos)
    y = np.array(y_pos)
    I = np.array(intensity)
    max_value = np.max(I); min_value = np.min(I)
    print('max_value='+str(max_value),'min_value='+str(min_value))
    rangeI = max(max_value**2, min_value**2)**(1/2)
    frm = plt.scatter(x,y,marker=",",c=I,cmap=cm.gist_heat,vmin=0,vmax=max_value,s=0.2)#vmin=-rangeI,vmax=rangeI,s=0.3)
    #set vmin=0 for positive-only, vmin=-rangeI for total J^2 plots
    plt.xlim(-7.5,7.5)
    plt.ylim(-7.5,7.5)
    plt.axes().set_aspect('equal')
    # plt.axis('off')
    plt.xticks([]); plt.yticks([])  #box but no tick marks
    plt.axhline(y=7.5,linewidth=3, color="gray")
    plt.axhline(y=-7.5,linewidth=3, color="gray")        # inc. width of x-axis and color it green
    plt.axvline(x=7.5,linewidth=3, color="gray")        # inc. width of x-axis and color it green
    plt.axvline(x=-7.5,linewidth=3, color="gray")        # inc. width of x-axis and color it green
    return frm


# def makeVideo():
#     '''create a video of the hotspot as the star rotates'''
#
#     num_periods = 1
#     x_pos,y_pos,intensity = mov1,mov2,mov3
#     fig = plt.figure()
#     movie_frames = []
#     N = int(len(intensity)/Nt) # N is the number of data points per time slice
#     print('creating video frames...')
#     for cycle in range(num_periods):
#         for t in range(Nt):
#             print('t = ' + str(t+1) + '/' + str(Nt),end='\r')
#             # separate the data by time slice and plot
#             newframe = plotHS(x_pos[N*t:N*(t+1)],y_pos[N*t:N*(t+1)],intensity[N*t:N*(t+1)])
#             movie_frames.append([newframe])
#     ani = animation.ArtistAnimation(fig, movie_frames, interval=150, blit=True, repeat_delay=1000)
#     print('saving video...')
#     videoname = filename + '.mp4'
#     ani.save(videoname)
#     print('complete.')
#     # plt.show()


# import matplotlib.gridspec as gridspec
# figgy = plt.figure(figsize=(8,8))
# gs = gridspec.GridSpec(2,1,height_ratios=[2,1.2])
# ax2 = plt.subplot(gs[0])
# ax1 = plt.subplot(gs[1])
# N = int(len(mov3)/Nt)
# # peaks = calcPeakData(column3)

def animate(i):
    print('creating frame ' + str(i) + ' ... ', end="\r")
    t = i%Nt
    ax1.clear()
    ax2.clear()

    #PLOT1: pulse profile
    ax1.plot(column1,column3)
    ax1.axvline(x=(t/Nt),color='red',linewidth=0.75)
    ax1.plot(t/Nt,column3[t],color='red',marker='o',markersize=2.5)
    #show peaks
    # for p in peaks: ax1.axvline(x=p/Nt,color='green',linewidth=0.5)
    # if len(peaks) == 2:
    #     t1 = peaks[0]; t2 = peaks[1]
    #     A = min(column3[t1],column3[t2])/max(column3[t1],column3[t2]) #ratio of the two peaks amplitudes
    #     phase_shift = abs((t1-t2)/Nt)
    #     #the real phase shift can't exceed 180º=0.5, if we wrap around
    #     if phase_shift > 0.5: phase_shift = 1-phase_shift
    #     # print('phase shift = ' + str(round(phase_shift,2)) + ', peak-to-peak = ' + str(round(peaktopeak,2)))
    #     ax1.set_title(r'$\Delta \phi$ = ' + str(round(phase_shift,3)) + '; $A_1/A_2$ = ' + str(round(A,3)),color='green', fontsize=14)

    #PLOT2: hot spot
    x = mov1[N*t:N*(t+1)]
    y = mov2[N*t:N*(t+1)]
    I = mov3[N*t:N*(t+1)]
    max_value = np.max(mov3); min_value = np.min(mov3)
    rangeI = max(max_value**2, min_value**2)**(1/2)
    ax2.scatter(x,y,c=I,cmap=cm.seismic,vmin=-rangeI,vmax=rangeI,s=0.2)

    # draw the boundary of the star
    radius = 7.125 # NB: Hard-coding! depends on compactness etc
    #max(np.max(column1),-1*np.min(column1),np.max(column2),-1*np.min(column2))
    circle1=plt.Circle((0,0),radius,color='black',fill=False)
    ax2.add_artist(circle1)
    # gcf() means Get Current Figure; gca() means Get Current Axis

    # ax2.set_title(r'$\theta_0 = 90$; $\alpha = 30$; $q = 2.0$')

    #formatting
    ax2.set(aspect=1)
    ax2.xaxis.set_tick_params(size=0)
    ax2.yaxis.set_tick_params(size=0)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax1.set_xlim([0, 1-(1/Nt)])
    ax1.set_xlabel('Phase'); ax1.set_ylabel('Total Flux (photons/cm^2/s)')


def fullVideo():
    num_cycles = 1
    ani = animation.FuncAnimation(
        figgy, animate, frames=num_cycles*Nt, interval=150, blit=False, repeat=False)
    writer = FFMpegWriter(fps=15) #, bitrate=1800)
    movfile = sourcefile + ".mp4"
    ani.save("FINALRUNS/" + movfile, writer=writer)
    print('\n video complete.')
    # plt.show()


# ###############################################################################

# fullVideo()
plotPulseProfile()
# plotHotSpot()
# temp_contours()
