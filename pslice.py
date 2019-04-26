import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm  #color maps
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import sys #this is for accessing command line input
import os; import glob #navigating directories
# global Nt #number of time steps one rotation of the star is divided into
# global Nen #number of energy bands the spectrum is divided into
global column1; global column2; global column3;
global sourcefile

Nt  = 32   #NB: One less than the input in definitions (default 31)
Nen = 64  #NB: Two less than the input in definitions (default 126)

'''
This code analyzes slices through the parameter space (pslice).
After a run, which produces a 1-parameter family of light-curves,
We can plot those curves together or make a movie cycling through them
so see how the profile changes as the single parameter is varied.
'''

'''
Nice Colors: http://ksrowell.com/blog-visualizing-data/2012/02/02/optimal-colors-for-graphs/
'''

###############################################################################

# # NB: Uncomment for spot type comparison subplots
# import matplotlib.gridspec as gridspec
# figgy = plt.figure()
# gs = gridspec.GridSpec(2,1,height_ratios=[2.2,1])
# ax1 = plt.subplot(gs[0])
# ax2 = plt.subplot(gs[1])
# filenames=[]


def plotcurves():
    path = 'BestFit/'
    for filename in glob.glob(os.path.join(path, '*.dat')):
        # print(filename)
        c1,c2,c3 = np.loadtxt(filename,unpack=True)
        x,y = c1,c3
        average_flux = np.sum(y)/len(y)
        x = np.append(x,c1+1)
        y = np.append(y,c3)
        # y *= 1/(average_flux) #normalize
        amp = round(np.max(y)/np.min(y),2)
        plt.plot(x,y,linewidth=1,label='Amplitude = '+str(amp))#,linestyle='dashed',color=colors[i])

        plt.ylabel('Flux')
        plt.tick_params(direction='in');
        plt.xlim(0,2); plt.ylim(0)
        plt.tick_params(labelbottom=False)
        # plt.yticks([0.75,1,1.25])
        plt.legend(loc=1)

        plt.savefig(str(filename)+'.png')
        plt.close()


def capCompare():  #FIGURE 4
    '''includes an option to plot the sum of the curves (Ytotal).
    This is useful for seeing how the contributions from the two hot spots
    sum to the total profile.'''

    path = 'FINAL_PLOTS/FIG4_SpotTypes_LC/300Hz'
    i=0
    deltaAB=np.zeros(Nt+1);
    deltaAC=np.zeros(Nt+1);
    # colors=['black', 'mediumblue', 'crimson']
    # colors=['#535154','#396AB1','#3E9651']
    # colors=['#396AB1','#3E9651','#DA7C30']
    colors=['#535154','#396AB1','#DA7C30']
    labels=['(A) uniform polar cap', '(B) uniform hotspot', '(C) realistic hotspot']
    plt.rc('font', size=12)
    for filename in glob.glob(os.path.join(path, '*.dat')):
        # print(filename)
        c1,c2,c3 = np.loadtxt(filename,unpack=True)
        x,y = c1,c3
        average_flux = np.sum(y)/len(y)
        x = np.append(x,1)
        y = np.append(y,c3[0])
        y *= 1/(average_flux) #normalize
        ymax = np.max(y)
        ymin = np.min(y)
        # plot the deltas (B-A) and (C-A)
        if i==0:
            deltaAB = np.add(deltaAB,-y)
            deltaAC = np.add(deltaAC,-y)
        if i==1: deltaAB = np.add(deltaAB,y)
        if i==2: deltaAC = np.add(deltaAC,y)
        ax1.plot(x,y,label=labels[i],linewidth=1.5,color=colors[i])#,linestyle='dashed',color=colors[i])
        i+=1

    # ax1.set_title('Evolution of the curve')
    ax2.set_xlabel('Phase');
    ax1.set_ylabel('Normalized Flux')
    ax1.tick_params(direction='in');
    ax1.set_xlim(0,1); ax1.set_ylim(0.9*ymin,ymax*1.35)
    ax1.tick_params(labelbottom=False)
    ax1.set_yticks([0.75,1,1.25])
    ax1.legend(loc=1)
    ax1.annotate(r'$f_s = 300$ Hz',(0.05,0.85),xycoords='axes fraction',
        bbox=dict(boxstyle="square, pad=0.5", fc="w",color='lightgray'))

    # second subplot shows the differences (B-A) and (C-A)
    ax2.plot(x,deltaAB,label='B-A',color=colors[1])
    ax2.plot(x,deltaAC,label='C-A',color=colors[2])

    ax2.set_xlim(0,1);
    ax2.set_ylim(-0.1,0.05);
    # ax2.set_ylim(-0.16,0.10);
    ax2.set_yticks([-0.05,0,0.05]);
    # ax2.set_yticks([-0.1,0,0.1]);
    ax2.tick_params(direction='in')
    ax2.set_ylabel(r'$\Delta F \, / \, \langle F \rangle$')
    # ax2.set_xticks([0,0.5,1])
    ax2.axhline(y=0,color='gray',linestyle='dashed',linewidth=0.75)
    ax2.legend(loc=4)

    plt.show()



def ampVSenergy():  #FIGURE 5
    '''plot the total amplitude of the light-curve (max-min)/average
       as a function of photon energy'''

    #NB: Nt=32, Nen=64
    path = '../FINAL_PLOTS_data/FIG5_SpotTypes_Amp/600Hz'
    i=0
    Temp = 0.35
    plt.rc('font', size=12)
    # colors=['#535154','#396AB1','#3E9651']
    # colors=['#396AB1','#3E9651','#DA7C30']
    colors=['#535154','#396AB1','#DA7C30']
    labels=['(A) uniform polar cap', '(B) uniform hotspot', '(C) realistic hotspot']
    for filename in glob.glob(os.path.join(path, '*.dat')):
        time_data,energy_data,flux_data = np.loadtxt(filename,unpack=True)
        Y_amplitudes=np.zeros(Nen)
        X_energies=np.zeros(Nen)
        for e in range(Nen):

            #Step1: compute the average flux over 1 period
            totalflux=0
            for phase in range(Nt):
                totalflux += flux_data[Nen*phase + e]
            avflux=totalflux/Nt
            #Step2: compute the total squared fractional amplitude
            totalSqDev=0
            for phase in range(Nt):
                # "RMS fractional amplitude" is defined as the square root
                # of the mean of the squared deviation from the average
                myflux = flux_data[Nen*phase + e]
                totalSqDev+= ((myflux-avflux)/avflux)**2
            #Step3: RMS
            Y_amplitudes[e] = (totalSqDev/Nt)**(0.5) #sqrt[mean fractional amplitude]
            X_energies[e] = energy_data[e]/Temp
        plt.plot(X_energies,Y_amplitudes,label=labels[i],linewidth=2,color=colors[i])
        i+=1
    plt.xscale('log');
    plt.xlabel(r'Photon Energy $(\epsilon/kT)$'); plt.ylabel('Pulse Amplitude')
    plt.annotate(r'$f_s = 600$ Hz',(0.6,0.86),xycoords='axes fraction',
        bbox=dict(boxstyle="square, pad=0.5", fc="w",color='lightgray'))
    plt.ylim(0)
    # plt.yticks([0,0.5,1,1.5,2])
    # plt.yticks([0,0.5,1])
    plt.tick_params(direction='in');
    plt.legend(loc=2)
    plt.show()


############################################################

def animate(t):
    print('creating frame ' + str(t) + ' ... ',end='\r')
    N = len(filenames)
    ax1.clear()
    c1,c2,c3 = np.loadtxt(filenames[t],unpack=True)
    x=c1; y=c3 #(c3-min(c3))/max(c3-min(c3))
    x = np.append(x,c1+1) #plot two periods
    y = np.append(y,c3)
    amp = round(np.max(y)/np.min(y),2)
    ax1.plot(x,y,label='Amplitude = '+str(amp))
    ax1.set_xlim([0, 2*(1-(1/Nt))])
    ax1.set_ylim(0,np.max(y)*1.1)
    ax1.set_xlabel('Phase'); ax1.set_ylabel('Total Flux')
    ax1.legend(loc=4)
    plt.title(filenames[t])


def lightcurvemovie():
    '''animate the evolution of the lightcurve through parameter space'''
    path = 'FINALRUNS/victoryT'
    for filename in glob.glob(os.path.join(path, '*.dat')):
        # print(filename)
        filenames.append(str(filename))
    N = len(filenames)
    ani = animation.FuncAnimation(
        figgy, animate, frames=N, interval=150, blit=False, repeat=True)
    writer = FFMpegWriter(fps=10) #, bitrate=1800)
    movfile = "victoryT.mp4"
    ani.save('FINALRUNS/'+movfile, writer=writer)
    print('\n video complete.')



def LCbreakdown():  #FIGURES 3 & 8
    '''includes an option to plot the sum of the curves (Ytotal).
    This is useful for seeing how the contributions from the two hot spots
    sum to the total profile.'''

    path = '../FINAL_PLOTS_data/FIG8_QDLightCurve/'
    i=0
    Ytotal=np.zeros(2*(Nt+1));
    # colors=['crimson','royalblue','purple']
    colors=['#CC2529','#396AB1','#6B4C9A'] #optimized red, blue, and purple
    # labels = ['North','South','Total']
    labels = ['Spot','Ring','Total']
    # plt.rc('font', size=12)

    # dirty method to find the over-all normalization first, then loop again to plot
    avflux=0
    for filename in glob.glob(os.path.join(path, '*.dat')):
        c1,c2,c3 = np.loadtxt(filename,unpack=True)
        avflux += np.sum(c3)/len(c3)

    for filename in glob.glob(os.path.join(path, '*.dat')):
        c1,c2,c3 = np.loadtxt(filename,unpack=True)
        x1 = np.append(c1,1);
        y1 = np.append(c3,c3[0]) #repeat first value no there's no space
        x = np.append(x1,x1+1) #plot two periods
        y = np.append(y1,y1)

        y*=(1/avflux) #normalize

        Ytotal = np.add(Ytotal,y)
        plt.plot(x,y,label=labels[i],linewidth=1.9,linestyle='dashed',color=colors[i])
        i+=1
    plt.plot(x,Ytotal,label='Total',color=colors[i],linewidth=2.2)
    plt.xlim([0,2])
    plt.ylim([0,1.6])#np.max(Ytotal)*1.3])
    # plt.title('title')
    plt.xlabel('Phase'); plt.ylabel('Normalized Flux')
    # plt.xticks([0,0.5,1])
    plt.tick_params(direction='in',pad=7); #pad to avoid 0.0's overlapping
    # plt.xticks([0,0.5,1])
    # plt.xticks([0,0.5,1,1.5,2])
    plt.axvline(x=1,color='gray',linestyle='dashed',linewidth=1)
    # plt.axes().set_aspect(0.6) #for FIG3, 1-period
    plt.legend(loc=1)
    # plt.annotate(r'$f_s = 100$ Hz',(0.05,0.9),xycoords='axes fraction',
        # bbox=dict(boxstyle="square, pad=0.5", fc="w",color='gray'))
    plt.show()



def LCcompare():  #FIGURE 2
    '''includes an option to plot the sum of the curves (Ytotal).
    This is useful for seeing how the contributions from the two hot spots
    sum to the total profile.'''

    path = 'FINAL_PLOTS/FIG2-3_DipoleLightCurves/'
    i=0
    labels=[r'$f_s = 10$ Hz',r'$f_s = 300$ Hz',r'$f_s = 600$ Hz']
    colors=['#948B3D','#DA7C30','#CC2529'] #yellow-orange-red
    # colors=['#191970','#4169e1','#3E9651'] #blue-green
    colors=['#922428','#3E9651','#7b68ee']
    # colors=['#222222','#4169e1','#3E9735']
    # colors=['#191970','#6B4C9A','#3E9651']
    # styles=['solid','dashdot','dashed']
    styles=['solid','solid','solid']
    # plt.rc('font', size=12)
    ymax=1
    for filename in glob.glob(os.path.join(path, '*.dat')):
        c1,c2,c3 = np.loadtxt(filename,unpack=True)
        x,y = c1,c3
        average_flux = np.sum(y)/len(y)
        x = np.append(x,c1+1) #plot two periods
        y = np.append(y,c3)
        y *= 1/(average_flux) #normalize
        ymax = np.max(y)
        # ymin = np.min(y)
        plt.plot(x,y,label=labels[i],color=colors[i],linewidth=2,linestyle=styles[i])
        i+=1
    plt.ylim(0,ymax*1.1)
    plt.xlim(0,2)
    plt.xticks([0,0.5,1,1.5,2])
    plt.axvline(x=1,color='gray',linestyle='dashed',linewidth=1)
    plt.tick_params(direction='in',pad=7)
    # plt.title('Skewness')
    plt.xlabel('Phase'); plt.ylabel('Normalized Flux')
    plt.legend(loc=2)
    plt.show()



###############################################################################

# lightcurvemovie()
# capCompare()
# LCbreakdown()
# LCcompare()
# plotcurves()
ampVSenergy()
