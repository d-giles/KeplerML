import numpy as np
np.set_printoptions(threshold='nan')
from scipy import stats
from multiprocessing import Pool,cpu_count

import sys

if sys.argv[1]:
    fl = sys.argv[1]
else:
    fl = raw_input("Input path: ")
print("Reading %s..."%fl)
lcs = np.load(fl)

def calc_outliers_pts(t, nf):
    # Is t really a necessary input? The answer is no, but eh
    posthreshold = np.mean(nf)+4*np.std(nf)
    negthreshold = np.mean(nf)-4*np.std(nf)
    
    numposoutliers,numnegoutliers,numout1s=0,0,0
    for j in range(len(nf)):
        # First checks if nf[j] is outside of 1 sigma
        if abs(np.mean(nf)-nf[j])>np.std(nf):
            numout1s += 1
            if nf[j]>posthreshold:
                numposoutliers += 1
            elif nf[j]<negthreshold:
                numnegoutliers += 1
    numoutliers=numposoutliers+numnegoutliers
    
    return numoutliers, numposoutliers, numnegoutliers, numout1s

def calc_slopes(t, nf, corrnf):

    # delta nf/delta t
    slopes=[(nf[j+1]-nf[j])/(t[j+1]-t[j]) for j in range(len(nf)-1)]
    
    #corrslopes removes the longterm linear trend (if any) and then looks at the slope
    corrslopes=[(corrnf[j+1]-corrnf[j])/(t[j+1]-t[j]) for j in range (len(corrnf)-1)]
    meanslope = np.mean(slopes)
    
    # by looking at where the 99th percentile is instead of just the largest number,
    # I think it avoids the extremes which might not be relevant (might be unreliable data)
    # Is the miniumum slope the most negative one, or the flattest one? Answer: Most negative
    maxslope=np.percentile(slopes,99)
    minslope=np.percentile(slopes,1)
    
    # Separating positive slopes and negative slopes
    # Should both include the 0 slope? I'd guess there wouldn't be a ton, but still...
    pslope=[slope for slope in slopes if slope>=0]
    nslope=[slope for slope in slopes if slope<=0]
    # Looking at the average (mean) positive and negative slopes
    if len(pslope)==0:
        meanpslope=0
    else:
        meanpslope=np.mean(pslope)
        
    if len(nslope)==0:
        meannslope=0
    else:
        meannslope=np.mean(nslope)
    
    # Quantifying the difference in shape.
    if meannslope==0:
        g_asymm = 10
    else:
        g_asymm=meanpslope / meannslope
        
    # Won't this be skewed by the fact that both pslope and nslope have all the 0's? Eh
    if len(nslope)==0:
        rough_g_asymm=10
    else:
        rough_g_asymm=len(pslope) / len(nslope)
    
    # meannslope is inherently negative, so this is the difference btw the 2
    diff_asymm=meanpslope + meannslope
    skewslope = stats.skew(slopes)
    absslopes=[abs(slope) for slope in slopes]
    meanabsslope=np.mean(absslopes)
    varabsslope=np.var(absslopes)
    varslope=np.var(slopes)
    
    #secder = Second Derivative
    # Reminder for self: the slope is "located" halfway between the flux and time points, 
    # so the delta t in the denominator is accounting for that.
    # secder = delta slopes/delta t, delta t = ((t_j-t_(j-1))+(t_(j+1)-t_j))/2
    #secder=[(slopes[j]-slopes[j-1])/((t[j+1]-t[j])/2+(t[j]-t[j-1])/2) for j in range(1, len(nf)-1)]
    #after algebraic simplification:
    secder=[2*(slopes[j]-slopes[j-1])/(t[j+1]-t[j-1]) for j in range(1, len(slopes)-1)]
    meansecder=np.mean(secder)
    
    #abssecder=[abs((slopes[j]-slopes[j-1])/((t[j+1]-t[j])/2+(t[j]-t[j-1])/2)) for j in range (1, len(slopes)-1)]
    # simplification:
    
    abssecder=np.abs(np.array(secder))
    absmeansecder=np.mean(abssecder)
    if len(pslope)==0:
        pslopestds=0
    else:
        pslopestds=np.std(pslope)
    if len(nslope)==0:
        nslopesstds=0
        stdratio=10
    else:
        nslopestds=np.std(nslope)
        stdratio=pslopestds/nslopestds
    
    sdstds=np.std(secder)
    meanstds=np.mean(secder)
    

    num_pspikes,num_nspikes,num_psdspikes,num_nsdspikes=0,0,0,0
    
    for slope in slopes:
        if slope>=meanpslope+3*pslopestds:
            num_pspikes+=1
        elif slope<=meanslope-3*nslopestds:
            num_nspikes+=1
    
    for sder in secder:
        if sder>=4*sdstds:
            num_psdspikes+=1
        elif sder<=-4*sdstds:
            num_nsdspikes+=1
    if nslopestds==0:
        stdratio=10
    else:
        stdratio = pslopestds / nslopestds
    
    # The ratio of postive slopes with a following postive slope to the total number of points.
    pstrendcount = 0
    for j,slope in enumerate(slopes[:-1]):
        if slope > 0 and slopes[j+1]>0:
            pstrendcount += 1
        
    pstrend=pstrendcount/len(slopes)

    slope_array = [maxslope, minslope, meanpslope, \
                   meannslope, g_asymm, rough_g_asymm, \
                   diff_asymm, skewslope, varabsslope, \
                   varslope, meanabsslope, absmeansecder, \
                   num_pspikes, num_nspikes, num_psdspikes, \
                   num_nsdspikes, stdratio, pstrend]

    return slopes, corrslopes, secder, slope_array

def calc_maxmin_periodics(t, nf, err):
    
    # This looks up the local maximums. Adds a peak if it's the largest within 10 points on either side.
    # Q: Is there a way to do this and take into account drastically different periodicity scales?
    
    naivemax,nmax_times = [],[]
    naivemins = []
    for j in range(len(nf)):
        if nf[j] == max(nf[max(j-10,0):min(j+10,len(nf)-1)]):
            naivemax.append(nf[j])
            nmax_times.append(t[j])
        elif nf[j] == min(nf[max(j-10,0):min(j+10,len(nf)-1)]):
            naivemins.append(nf[j])
    len_nmax=len(naivemax) #F33
    len_nmin=len(naivemins) #F34    
    if len(naivemax)>2:
        mautocorrcoef = np.corrcoef(naivemax[:-1], naivemax[1:])[0][1] #F35
    else:
        mautocorrcoef = 0
    """peak to peak slopes"""
    ppslopes = [abs((naivemax[j+1]-naivemax[j])/(nmax_times[j+1]-nmax_times[j])) \
                for j in range(len(naivemax)-1)]
    if len(ppslopes)==0:
        ptpslopes = 0
    else:
        ptpslopes=np.mean(ppslopes) #F36

    maxdiff=[nmax_times[j+1]-nmax_times[j] for j in range(len(naivemax)-1)]

    if len(maxdiff)==0:
        periodicity=0
        periodicityr=0
        naiveperiod=0
    else:
        periodicity=np.std(maxdiff)/np.mean(maxdiff) #F37
        periodicityr=np.sum(abs(maxdiff-np.mean(maxdiff)))/np.mean(maxdiff) #F38
        naiveperiod=np.mean(maxdiff) #F39
    if len(naivemax)==0:
        maxvars=0
        maxvarsr=0
    else:
        maxvars = np.std(naivemax)/np.mean(naivemax) #F40
        maxvarsr = np.sum(abs(naivemax-np.mean(naivemax)))/np.mean(naivemax) #F41

    emin = naivemins[::2] # even indice minimums
    omin = naivemins[1::2] # odd indice minimums
    meanemin = np.mean(emin)
    if len(omin)==0:
        meanomin=0
    else:
        meanomin = np.mean(omin)
    oeratio = meanomin/meanemin #F42

    peaktopeak_array = [len_nmax, len_nmin, mautocorrcoef,\
                        ptpslopes, periodicity, periodicityr, \
                        naiveperiod, maxvars, maxvarsr, oeratio]

    return peaktopeak_array, naivemax, naivemins

def featureCalculation(lc):
    t = lc[0]
    nf = lc[1]
    err = lc[2]
    # t = time
    # err = error
    # nf = normalized flux.
    
    longtermtrend = np.polyfit(t, nf, 1)[0] # Feature 1 (Abbr. F1) overall slope
    yoff = np.polyfit(t, nf, 1)[1] # Not a feature? y-intercept of linear fit
    meanmedrat = np.mean(nf) / np.median(nf) # F2
    skews = stats.skew(nf) # F3
    varss = np.var(nf) # F4
    coeffvar = np.std(nf)/np.mean(nf) #F5
    stds = np.std(nf) #F6

    corrnf = nf - longtermtrend*t - yoff #this removes any linear slope to lc so you can look at just troughs - is this a sign err tho?
    # D: I don't think there's a sign error

    # Features 7 to 10
    numoutliers, numposoutliers, numnegoutliers, numout1s = calc_outliers_pts(t, nf)
    kurt = stats.kurtosis(nf)

    mad = np.median([abs(nf[j]-np.median(nf)) for j in range(len(nf))])

    # slopes array contains features 13-30
    slopes, corrslopes, secder, slopes_array = calc_slopes(t, nf, corrnf) 
    maxslope = slopes_array[0]
    minslope = slopes_array[1]
    meanpslope  = slopes_array[2]
    meannslope  = slopes_array[3]
    g_asymm = slopes_array[4]
    rough_g_asymm  = slopes_array[5]
    diff_asymm  = slopes_array[6]
    skewslope  = slopes_array[7]
    varabsslope  = slopes_array[8]
    varslope  = slopes_array[9]
    meanabsslope  = slopes_array[10]
    absmeansecder = slopes_array[11]
    num_pspikes = slopes_array[12]
    num_nspikes  = slopes_array[13]
    num_psdspikes = slopes_array[14]
    num_nsdspikes = slopes_array[15]
    stdratio = slopes_array[16]
    pstrend = slopes_array[17]

    # Checks if the flux crosses the zero line.
    zcrossind= [j for j in range(len(nf)-1) if corrnf[j]*corrnf[j+1]<0]
    num_zcross = len(zcrossind) #F31

    plusminus=[j for j in range(1,len(slopes)) if (slopes[j]<0)&(slopes[j-1]>0)]
    num_pm = len(plusminus)

    # peak to peak array contains features 33 - 42
    peaktopeak_array, naivemax, naivemins = calc_maxmin_periodics(t, nf, err)
    len_nmax=peaktopeak_array[0]
    len_nmin=peaktopeak_array[1]
    mautocorrcoef=peaktopeak_array[2]
    ptpslopes=peaktopeak_array[3]
    periodicity=peaktopeak_array[4]
    periodicityr=peaktopeak_array[5]
    naiveperiod=peaktopeak_array[6]
    maxvars=peaktopeak_array[7]
    maxvarsr=peaktopeak_array[8]
    oeratio=peaktopeak_array[9]

    # amp here is actually amp_2 in revantese
    # 2x the amplitude (peak-to-peak really), the 1st percentile will be negative, so it's really adding magnitudes
    amp = np.percentile(nf,99)-np.percentile(nf,1) #F43
    normamp = amp / np.mean(nf) #this should prob go, since flux is norm'd #F44

    # ratio of points within 10% of middle to total number of points 
    mbp = len([nf[j] for j in range(len(nf))\
               if (nf[j] < (np.median(nf) + 0.1*amp)) \
               & (nf[j] > (np.median(nf)-0.1*amp))]) / len(nf) #F45
    
    f595 = np.percentile(nf,95)-np.percentile(nf,5)
    f1090 =np.percentile(nf,90)-np.percentile(nf,10)
    f1782 =np.percentile(nf, 82)-np.percentile(nf, 17)
    f2575 =np.percentile(nf, 75)-np.percentile(nf, 25)
    f3267 =np.percentile(nf, 67)-np.percentile(nf, 32)
    f4060 =np.percentile(nf, 60)-np.percentile(nf, 40)
    mid20 =f4060/f595 #F46
    mid35 =f3267/f595 #F47
    mid50 =f2575/f595 #F48
    mid65 =f1782/f595 #F49
    mid80 =f1090/f595 #F50 

    percentamp = max([(nf[j]-np.median(nf)) / np.median(nf) for j in range(len(nf))]) #F51
    magratio = (max(nf)-np.median(nf)) / amp #F52

    autocorrcoef = np.corrcoef(nf[:-1], nf[1:])[0][1] #F54

    sautocorrcoef = np.corrcoef(slopes[:-1], slopes[1:])[0][1] #F55
    
    #measures the slope before and after the maximums
    flatness = [np.mean(slopes[max(0,j-6):min(max(0,j-1), len(slopes)-1):1])\
                - np.mean(slopes[max(0,j):min(j+4, len(slopes)-1):1])\
                for j in range(6,len(slopes)-6) \
                if nf[j] in naivemax]
    
    if len(flatness)==0: flatmean=0
    else: flatmean = np.nansum(flatness)/len(flatness) #F55
    
    #measures the slope before and after the minimums
    # trying flatness w slopes and nf rather than "corr" vals, despite orig def in RN's program
    tflatness = [-np.mean(slopes[max(0,j-6):min(max(j-1,0),len(slopes)-1):1])\
                 + np.mean(slopes[j:min(j+4,len(slopes)-1):1])\
                 for j in range(6,len(slopes)-6)\
                 if nf[j] in naivemins] 
    
    # tflatness for mins, flatness for maxes
    if len(tflatness)==0: 
        tflatmean=0
    else: 
        tflatmean = np.nansum(tflatness) / len(tflatness) #F56

    roundness = [np.mean(secder[max(0,j-6):j:1])\
                  +np.mean(secder[j:min(j+6,len(secder)-1):1])\
                  for j in range(6,len(secder)-6)\
                  if nf[j] in naivemax]
    
    if len(roundness)==0: 
        roundmean=0
    else: 
        roundmean = np.nansum(roundness) / len(roundness) #F57
        
    troundness = [np.mean(secder[max(0,j-6):j])\
                  +np.mean(secder[j:min(j+6,len(secder)-1)])\
                  for j in range(6,len(secder)-6)\
                  if nf[j] in naivemins]
    
    if len(troundness)==0:
        troundmean=0
    else:
        troundmean = np.nansum(troundness)/len(troundness) #F58
        
    if troundmean==0 and roundmean==0: 
        roundrat=1
    elif troundmean==0: 
        roundrat=10
    else: 
        roundrat = roundmean / troundmean #F59
        
    if flatmean==0 and tflatmean==0: 
        flatrat=1
    elif tflatmean==0: 
        flatrat=10
    else: 
        flatrat = flatmean / tflatmean #F60
        
    ndata = [longtermtrend, meanmedrat, skews, varss, coeffvar, stds, numoutliers, numnegoutliers, numposoutliers, numout1s, kurt, mad, maxslope, minslope, meanpslope, meannslope, g_asymm, rough_g_asymm, diff_asymm, skewslope, varabsslope, varslope, meanabsslope, absmeansecder, num_pspikes, num_nspikes, num_psdspikes, num_nsdspikes,stdratio, pstrend, num_zcross, num_pm, len_nmax, len_nmin, mautocorrcoef, ptpslopes, periodicity, periodicityr, naiveperiod, maxvars, maxvarsr, oeratio, amp, normamp,mbp, mid20, mid35, mid50, mid65, mid80, percentamp, magratio, sautocorrcoef, autocorrcoef, flatmean, tflatmean, roundmean, troundmean, roundrat, flatrat]
    return ndata

fdata = []

numCpus = cpu_count()
print("Using %s cpus to calculate features..."%numCpus)
p = Pool(numCpus)
fdata = p.map(featureCalculation,lcs)
p.close()
p.join()
fdata = np.array(fdata)
if sys.argv[2]:
    of = sys.argv[2]
else:
    of = raw_input('Output path: ')
print("Saving output to %s"%of)
np.save(of,fdata)
print("Done.")