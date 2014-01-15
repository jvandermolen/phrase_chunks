import pylab as pl
import scikits.audiolab as al
import numpy as np
import subprocess
import os
import errno

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def str2seconds(time):
    hours = int(time[:2])
    minutes = int(time[3:5])
    seconds = int(time[6:8])
    return seconds + 60*minutes + 3600*hours

def getLength(filename):
    result = subprocess.Popen(['ffmpeg', '-i', filename], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    result = [x for x in result.stdout.readlines() if "Duration" in x][0]
    two_points = result.index(':')
    return str2seconds(result[two_points + 2 : two_points + 13])

def split(wavFile, start=-1, end=-1):

    signal, fs, enc = al.wavread(wavFile)
    
    if start != -1 and end != -1:
        if start < end and start >= 0 and end <= len(signal)/fs:
            start *= fs
            end *= fs
            signal = signal[int(start):int(end)]
        else:
            raise NameError('error on start or end inputs')

    Pxx, freqs, bins, im = pl.specgram(signal, Fs=fs)

    mu = np.mean(Pxx, dtype=np.float64) #energy average
    means = np.mean(Pxx, 0, np.float64) #means through time
    oz = [np.int(x > mu) for x in means] #one if mean greater than mu, zero otherwise
    i = [(ind, oz.index(1,ind)-ind if oz[ind:].count(1)>0 else 0) for ind in range(len(oz)) if oz[ind]==0 and ((ind>0 and oz[ind-1]==1) or ind == 0)] #list of (index, number) pairs, where index is the position of a run of zeros and number is the length of the run

    minSilence = 0.5 #seconds
    cuts = [bins[ind+int(number/2)] for ind, number in i if ind > 0 and number > minSilence/bins[0]] #list of cuts put in the middle of a zero run
    cuts.append(getLength(wavFile))

    return cuts

def adjustCuts(cuts, minPhrase, maxPhrase, wavFile):
    durations = [cuts[i] - (cuts[i-1] if i > 0 else 0) for i in range(len(cuts))]
    minDuration = min(durations)
    maxDuration = max(durations)
    if minDuration >= minPhrase and maxDuration <= maxPhrase:
        return cuts
    else:
        if minDuration < minPhrase:
            iMinDuration = durations.index(minDuration)
            if (durations[(iMinDuration + 1) % len(durations)] < durations[iMinDuration - 1] and iMinDuration < len(durations)-1) or iMinDuration == 0:
                cuts.pop(iMinDuration)
            else:
                cuts.pop(iMinDuration-1)
            return adjustCuts(cuts, minPhrase, maxPhrase, wavFile)
        if maxDuration > maxPhrase:
            iMaxDuration = durations.index(maxDuration)
            start = cuts[iMaxDuration-1] if iMaxDuration > 0 else 0
            end = cuts[iMaxDuration]
            newCuts = split(wavFile, start, end)
            cuts = cuts[:iMaxDuration].extend(list(np.array(newCuts)+start)).extend(cuts[iMaxDuration+1:])
            return adjustCuts(cuts, minPhrase, maxPhrase, wavFile)

if __name__ == '__main__':
    import sys, getopt
    #name of directories with a '/' at the end
    args, video_src = getopt.getopt(sys.argv[1:], 'i:o:w:', ['frpersec='])
    args = dict(args)
    inDir = args.get('-i', './16khz/')
    outDir = args.get('-o', './chunks/')
    wavFile = args.get('-w', 'test.wav')
    frSec = int(args.get('--frpersec', '100')) #frames per second

    mkdir(outDir)
"""
    dirs = os.listdir(inDir)
    
    for wavFile in dirs:

        if wavFile[-4:] != '.wav':
            continue
"""
    dirName = wavFile.replace('.wav','')

    if os.path.isdir(outDir + dirName):
        continue
    mkdir(outDir + dirName)

    cuts = split(inDir + wavFile)
    minPhrase = 5 #seconds
    maxPhrase = 30 #seconds
    cuts = adjustCuts(cuts, minPhrase, maxPhrase, wavFile) #adjust the cuts to have chunks with a minimum length of minPhrase seconds and maximum length of maxPrhase seconds

    ctlFile = outDir + dirName + '/' + dirName + '.ctl'
    transFile = outDir + dirName + '/' + dirName + '.trans'   

    fctl = open(ctlFile, 'w')
    ftrans = open(transFile, 'w')

    for i in range(len(cuts)):
        seekTime = cuts[i-1] if i > 0 else 0
        duration = cuts[i] - seekTime

        #add line to control file
        idctlline = dirName + '_' + str(int(np.round(seekTime*frSec))) + '_' + str(int(np.round(cuts[i]*frSec)))
        ctlline = idctlline.replace('_', ' ') + ' ' + idctlline + '\n'
        fctl.write(ctlline)

        #add line to transcription file
        transline = '<s>  </s> (' + idctlline + ')' + '\n'
        ftrans.write(transline) 

        outFile = outDir + dirName + '/' + idctlline + '.wav'
        ffmpegCommand = ['ffmpeg', '-i', inDir + wavFile, '-ss', str(seekTime), '-t', str(duration), '-acodec', 'copy', outFile]
        subprocess.call(ffmpegCommand)

    fctl.close()
    ftrans.close()

# Pxx is the segments x freqs array of instantaneous power, freqs is
# the frequency vector, bins are the centers of the time bins in which
# the power is computed, and im is the matplotlib.image.AxesImage
# instance

"""
ax1 = subplot(211)
plot(t, x)
subplot(212, sharex=ax1)           
Pxx, freqs, bins, im = specgram(x, NFFT=NFFT, Fs=Fs, noverlap=900, 
                                cmap=cm.gist_heat)
show()
"""
