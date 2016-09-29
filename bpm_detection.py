import wave, array, math, time, argparse, sys
import numpy, pywt
from scipy import signal
import pdb
from pydub import AudioSegment
import matplotlib.pyplot as plt
from pydub import AudioSegment
from mutagen.id3 import ID3, TPE1,TBPM,ID3NoHeaderError

def read_wav(filename):

    #open file, get metadata for audio
    try:
        wsound = AudioSegment.from_mp3(filename)
        wsound.export(filename, format="wav")

        wave_file = wave.open(filename,'rb')
	
    except IOError:
        print ("File Not Found")
        return


    nsamps = wave_file.getnframes();
    assert(nsamps > 0);

    frame_per_second = wave_file.getframerate()
    assert(frame_per_second > 0)

    # read entire file and make into an array
    samps = list(array.array('i',wave_file.readframes(nsamps)))
    #print 'Read', nsamps,'samples from', filename
    try:
        assert(nsamps == len(samps))
    except AssertionError :
        print  (nsamps, "not equal to", len(samps))
    
    return samps, frame_per_second
    

def no_audio_data():
    print ("No audio data for sample, skipping...")
    return None, None


    
# simple peak detection
def peak_detect(data):
    max_val = numpy.amax(abs(data)) 
    peak_ndx = numpy.where(data==max_val)
    if len(peak_ndx[0]) == 0: #if nothing found then the max must be negative
        peak_ndx = numpy.where(data==-max_val)
    return peak_ndx
    
def bpm_detector(data,frame_per_second):
    cA = [] 
    cD = []
    correl = []
    cD_sum = []
    levels = 4
    max_decimation = 2**(levels-1);
    min_ndx = 60./ 220 * (frame_per_second/max_decimation)
    max_ndx = 60./ 40 * (frame_per_second/max_decimation)
    
    for loop in range(0,int(levels)):
        cD = []
        # 1) DWT
        if loop == 0:
            [cA,cD] = pywt.dwt(data,'db4');
            cD_minlen = len(cD)/max_decimation+1;
            cD_sum = numpy.zeros(cD_minlen);
        else:
            [cA,cD] = pywt.dwt(cA,'db4');
        # 2) Filter
        cD = signal.lfilter([0.01],[1 -0.99],cD);

        # 4) Subtractargs.filename out the mean.

        # 5) Decimate for reconstruction later.
        cD = abs(cD[::(2**(levels-loop-1))]);
        cD = cD - numpy.mean(cD);
        # 6) Recombine the signal before ACF
        #    essentially, each level I concatenate 
        #    the detail coefs (i.e. the HPF values)
        #    to the beginning of the array
        cD_sum = cD[0:cD_minlen] + cD_sum;

    if [b for b in cA if b != 0.0] == []:
        return no_audio_data()
    # adding in the approximate data as well...    
    cA = signal.lfilter([0.01],[1 -0.99],cA);
    cA = abs(cA);
    cA = cA - numpy.mean(cA);
    cD_sum = cA[0:cD_minlen] + cD_sum;
    
    # ACF
    correl = numpy.correlate(cD_sum,cD_sum,'full') 
    
    midpoint = len(correl) / 2
    correl_midpoint_tmp = correl[midpoint:]
    peak_ndx = peak_detect(correl_midpoint_tmp[min_ndx:max_ndx]);
    if len(peak_ndx) > 1:
        return no_audio_data()
        
    peak_ndx_adjusted = peak_ndx[0]+min_ndx;
    bpm = 60./ peak_ndx_adjusted * (frame_per_second/max_decimation)
    print (bpm)
    return bpm,correl
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .wav file to determine the Beats Per Minute.')
    parser.add_argument('--filename', required=True,
                   help='.wav file for processing')
    parser.add_argument('--window', type=float, default=3,
                   help='size of the the window (seconds) that will be scanned to determine the bpm.  Typically less than 10 seconds. [3]')

    args = parser.parse_args()
    samps,frame_per_second = read_wav(args.filename)
    
    data = []
    correl=[]
    bpm = 0
    n=0;
    nsamps = len(samps)
    window_samps = int(args.window*frame_per_second)         
    samps_ndx = 0;  #first sample in window_ndx 
    max_window_ndx = nsamps / window_samps;
    bpms = numpy.zeros(max_window_ndx)

    #iterate through all windows
    for window_ndx in range(0,int(max_window_ndx)):

        #get a new set of samples
        data = samps[samps_ndx:samps_ndx+window_samps]
        if not ((len(data) % window_samps) == 0):
            raise AssertionError( str(len(data) ) ) 
        
        bpm, correl_temp = bpm_detector(data,frame_per_second)
        if bpm == None:
            continue
        bpms[window_ndx] = bpm
        correl = correl_temp
        
        #iterate at the end of the loop
        samps_ndx = samps_ndx+window_samps;
        n=n+1; #counter for debug...

    bpm = numpy.median(bpms)
    print ('Completed.  Estimated Beats Per Minute:', bpm)
    
    #Adding bpm to the mp3 file back
    try:
       audio = ID3(args.filename)
    except ID3NoHeaderError:
       audio = ID3()
       audio.add(TBPM(encoding=3,text=str(bpm)))
       audio.save(args.filename)
    
    n = range(0,int(len(correl)))
    plt.plot(n,abs(correl)); 
    plt.show(False); #plot non-blocking
    time.sleep(10);
    plt.close()
