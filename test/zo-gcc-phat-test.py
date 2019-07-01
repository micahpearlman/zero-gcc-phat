'''
Implementation of GCC PHAT in Python
See: https://github.com/respeaker/mic_array/blob/master/gcc_phat.py
'''

#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'test'))
    print(os.getcwd())
except:
    pass

#%% Define Python GCC Phat
import numpy as np


def gcc_phat(sig, refsig, fs=1, max_tau=None, interp=16):
    '''
    This function computes the offset between the signal sig and the reference signal refsig
    using the Generalized Cross Correlation - Phase Transform (GCC-PHAT)method.
    '''
    
    # make sure the length for the FFT is larger or equal than len(sig) + len(refsig)
    n = sig.shape[0] + refsig.shape[0]

    # Generalized Cross Correlation Phase Transform
    SIG = np.fft.rfft(sig, n=n)    
    REFSIG = np.fft.rfft(refsig, n=n)    
    REFSIG_CONJ = np.conj(REFSIG)
    
    R = SIG * REFSIG_CONJ
    IR = R / np.abs(R)
    cc = np.fft.irfft(IR, n=(interp * n))

    max_shift = int(interp * n / 2)
    if max_tau:
        max_shift = np.minimum(int(interp * fs * max_tau), max_shift)

    cc_mn = cc[-max_shift:]
    cc_mx = cc[:max_shift+1]
    cc = np.concatenate((cc_mn, cc_mx))

    # find max cross correlation index
    shift = np.argmax(np.abs(cc)) - max_shift

    tau = shift / float(interp * fs)
    
    return tau, cc


#%% Test GCC PHAT

refsig = np.linspace(1, 10, 10)

for i in range(0, 10):
    sig = np.concatenate((np.linspace(0, 0, i), refsig, np.linspace(0, 0, 10 - i)))
    offset, _ = gcc_phat(sig, refsig)
    print(offset)





#%%
