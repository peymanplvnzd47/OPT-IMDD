
"""
Created on Mon Jun  5 17:21:11 2023

@author: Peyman, Ale
"""



"""BER_OPTIMDD_ER.ipynb




"""

"""# IMPORT LIBRARIES"""

import logging as logg
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.constants as const
import numpy as np
from numpy.fft import fft, fftfreq, ifft
from numpy.random import normal
from scipy.ndimage.filters import gaussian_filter
import scipy.constants as const
from scipy.special import erf, erfc
from tqdm.notebook import tqdm
import math
import scipy as sp
from PIL import Image, ImageTk
import tkinter as tk
from PIL import Image, ImageTk
from numba import njit, prange





print('HELLO MY OPT-IMDD')



"""# FUNCTIONS"""





def bitarray2dec(in_bitarray):
    """
    Converts an input NumPy array of bits (0 and 1) to a decimal integer.
    """

    number = 0

    for i in range(len(in_bitarray)):
        number = number + in_bitarray[i] * pow(2, len(in_bitarray) - 1 - i)

    return number


def GrayMapping(M, constType):
    """
    Gray Mapping for digital modulations.
  
    """
    L = int(M - 1)
    const = np.arange(-L, L + 1, 2)

    return const






def modulateGray(bits, M, constType):
    """
    Modulate bit sequences to constellation symbol sequences (w/ Gray mapping).
    """
    bitsSymb = int(np.log2(M))
    const = GrayMapping(M, constType)
    symb = bits.reshape(-1, bitsSymb).T
    symbInd = bitarray2dec(symb)
    return const[symbInd]


def pnorm(x):
    """
    Normalize the average power of each componennt of x.
    """
    return x / np.sqrt(np.mean(x * np.conj(x)).real)


def decimal2bitarray(number, bit_width):
    """
    Converts a positive integer to NumPy array of the specified size containing bits (0 and 1). This version is slightly
    quicker that dec2bitarray but only work for one integer.
  
    """
    result = np.zeros(bit_width, np.int8)
    i = 1
    pox = 0
    while i <= number:
        if i & number:
            result[bit_width - pox - 1] = 1
        i <<= 1
        pox += 1
    return result


def dec2bitarray(in_number, bit_width):
    """
    Converts a positive integer or an array-like of positive integers to NumPy array of the specified size containing
    bits (0 and 1).
    """

    if isinstance(in_number, (np.integer, int)):
        return decimal2bitarray(in_number, bit_width).copy()
    result = np.zeros(bit_width * len(in_number), np.int8)
    for pox, number in enumerate(in_number):
        result[pox * bit_width:(pox + 1) * bit_width] = decimal2bitarray(number, bit_width).copy()
    return result



def upsample(x, n):
    """
    Upsample the input array by a factor of n
    Adds n-1 zeros between consecutive samples of x
   
    """
    y = np.empty(len(x) * n, dtype=float)
    y[0::n] = x
    zero_array = np.zeros(len(x), dtype=float)
    for i in range(1, n):
        y[i::n] = zero_array

    return y


def rcosfilter(N, alpha, Ts, Fs):
    """
    Generates a raised cosine (RC) filter (FIR) impulse response.
    Parameters
    ----------
    N : int
        Length of the filter in samples.
    alpha : float
        Roll off factor (Valid values are [0, 1]).
    Ts : float
        Symbol period in seconds.
    Fs : float
        Sampling Rate in Hz.
    Returns
    -------
    time_idx : 1-D ndarray (float)
        Array containing the time indices, in seconds, for the impulse response.
    h_rc : 1-D ndarray (float)
        Impulse response of the raised cosine filter.
    """

    T_delta = 1 / float(Fs)
    time_idx = ((np.arange(N) - N / 2)) * T_delta
    sample_num = np.arange(N)
    h_rc = np.zeros(N, dtype=float)

    for x in sample_num:
        t = (x - N / 2) * T_delta
        if t == 0.0:
            h_rc[x] = 1.0
        elif alpha != 0 and t == Ts / (2 * alpha):
            h_rc[x] = (np.pi / 4) * (np.sin(np.pi * t / Ts) / (np.pi * t / Ts))
        elif alpha != 0 and t == -Ts / (2 * alpha):
            h_rc[x] = (np.pi / 4) * (np.sin(np.pi * t / Ts) / (np.pi * t / Ts))
        else:
            h_rc[x] = (np.sin(np.pi * t / Ts) / (np.pi * t / Ts)) * \
                      (np.cos(np.pi * alpha * t / Ts) / (1 - (((2 * alpha * t) / Ts) * ((2 * alpha * t) / Ts))))

    return time_idx, h_rc


def rrcosfilter(N, alpha, Ts, Fs):
    """
    Generates a root raised cosine (RRC) filter (FIR) impulse response.
    Parameters
    ----------
    N : int
        Length of the filter in samples.
    alpha : float
        Roll off factor (Valid values are [0, 1]).
    Ts : float
        Symbol period in seconds.
    Fs : float
        Sampling Rate in Hz.
    Returns
    ---------
    time_idx : 1-D ndarray of floats
        Array containing the time indices, in seconds, for
        the impulse response.
    h_rrc : 1-D ndarray of floats
        Impulse response of the root raised cosine filter.
    """

    T_delta = 1 / float(Fs)
    time_idx = ((np.arange(N) - N / 2)) * T_delta
    sample_num = np.arange(N)
    h_rrc = np.zeros(N, dtype=float)

    for x in sample_num:
        t = (x - N / 2) * T_delta
        if t == 0.0:
            h_rrc[x] = 1.0 - alpha + (4 * alpha / np.pi)
        elif alpha != 0 and t == Ts / (4 * alpha):
            h_rrc[x] = (alpha / np.sqrt(2)) * (((1 + 2 / np.pi) * \
                                                (np.sin(np.pi / (4 * alpha)))) + (
                                                           (1 - 2 / np.pi) * (np.cos(np.pi / (4 * alpha)))))
        elif alpha != 0 and t == -Ts / (4 * alpha):
            h_rrc[x] = (alpha / np.sqrt(2)) * (((1 + 2 / np.pi) * \
                                                (np.sin(np.pi / (4 * alpha)))) + (
                                                           (1 - 2 / np.pi) * (np.cos(np.pi / (4 * alpha)))))
        else:
            h_rrc[x] = (np.sin(np.pi * t * (1 - alpha) / Ts) + \
                        4 * alpha * (t / Ts) * np.cos(np.pi * t * (1 + alpha) / Ts)) / \
                       (np.pi * t * (1 - (4 * alpha * t / Ts) * (4 * alpha * t / Ts)) / Ts)

    return time_idx, h_rrc



def srrcosfilter(N, alpha, Ts, Fs):
    """
    Generates a square-root raised cosine (SRRC) filter (FIR) impulse response.
    Parameters
    ----------
    N : int
        Length of the filter in samples.
    alpha : float
        Roll off factor (Valid values are [0, 1]).
    Ts : float
        Symbol period in seconds.
    Fs : float
        Sampling Rate in Hz.
    Returns
    -------
    time_idx : 1-D ndarray (float)
        Array containing the time indices, in seconds, for the impulse response.
    h_srrc : 1-D ndarray (float)
        Impulse response of the square-root raised cosine filter.
    """
    T_delta = 1 / float(Fs)
    time_idx = ((np.arange(N) - N / 2)) * T_delta
    sample_num = np.arange(N)
    h_srrc = np.zeros(N, dtype=float)

    for x in sample_num:
        t = (x - N / 2) * T_delta
        if t == 0.0:
            h_srrc[x] = (1 - alpha + (4 * alpha / np.pi)) / np.sqrt(Ts)
        elif alpha != 0 and t == Ts / (4 * alpha):
            h_srrc[x] = ((alpha * (np.pi + 2)) / np.sqrt(2 * Ts)) * \
                        ((np.sin(np.pi * (1 + alpha) * Ts / (4 * alpha)) + \
                          (np.cos(np.pi * (1 - alpha) * Ts / (4 * alpha))) / (4 * alpha * (1 + alpha) / Ts)))
        elif alpha != 0 and t == -Ts / (4 * alpha):
            h_srrc[x] = ((alpha * (np.pi + 2)) / np.sqrt(2 * Ts)) * \
                        ((np.sin(np.pi * (1 - alpha) * Ts / (4 * alpha)) + \
                          (np.cos(np.pi * (1 + alpha) * Ts / (4 * alpha))) / (4 * alpha * (1 - alpha) / Ts)))
        else:
            h_srrc[x] = (((np.sin(np.pi * (1 - alpha) * t / Ts) + \
                           (4 * alpha * t / Ts) * np.cos(np.pi * (1 + alpha) * t / Ts)) / \
                          (np.pi * t * (1 - (4 * alpha * t / Ts) * (4 * alpha * t / Ts)) / Ts)) / \
                         np.sqrt(Ts))

    return time_idx, h_srrc


def pulseShape(pulseType, SpS=2, N=512, alpha=0.8, Ts=1):
    """
    Generate a pulse shaping filter.
    Parameters
    ----------
    pulseType : string ('rect','nrz','rrc')
        type of pulse shaping filter.
    SpS : int, optional
        Number of samples per symbol of input signal. The default is 2.
    N : int, optional
        Number of filter coefficients. The default is 1024.
    alpha : float, optional
        Rolloff of RRC filter. The default is 0.1.
    Ts : float, optional
        Symbol period in seconds. The default is 1.
    Returns
    -------
    filterCoeffs : np.array
        Array of filter coefficients (normalized).
    """
    fa = (1 / Ts) * SpS

    t = np.linspace(-2, 2, SpS)
    Te = 1

    if pulseType == "rect":
        filterCoeffs = np.concatenate(
            (np.zeros(int(SpS / 2)), np.ones(SpS), np.zeros(int(SpS / 2)))
        )
    elif pulseType == "nrz":
        filterCoeffs = np.convolve(
            np.ones(SpS),
            2 / (np.sqrt(np.pi) * Te) * np.exp(-(t ** 2) / Te),
            mode="full",
        )
    elif pulseType == "rrc":
        tindex, filterCoeffs = rrcosfilter(N, alpha, Ts, fa)
    elif pulseType == "rc":
        tindex, filterCoeffs = rcosfilter(N, alpha, Ts, fa)
    elif pulseType == "srrc":
        tindex, filterCoeffs = srrcosfilter(N, alpha, Ts, fa)

    filterCoeffs = filterCoeffs / np.sqrt(np.sum(filterCoeffs ** 2))

    return filterCoeffs


def firFilter(h, x):
    """
    Perform FIR filtering and compensate for filter delay.
    Parameters
    ----------
    h : np.array
        Coefficients of the FIR filter (impulse response, symmetric).
    x : np.array
        Input signal.
    prec: cp.dtype
        Size of the complex representation.
    Returns
    -------
    y : np.array
        Output (filtered) signal.
    """
    try:
        x.shape[1]
    except IndexError:
        x = x.reshape(len(x), 1)
    y = x.copy()
    nModes = x.shape[1]

    for n in range(nModes):
        y[:, n] = np.convolve(x[:, n], h, mode="same")
    if y.shape[1] == 1:
        y = y[:, 0]
    return y


def mzm(Ai, u, Vπ, Vb):
    """
    Optical Mach-Zehnder Modulator (MZM).
    Parameters
    ----------
    Ai : scalar or np.array
        Amplitude of the optical field at the input of the MZM.
    u : np.array
        Electrical driving signal.
    Vπ : scalar
        MZM's Vπ voltage.
    Vb : scalar
        MZM's bias voltage.
    Returns
    -------
    Ao : np.array
        Modulated optical field at the output of the MZM.
    """
    pi = np.pi
    return Ai * np.cos(0.5 / Vπ * (u + Vb) * pi)


def linFiberCh(Ei, L, alpha, D, Fc, Fs):
    """
    Simulate signal propagation through a linear fiber channel.
    Parameters
    ----------
    Ei : np.array
        Input optical field.
    L : real scalar
        Length of the fiber.
    alpha : real scalar
        Fiber's attenuation coefficient in dB/km.
    D : real scalar
        Fiber's chromatic dispersion (2nd order) coefficient in ps/nm/km.
    Fc : real scalar
        Optical carrier frequency in Hz.
    Fs : real scalar
        Sampling rate of the simulation.
    Returns
    -------
    Eo : np.array
        Optical field at the output of the fiber.
    """
    # c  = 299792458   # speed of light [m/s](vacuum)
    c_kms = const.c / 1e3
    λ = c_kms / Fc
    α = alpha / (10 * np.log10(np.exp(1)))
    β2 = -(D * λ ** 2) / (2 * np.pi * c_kms)

    Nfft = len(Ei)

    ω = 2 * np.pi * Fs * fftfreq(Nfft)
    ω = ω.reshape(ω.size, 1)

    Nmodes = 1
    Ei = Ei.reshape(Ei.size, Nmodes)
    ω = np.tile(ω, (1, Nmodes))
    Eo = ifft(
        fft(Ei, axis=0) * np.exp(-α * L + 1j * (β2 / 2) * (ω ** 2) * L), axis=0
    )

    Eo = Eo.reshape(
        Eo.size,
    )

    return Eo


def eyediagramG(sigIn, Nsamples, SpS, n=3, ptype="fast", plotlabel=None):
    """
    Plot the eye diagram of a modulated signal waveform.
    :param Nsamples: number of samples to be plotted
    :param SpS: samples per symbol
    :param n: number of symbol periods
    :param plotlabel: label for the plot legend
    :return: A tuple containing the values A, B, C, D
    """
    sig = sigIn.copy()

    if not plotlabel:
        plotlabel = " "

 
    #plotlabel_ = plotlabel 
    y = sig[:Nsamples].real
    x = np.arange(0, y.size, 1) % (n * SpS)          
    plt.figure()
    y[x == n * SpS] = np.nan
    y[x == 0] = np.nan

    
  

    return x / SpS, y, min(x / SpS), max(x / SpS)


def awgn(sig, snr, Fs=1, B=1):
    """
    Implement an AWGN channel.
    Parameters
    ----------
    sig : np.array
        Input signal.
    snr : scalar
        Signal-to-noise ratio in dB.
    Fs : real scalar
        Sampling frequency. The default is 1.
    B : real scalar
        Signal bandwidth. The default is 1.
    Returns
    -------
    np.array
        Input signal plus noise.
    """
    snr_lin = 10 ** (snr / 10)
    noiseVar = signal_power(sig) / snr_lin
    σ = np.sqrt((Fs / B) * noiseVar)
    noise = normal(0, σ, sig.shape) + 1j * normal(0, σ, sig.shape)
    noise = 1 / np.sqrt(2) * noise

    return sig + noise


def sigPow(x):
    """
    Calculate the average power of x.
    Parameters
    ----------
    x : np.array
        Signal.
    Returns
    -------
    scalar
        Average power of x: P = Nmodes*mean(abs(x)**2).
    """
    return np.mean(np.abs(x) ** 2)


def signal_power(x):
    """
    Calculate the total average power of x.
    Parameters
    ----------
    x : np.array
        Signal.
    Returns
    -------
    scalar
        Total average power of x: P = Nmodes*mean(abs(x)**2).
    """
    try:
        Nmodes = x.shape[1]
    except IndexError:
        Nmodes = 1

    return Nmodes * sigPow(x)


def lowPassFIR(fc, fa, N, typeF="rect"):
    """
    Calculate FIR coefficients of a lowpass filter.
    Parameters
    ----------
    fc : float
        Cutoff frequency.
    fa : float
        Sampling frequency.
    N : int
        Number of filter coefficients.
    typeF : string, optional
        Type of response ('rect', 'gauss'). The default is "rect".
    Returns
    -------
    h : np.array
        Filter coefficients.
    """
    fu = fc / fa
    d = (N - 1) / 2
    n = np.arange(0, N)

    # calculate filter coefficients
    if typeF == "rect":
        h = (2 * fu) * np.sinc(2 * fu * (n - d))
    elif typeF == "gauss":
        h = (
                np.sqrt(2 * np.pi / np.log(2))
                * fu
                * np.exp(-(2 / np.log(2)) * (np.pi * fu * (n - d)) ** 2)
        )
    return h





def minEuclid(symb, const):
    """
    Find minimum Euclidean distance.
    Find closest constellation symbol w.r.t the Euclidean distance 
    """
    ind = np.zeros(symb.shape, dtype=np.int64)
    for ii in prange(len(symb)):
        ind[ii] = np.abs(symb[ii] - const).argmin()
    return ind


def Qfunc(x):
  
    return 0.5 - 0.5 * erf(x / np.sqrt(2))





def theoryBER(M, EbN0, constType):
    """
    Theoretical (approx.) bit error probability for PAM in AWGN channel.
    Parameters
    ----------
    M : int
        Modulation order.
    EbN0 : scalar
        Signal-to-noise ratio (SNR) per bit in dB.
    constType : string
        Modulation type: 'pam'
    Returns
    -------
    Pb : scalar
        Theoretical probability of bit error.
    """
    EbN0lin = 10 ** (EbN0 / 10)
    k = np.log2(M)
    if constType == "pam":
        Ps = (2 * (M - 1) / M) * Qfunc(
            np.sqrt(6 * np.log2(M) / (M ** 2 - 1) * EbN0lin)
        )
        Pb = Ps / k
    return Pb, Ps


class parameters:
    """
    Basic class to be used as a struct of parameters
    """
    pass

def demap(indSymb, bitMap):
    """
    Contellation symbol index to bit sequence demapping.
    Parameters
    ----------
    indSymb : np.array of ints
        Indexes of received symbol sequence.
    bitMap : (M, log2(M)) np.array
        bit-to-symbol mapping.
    Returns
    -------
    decBits : np.array
        Sequence of demapped bits.
    """
    M = bitMap.shape[0]
    b = int(np.log2(M))

    decBits = np.zeros(len(indSymb) * b, dtype="int")

    for i in prange(len(indSymb)):
        decBits[i * b : i * b + b] = bitMap[indSymb[i], :]
    return decBits

def demodulateGray(symb, M, constType):
    """
    Demodulate symbol sequences to bit sequences (w/ Gray mapping).
    Hard demodulation is based on minimum Euclidean distance.
    Parameters
    ----------
    symb : array of complex constellation symbols
        sequence of constellation symbols to be demodulated.
    M : int
        order of the modulation format.
    constType : string
         'pam' 
    Returns
    -------
    array of ints
        sequence of demodulated bits.
    """
 
    const = GrayMapping(M, constType)

    # get bit to symbol mapping
    indMap = minEuclid(const, const)
    bitMap = dec2bitarray(indMap, int(np.log2(M)))
    b = int(np.log2(M))
    bitMap = bitMap.reshape(-1, b)

    # demodulate received symbol sequence
    indrx = minEuclid(symb, const)

    return demap(indrx, bitMap)




def photodiode_V2(E, paramPD=None):
    """
    Pin photodiode (PD).
    Parameters
    ----------
    E : np.array
        Input optical field.
    paramPD : parameter object (struct), optional
        Parameters of the photodiode.
    paramPD.R: photodiode responsivity [A/W][default: 1 A/W]
    paramPD.IRND: Input-Refrred Noise Density [A/sqrt(Hz)] [Typical value: [10-20] pA/sqrt(Hz)]
    paramPD.RL: impedance load [Ω] [default: 50Ω]
    paramPD.B bandwidth [Hz][default: 30e9 Hz]
    paramPD.Fs: sampling frequency [Hz] [default: 60e9 Hz]
    paramPD.fType: frequency response type [default: 'rect']
    paramPD.N: number of the frequency resp. filter taps. [default: 8001]
    paramPD.ideal: ideal PD?(i.e. no noise, no frequency resp.) [default: True]
    paramPD.G_TIA: TIA Gain [V/A][default: 1 V/A]
    Returns
    -------
    ipd : np.array
          photocurrent.
    """
    q = const.value("elementary charge")
    kB = const.value("Boltzmann constant")
    if paramPD is None:
      
        paramPD = []
    # check input parameters
    R = getattr(paramPD, "R", 1)
    B = getattr(paramPD, "B", 30e9)
    Fs = getattr(paramPD, "Fs", 60e9)
    N = getattr(paramPD, "N", 8000)
    Id = getattr(paramPD, "Id", 5e-9)
    fType = getattr(paramPD, "fType", "rect")
    ideal = getattr(paramPD, "ideal", True)
    Tc = getattr(paramPD, "Tc", 25)
    kB = const.value("Boltzmann constant")
    RL = getattr(paramPD, "RL", 50)    
    G_TIA = getattr(paramPD, "G_T", 1)

    IRND = getattr(paramPD, "IRND", 10e-12)

    assert R > 0, "PD responsivity should be a positive scalar"
    assert (
        Fs >= 2 * B
    ), "Sampling frequency Fs needs to be at least twice of B."

    ipd = R * E    * np.conj(E)  # ideal fotodetected current

    if not (ideal):

        Pin = (np.abs(E) ** 2).mean()

        # Shot noise
        σ2_s = 2 * q * (R * Pin + Id) * B  # shot noise variance

        #TIA noise
        σ2_I = IRND**2 * B  # TIA noise variance

        # thermal noise
        T = Tc + 273.15  # temperature in Kelvin
        σ2_T = 4 * kB * T * B / RL  # thermal noise variance

        # add noise sources to the p-i-n receiver
        I_tia = normal(0, np.sqrt(σ2_I), ipd.size)
        It =  normal(0, np.sqrt(Fs * (σ2_T / (2 * B))), ipd.size)
        Is = normal(0, np.sqrt(σ2_s), ipd.size)

        ipd +=  I_tia  + Is # +  It

        # lowpass filtering
        h = lowPassFIR(B, Fs, N, typeF=fType)

        ipd = firFilter(h, ipd) 
        

    return ipd.real * G_TIA




def calculate_BER(SYMV_values, powerValues):
    # simulation parameters
    SpS = 4            # Samples per symbol
    M = 2               # order of the modulation format
    Rs = 40e9           # Symbol rate (for the PAM-2 case, Rs = Rb)
    Tsymb = 1/Rs        # Symbol period in seconds
    Fs = 1/(Tsymb/SpS)  # Signal sampling frequency (samples/second)
    Ts = 1/Fs           # Sampling period

    # MZM parameters
    Vπ = 2
    Vb = -Vπ/2

    # typical rc pulse
    pulse = pulseShape('rc', SpS)
    pulse = pulse/max(abs(pulse))

    # Number of bits
    Num_bits = 1e5

    BER = np.zeros((powerValues.shape[0], len(SYMV_values)))  # 2D array to store BER values for each SYMV

    constel = GrayMapping(M,'pam')  # get PAM constellation
    Es = signal_power(constel)  # calculate the average energy per symbol of the PAM constellation

    Beq = Rs    
    discard = 100
    IRND_c = 10e-12
    R = 1

    for indSYMV, SYMV in enumerate(SYMV_values):
        for indPi, Pi_dBm in enumerate(tqdm(powerValues)):
            Pi = 10**((Pi_dBm)/10)*1e-3  # optical signal power in W at the MZM input

            # generate pseudo-random bit sequence
            bitsTx = np.random.randint(2, size=int(np.log2(M)*Num_bits))
            n = np.arange(0, bitsTx.size)

            # generate pam-m modulated symbol sequence
            symbTx = modulateGray(bitsTx, M, 'pam')    
            # symbTx = pnorm(symbTx) # power normalization

            # upsampling
            symbolsUp = upsample(symbTx, SpS)

            # pulse formatting
            sigTx = firFilter(pulse, symbolsUp)

            sigTx = pnorm(sigTx) # power normalization

            # optical modulation
            original_vector = sigTx

            new_min = -SYMV
            new_max = SYMV

            # OMA = (ER-1)/(ER+1)
            mapped_vector = (original_vector - np.min(original_vector)) * (new_max - new_min) / (np.max(original_vector) - np.min(original_vector)) + new_min
            sigTx = mapped_vector

            Ai = np.sqrt(Pi)*np.ones(sigTx.size)
            sigTxo = mzm(Ai, 0.25*sigTx, Vπ, Vb)

            ER = np.max(sigTxo)/np.min(sigTxo)

            Prx = signal_power(sigTxo)

            SER_TIA = SER_TIA_PIN(M, IRND_c, Prx, ER, Beq, R)

            # pin receiver
            paramPD = parameters()
            paramPD.ideal = False
            paramPD.B = Rs
            paramPD.Fs = Fs

            # I_Rx = photodiode(sigTxo.real, paramPD)
            I_Rx = photodiode_V2(sigTxo.real, paramPD) # v2 PHOTODIODE WITH TIA + SHOT + THERMAL NOISE 

            I_Rx = I_Rx/np.std(I_Rx)

            # capture samples in the middle of signaling intervals
            symbRx = I_Rx[0::SpS]

            # subtract DC level and normalize power
            symbRx = symbRx - symbRx.mean()
            symbRx = pnorm(symbRx)

            snr = signal_power(symbRx)/(2*signal_power(symbRx-symbTx))
            EbN0 = 10*np.log10(snr/np.log2(M))

            # demodulate symbols to bits with minimum Euclidean distance 
            bitsRx = demodulateGray(np.sqrt(Es)*symbRx, M, 'pam')

            err = np.logical_xor(bitsRx[discard:bitsRx.size-discard], bitsTx[discard:bitsTx.size-discard])
            BER[indPi, indSYMV] = np.mean(err)

    return BER





def edfa(Ei, Fs, G=20, NF=4.5, Fc=193.1e12):
    """
    Implement simple EDFA model.

    Parameters
    ----------
    Ei : np.array
        Input signal field.
    Fs : scalar
        Sampling frequency in Hz.
    G : scalar, optional
        Amplifier gain in dB. The default is 20.
    NF : scalar, optional
        EDFA noise figure in dB. The default is 4.5.
    Fc : scalar, optional
        Central optical frequency. The default is 193.1e12.

    Returns
    -------
    Eo : np.array
        Amplified noisy optical signal.

    """
    assert G > 0, "EDFA gain should be a positive scalar"
    assert NF >= 3, "The minimal EDFA noise figure is 3 dB"

    NF_lin = 10 ** (NF / 10)
    G_lin = 10 ** (G / 10)
    nsp = (G_lin * NF_lin - 1) / (2 * (G_lin - 1))

  

    N_ase = (G_lin - 1) * nsp * const.h * Fc
    p_noise = N_ase * Fs

    noise = normal(0, np.sqrt(p_noise / 2), Ei.shape)
    return Ei * np.sqrt(G_lin) + noise




def theoryBERc(M, EbN0, constType):
    EbN0lin = 10 ** (EbN0 / 10)
    k = np.log2(M)
  
    Ps = (2 * (M - 1) / M) * Qfunc(
        np.sqrt(6 * np.log2(M) / (M**2 - 1) * EbN0lin)
    )
    Pb = Ps / k
    return Pb



def lms_equalizer(I_Tx, I_Rx, num_taps, step_size, num_iterations):
    """
    LMS Equalizer Function
    
    Applies the Least Mean Squares (LMS) algorithm for equalization.
    
    Parameters:
        I_Rx (numpy.ndarray): Received signal.
        I_Tx (numpy.ndarray): Transmitted signal.
        num_taps (int): Number of taps for the equalizer.
        step_size (float): Step size or learning rate for the LMS update rule.
        num_iterations (int): Number of iterations for the LMS algorithm.
        
    Returns:
        equalized_output (numpy.ndarray): Equalized output signal.
    """
    # Initialization
    tap_weights = np.zeros(num_taps)
    equalized_output = np.zeros_like(I_Tx)

    # Apply LMS algorithm
    for iteration in range(num_iterations):
        # Iterate through each symbol
        for i in range(num_taps, len(I_Tx)):
            # Extract the current window of the signal
            window = I_Tx[i - num_taps : i]
            
            # Compute the equalized output for the current symbol
            equalized_output[i] = np.dot(window, tap_weights)
            
            # Update the tap weights using the LMS update rule
            error = I_Rx[i] - equalized_output[i]
            tap_weights += step_size * error * window

    return equalized_output








