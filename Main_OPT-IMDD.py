
"""
Created on Mon Jun  5 20:43:28 2023

@author: Peyman, Ale
"""

import tkinter as tk
from tkinter import ttk
import tkinter as tk
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from optimdd import *







########################################
'''
 GUI

'''

#########################################






import tkinter as tk
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk




# Define global variables to store the input parameters and generated signals

global I_Rx_ideal, sigTxo_c, sigCh_c



SpS = None
M = None
Rs = None
I_Rx_ideal = None
ERdB = None
Symbol_ER  = None

def get_parameters():
    global SpS, M, Rs, L, CD, alpha, SYMV, ERdB, Pi_dBm, sigTxo, sigCh, SYMV_values, powerValues, photodiode_choice, paramPD
    global IRND_inp, Resp, NF, Symbol_ER, Fs, Gain, PD_BW
    
    Pi_dBm = float(entry_Pi_dBm.get())
    SpS = int(entry_SpS.get())
    M = int(entry_M.get())
    Rs = float(entry_Rs.get())
    L = float(entry_L.get())
    CD = float(entry_CD.get())
    alpha = float(entry_alpha.get())
    SYMV = float(entry_SYMV.get())
    
    IRND_inp = float(entry_IRND_inp.get())
    Resp = float(entry_Resp.get())
    NF = float(entry_NF.get())
    Gain = float(entry_Gain.get())
    
    PD_BW = float(entry_PD_BW.get())
    
    

    
    
    


def generate_signals():
    global I_Rx_ideal, I_Rx, sigTx, CD, alpha, SYMV, ERdB, Pi_dBm, sigTxo, sigCh, SYMV_values, powerValues, photodiode_choice, paramPD, IRND_inp
    global Resp, NF, G, Symbol_ER , err, bitsRx, constel, snr, EbN0, Fs, Rs, M, pulse, bit_ER_sim, err_hist
    global sigCh_c, sigTxo_c, PD_BW
    # Check if the necessary parameters have been provided
    if SpS is None or M is None or Rs is None:
        print("Please enter values for SpS, M, and Rs.")
        return
    
    # simulation parameters
    Tsymb = 1 / Rs       # Symbol period in seconds
    Fs = SpS * Rs       # Signal sampling frequency (samples/second)
    Ts = 1 / Fs         # Sampling period

    # MZM parameters
    Vπ = 2
    Vb = -Vπ/2
    #Pi_dBm = 20  # laser optical power at the input of the MZM in dBm
    Pi = 10**(Pi_dBm/10) * 1e-3  # convert from dBm to W
    
    constel = GrayMapping(M,'pam') # get PAM constellation
    Es = signal_power(constel) # calculate the average energy per symbol of the PAM constellation
    


    # generate pseudo-random bit sequence
    bitsTx = np.random.randint(2, size=int(np.log2(M) * 1e4))

    # generate ook modulated symbol sequence
    symbTx = modulateGray(bitsTx, M, 'pam')


    # upsampling
    symbolsUp = upsample(symbTx, SpS)

    # Typical RC pulse
    pulse = pulseShape('nrz', SpS)
    pulse = pulse / max(abs(pulse))

    # pulse shaping
    sigTx = firFilter(pulse, symbolsUp)
    
    sigTxo_c = sigTx
    
    sigTx = pnorm(sigTx)
    original_vector = sigTx
    new_min = -SYMV
    new_max = SYMV
    mapped_vector = (original_vector - np.min(original_vector)) * (new_max - new_min) / (np.max(original_vector) - np.min(original_vector)) + new_min
    sigTx = mapped_vector
    
    
    


    # optical modulation
    Ai = np.sqrt(Pi) * np.ones(sigTx.size)
    sigTxo = mzm(Ai, 0.25 * sigTx, Vπ, Vb)


 
    Fc = 193.1e12  # central optical frequency [Hz]
    
    sigCh = linFiberCh(sigTxo, L, alpha, CD, Fc, Fs)
    
    
    
    
    # EDFA
    # G = alpha*L    # edfa gain
    G = Gain
    NF = 4.5   # edfa noise figure
    sigCh = edfa(sigCh, Fs, G, NF, Fc)
    
    

    
    # ideal photodiode (noiseless, no bandwidth limitation)
    paramPD = parameters()
    paramPD.ideal = True
    
    paramPD.B = Rs
    paramPD.Fs = Fs
    paramPD.R = Resp
    

    I_Rx_ideal = photodiode_V2(sigTxo.real, paramPD)
    
    
    

    
    
    
    ExR = np.max(I_Rx_ideal)/np.min(I_Rx_ideal)
    
    ERdB = 10*math.log10(ExR)
    
    
    # noisy photodiode (thermal noise + shot noise + bandwidth limitation)
    paramPD = parameters()
    paramPD.ideal = False#False # False #photodiode_choice#True#False
    
    
    paramPD.B = PD_BW
    paramPD.Fs = Fs
    paramPD.IRND = IRND_inp * 1e-12
    paramPD.R = Resp
   # paramPD.fType = 'rc'
    
  
    
 

    I_Rx = photodiode_V2(sigCh, paramPD)


 
    
    Prx = signal_power(I_Rx )


    I_Rxb = I_Rx/np.std(I_Rx)
    sigCh_c = I_Rxb - I_Rxb.mean()
    


    symbRx = I_Rxb[0::SpS]
    
    
    # Downsampling with optimal instances
    symbol_count = len(I_Rx) // SpS
    symbRx_hist = np.zeros((SpS, symbol_count), dtype=I_Rx.dtype)
    err_hist = np.zeros((SpS, 1), dtype=I_Rx.dtype)
    
    for i in range(symbol_count):
        start_idx = i * SpS
        end_idx = start_idx + SpS
        symbRx_hist[:, i] = I_Rx[start_idx:end_idx]
        
    for c in range(SpS):
      symbRx1 = symbRx_hist[c,:] - symbRx_hist[c,:].mean()
      symbRx1 = pnorm(symbRx1)       
      bitsRx1 = demodulateGray(np.sqrt(Es)*symbRx1, M, 'pam')  
      discard = 100
      err1 = np.logical_xor(bitsRx1[discard:bitsRx1.size-discard], bitsTx[discard:bitsTx.size-discard])
      err_hist[c,0] = np.mean(err1)
        
    

   
    symbRx = symbRx - symbRx.mean()
    symbRx = pnorm(symbRx)
    snr = signal_power(symbRx)/(2*signal_power(symbRx-symbTx))
    EbN0 = 10*np.log10(snr/np.log2(M))
    
    bitsRx = demodulateGray(np.sqrt(Es)*symbRx, M, 'pam')
    discard = 100
    err = np.logical_xor(bitsRx[discard:bitsRx.size-discard], bitsTx[discard:bitsTx.size-discard])
    bit_ER_sim = np.mean(err)
    
    

        

def generate_signals_vec(L_vect):
    global I_Rx_ideal, I_Rx, sigTx, CD, alpha, SYMV, ERdB, Pi_dBm, sigTxo, sigCh, SYMV_values, powerValues, photodiode_choice, paramPD, IRND_inp
    global Resp, NF, G, Symbol_ER , err, bitsRx, constel, snr, EbN0, Fs, Rs, M, pulse, Symbol_ER_sim, err_hist
    global opt_sample_instances
    # Check if the necessary parameters have been provided
    if SpS is None or M is None or Rs is None:
        print("Please enter values for SpS, M, and Rs.")
        return
    
    # simulation parameters
    Tsymb = 1 / Rs       # Symbol period in seconds
    Fs = SpS * Rs       # Signal sampling frequency (samples/second)
    Ts = 1 / Fs         # Sampling period

    # MZM parameters
    Vπ = 2
    Vb = -Vπ/2
    #Pi_dBm = 20  # laser optical power at the input of the MZM in dBm
    Pi = 10**(Pi_dBm/10) * 1e-3  # convert from dBm to W
    
    constel = GrayMapping(M,'pam') # get PAM constellation
    Es = signal_power(constel) # calculate the average energy per symbol of the PAM constellation
    


    # generate pseudo-random bit sequence
    bitsTx = np.random.randint(2, size=int(np.log2(M) * 1e4))

    # generate ook modulated symbol sequence
    symbTx = modulateGray(bitsTx, M, 'pam')


    # upsampling
    symbolsUp = upsample(symbTx, SpS)

    # typical RC pulse
    pulse = pulseShape('rc', SpS)
    pulse = pulse / max(abs(pulse))

    # pulse shaping
    sigTx = firFilter(pulse, symbolsUp)
    
    sigTx = pnorm(sigTx)
    original_vector = sigTx
    new_min = -SYMV
    new_max = SYMV
    mapped_vector = (original_vector - np.min(original_vector)) * (new_max - new_min) / (np.max(original_vector) - np.min(original_vector)) + new_min
    sigTx = mapped_vector
    


    # optical modulation
    Ai = np.sqrt(Pi) * np.ones(sigTx.size)
    sigTxo = mzm(Ai, 0.25 * sigTx, Vπ, Vb)
    


 
    Fc = 193.1e12  # central optical frequency [Hz]
    
    sigCh = linFiberCh(sigTxo, L_vect, alpha, CD, Fc, Fs)
    
    
    # EDFA
    # G = alpha*L    # edfa gain
    # G = Gain
    # NF = 4.5   # edfa noise figure
    sigCh = edfa(sigCh, Fs, G, NF, Fc)
    
    
    # ideal photodiode (noiseless, no bandwidth limitation)
    paramPD = parameters()
    paramPD.ideal = True
    paramPD.R = Resp
    I_Rx_ideal = photodiode_V2(sigTxo.real, paramPD)
    
    
    
    ExR = np.max(I_Rx_ideal)/np.min(I_Rx_ideal)
    
    ERdB = 10*math.log10(ExR)
    
    
    # noisy photodiode (thermal noise + shot noise + bandwidth limitation)
    paramPD = parameters()
    paramPD.ideal = False#False # False #photodiode_choice#True#False
    paramPD.B = Rs
    paramPD.Fs = Fs
    paramPD.IRND = IRND_inp * 1e-12
    paramPD.R = Resp
   # paramPD.fType = 'rc'
    
  
    
 

    I_Rx = photodiode_V2(sigCh, paramPD)
    
    Prx = signal_power(sigCh)
    Prx_dbm = (10*np.log10(signal_power(Prx)/1e-3))

   
    I_Rxb = I_Rx/np.std(I_Rx)

    symbRx = I_Rxb[0::SpS]
    
    
    # Downsampling with optimal instances
    symbol_count = len(I_Rx) // SpS
    symbRx_hist = np.zeros((SpS, symbol_count), dtype=I_Rx.dtype)
    err_hist = np.zeros((SpS, 1), dtype=I_Rx.dtype)   
    
    for i in range(symbol_count):
        start_idx = i * SpS
        end_idx = start_idx + SpS
        symbRx_hist[:, i] = I_Rx[start_idx:end_idx]
        
    for c in range(SpS):
      symbRx1 = symbRx_hist[c,:] - symbRx_hist[c,:].mean()
      symbRx1 = pnorm(symbRx1)       
      bitsRx1 = demodulateGray(np.sqrt(Es)*symbRx1, M, 'pam')  
      discard = 100
      err1 = np.logical_xor(bitsRx1[discard:bitsRx1.size-discard], bitsTx[discard:bitsTx.size-discard])
      err_hist[c,0] = np.mean(err1)
    
    
    opt_sample_instances = np.min(err_hist)
    
    
    
    
    
    
    symbRx = symbRx - symbRx.mean()
    symbRx = pnorm(symbRx)
    snr = signal_power(symbRx)/(2*signal_power(symbRx-symbTx))
    EbN0 = 10*np.log10(snr/np.log2(M))
    
    bitsRx = demodulateGray(np.sqrt(Es)*symbRx, M, 'pam')
    discard = 100
    err = np.logical_xor(bitsRx[discard:bitsRx.size-discard], bitsTx[discard:bitsTx.size-discard])
    bit_ER_sim = np.mean(err)
    
    return opt_sample_instances, Prx_dbm
    
    

                
        
        
    
    
    
    
    
    
    

def plot_eyediagram():
    plt.subplots_adjust(hspace=0.5)
    # Check if I_Rx_ideal has been generated
    if I_Rx_ideal is None:
        print("Please generate the signals first.")
        return
    
    # Clear the previous plot, if any
    ax1.cla()
    ax2.cla()
    ax3.cla()
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    
    # Call the eyediagram function with the updated parameters
    discard = 100
    A, B, C, D = eyediagramG(sigTx[discard:-discard], I_Rx.size-2*discard, SpS)
    A1, B1, C1, D1 = eyediagramG(I_Rx_ideal[discard:-discard], I_Rx.size-2*discard, SpS)
    A2, B2, C2, D2 = eyediagramG(I_Rx[discard:-discard], I_Rx.size-2*discard, SpS)
    # Plot in the first section
  
    ax1.set_xlim(-3*Rs,3*Rs);
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('PSD of the Optical Signal')
    ax1.psd(np.abs(sigTxo)**2, Fs=Fs, NFFT = 16*1024, sides='twosided', label = 'Optical signal spectrum')
    
    # Plot in the second section
   
    ax2.plot(A1, B1*1000 , color="blue", alpha=0.9)
    ax2.set_title('Eyediagram - Ideal Fiber')
    ax2.set_xlabel('Symbol Time')
    ax2.set_ylabel('Current (mA)')
    ax2.set_xlim(C1, D1)
    
    # Plot in the third section
    ax3.plot(A2, B2*1000, color="red", alpha=0.9)
    ax3.set_title('Eyediagram - Non-Ideal Fiber')
    ax3.set_ylabel('Current (mA)')
    ax3.set_xlim(C2, D2)
    
    # Refresh the canvas
    canvas.draw_idle()
    update_output_text()
 
 

import numpy as np
import matplotlib.pyplot as plt

def BER_plotter(SpS, M, SYMV_values):
    global ERdB
    Rs = 40e9
    Tsymb = 1/Rs
    Fs = 1/(Tsymb/SpS)
    Ts = 1/Fs
    Vπ = 2
    Vb = -Vπ/2
    pulse = pulseShape('rc', SpS)
    pulse = pulse/max(abs(pulse))
    Num_bits = 1e5
    SYMV_values = np.array(SYMV_values)
    powerValues = np.arange(-30, -5)
    BER = np.zeros((powerValues.shape[0], len(SYMV_values)))
    Pb = np.zeros(powerValues.shape)
    constel = GrayMapping(M, 'pam')
    Es = signal_power(constel)
    Beq = Rs
    discard = 100
    IRND_c = 10e-12
    R = 1

    for indSYMV, SYMV in enumerate(SYMV_values):
        for indPi, Pi_dBm in enumerate(powerValues):
            Pi = 10**((Pi_dBm)/10)*1e-3

            bitsTx = np.random.randint(2, size=int(np.log2(M)*Num_bits))
            n = np.arange(0, bitsTx.size)
            symbTx = modulateGray(bitsTx, M, 'pam')
            symbolsUp = upsample(symbTx, SpS)
            sigTx = firFilter(pulse, symbolsUp)
            sigTx = pnorm(sigTx)

            original_vector = sigTx
            new_min = -SYMV
            new_max = SYMV
            mapped_vector = (original_vector - np.min(original_vector)) * (new_max - new_min) / (np.max(original_vector) - np.min(original_vector)) + new_min
            sigTx = mapped_vector

            Ai = np.sqrt(Pi)*np.ones(sigTx.size)
            sigTxo = mzm(Ai, 0.25*sigTx, Vπ, Vb)
            
            ExR = np.max(sigTxo)/np.min(sigTxo)
            ERdB = 10*math.log10(ExR)
            
            ER = np.max(sigTxo)/np.min(sigTxo)
            Prx = signal_power(sigTxo)

            paramPD = parameters()
            paramPD.ideal = False
            paramPD.B = Rs
            paramPD.Fs = Fs
            I_Rx = photodiode_V2(sigTxo.real, paramPD)
            
           
            
            
            I_Rx = I_Rx/np.std(I_Rx)
            
            

            symbRx = I_Rx[0::SpS]
       
            symbRx = symbRx - symbRx.mean()
            symbRx = pnorm(symbRx)
            snr = signal_power(symbRx)/(2*signal_power(symbRx-symbTx))
            EbN0 = 10*np.log10(snr/np.log2(M))
            bitsRx = demodulateGray(np.sqrt(Es)*symbRx, M, 'pam')
            err = np.logical_xor(bitsRx[discard:bitsRx.size-discard], bitsTx[discard:bitsTx.size-discard])
            BER[indPi, indSYMV] = np.mean(err)
            Pb[indPi] = theoryBER(M, EbN0, 'pam')[0]

    return powerValues, np.log10(BER)











def open_analysis_window():
    global SpS, M, SYMV_values
    analysis_window = tk.Toplevel(root)
    analysis_window.title("BER Analysis Window")
    def get_parameters_ber():
        global SpS, M, SYMV_values, log_BER
        
        SpS = int(entry_SpS.get())
        M = int(entry_M.get())
        SYMV_values = [float(entry_SYMV_values.get())]
      
        
      
        
    
    
    label_SpS = tk.Label(analysis_window, text="SpS:")
    label_SpS.pack()
    entry_SpS = tk.Entry(analysis_window)
    entry_SpS.pack()
    
    
    
    label_SYMV_values = tk.Label(analysis_window, text="V_pp/2:")
    label_SYMV_values.pack()
    entry_SYMV_values = tk.Entry(analysis_window)
    entry_SYMV_values.pack()
    
    label_M = tk.Label(analysis_window, text="PAM-M:")
    label_M.pack()
    entry_M = tk.Entry(analysis_window)
    entry_M.pack()
    


    params_button_ber = tk.Button(analysis_window, text="Get Parameters", command=get_parameters_ber)
    params_button_ber.pack()

  
    
    plt.style.use("dark_background")
    
    # Create subplots for each section of the plot
    fig, (ax11) = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    
    # Create a matplotlib canvas to display the plot inside the plot frame
    canvas = FigureCanvasTkAgg(fig, master=analysis_window)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    # Add the navigation toolbar
    toolbar = NavigationToolbar2Tk(canvas, analysis_window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    
    
    
    
    
    def BER_plotter_in():
        global ERdB
        
    
        # Set the dark_background style for the plot
      
        plt.subplots_adjust(hspace=0.5)
        # Check if I_Rx_ideal has been generated
        
     
        
        # Clear the previous plot, if any
        # ax11.cla()
        ax11.grid(True)
        
    

        powerValues, log_BER = BER_plotter(SpS, M, SYMV_values)
    
        # Plot the results
        for indSYMV, SYMV in enumerate(SYMV_values):
            ax11.plot(powerValues, log_BER[:, indSYMV],'-o' ,alpha=0.9, label='ER(dB)'+str(np.round(ERdB,2)))
        ax11.set_xlabel('Input Power (dBm)')
        ax11.set_ylabel('Log10 BER')
        ax11.set_title('BER Analysis')
        ax11.legend(loc='upper right')
    
     
        canvas.draw_idle()
        update_output_text()
        
        
    params_button_ber = tk.Button(analysis_window, text="Plot BER", command=BER_plotter_in)
    params_button_ber.pack()
    
    

    
    
 



def open_analysis_eq():
    global num_taps, step_size, num_iterations, I_Rx_e, I_Tx_e, sigCh, sigTxo
    analysis_window = tk.Toplevel(root)
    analysis_window.title("Equalizer Analysis Window")
    def get_parameters_ber():
        global num_taps, step_size, num_iterations, I_Rx_e, I_Rx_e
        
        num_taps = int(entry_num_taps.get())
        step_size = float(entry_step_size.get())
        num_iterations = int(entry_num_iterations.get())

      
        
      
        
    
    
    label_num_taps = tk.Label(analysis_window, text="num_taps:")
    label_num_taps.pack()
    entry_num_taps = tk.Entry(analysis_window)
    entry_num_taps.pack()
    
    
    
    label_step_size = tk.Label(analysis_window, text="step_size")
    label_step_size.pack()
    entry_step_size = tk.Entry(analysis_window)
    entry_step_size.pack()
    
    label_num_iterations = tk.Label(analysis_window, text="num_iterations")
    label_num_iterations.pack()
    entry_num_iterations = tk.Entry(analysis_window)
    entry_num_iterations.pack()
    


    params_button_ber = tk.Button(analysis_window, text="Get Parameters", command=get_parameters_ber)
    params_button_ber.pack()
     
    
    plt.style.use("dark_background")
    # Create subplots for each section of the plot
    fig, (ax11q, ax12q) = plt.subplots(nrows=2, ncols=1, figsize=(5, 5))
    
 
  
    
    # Create a matplotlib canvas to display the plot inside the plot frame
    canvas = FigureCanvasTkAgg(fig, master=analysis_window)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    # Add the navigation toolbar
    toolbar = NavigationToolbar2Tk(canvas, analysis_window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    
    def EQ_plotter_in():
        global num_taps, step_size, num_iterations, I_Rx_e, I_Tx_e, sigTxo_c, sigCh_c
        
        ax11q.cla()
        ax11q.grid(True)
        ax12q.cla()
        ax12q.grid(True)
        
    
        # Set the dark_background style for the plot
        plt.subplots_adjust(hspace=0.5)
        I_Tx_e =  sigTxo_c.real
        I_Rx_e =  sigCh_c.real 
        
        # Apply the LMS equalizer
        equalized_output = lms_equalizer(I_Tx_e, I_Rx_e, num_taps, step_size, num_iterations)
        discard = 400
        
        A1q, B1q, C1q, D1q = eyediagramG(I_Rx_e[discard:-discard], I_Rx_e.size-2*discard, SpS)
        A2q, B2q, C2q, D2q = eyediagramG(equalized_output[discard:-discard], I_Rx_e.size-2*discard, SpS)
       
        ax11q.plot(A1q, B1q, color="blue", alpha=0.9)
        ax11q.set_title('Before equalizer')
        ax11q.set_xlabel('Symbol Time')
        ax11q.set_ylabel('Current (mA)')
        ax11q.set_xlim(C1q, D1q)
        
        ax12q.plot(A2q, B2q, color="blue", alpha=0.9)
        ax12q.set_title('After equalizer')
        ax12q.set_xlabel('Symbol Time')
        ax12q.set_ylabel('Current (mA)')
        ax12q.set_xlim(C1q, D1q)
        
        
        canvas.draw_idle()
        update_output_text()
        
        
    params_button_ber = tk.Button(analysis_window, text="Adaptive Equalizer", command=EQ_plotter_in)
    params_button_ber.pack()
    

    
    
    
    
    
    
    


def open_analysis_vec():
    global num_taps, step_size, I_Rx_e, I_Tx_e, startm, endm, L_vector
    analysis_window = tk.Toplevel(root)
    analysis_window.title("Distance")
    def get_parameters_ber():
        global num_taps, step_size, I_Rx_e, I_Tx_e, startm, endm, L_vector
        
        startm = int(entry_startm.get())
        endm = int(entry_endm.get())
        step_size = float(entry_step_size.get())
        L_vector = np.arange(startm, endm, step_size)

       

      
        
      
        
    
    
    label_startm = tk.Label(analysis_window, text="starting distance (km):")
    label_startm.pack()
    entry_startm = tk.Entry(analysis_window)
    entry_startm.pack()
    
    label_endm = tk.Label(analysis_window, text="ending distnace (km):")
    label_endm.pack()
    entry_endm = tk.Entry(analysis_window)
    entry_endm.pack()
    
  

    
    label_step_size = tk.Label(analysis_window, text="step_size (km)")
    label_step_size.pack()
    entry_step_size = tk.Entry(analysis_window)
    entry_step_size.pack()
    
  


    params_button_ber = tk.Button(analysis_window, text="Get Parameters", command=get_parameters_ber)
    params_button_ber.pack()
    
    
    
     
    
    plt.style.use("dark_background")
    # Create subplots for each section of the plot
    fig, (ax11q, ax12q) = plt.subplots(nrows=2, ncols=1, figsize=(5, 5))
    
 
  
    
    # Create a matplotlib canvas to display the plot inside the plot frame
    canvas = FigureCanvasTkAgg(fig, master=analysis_window)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    # Add the navigation toolbar
    toolbar = NavigationToolbar2Tk(canvas, analysis_window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    
    
    
    def EQ_plotter_in():
        global num_taps, step_size, L_vector
        
        bit_er_sim_array = []
        prx_array = []
        

        for value in L_vector:
            bit_ER_sim1 , Prx1 = generate_signals_vec(value)
            
            bit_er_sim_array.append(bit_ER_sim1)
            prx_array.append(Prx1)
            
            
                
        ax11q.cla()
        ax11q.grid(True)
        ax12q.cla()
        ax12q.grid(True)
        
    
        # Set the dark_background style for the plot
        plt.subplots_adjust(hspace=0.5)
      
        
      
        
       
        ax11q.plot(L_vector, np.log10(bit_er_sim_array),'-o',color="green", alpha=0.9)
        ax11q.set_title('Log(BER)')
        ax11q.set_xlabel('Distance (km)')
        ax11q.set_ylabel('Log(BER)')
        ax11q.set_xlim(startm,endm)
        
    
        
        ax12q.plot(L_vector, prx_array, '-o',color="red",alpha=0.9)
        ax12q.set_title('Received power')
        ax12q.set_xlabel('Distance (km)')
        ax12q.set_ylabel('P_Rx (dBm)')
        ax12q.set_xlim(startm,endm)
    
        
        
        canvas.draw_idle()
        update_output_text()
        
        
    params_button_ber = tk.Button(analysis_window, text="Plot", command=EQ_plotter_in)
    params_button_ber.pack()
    
    
    
     






# Create the GUI
root = tk.Tk()
root.title("Eyediagram Plotter")

root.tk.call("source", "azure.tcl")
root.tk.call("set_theme", "light")


plot_frame = tk.Frame(root)
plot_frame.grid(row=0, column=0, padx=(20, 0), pady=20, sticky="nsew")


text_frame = tk.Frame(root)
text_frame.grid(row=0, column=2, padx=20, pady=15, sticky="nsew")

text_frameT = tk.Frame(root)

text_frameT.grid(row=0, column=3, padx=20, pady=15, sticky="nsew")







output_text = tk.Text(text_frameT, height=30, width=30)
output_text.grid(row=0, column=1, padx=2, pady=5)
output_text.configure(state="disabled")  # Disable editing

# Set a fixed size for the widget



input_frame = tk.Frame(root)
input_frame.grid(row=0, column=1, padx=20, pady=5, sticky="nsew")




label_frame = tk.LabelFrame(input_frame, text="MZM Parameters", padx=10, pady=10, bg="white")
label_frame.pack()




label_Pi_dBm = tk.Label(label_frame, text="Laser optical power at the input of the MZM in dBm")
label_Pi_dBm.pack()
entry_Pi_dBm = tk.Entry(label_frame)
entry_Pi_dBm.pack()






#Vπ = 2
#Vb = -Vπ/2

label_Vπ = tk.Label(label_frame, text="Vπ (Volt):")
label_Vπ.pack()
entry_Vπ = tk.Entry(label_frame)
entry_Vπ.pack()

# Set the default value
default_valueVpi = "2"
entry_Vπ.insert(tk.END, default_valueVpi)
#entry_Vπ.configure(state="disabled")



label_Vb = tk.Label(label_frame, text="Bias Voltage (Volt):")
label_Vb.pack()

entry_Vb = tk.Entry(label_frame)
entry_Vb.pack()

# Set the default value
default_valueVb = "-1"
entry_Vb.insert(tk.END, default_valueVb)
#entry_Vb.configure(state="disabled")


label_frame_s = tk.LabelFrame(input_frame, text="TX signal Parameters", padx=10, pady=10, bg="white")
label_frame_s.pack()

# Create input fields for SpS, M, and Rs
label_SpS = tk.Label(label_frame_s, text="SpS:")
label_SpS.pack()
entry_SpS = tk.Entry(label_frame_s)
entry_SpS.pack()

label_M = tk.Label(label_frame_s, text="PAM-M:")
label_M.pack()
entry_M = tk.Entry(label_frame_s)
entry_M.pack()

label_Rs = tk.Label(label_frame_s, text="Rs Symbol rate:")
label_Rs.pack()
entry_Rs = tk.Entry(label_frame_s)
entry_Rs.pack()

# Set the default value
default_value = "10e9"
entry_Rs.insert(tk.END, default_value)


label_frame_fb= tk.LabelFrame(input_frame, text="Fiber Parameters", padx=10, pady=10, bg="white")
label_frame_fb.pack()



label_L = tk.Label(label_frame_fb, text="Distance [km]:")
label_L.pack()
entry_L = tk.Entry(label_frame_fb)
entry_L.pack()


label_CD = tk.Label(label_frame_fb, text="Chromatic Dispersion [ps/nm/km]:")
label_CD.pack()
entry_CD = tk.Entry(label_frame_fb)
entry_CD.pack()


label_alpha = tk.Label(label_frame_fb, text="Loss [dB/km]:")
label_alpha.pack()
entry_alpha = tk.Entry(label_frame_fb)
entry_alpha.pack()




label_SYMV = tk.Label(label_frame_s, text="V_pp/2 (0<V_pp<4):")
label_SYMV.pack()
entry_SYMV = tk.Entry(label_frame_s)
entry_SYMV.pack()


label_frame_pd= tk.LabelFrame(input_frame, text="Photodiode Parameters", padx=10, pady=10, bg="white")
label_frame_pd.pack()



label_IRND_inp = tk.Label(label_frame_pd, text="IRND[pA/sqrt(Hz)]:")
label_IRND_inp.pack()
entry_IRND_inp = tk.Entry(label_frame_pd)
entry_IRND_inp.pack()

# Set the default value
default_valueIRND = "10"
entry_IRND_inp.insert(tk.END, default_valueIRND)



label_Resp = tk.Label(label_frame_pd, text="Responsivity [A/W][default: 1 A/W]")
label_Resp.pack()
entry_Resp = tk.Entry(label_frame_pd)
entry_Resp.pack()

# Set the default value
default_valueR = "1"
entry_Resp.insert(tk.END, default_valueR)






label_PD_BW= tk.Label(label_frame_pd, text="Lowpass Filter BW: Rs]")
label_PD_BW.pack()
entry_PD_BW = tk.Entry(label_frame_pd)
entry_PD_BW.pack()

# Set the default value
default_PD_BWr = '10e9'
entry_PD_BW.insert(tk.END, default_PD_BWr)



label_frame_edf = tk.LabelFrame(text_frame, text="EDFA Parameters", padx=10, pady=10, bg="white")
label_frame_edf.grid(row=2, column=2, padx=2, pady=5)
#label_frame_edf.pack()

label_NF = tk.Label(label_frame_edf, text="Noise Figure(dB)")
label_NF.pack()
entry_NF = tk.Entry(label_frame_edf)
entry_NF.pack()

# Set the default value
default_valueRe = "4.5"
entry_NF.insert(tk.END, default_valueRe)



label_Gain = tk.Label(label_frame_edf, text="EDFA Gain(dB)>0")
label_Gain.pack()
entry_Gain = tk.Entry(label_frame_edf)
entry_Gain.pack()





# Set the dark_background style for the plot
plt.style.use("dark_background")

# Create subplots for each section of the plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6, 7))

# Create a matplotlib canvas to display the plot inside the plot frame
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Add the navigation toolbar
toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)






# Create the get parameters button inside the input frame
params_button = tk.Button(text_frame, text="Get Parameters", command=get_parameters)
params_button.grid(row=3, column=2, padx=2, pady=5)
#params_button.pack()

# Create the generate signals button inside the input frame
gen_signals_button = tk.Button(text_frame, text="Generate Signals", command=generate_signals)
gen_signals_button.grid(row=4, column=2, padx=2, pady=5)
#gen_signals_button.pack()

# Create the plot button inside the input frame
plot_button = tk.Button(text_frame, text="Start simulation", command=plot_eyediagram)
plot_button.grid(row=5, column=2, padx=2, pady=5)
#plot_button.pack()









output_text.insert(tk.END, "Eyediagram plotted successfully.\n")
def update_output_text():
    output_text.configure(state="normal")  # Enable editing
    output_text.delete("1.0", tk.END)  # Clear existing text
    output_text.insert(tk.END, "The extension ratio in dB is: " + str(np.round(ERdB,3)) + "\n")  # Insert new text
    output_text.insert(tk.END, 'Average power of the modulated \n optical\n signal [mW]: %.3f mW\n'%(signal_power(sigTxo)/1e-3) + "\n")
    output_text.insert(tk.END,'Average power of the modulated optical\n signal at TX [dBm]:\n %.3f dBm\n'%(10*np.log10(signal_power(sigTxo)/1e-3))+ "\n")
    output_text.insert(tk.END,'Average power of the optical signal at\n RX [dBm](before photodiode):\n %.3f dBm\n'%(10*np.log10(signal_power(sigCh)/1e-3))+ "\n")
    output_text.insert(tk.END,'Average power of the optical signal at\n RX [dBm](after photodiode):\n %.3f dBm\n'%(10*np.log10(signal_power(I_Rx)/1e-3))+ "\n")
    #output_text.insert(tk.END, 'Bit error rate:\n {:.3e} \n'.format(Symbol_ER) + "\n")
    output_text.insert(tk.END, 'Best acheived BER:\n {:.3e} \n'.format(np.min(err_hist)) + "\n")
    output_text.insert(tk.END, 'Bit error rate of each samples(Sim):\n' + str(err_hist) + "\n") 

   # output_text.configure(state="disabled")  # Disable editing
    if SYMV>4 or SYMV<0:
        
         output_text.tag_config("red", foreground="red")
         output_text.insert(tk.END, "Error!!! Wrong value for SYMV: This value must be between 0 and 4 based on the design of the Mach-Zehnder modulator" + "\n")
        # output_text.insert(tk.END, 'Average power of the modulated optical signal [mW]: %.3f mW'%(signal_power(sigTxo)/1e-3) + "\n")
       




import tkinter as tk
from PIL import ImageTk, Image

def help_window():
    h_window = tk.Toplevel(root)
    h_window.title("Help")
    
    # Create a label to display the text
    text1 = """Welcome to the OPTICAL IMDD SIMULATOR, a powerful tool to evaluate next-generation optical fiber transceivers for high-speed data centers.Assess Bit Error Rate
    (BER), Power Spectral Analysis, plot Eye Diagrams, analyze the impact of Equalizer adoption,and explore BER trends with distance and input power.Optimize your optical interconnections with confidence using our simulator."""

    text_label = tk.Label(h_window, text=text1, anchor='center')
    text_label.pack()
    
  
    
    # Load and convert the image to a supported format (e.g., PNG)
    image = Image.open("E:\Polito Courses\Semester 2\Project\Phase 1\OPT-IMDD_VERSION 1 presentation\ER_vpp.jpg")
    image = image.convert("RGBA")  # Convert to RGBA format

    
    # Save the image as a temporary PNG file
    temp_image_path = "temp_image.png"
    image.save(temp_image_path)
    
 
    
    # Create a Tkinter-compatible photo image from the saved PNG
    # photo = ImageTk.PhotoImage(file=temp_image_path)

    
    # # Create a label to display the image
    # image_label = tk.Label(h_window, image=photo)
    # image_label.pack()
    
    # # Keep a reference to the image to prevent it from being garbage collected
    # image_label.image = photo
    
    text2 = """For detailed instructions on how to use this simulator, please refer to the user manual document."""
    
    text_label = tk.Label(h_window, text=text2, anchor='center')
    text_label.pack()
    
    text1 = """COPYRIGHT

    This simulator has been developed by Peyman Pahlevanzadeh and Leo Alessandra at Politecnico di Torino. All rights reserved.
    
    © June 2023
    """
    text_label = tk.Label(h_window, text=text1, anchor='center')
    text_label.pack()    
   
    
    # Run the GUI main loop
    root.mainloop()
    








 
# Create the menu bar
menu_bar = tk.Menu(root)

# Create the File menu and its options
file_menu = tk.Menu(menu_bar, tearoff=0)
file_menu.add_command(label="Open")
file_menu.add_command(label="Save")
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

# Create the Edit menu and its options
edit_menu = tk.Menu(menu_bar, tearoff=0)
edit_menu.add_command(label="Cut")
edit_menu.add_command(label="Copy")
edit_menu.add_command(label="Paste")

# Create the View menu and its options
view_menu = tk.Menu(menu_bar, tearoff=0)
view_menu.add_command(label="Zoom In")
view_menu.add_command(label="Zoom Out")

# Create the Help menu and its options
help_menu = tk.Menu(menu_bar, tearoff=0)
help_menu.add_command(label="About")




# Create the Help menu and its options
BER_menu = tk.Menu(menu_bar, tearoff=0)
BER_menu.add_command(label="BER Analysis",command=open_analysis_window)



# Create the Help menu and its options
Perdist_menue = tk.Menu(menu_bar, tearoff=0)
Perdist_menue.add_command(label="Vector",command=open_analysis_vec)



# Create the Help menu and its options
EQ_menue = tk.Menu(menu_bar, tearoff=0)
EQ_menue.add_command(label="Equalizer",command=open_analysis_eq)






# Create the Help menu and its options
help_menu = tk.Menu(menu_bar, tearoff=0)
help_menu.add_command(label="Help",command=help_window)



# Add the menus to the menu bar
#menu_bar.add_cascade(label="File", menu=file_menu)
#menu_bar.add_cascade(label="Edit", menu=edit_menu)
#menu_bar.add_cascade(label="View", menu=view_menu)
menu_bar.add_cascade(label="Help", menu=help_menu)
menu_bar.add_cascade(label="BER", menu=BER_menu)
menu_bar.add_cascade(label="Equalizer", menu=EQ_menue)

menu_bar.add_cascade(label="VEC", menu=Perdist_menue)




# Configure the menu bar as the menu for the root window
root.config(menu=menu_bar)








# Run the GUI main loop
root.mainloop()
