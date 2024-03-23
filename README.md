# Project Overview

The objective of this project is to develop next-generation optical fiber transceivers using PAM-4 modulation to achieve a desired 400 Gbit/s data rate in data centers. The proposed product is a software-based digital twin of the transmission system, serving as a numerical simulator in the time-domain to replicate and predict the performance of the real physical system.

## Components

### Mach-Zehnder Modulator (MZM)

Implemented to enable a communication system based on transmitting and receiving bits. The MZM modulator facilitates manipulation of the extinction ratio, which directly affects its input and output characteristics.

### Real-world Elements

Incorporated elements such as optical fiber with characteristics like attenuation and chromatic dispersion. At the receiver side, a photodiode accounts for shot noise and TIA (Transimpedance Amplifier) noise, along with an EDFA (Erbium-Doped Fiber Amplifier) with its own noise characteristics.

## Workflow

The project initiates with transmitting bits at the transmitter side, generated using the `np.random.randint()` function. Subsequently, the `modulateGray` function is utilized to modulate the transmitted bits using a Gray coding scheme. This function requires supporting functions such as `GrayMapping`, `bitarray2dec`, and `GrayCode`.

- **GrayCode:** Generates a Gray code, a binary numeral system where successive values differ by only one bit. It returns a list of binary strings representing the Gray code.

- **bitarray2dec:** Converts an array of bits (0 and 1) to a decimal integer.

- **GrayMapping:** Implements Gray Mapping for digital modulations. It calculates the required number of bits for the constellation symbols based on the modulation order, generates the Gray code for those bits, and assigns constellation symbols based on the code and the specified PAM constellation.

## Signal Processing

The signal processing pipeline involves:

- **Upsampling Filter (`upsample`):** Performs upsampling on the input signal array by a factor of `n`, which equals samples per symbol. This operation inserts `n-1` zeros between consecutive samples of the input signal.

- **Pulse Shaping (`pulseShape`):** Creates a pulse shaping function for a typical Raised Cosine (RC) pulse. Users can specify different types of filters (e.g., 'rect', 'nrz', 'rrc', 'rc', and 'srrc') by calling respective functions (`rrcosfilter()`, `rcosfilter()`, or `srrcosfilter()`), along with parameters such as samples per symbol, number of filter coefficients, roll-off factor, and symbol period in seconds.

This project aims to provide a comprehensive toolset for simulating and analyzing next-generation optical fiber transceivers, catering to the demands of high-speed data transmission in modern data center environments.



## Bit Error Rate (BER) Calculation

To analyze the performance, the Bit Error Rate (BER) is calculated by comparing the received bits (`bitsRx`) with the transmitted bits (`bitsTx`). The code computes the BER by performing a bitwise XOR (exclusive OR) operation between the received and transmitted bits. 

A `discard` variable is used to exclude a certain number of bits at the beginning and end of the sequences, which may be affected by initialization or synchronization issues. The resulting vector indicates the positions where the transmitted bits differ from the received bits. 

The mean of the resulting error vector (`err`) is computed using `np.mean(err)`, representing the ratio of bit errors to the total number of compared bits.

For a more accurate BER calculation, an optimal time sampling instance is determined through a search process.

## Erbium-Doped Fiber Amplifier (EDFA)

The implemented EDFA function calculates the noise spectral power (`nsp`) using the formula `(Glin·NFlin−1)/(2·(Glin−1))`, where `Glin` represents the linear gain and `NFlin` denotes the linear noise figure.

The Amplified Spontaneous Emission (ASE) noise power (`Nase`) is determined using the formula `(Glin−1)·nsp·h·Fc`, where `h` represents the Planck constant.

The total noise power (`pnoise`) is computed by multiplying `Nase` with the sampling frequency (`Fs`). 

At the end, a random noise sample with zero mean and a standard deviation determined by `ppnoise` is generated. The function returns the amplified and noisy optical signal `Ei ·√(Glin + noise)`.

