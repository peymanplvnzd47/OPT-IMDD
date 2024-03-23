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
