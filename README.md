# AAE6102-Assignment-1
A technical report for GNSS Software-Defined Receiver (SDR) signal processing. An open-source SDR program is used in this project [1], which can be found at https://www.ngs.noaa.gov/gps-toolbox/GPS_VT_SDR.htm.

## Usage

To execute the computation program, run the `SDR_main.m` script.

### Configuration for Different Scenarios

The script includes two configurations for different input data scenarios: **Open Sky** and **Urban**. These configurations are controlled by the following lines in the code:

- **Open Sky Configuration**:
  ```matlab
  %[file, signal, acq, track, solu, dyna, cmn] = initParametersOpenSky();

- **Urban Configuration**:
  ```matlab
  %[file, signal, acq, track, solu, dyna, cmn] = initParametersUrban();

By commenting out one of them, the other one can be executed.

### Initial input parameters

The initial input parameters, including file names and paths, can be modified in the `initParametersOpenSky.m` and 'initParametersUrban.m' scripts.






## Example of results
### Acquisition

The acquisition was performed using a parallel code phase search algorithm, which efficiently searches for satellite signals in both the frequency and code delay domains. The algorithm estimates the Signal-to-Noise Ratio (SNR), Doppler frequency, and code delay for each detected satellite. The initial acquisition results for the detected satellites are summarized in the table below. The table includes the PRN (Pseudo-Random Noise) number, SNR (signal-to-noise ratios), doppler frequency, code delay, and fine frequency for each satellite. Notably, PRN 16 exhibited the highest SNR (32.42 dB), indicating a strong signal, while PRN 3 and PRN 4 showed moderate SNRs of 20.62 dB and 21.15 dB, respectively. The Doppler frequencies ranged from -3,000 Hz to 3,500 Hz, reflecting the relative motion between the satellites and the receiver.

| PRN | SNR (dB) | Doppler (Hz) | Code delay (chips) | Fine frequency (Hz)|
|-----|----------|--------------|--------------------|--------------------|
|3    |20.6226   |1,000         |2,865               |4,580,985           |
|4    |21.1472   |-3,000        |12,134              |4,576,910           |
|16   |32.4237   |0             |26,007              |4,579,755           |
|22   |26.3315   |1,500         |2,899               |4,581,565           |
|26   |31.5943   |2,000         |247                 |4,581,915           |
|27   |28.2166   |-3,000        |49,185              |4,576,770           |
|31   |29.0355   |1,000         |39,257              |4,581,060           |
|32   |26.1496   |3,500         |20,787              |4,583,350           |

### Tracking

### Navigation Data Decoding

### Position and Velocity Estimation

### Kalman Filter-Based Positioning


## References
[1] B. Xu and L.-T. Hsu, ‘Open-source MATLAB code for GPS vector tracking on a software-defined receiver’, GPS Solut, vol. 23, no. 2, p. 46, Apr. 2019, doi: 10.1007/s10291-019-0839-x.
