# GNSS Software-Defined Receiver (SDR) signal processing
A technical report for GNSS Software-Defined Receiver (SDR) signal processing. An open-source SDR program is used in this project [1], which can be found at https://www.ngs.noaa.gov/gps-toolbox/GPS_VT_SDR.htm.

## Usage

To execute the computation program, run the `SDR_main.m` script.

### Configuration for Different Scenarios

The script includes two configurations for different input data scenarios: **Open Sky** and **Urban**. The following lines in the code control these configurations:

- **Open Sky Configuration**:
  ```matlab
  %[file, signal, acq, track, solu, dyna, cmn] = initParametersOpenSky();

- **Urban Configuration**:
  ```matlab
  %[file, signal, acq, track, solu, dyna, cmn] = initParametersUrban();

By commenting out one of them, the other one can be executed.

### Initial input parameters

The initial input parameters, including file names and paths, can be modified in the `initParametersOpenSky.m` and 'initParametersUrban.m' scripts.

### Post-Processing

To generate visualizations in the  `Figures` folder, simply run the `PostProcessing.m` script in the root directory.







## Example of results
### Acquisition

The acquisition was performed using a parallel code phase search algorithm, which efficiently searches for satellite signals in both the frequency and code delay domains. The initial acquisition results for the detected satellites are summarized in the tables below. The tables include the PRN (Pseudo-Random Noise) number, SNR (signal-to-noise ratios), doppler frequency, code delay, and fine frequency for each satellite. These results can be found in the output file named `Acquired_***_***.mat`.

| PRN (open sky)  | SNR (dB) | Doppler (Hz) | Code delay (chips) | Fine frequency (Hz)|
|-----------------|----------|--------------|--------------------|--------------------|
|3                |20.6226   |1,000         |2,865               |4,580,985           |
|4                |21.1472   |-3,000        |12,134              |4,576,910           |
|16               |32.4237   |0             |26,007              |4,579,755           |
|22               |26.3315   |1,500         |2,899               |4,581,565           |
|26               |31.5943   |2,000         |247                 |4,581,915           |
|27               |28.2166   |-3,000        |49,185              |4,576,770           |
|31               |29.0355   |1,000         |39,257              |4,581,060           |
|32               |26.1496   |3,500         |20,787              |4,583,350           |

| PRN (urban)  | SNR (dB) | Doppler (Hz) | Code delay (chips) | Fine frequency (Hz)|
|--------------|----------|--------------|--------------------|--------------------|
|1             |42.2098   |1,000         |22,672              |1,205               |
|3             |29.7096   |4,500         |829                 |4,282               |
|7             |20.4342   |500           |10,810              |370                 |
|11            |26.1329   |500           |24,845              |410                 |
|18            |21.7993   |-500          |15,421              |-325                |


### Tracking

In the source code, the in-phase component and quadrature component of the early correlator prompt correlator and late correlator are expressed by $E_i$, $E_q$, $P_i$, $P_q$, $L_i$, $L_q$, respectively. Then, the three correlators are calculated as:

$$E = \sqrt{E_i^2 + E_q^2}$$

$$P = \sqrt{P_i^2 + P_q^2}$$

$$L = \sqrt{L_i^2 + L_q^2}$$

They correspond to time delays of -0.5, 0, and 0.5 chips, respectively. For both open sky and urban cases, the correlation plots of satellite 3 at four time points are shown below.

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/correlation_analysis_plotOpensky.png)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/correlation_analysis_plotUrban.png)

### Navigation Data Decoding

Navigation data decoding aims to obtain pseudorange, GPS time, ephemeris, almanac, and klobuchar information, which are included in the output file named `eph_***_***.mat`. Again, use satellite 3 for both cases as an example; some decoding results are listed below.

| Parameters            | Open sky    | Urban    |
|-----------------------|-------------|----------|
|Time of week (s)       |390,108      |449,358   |
|Week number            |2,179        |2,056     |
|Time of clock(s)       |296,000      |453,600   |
|Time of ephemeris (s)  |396,000      |453,600   |
|...|







### Position and Velocity Estimation (WLS)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Open%20sky%20latitude.png)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Opensky%20longitude.png)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Opensky%20velocity.png)


![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Urban%20latitude.png)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Urban%20longitude.png)

![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Urban%20velocity.png)


![image](https://github.com/ZimoJupiter/AAE6102-Assignment-1/blob/main/Figures/Map.png)



### Kalman Filter-Based Positioning





## References
[1] B. Xu and L.-T. Hsu, ‘Open-source MATLAB code for GPS vector tracking on a software-defined receiver’, GPS Solut, vol. 23, no. 2, p. 46, Apr. 2019, doi: 10.1007/s10291-019-0839-x.
