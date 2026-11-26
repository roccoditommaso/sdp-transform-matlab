# SDP Transform for MATLAB

The Stockwell–Ditommaso–Ponzo (SDP) transform evaluates the Stockwell transform on a rational frequency grid. With the sampling frequency, grid parameter (conductor), and analysis band held fixed, the frequency nodes do not change with record length. This makes it easier to compare time–frequency maps across signals of different durations without interpolating onto a shared frequency grid.

## Quick start

Download the release and add its folder to the MATLAB path, or open that folder in MATLAB. Run the included worked example:

```matlab
sdpt_example
```

To obtain the transform and its metadata directly:

```matlab
S = load('sdpt_example_data.mat');
[M, f, t, info] = sdpt(S.x, 0.8, 12, S.fs);
```

The last call uses the illustrative frequency band 0.8–12 Hz. The outputs include the time–frequency map, frequency and time axes, and metadata such as the extracted ridge and frequency-node labels. The example data are synthetic. The FFT-based convolution is circular: use care when interpreting transients near the ends of the record, especially at low frequencies.

## Files

- `sdpt.m`: MATLAB implementation.
- `sdpt_example.m` and `sdpt_example_data.mat`: worked example and synthetic data.
The illustrated guide is available with the original [File Exchange package](https://it.mathworks.com/matlabcentral/fileexchange/184677-sdp-transform-stockwell-time-frequency-analysis-on-a-record).

The software requires MATLAB. See the [MATLAB Central File Exchange listing](https://it.mathworks.com/matlabcentral/fileexchange/184677-sdp-transform-stockwell-time-frequency-analysis-on-a-record) for the original distribution.

## Citation

Rocco Ditommaso and Felice Carlo Ponzo (2026). *Time–frequency representation of nonlinear structural response through rationally sampled Stockwell analysis*. Mechanical Systems and Signal Processing, 260, 114920. [doi:10.1016/j.ymssp.2026.114920](https://doi.org/10.1016/j.ymssp.2026.114920).

See also R. Ditommaso, M. Mucciarelli, and F. C. Ponzo (2012). *Analysis of non-stationary structural systems by using a band-variable filter*. Bulletin of Earthquake Engineering, 10, 895–911. [doi:10.1007/s10518-012-9338-y](https://doi.org/10.1007/s10518-012-9338-y).

## Authors and funding

Rocco Ditommaso and Felice Carlo Ponzo, Department of Engineering, University of Basilicata, Potenza, Italy. This work was supported by the Italian Civil Protection Department through the DPC–ReLUIS 2024–2026 project, WP6 (Monitoring and satellite data).

## License

The software and bundled example are distributed under the [BSD 3-Clause License](LICENSE.txt). The non-endorsement clause requires specific prior written permission to use the authors’ names to endorse or promote derived products. This does not preclude scholarly citation or the attribution required by the license. The paper is linked above rather than bundled here.
