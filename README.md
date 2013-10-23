##Spiky
Spiky is a Matlab GUI wrapper for the Chronux spike sorting library (http://www.chronux.org),
with a number of extra features for improved sorting, visualization and analysis of spike
trains and continuous data.

###Requirements
Matlab 6.5 or higher on Linux, Windows or Mac.

##Features
- GUI with zoom, scroll and pan abilities
- 1-click automatic spike sorting*
- Visualization and validation tools
- Basic channel operators
- Bandpass filtering
- Custom math operators
- Analysis tools for spike trains and continuous data (see Plug-ins below)
- Merging of multiple data files
- Batch processing
- Command line interface
- Extendible through scripting

###Install
Synchronize the master branch with `git`, or download the latest '.zip' archive (see ZIP
button above) and unpack to a local location. Add all sub-folders to the Matlab path. Type
`spiky` to start.

###Documentation
https://github.com/pmknutsen/spiky/wiki

###Issue/Bug Tracker
https://github.com/pmknutsen/spiky/issues

##Data Formats
Import:
- Matlab Data Acquisition Toolbox (.daq)

Export:
- Matlab (.mat)

##Analysis Plug-Ins
Continuous:
- Cross correlation
- Event triggered average
- Histogram
- Spectral coherence*
- Power spectral density*
- Spectrogram*

Discrete:
- Cross correlations
- Peristimulus time histogram (PSTH)
- Spike triggered average (STA)
- Temporal distribution

* via Chronux

