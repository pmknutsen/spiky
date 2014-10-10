##Spiky
Spiky is a Matlab application for the analysis of electrophysiological and behavioral data.

Spiky includes the Chronux library for spike sorting and spectral analysis (http://www.chronux.org),
along with extra features for convenient sorting, data visualization and analysis of discrete
and continuous data.

##Features
- 1-click automatic spike sorting
- Visualization and validation tools
- Basic channel operators
- Per-channel bandpass filtering
- Custom math operators
- Zoom, scroll and pan data
- Analysis of discrete and continuous data
- Merging across data files
- Batch processing
- Matlab API
- Extendible through scripting

###Requirements
Matlab 6.5 (recommended) or higher on Linux, Windows or Mac.

###Install
Synchronize the master branch with `git`, or download the latest '.zip' archive (see ZIP
button above) and unpack to a local location. Add all sub-folders to the Matlab path. Type
`spiky` to start.

###API
Spiky is accessible from the command line or your own functions and scripts via a global
structure containing function handles to all Spiky routines.

###Documentation
https://github.com/pmknutsen/spiky/wiki

###Issue/Bug Tracker
https://github.com/pmknutsen/spiky/issues

##Data Formats
Import:
- Matlab Data Acquisition Toolbox (.daq)
- Hierarchical Data Format (.hd5)
- Matlab vectors (.mat)

Export:
- Matlab (.mat)

##Analysis Routines
Continuous:
- Cross correlation
- Event triggered average
- Histogram
- Spectral coherence
- Power spectral density
- Spectrogram

Discrete:
- Cross correlations
- Peristimulus time histogram (PSTH)
- Spike triggered average (STA)
- Temporal distribution
