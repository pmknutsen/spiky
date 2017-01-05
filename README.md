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
- Customizable extensions

Screenshot 1

![Screenshot:](https://github.com/pmknutsen/spiky/blob/master/themes/spiky_theme.png "Screenshot - Spiky theme")

###Requirements
Matlab 6.5 or higher (preferred) on Linux, Windows or Mac.

###Install
Synchronize the master branch with `git`, or download the latest '.zip' archive (see ZIP
button above) and unpack to a local location. Add all sub-folders to the Matlab path. Type
`spiky` to start.

###API
Spiky is accessible from the command line or your own functions and extensions via a global
structure containing function handles to all Spiky routines.

###Documentation
https://github.com/pmknutsen/spiky/wiki

###Issue/Bug Tracker
https://github.com/pmknutsen/spiky/issues

##Data Formats
Import:
- Matlab Data Acquisition Toolbox (.daq)
- Open ePhys (.openephys)
- Axona (.set)
- Alpha Omega (.mat)
- Hierarchical Data Format (.hd5; requires R2012b or higher)
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

##Extensions
Additional code can be added to extend Spiky functionality. Extensions can be placed in
the `/extensions` folder, or their filenames prefixed with `spiky_` and be located
anywhere in the Matlab search path. See the file `/extensions/Sample_Extension.m` for
details.




