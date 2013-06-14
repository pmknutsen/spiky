##Analysis
Analysis functions are grouped as functions that process either continuous (analog) or discrete
(digital) signals. These are placed in the `/continuous` or `/discrete` folders, respectively.

You may add your own analysis functions by copying these into one of these two directories.
These functions can then be accessed via the menu from within Spiky.

Analysis functions should adhere to the following formatting.


##Output arguments
Analysis functions can return 0 or 1 output arguments.

Examples:
`S = MYAnalysis(FV)`
where FV is the Spiky data structure and S is a structure that minimally contains the following
fields.

`PRE`
A [1xT] or [NxT] matrix containing signal data.

`PRE_KHz`
A [1x1] scalar containing a sampling rate in kHz

`PRE_TimeBegin`
A [1x1] scalar containing the start time of the signal in seconds.

`PRE_TimeEnd`
A [1x1] scalar containing the end time of the signal in seconds.

PRE is a prefix that can be any string.

Spiky will validate the fields in S and, if these are valid, insert the new signal into its
data structure and display it in the GUI.

Optional fields:
`PRE_Unit`
String describing the signal unit (e.g. Hz).

`PRE_Scale`
A vector describing the scale of the signal when signal is a 2D matrix.


##Running Spiky sub-routines
You can access any of the Spiky sub-routines from within an analysis function by loading the
global variable `Spiky` into the local workspace:

`global Spiky`

You can then run any sub-routine with the syntax:

`Spiky.SubRoutine()`


##Accessing the Spiky data structure
Load the current Spiky data structure with the command:

`[FV,hWin] = Spiky.GetStruct();`
