# Signal-Processing-for-Surface-EMG-Data

Contributors
```
Dhruv Gupta
Elena Peng
```
Goals
```
To determine muscle onset time and muscle activation amplitude in surface EMG signals 
obtained from quadriceps muscles when performing isometric straight leg raise tasks.
```
Input of the Code
```
(1) Input the data files needed for onset time and activation amplitude analysis. 
    - A testing trial file
    - Maximum voluntary contraction (MVC) files
    - A resting trial during MVC file

(2) Input the parameters for onset time calculation.
    - RMS window (ms)
    - Smoothing window for RMS (ms)
    - Threshold of standard deviation to define onset (SD)
```
Output of the Code
```
(1) Output of onset time (ms) of the 4 muscles (Rectus femoris, Vastus medialis oblique, 
Vastus lateralis, and Vastus medialis) and onset time (ms) of force.

(2) Output of activation amplitudes normalized to the MVC (%) for the 4 muscles.

(3) Output plots to visualize the onset time marked with a red circle for the 4 muscle 
and the force channel. 
```
