# features

**IDL code for identifying and tracking features of interest in video microscopy data**

IDL is the Interactive Data Language, and is a product of
[Exelis Visual Information Solutions](http://www.exelisvis.com)

The code in **features** is licensed under the GPLv3.

## What it does

1. **Crocker-Grier Particle Tracking**
    * **featuretool**: GUI-based interactive implementation of the Crocker-Grier algorithm
    * **bpass**: Spatial bandpass filter for minimizing noise and background variations
in images before attempting feature identification.
    * **feature**: Find and characterize circular features within an image.
NOTE: the updated local-maximum detection algorithm may require revisions
to existing particle-tracking scripts.
    * **fastfeature**: Thresholding algorithm for identifying features in an image.
    * **track**: Link features in consecutive images into time-resolved trajectories.
    * **sg_lmax**: Locate local maxima in an image using Savitzky-Golay filters.
2. **Feature Tracking Based on the Circletransform Algorithm**
    * **circletransform**: Transform ring-like features in an image into simply peaked
features suitable for tracking with CG algorithms.
    * **ctfeature**: Uses circletransform and fastfeature to identify and characterize
ring-like features in an image.

### Related Projects

1. [Particle Tracking in IDL](http://www.physics.emory.edu/faculty/weeks//idl/index.html) by Eric Weeks and John Crocker.
2. [Trackpy](http://soft-matter.github.io/trackpy) Particle tracking in Python.

### References
1. J. C. Crocker and D. G. Grier, "Methods of digital video microscopy for colloidal studies,"
_Journal of Colloid and Interface Science_ **179**, 298-310 (1996).

2. F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao, L. Dixon and D. G. Grier, 
"Flow visualization and flow cytometry with holographic video microscopy,"
_Optics Express_ **17**, 13071-13079 (2009).
