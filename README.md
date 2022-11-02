# Pattern_separation_toolbox
A set of functions allowing for easy access of pattern separation in spike train ensembles. In particular it includes the new information-theoretic measures described in 'Robust and consistent measures of pattern separation based on information theory and demonstrated in the dentate gyrus' (Bird, Cuntz, &amp; Jedlicka, 2022).

The toolbox is designed to be simple to use and is centred on the function 'analyse_pattern_separation'. Calling this function on an input spike train pattern and an output spike train pattern will compute common classical measures of pattern separation, such as increased cosine distances or reduced correlations, as well as the new information theoretic measures introduced in the above paper, such as sparsity-weighted mutual information and relative redundancy reduction. 

## Installation

The pattern separation toolbox is available in Matlab (currently R2022a). To install the toolbox, simply unzip the folder in an appropriate directory. Running the script 'start_pattern_separation' will add the necessary subfolders to the path and raise warnings if any dependencies cannot be found.

## Dependencies

The pattern separation toolbox relies on the following Matlab toolboxes:

  **a)** **Necessary dependencies**. These toolboxes are necessary to compute the parameter-free information theoretic measures of pattern separation and many functions in the toolbox will not work without them.
 
 i) *Optimization Toolbox*
 
  ii) *Global Optimization Toolbox*

 **b)** **Desirable dependencies**. These toolboxes are necessary to access the full functionality of the pattern separation toolbox, but it is possible to avoid their usage.
  
  i) *Parallel Computing Toolbox*. This toolbox is necessary to carry out information theoretic calculations in parallel and will greatly speed up some computations. If this toolbox is not available, the name-value argument 'parallel' in the function 'analyse_pattern_separation' can be set to 'false' to allow for serial computations instead.
 
 ii) *Signal Processing Toolbox*. This toolbox is necessary to convert raw voltage traces to spike times. If this is not needed, it is possible to ignore this toolbox.
  
  iii) *Curve Fitting Toolbox*. This toolbox is necessary to use the modified Kozachenko-Leonenko estimators for mutual information and redundancy reduction. If this is not needed, it is possible to ignore this toolbox.


## The tutorial

The tutorial is included in the folder 'tutorials' as a Matlab live script (.mlx) divided into sections. The tutorial is designed to demonstrate the key functionality of the pattern separation toolbox.
