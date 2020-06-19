This dataset is associated with the paper by Connolly, Mackenzie and Gorash wherein an 
 implementation of hyperelasticity in principal stretches is developed.


This dataset contains Fortran (.f90) programs and Fortran (.for) UMAT subroutines for Abaqus, for
 numerical computation of stress and elasticity tensors for hyperelastic constitutitive models
 defined in terms of Cauchy-Green invariants and equivalently in principal stretches.


Also provided are the results of a study to find an approximate optimal tolerance value for the
 proposed algorithm to avoid numerical errors when principal stretches are numerically similar or
 equal.


******************************* COMPATIBILITY ***************************************************
These codes been tested using Abaqus/Standard v2016, Microsoft Visual Studio 2010 and Intel 
 Parallel Studio XE 2013. No guarantee is made for their functionality, in these or any other
 software combinations. However, any issues or bugs may be reported by email to:
 stephen.connolly@strath.ac.uk. The author also appreciates any generic feedback on programming
 methods, code improvements and collaboration opportunities.
The provided excel results files were created using Microsoft Excel 2013 (15.0.5093.1000) 
 MSO (15.0.5085.1000) 64-bit
*************************************************************************************************


The dataset is organised by folders for the 4 constitutive models used in the numerical 
 studies in the associated paper. These folders contain programs to compute and compare the
 error of the developed principal stretch implementation to a well-established Cauchy-Green 
 invariant implementation. These are provided in the reference and spatial configurations.
 Additionally, two UMAT subroutines are provided for each of the constitutive models, for 
 both principal strethc and Cauchy-Green invariant implementations.


The dataset also contains templates for principal stretch and Cauchy-Green invariant
 implementations for UMAT user-subroutines and Fortran programs. These codes contain
 annotations where user input is required.
- Use of "Ctrl+F" of "user input" is recommended


In the programs and subroutines, a tolerance value (tol1=1d-6) is found to be approximately
 optimal. However, the tolerance value may be modified to compare: no L'Hopital's rule; 
 L'Hopital's rule and a other tolerance values.


These codes have been developed for academic use. The author therefore appreciates academic
 credit if they are used directly or modified in the production of publications, communications 
 or other intellectual property. References should be made to the associated journal article and
 this dataset.


Author information (correct as of 02-04-2018):
Stephen Connolly, MEng
PhD Student at the University of Strathclyde, Glasgow
Graduation year: 2020 (expected)
email: stephen.connolly@strath.ac.uk