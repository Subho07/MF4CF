# Model free 4-component scattering power decomposition for full polarimetric SAR data
[![DOI](https://img.shields.io/badge/DOI-https%3A%2F%2Fdoi.org%2F10.36227%2Ftechrxiv.13298033.v1-blue)](https://doi.org/10.36227/techrxiv.13298033.v1)

## General information
Target decomposition methods of polarimetric Synthetic Aperture Radar (PolSAR) data explain scattering information from a target. In MF4CF, the non-conventional 3D Barakat degree of polarization is used to obtain the scattered electromagnetic wave's polarization state. 
The degree of polarization is used to obtain the even-bounce, odd-bounce, and diffused scattering power components. Along with this, a measure of target scattering asymmetry is calculated, which is then suitably utilized to obtain the helicity power. 
All the power components are roll-invariant, non-negative and unambiguous.

## Function description

<p align="center">
  <img src="flow_chart_MF4CF.png" alt="Opening the plugin"/>
</p>

## Up and running
This is a `MATLAB` based code. Users will have to select `T3` matrix only in `PolSARpro` format. The default window size is set to 7, which is defined under variable `wsi`. Please change this window size as per your requirement. 

## References
Dey, Subhadip; Bhattacharya, Avik; Frery, Alejandro C.; López-Martínez, Carlos (2020): A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.13298033.v1 
