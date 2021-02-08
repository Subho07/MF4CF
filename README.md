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

## Algorithm

<img src="https://latex.codecogs.com/gif.latex?\\&space;\textbf{{\color{Blue}&space;Input:}}&space;~\mathbf{T}&space;\\&space;\textbf{{\color{Blue}&space;Result:}}&space;~P_{c},&space;P_{v},&space;P_{s},&space;P_{d}\\&space;\textbf{{\color{Blue}&space;Kennaugh&space;matrix:}}&space;~\mathbf{K}\\&space;\textbf{{\color{Blue}&space;Scattering&space;Descriptors:}}\\&space;\\&space;m_{\text{FP}}&space;\leftarrow&space;\sqrt{1&space;-&space;\dfrac{27\left|\mathbf{T}\right|}{\text{tr}^{3}(\mathbf{T})}}&space;\;\\&space;\\&space;\theta_{\text{FP}}&space;\leftarrow&space;\tan^{-1}\dfrac{4m_{\text{FP}}K_{11}K_{44}}&space;{K_{44}^2&space;-&space;(1&space;&plus;&space;4m_{\text{FP}}^2)K_{11}^2}\;\\&space;\\&space;\tau_{\text{FP}}&space;\leftarrow&space;\tan^{-1}\dfrac{|K_{14}|}{K_{11}}\;\\&space;\textbf{{\color{Blue}&space;Scattering&space;Powers:}}\\&space;\\&space;P_{c}&space;\leftarrow&space;2&space;m_{\text{FP}}&space;K_{11}&space;\sin&space;2\tau_{\text{FP}}\;\\&space;P_{v}&space;\leftarrow&space;2(1&space;-&space;m_{\text{FP}})K_{11}\;\\&space;\\&space;\textbf{{\color{Blue}&space;Residue:}}&space;P_{r}\leftarrow&space;2&space;m_{\text{FP}}&space;K_{11}&space;(1&space;-&space;\sin&space;2\tau_{\text{FP}})\;\\&space;\\&space;P_{s}&space;\leftarrow&space;\dfrac{P_{r}}{2}(1&space;&plus;&space;\sin&space;2\theta_{\text{FP}})\;\\&space;\\&space;P_{d}&space;\leftarrow&space;\dfrac{P_{r}}{2}(1&space;-&space;\sin&space;2\theta_{\text{FP}})\;\\" title="\\ \textbf{{\color{Blue} Input:}} ~\mathbf{T} \\ \textbf{{\color{Blue} Result:}} ~P_{c}, P_{v}, P_{s}, P_{d}\\ \textbf{{\color{Blue} Kennaugh matrix:}} ~\mathbf{K}\\ \textbf{{\color{Blue} Scattering Descriptors:}}\\ \\ m_{\text{FP}} \leftarrow \sqrt{1 - \dfrac{27\left|\mathbf{T}\right|}{\text{tr}^{3}(\mathbf{T})}} \;\\ \\ \theta_{\text{FP}} \leftarrow \tan^{-1}\dfrac{4m_{\text{FP}}K_{11}K_{44}} {K_{44}^2 - (1 + 4m_{\text{FP}}^2)K_{11}^2}\;\\ \\ \tau_{\text{FP}} \leftarrow \tan^{-1}\dfrac{|K_{14}|}{K_{11}}\;\\ \textbf{{\color{Blue} Scattering Powers:}}\\ \\ P_{c} \leftarrow 2 m_{\text{FP}} K_{11} \sin 2\tau_{\text{FP}}\;\\ P_{v} \leftarrow 2(1 - m_{\text{FP}})K_{11}\;\\ \\ \textbf{{\color{Blue} Residue:}} P_{r}\leftarrow 2 m_{\text{FP}} K_{11} (1 - \sin 2\tau_{\text{FP}})\;\\ \\ P_{s} \leftarrow \dfrac{P_{r}}{2}(1 + \sin 2\theta_{\text{FP}})\;\\ \\ P_{d} \leftarrow \dfrac{P_{r}}{2}(1 - \sin 2\theta_{\text{FP}})\;\\" />


## Up and running
This is a `MATLAB` based code. Users will have to select `T3` matrix in `PolSARpro` format. The default window size is set to 7, which is defined under variable `wsi`. Please change this window size as per your requirement. 

## References
Dey, Subhadip; Bhattacharya, Avik; Frery, Alejandro C.; López-Martínez, Carlos (2020): A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.13298033.v1 
