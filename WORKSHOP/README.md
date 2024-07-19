### WCCM2024 workshop
Source code, model database, and FEA input files for the  
WCCM2024 workshop on _Automated Model Discovery_,  
Vancouver, July 21st 2024.

Files:
- **BLOCKcae.cae** - abaqus model database used to generate finite element model of 10x10x10 hexahedral element block undergoing twist and pull loading
- **UANIuniversal.f** - universal material subroutine
- **UNIVERSAL_PARAM_TYPES.inc** - universal material subroutine parameter table type definition file
- **BLOCK10-isotropic.inp** - FEA input file twist pull block simulation, incorporating simple build-in polynomial isotropic hyperelastic strain energy function
- **BLOCK10-brain.inp** - FEA input file twist pull block simulation, leveraging universal material subroutine to model constitutive neural network-discovered gray matter material model
- **BLOCK10-arterial.inp** - FEA input file twist pull block simulation, leveraging universal material subroutine to model constitutive neural network-discovered arterial adventitia tissue material model
- **BLOCK10-myocard.inp** - FEA input file twist pull block simulation, leveraging universal material subroutine to model constitutive neural network-discovered myocardial tissue material model

  

When using, please cite  
**"A universal material model subroutine for soft matter systems",  
M. Peirlinck, J.A. Hurtado, M.K. Rausch, A. Buganza Tepole, E. Kuhl,  
Engineering with Computers, 2024**,
[doi](https://doi.org/10.48550/arXiv.2404.13144).
