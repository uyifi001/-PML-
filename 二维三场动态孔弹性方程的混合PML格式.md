# A hybrid PML formulation for the 2D three-field dynamic poroelastic equations

Hernán Mella \(^{a,*}\) , Esteban Sáez \(^{b,c}\) , Joaquín Mura \(^{d}\)

\(^{a}\) School of Electrical Engineering, Pontificia Universidad Católica de Valparaíso, Valparaíso, Chile  
\(^{b}\) Department of Structural and Geotechnical Engineering, Pontificia Universidad Católica de Chile, Santiago, Chile  
\(^{c}\) Center for Integrated Natural Disaster Management CONICYT/FONDAP/15110017, Santiago, Chile  
\(^{d}\) Department of Mechanical Engineering, Universidad Técnica Federico Santa María, Santiago, Chile  

Received 25 April 2023; received in revised form 14 August 2023; accepted 15 August 2023  
Available online 1 September 2023  

Dataset link: https://github.com/hmella/POROUS-HYBRID-PML  

## Abstract

Simulation of wave propagation in poroelastic half-spaces presents a common challenge in fields like geomechanics and biomechanics, requiring Absorbing Boundary Conditions (ABCs) at the semi-infinite space boundaries. Perfectly Matched Layers (PML) are a popular choice due to their excellent wave absorption properties. However, PML implementation can lead to problems with unknown stresses or strains, time convolutions, or PDE systems with Auxiliary Differential Equations (ADEs), which increases computational complexity and resource consumption.

This article introduces a novel hybrid PML formulation for arbitrary poroelastic domains. Instead of using ADEs, this formulation utilizes time-history variables to reduce the number of unknowns and mathematical operations. The modification of the PDEs to account for the PML is limited to the PML domain only, resulting in smaller matrices while maintaining the governing equations in the interior domain and preserving the temporal structure of the problem. The hybrid approach introduces three scalar variables localized within the PML domain.

The proposed formulation was tested in three numerical experiments in geophysics using realistic parameters for soft sites with free surfaces. The results were compared with numerical solutions from extended domains and simpler ABCs, such as paraxial approximation, demonstrating the accuracy, efficiency, and precision of the proposed method. The article also discusses the applicability of the method to complex media and its extension to the Multiaxial PML formulation.

The codes used for the simulations are available for download from https://github.com/hmella/POROUS-HYBRID-PML.  
© 2023 Elsevier B.V. All rights reserved.

**Keywords:** Perfectly Matched Layers; Poroelastic wave propagation; Absorbing Boundary Condition; Three-field Biot's equations

---

## 1. Introduction

Fluid-saturated porous media are a common occurrence in nature. For instance, soils and rocks are often saturated with water in practical cases, while living tissues are saturated with blood and air. In both cases, if the solid skeleton's displacements and strains are relatively small, linear elasticity provides an accurate representation of the underlying dynamics. Additionally, when loads are applied quickly and inertial forces play a significant role, a proper modeling strategy for wave propagation in poroelastic media is necessary. Biot's poroelastic theory can be employed to describe wave propagation in poroelastic media, for instance in problems of ultrasound imaging of the human body or geophysical applications involving seismic wave propagation. Nevertheless, as noted by Zienkiewicz et al. [1], there are many applications where Biot's theory is not strictly necessary. For instance, in [2], the authors show that if the phenomena are relatively slow speed (which is characterized by medium properties such as permeability and fluid compressibility), compressive waves in the fluid can be neglected, leading to a simplified form in terms only of the pore pressure and solid displacement. The main challenge in these types of problems is properly handling outgoing waves. In the directions where outgoing waves travel, finite energy considerations lead to the so-called "radiation conditions" towards infinity, which are used by Integral Equations or Boundary-Element methods to determine Green kernels and solve the problem rigorously. However, these conditions are often difficult to calculate and are limited to homogeneous and isotropic material properties at infinity [3]. An alternative solution is to use a foam-like subdomain to confine the region of interest, creating virtual windows in space to focus computational efforts on a specific area of the problem. The subdomain must effectively absorb the outgoing waves from the virtual window.

Numerous numerical methods have been proposed as Absorbing Boundary Conditions (ABCs). Local ABCs are often used for dry elastic problems or single-phase media due to their ease of implementation and local character in both time and space [4]. However, fluid-saturated porous media or two-phase media presents a different challenge due to an interaction between the solid skeleton and fluid flow, which depends on the loading rate. According to Biot's theory, high-frequency loading generates two dilatational waves and one shear wave. When the porous media has low permeability and the loading is within the low frequency range, the fluid's relative motion with respect to the soil is negligible, and viscous coupling dampens out the second dilatational wave [5,6]. In this case, the fluid-saturated porous media behaves like a single-phase medium, where only one dilatational and one shear wave propagate.

In recent decades, the Perfectly Matched Layer (PML) has gained popularity as an ABC due to its excellent energy-absorbing properties. PML was first developed by Berenger [7] in the context of electromagnetism. Although the initial development was for Maxwell's equations, its use as an ABC for acoustic [8], elastic [9], and poroelastic [10] domains was later extended. Since then, the technique has been widely used to simulate the propagation of elastic [11-15] and poroelastic waves [16-19] and new and novel formulations have been introduced. For instance, Wang et al. developed an unsplit-field PML formulation by perturbing the elastodynamics equation using convolutional operators [11]. Later, Basu & Chopra introduced an unsplit-field displacement-based time-domain formulation, using time integrals for the strain and stress tensors instead of convolutions [12]. Inspired by Basu & Chopra's approach, Kucukcoban et al. presented a time-domain mixed formulation, incorporating displacement and stress history as unknowns. Their objective was to overcome the estimation of time integrals while offering a simpler scheme [13]. Based on this work, the same authors further improved the approach by proposing a hybrid PML formulation, where the mixed variables are exclusively defined within the PML layer [14]. Similarly, Zhou and colleagues introduced a novel unsplit-field scheme that employed a single vector auxiliary variable confined to the PML layer together with a first-order system of Auxiliary Differential Equations (ADEs) [15].

A similar development has been made in poroelastodynamics starting by the pioneer work made by Zheng et al. who introduced a split-field poroelastic PML formulation [10]. Influenced by the work of Wang et al. in elastodynamics, Song et al. extended their unsplitted formulation introducing a perturbed set of PDEs using convolutions [16]. Martin et al. introduced later an unsplit-field formulation of the convolutional PML method for poroelasticity [17], which is an extension of their previous work in elasticity [20]. More recently, He et al. introduced an unsplit-field second-order poroelastic formulation of the PML method by using an additional system of ADEs. Overall, in either the elastic or poroelastic case, the PML formulations can be divided into split-field and unsplit-field approaches, both of which have drawbacks. Split-field formulations often result in mixed problems where stresses or strains are unknowns, increasing the computational cost of solving the problem [10,21-23]. Unsplit-field formulations typically require the estimation of convolutions or solving ADEs [17,18], which can also be expensive due to the increased number of mathematical operations or the introduction of auxiliary variables. Additionally, little attention has been paid to simulating poroelastic waves in arbitrary domains with realistic subsoil properties.

In this article, we propose a hybrid formulation of the PML method for the second order three-field Biot's equations to address the previously mentioned limitations. The proposed formulation maintain the second-order in time structure of the original equations, making it compatible with most time integration schemes. Furthermore, only three additional scalar variables are introduced and no splitted variables, convolutions, or ADEs are needed, making the proposal easier to implement and less computationally expensive than previous developments [17-19]. The hybrid formulation modifies Biot's equations only in the PML region, further improving the computational cost savings. Our method perform well under challenging conditions, such as free-surface wave propagation, transitions between water and air-filled soft media, and complex geometries, making them suitable for simulating realistic media.

This work is inspired on previous contributions made by Kucukcoban et al. [13,14] and Fathi et al. [24] on the application of the PML method in the elastic wave propagation phenomenon, and can be considered as an extension of these previous work to the poroelastic case. Moreover, the notation adopted in this document follows the notation introduced by these authors.

---

## 2. Poro-elastodynamic equations

The three-field model proposed by Biot [1,2,25] considers a domain \(\Omega \subseteq \mathbb{R}^2\) where the solid displacement \(\mathbf{u}(\mathbf{x},t)\), the displacement of the fluid phase relative to the solid \(\mathbf{w}(\mathbf{x},t)\), and the pore pressure in the fluid \(p(\mathbf{x},t)\) interacts at any position \(\mathbf{x}\) and time \(t\) such that \((\mathbf{x},t)\in \Omega \times T\) with \(T = (0,\infty)\), according to:

\[
\begin{array}{r}
\rho \ddot{\mathbf{u}} +\rho_f\ddot{\mathbf{w}} = \nabla \cdot \boldsymbol{\sigma} \qquad \text{in } \Omega \times T \\
\rho_f\ddot{\mathbf{u}} +\rho_w\ddot{\mathbf{w}} +\frac{\eta}{\kappa}\dot{\mathbf{w}} = -\nabla p \qquad \text{in } \Omega \times T \\
-\dot{p} = \nabla \cdot \{M(\alpha \dot{\mathbf{u}} +\dot{\mathbf{w}})\} \qquad \text{in } \Omega \times T 
\end{array} \quad (1)
\]

where the effective density \(\rho = \rho_{s}(1 - \phi) + \rho_{f}\phi\) depends on the solid \((\rho_{s})\) and fluid densities \((\rho_{f})\), as well as the porosity \(\phi\). The fluid density \(\rho_{w} = \tau \rho_{f} / \phi\) depends on the tortuosity \(\tau\), the dynamic viscosity of the fluid \(\eta\), and the saturated permeability of the porous media \(\kappa\). The parameter \(\alpha\) represents the Biot-Willis coefficient and \(M\) the fluid-solid coupling bulk modulus, defined as

\[
\alpha = 1 - \frac{K_b}{K_s},\quad M = \left(\frac{\phi}{K_f} +\frac{\alpha - \phi}{K_s}\right)^{-1} \quad (2)
\]

where \(K_{b},K_{s},K_{f}\) are the bulk moduli of the dry porous skeleton, solid, and fluid, respectively.

The constitutive law for the poroelastic media in (1) is given by

\[
\begin{array}{l}
\boldsymbol{\sigma}(\mathbf{u},p) = \mathbf{C}\mathbf{e}(\mathbf{u}) - \alpha p\mathbf{I}\\
\mathbf{e}(\mathbf{u}) = \frac{1}{2}\left\{\nabla \mathbf{u} + (\nabla \mathbf{u})^T\right\}
\end{array} \quad (3)
\]

where \(C\) is the fourth-order elastic tensor with components \(C_{ijkl} = \lambda_{b}\delta_{ij}\delta_{kl} + \mu_{b}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\) (with \(\delta_{ij}\) the Kronecker delta), \(\mathbf{e}\) the linear strain tensor, and \(\mathbf{I}\) the identity tensor.

The domain \(\Omega\) is intended to be an unbounded semi-space, with the top boundary consisting of two parts: \(\Gamma = \Gamma_{N}\cup \Gamma_{g}\). The free surface, denoted by \(\Gamma_{N}\), is the region where tractions must vanish [26]. On the other hand, \(\Gamma_{g}\) is the part of the top boundary where an external load \(\mathbf{g}\) is being applied. Therefore, the boundary conditions of the problem are:

\[
\begin{array}{rl}
\boldsymbol{\sigma}(\mathbf{u},p)\cdot \mathbf{n} = \mathbf{g} & \text{in } \Gamma_g\times T\\
\boldsymbol{\sigma}(\mathbf{u},p)\cdot \mathbf{n} = \mathbf{0} & \text{in } \Gamma_N\times T\\
p = 0 & \text{in } \Gamma_N\times T
\end{array} \quad (4)
\]

The variables \(\mathbf{u}\), \(\mathbf{w}\), and \(p\) must vanish as the distance from the source tends to infinity.

---

## 3. Derivation of PML formulas

The region of interest where boundary-reflected waves are not desired (namely Regular Domain) \(\Omega^{\mathrm{RD}}\) is surrounded and truncated by a thin Perfectly Matched Layer \(\Omega^{\mathrm{PML}}\) such that \(\Omega = \Omega^{\mathrm{RD}}\cup \Omega^{\mathrm{PML}}\) (see Fig. 1(b)). Thus, \(\Omega^{\mathrm{PML}}\) is defined as an extension of \(\Omega^{\mathrm{RD}}\) where the outgoing waves are attenuated by a complex coordinate stretching (see Fig. 1).

### 3.1. Complex-coordinate stretching

A complex-coordinate stretching applied to Biot's Eqs. (1) leads to a modified set of equations within \(\Omega\). Thus, the complex-coordinate stretching follows:

\[
r \longmapsto \int_{0}^{r}\epsilon_{r}(r^{\prime},s) \, dr^{\prime} \quad (5)
\]

where \(r\) denotes the spatial coordinate being transformed (namely \(x\) or \(y\) in the two-dimensional case), \(s\) the dual variable of the time in the Laplace domain, and \(\epsilon_{r}\) a complex-coordinate stretching function along the coordinate \(r\). The coordinate transformation given in (5) implies:

\[
\frac{\partial}{\partial r} \longmapsto \frac{1}{\epsilon_r(r)}\frac{\partial}{\partial r} \quad (6)
\]

which is the fundamental relation used to transform the governing equations.

The function \(\epsilon_{r}\) in (7) can adopt diverse forms to achieve different absorption properties. For instance, Kuzuoglu and Mitra [27] introduced the Convolutional Frequency Shifted (CFS) PML method by defining a frequency-dependent stretching function to improve the absorption efficiency at different frequencies and also to improve the time stability. Later, Correia and Jin [28] introduced the higher-order PML with a better absorption rate than CFS-PML. However, both approaches make the real and imaginary parts of \(\epsilon_{r}\) frequency-dependent, leading to convolution terms in the PML formulation. Later, Meza-Fajardo and Papageorgiou developed the Multiaxial PML (M-PML) method [29] and introduced multi-directional attenuation functions with better time stability properties. Recently, François et al. [30] proposed a non-convolutional version of the CFS-PML method in elastodynamics by introducing auxiliary variables. A standard selection for \(\epsilon_{r}\) would be:

\[
\epsilon_{r}(r,s) = \alpha_{r}(r) + \frac{\beta_{r}(r)}{s}, \quad (7)
\]

with \(\alpha_{r}\) and \(\beta_{r}\) denoting real-valued scaling and attenuation functions, respectively.

The real part of \(\epsilon_{r}\) scales the spatial coordinate \(r\) while the imaginary part is responsible for the amplitude decay of the waves entering the PML region (see Fig. 1). To avoid modifying the propagating waves inside \(\Omega^{\mathrm{RD}}\) and ensure the wave attenuation inside \(\Omega^{\mathrm{PML}}\), the following conditions must be fulfilled: \(\alpha_{r}(r)\) and \(\beta_{r}(r)\) are constant inside \(\Omega^{\mathrm{RD}}\) and take values of 1 and 0 respectively, and both functions increases monotonically with \(r\) within \(\Omega^{\mathrm{PML}}\). Several ways of defining \(\alpha_{r}\) and \(\beta_{r}\) have been proposed in the literature in different contexts [7,13,18,31], but in this investigation, polynomial profiles were chosen according to:

\[
\alpha_{r}(r) = \left\{ \begin{array}{ll}
1 & \mathrm{if~}0\leq r\leq r_{0}\\
1 + \alpha_{0}\left\{\frac{(r - r_{0})n_{r}}{L_{\mathrm{PML}}}\right\}^{m} & \mathrm{if~}r_{0}\leq r\leq r_{t}
\end{array} \right. \quad (8a)
\]

\[
\beta_{r}(r) = \left\{ \begin{array}{ll}
0 & \mathrm{if~}0\leq r\leq r_{0}\\
\beta_{0}\left\{\frac{(r - r_{0})n_{r}}{L_{\mathrm{PML}}}\right\}^{m} & \mathrm{if~}r_{0}\leq r\leq r_{t}
\end{array} \right. \quad (8b)
\]

where \(n_{r}\) denotes the \(r\)th component of the outward normal to the interface between \(\Omega^{\mathrm{RD}}\) and \(\Omega^{\mathrm{PML}}\), \(L_{\mathrm{PML}}\) the width of the PML layer, and \(m\) the order of the attenuation profiles, \(r_0\) and \(r_t\) represent the start and end of the absorbing layer (see Fig. 1). The constants \(\alpha_0\) and \(\beta_0\) define the absorption rate in \(\Omega^{\mathrm{PML}}\), and can be chosen as [13]:

\[
\alpha_0 = \frac{(m + 1)b}{2L_{\mathrm{PML}}}\log \left(\frac{1}{|R|}\right),\qquad \beta_0 = \frac{(m + 1)c_p}{2L_{\mathrm{PML}}}\log \left(\frac{1}{|R|}\right) \quad (9)
\]

where \(b\) denotes a characteristic length (e.g., the width of the distributed load applied in the Neumann boundary or the cell size of the finite element mesh), \(c_p\) the propagation velocity of the fastest wave, and \(R\) a reflection coefficient.

### 3.2. PML derivation in the Laplace domain

The complex-coordinate stretching is enforced by introducing the relation (6) into the Laplace-transformed motion equations. Applying the Laplace transform to Biot's equations leads:

\[
\begin{array}{r l}
s^{2}\rho \hat{\mathbf{u}} +s^{2}\rho_{f}\hat{\mathbf{w}} = \nabla \cdot \hat{\boldsymbol{\sigma}} \qquad \text{(linear momentum conservation)}\\
s^{2}\rho_{f}\hat{\mathbf{u}} +s^{2}\rho_{w}\hat{\mathbf{w}} +s\frac{\eta}{\kappa}\hat{\mathbf{w}} = -\nabla \hat{p}\\
-s\hat{p} = \nabla \cdot \left\{M(\alpha s\hat{\mathbf{u}} +s\hat{\mathbf{w}})\right\}\\
\hat{\boldsymbol{\sigma}} = \mathbf{C}\hat{\boldsymbol{\epsilon}} -\alpha \hat{p}\mathbf{I}\\
\hat{\boldsymbol{\epsilon}} = \frac{1}{2}\left\{\nabla \hat{\mathbf{u}} +(\nabla \hat{\mathbf{u}})^{T}\right\}
\end{array} \quad (10)
\]

where \(\hat{f}\) denotes the variable \(f\) in the Laplace domain.

#### 3.2.1. Linear momentum conservation

Replacing (6) into (10a) for the plane strain case gives:

\[
\begin{array}{l}
s^{2}\rho \hat{u}_{x} + s^{2}\rho_{f}\hat{w}_{x} = \frac{1}{\epsilon_{x}}\frac{\partial\sigma_{xx}}{\partial x} +\frac{1}{\epsilon_{y}}\frac{\partial\sigma_{xy}}{\partial y} -\frac{1}{\epsilon_{x}}\frac{\partial\alpha\hat{p}}{\partial x}\\
s^{2}\rho \hat{u}_{y} + s^{2}\rho_{f}\hat{w}_{y} = \frac{1}{\epsilon_{x}}\frac{\partial\sigma_{yx}}{\partial x} +\frac{1}{\epsilon_{y}}\frac{\partial\sigma_{yy}}{\partial y} -\frac{1}{\epsilon_{y}}\frac{\partial\alpha\hat{p}}{\partial y}
\end{array} \quad (11)
\]

which after multiplying both equations by \(\epsilon_{x}\epsilon_{y}\), introducing the variables \(a = \alpha_{x}\alpha_{y}\), \(b = \alpha_{x}\beta_{y} + \alpha_{y}\beta_{x}\), \(c = \beta_{x}\beta_{y}\), and rearranging terms can be rewritten as:

\[
(s^{2}a + sb + c)(\rho \hat{\mathbf{u}} +\rho_{f}\hat{\mathbf{w}}) = \nabla \cdot \left\{\boldsymbol{\sigma}(\mathbf{u},p)\left(\boldsymbol{\Lambda}_{e} + \frac{1}{s}\boldsymbol{\Lambda}_{p}\right)\right\} \quad (12)
\]

where the tensors \(\boldsymbol{\Lambda}_{e}\) and \(\boldsymbol{\Lambda}_{p}\) are defined as

\[
\begin{bmatrix} \alpha_{y} & 0 \\ 0 & \alpha_{x} \end{bmatrix} + \frac{1}{s}\begin{bmatrix} \beta_{y} & 0 \\ 0 & \beta_{x} \end{bmatrix} = \boldsymbol{\Lambda}_{e} + \frac{1}{s}\boldsymbol{\Lambda}_{p} \quad (13)
\]

#### 3.2.2. Darcy equation

Proceeding similarly with (10b) results in:

\[
\begin{array}{l}
s^{2}\rho_{f}\hat{u}_{x} + s^{2}\rho_{w}\hat{w}_{x} + s\frac{\eta}{\kappa}\hat{w}_{x} = -\frac{1}{\epsilon_{x}}\frac{\partial\hat{p}}{\partial x}\\
s^{2}\rho_{f}\hat{u}_{y} + s^{2}\rho_{w}\hat{w}_{y} + s\frac{\eta}{\kappa}\hat{w}_{y} = -\frac{1}{\epsilon_{y}}\frac{\partial\hat{p}}{\partial y}
\end{array} \quad (14)
\]

Multiplying by \(\epsilon_{x}\epsilon_{y}\) and rearranging yields:

\[
(\tilde{\boldsymbol{\Lambda}}_{e}s + \tilde{\boldsymbol{\Lambda}}_{p})\left(s\rho_{f}\hat{\mathbf{u}} +s\rho_{w}\hat{\mathbf{w}} +\frac{\eta}{\kappa}\hat{\mathbf{w}}\right) = -\nabla \hat{p}, \quad (15)
\]

where the tensors \(\tilde{\boldsymbol{\Lambda}}_{e}\) and \(\tilde{\boldsymbol{\Lambda}}_{p}\) are defined as:

\[
\begin{bmatrix} \alpha_{x} & 0 \\ 0 & \alpha_{y} \end{bmatrix} + \frac{1}{s}\begin{bmatrix} \beta_{x} & 0 \\ 0 & \beta_{y} \end{bmatrix} = \tilde{\boldsymbol{\Lambda}}_{e} + \frac{1}{s}\tilde{\boldsymbol{\Lambda}}_{p}, \quad (16)
\]

The tensors \(\tilde{\boldsymbol{\Lambda}}_{e}\) and \(\tilde{\boldsymbol{\Lambda}}_{p}\) are similar to \(\boldsymbol{\Lambda}_{e}\) and \(\boldsymbol{\Lambda}_{p}\) (respectively) but have their diagonal reversed.

#### 3.2.3. Constitutive relations

Finally, multiplying the Eqs. (10c) and (10e) by \(\epsilon_{x}\epsilon_{y}\), and rearranging terms, we obtain the following PML equations in the Laplace domain:

\[
\begin{array}{l}
s a\hat{\mathbf{e}} +b\hat{\mathbf{e}} +\frac{1}{s} c\hat{\mathbf{e}} = \frac{1}{2} s\left\{\nabla \hat{\mathbf{u}}\boldsymbol{\Lambda}_{e} + (\nabla \hat{\mathbf{u}}\boldsymbol{\Lambda}_{e})^{T}\right\} +\frac{1}{2}\left\{\nabla \hat{\mathbf{u}}\boldsymbol{\Lambda}_{p} + (\nabla \hat{\mathbf{u}}\boldsymbol{\Lambda}_{p})^{T}\right\}\\
-\left(s a\hat{p} +b\hat{p} +\frac{1}{s} c\hat{p}\right) = \nabla \cdot \left\{M(\boldsymbol{\Lambda}_{e}s + \boldsymbol{\Lambda}_{p})(\alpha \hat{\mathbf{u}} +\hat{\mathbf{w}})\right\}
\end{array} \quad (17)
\]

### 3.3. Time-domain formulation of the PML equations

Applying the inverse Laplace transform to the Eqs. (12), (15), (17a), and (17b) we obtain the time domain PML formulation of the Biot's equations:

\[
\begin{array}{l}
\rho (a\ddot{\mathbf{u}} + b\dot{\mathbf{u}} + c\mathbf{u}) + \rho_{f}(a\ddot{\mathbf{w}} + b\dot{\mathbf{w}} + c\mathbf{w}) = \nabla \cdot \left(\boldsymbol{\sigma} \boldsymbol{\Lambda}_{e} + \int_{0}^{t}\boldsymbol{\sigma} \, d\tau \, \boldsymbol{\Lambda}_{p}\right)\\
a\dot{\mathbf{e}} + b\mathbf{e} + c\int_{0}^{t}\mathbf{e} \, d\tau = \frac{1}{2}\left\{\nabla \dot{\mathbf{u}}\boldsymbol{\Lambda}_{e} + (\nabla \dot{\mathbf{u}}\boldsymbol{\Lambda}_{e})^{T} + \nabla \mathbf{u}\boldsymbol{\Lambda}_{p} + (\nabla \mathbf{u}\boldsymbol{\Lambda}_{p})^{T}\right\}\\
-\nabla p = \left(\tilde{\boldsymbol{\Lambda}}_{e}\frac{\partial}{\partial t} +\tilde{\boldsymbol{\Lambda}}_{p}\right)\left(\rho_{f}\dot{\mathbf{u}} +\rho_{w}\dot{\mathbf{w}} +\frac{\eta}{\kappa}\mathbf{w}\right)\\
-\left(a\dot{p} + b p + c\int_{0}^{t}p \, d\tau\right) = \nabla \cdot \left\{M\left(\boldsymbol{\Lambda}_{e}\frac{\partial}{\partial t} +\boldsymbol{\Lambda}_{p}\right)(\alpha \mathbf{u} + \mathbf{w})\right\}
\end{array} \quad (18)
\]

To avoid using auxiliary differential equations or the discrete evaluation of time integrals, we introduce the auxiliary memory variables \(\mathbf{S}(\mathbf{x},t)\), \(\mathbf{E}(\mathbf{x},t)\), and \(\pi(\mathbf{x},t)\) for the stress, strain, and pressure (respectively), defined as:

\[
\mathbf{S}(\mathbf{x},t) = \int_{0}^{t}\mathbf{C}\mathbf{e}(\mathbf{x},\tau) \, d\tau ,\quad \mathbf{E}(\mathbf{x},t) = \int_{0}^{t}\mathbf{e}(\mathbf{x},\tau) \, d\tau ,\quad \pi(\mathbf{x},t) = \int_{0}^{t}p(\mathbf{x},\tau) \, d\tau , \quad (19)
\]

Consequently

\[
\begin{array}{l}
\dot{\mathbf{S}}(\mathbf{x},t) = \mathbf{C}\mathbf{e}(\mathbf{x},t),\quad \ddot{\mathbf{S}}(\mathbf{x},t) = \mathbf{C}\dot{\mathbf{e}}(\mathbf{x},t),\\
\dot{\mathbf{E}}(\mathbf{x},t) = \mathbf{e}(\mathbf{x},t),\quad \ddot{\mathbf{E}}(\mathbf{x},t) = \dot{\mathbf{e}}(\mathbf{x},t),\\
\dot{\pi}(\mathbf{x},t) = p(\mathbf{x},t),\quad \ddot{\pi}(\mathbf{x},t) = \dot{p}(\mathbf{x},t).
\end{array} \quad (20)
\]

For the sake of simplicity, we introduce the following definitions:

\[
\begin{array}{l}
\mathcal{I}f = a\ddot{f} +b\dot{f} +cf\\
\boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S},\pi) = (\dot{\mathbf{S}} -\alpha \dot{\pi}\mathbf{I})\boldsymbol{\Lambda}_{e} + (\mathbf{S} - \alpha \pi \mathbf{I})\boldsymbol{\Lambda}_{p}
\end{array} \quad (21)
\]

where \(\mathcal{I}\) denotes an operator that acts on any scalar, vector, or tensor function \(f\) and \(\boldsymbol{\sigma}^{\mathrm{PML}}\) is the PML stress tensor. Thus, replacing the relations given in (20) into (18) and introducing the previous definitions, the time-domain PML formulation becomes: find \(\mathbf{u}\), \(\mathbf{w}\), \(\pi\), and \(\mathbf{S}\) satisfying:

\[
\begin{array}{l}
\rho \mathcal{I}\mathbf{u} + \rho_{f}\mathcal{I}\mathbf{w} = \nabla \cdot \boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S},\pi) \qquad \mathrm{in} \ \Omega \times T \\
\mathcal{D}(\mathcal{I}\mathbf{S}) = \frac{1}{2}\left\{\nabla \dot{\mathbf{u}}\boldsymbol{\Lambda}_{e} + \boldsymbol{\Lambda}_{e}(\nabla \dot{\mathbf{u}})^{T} + \nabla \mathbf{u}\boldsymbol{\Lambda}_{p} + \boldsymbol{\Lambda}_{p}(\nabla \mathbf{u})^{T}\right\} \qquad \mathrm{in} \ \Omega \times T \\
-\nabla \dot{\pi} = \left(\tilde{\boldsymbol{\Lambda}}_{e}\frac{\partial}{\partial t} + \tilde{\boldsymbol{\Lambda}}_{p}\right) \left( \rho_{f} \dot{\mathbf{u}} + \rho_{w} \dot{\mathbf{w}} + \frac{\eta}{\kappa} \mathbf{w} \right) \qquad \mathrm{in} \ \Omega \times T \\
-\mathcal{I} \pi = \nabla \cdot \left\{ M \left( \boldsymbol{\Lambda}_{e} \frac{\partial}{\partial t} + \boldsymbol{\Lambda}_{p} \right) (\alpha \mathbf{u} + \mathbf{w}) \right\} \qquad \mathrm{in} \ \Omega \times T
\end{array} \quad (22)
\]

where \(\mathcal{D}\) denotes the compliance operator, which takes the stress tensor as argument and returns the strain tensor (\(\mathbf{e} = \mathcal{D}(\boldsymbol{\sigma})\)).

**Remark.** The tensor \(\mathbf{S}\) is symmetric and, therefore, only three additional scalar functions are introduced as unknowns by the time-domain PML formulation.

---

## 4. Hybrid formulation

With the time-domain formulation given in (22), the vector functions \(\mathbf{u}\), \(\mathbf{w}\), the scalars \(p\), \(\pi\), and the symmetric tensor field \(\mathbf{S}\) must be solved simultaneously on the whole domain \(\Omega = \Omega^{\mathrm{RD}} \cup \Omega^{\mathrm{PML}} \subset \mathbb{R}^2\), which may lead to computationally expensive problems. Therefore, we define a hybrid formulation where the problem given in (22) is split into two sub-problems defined separately on \(\Omega^{\mathrm{RD}}\) and \(\Omega^{\mathrm{PML}}\) but coupled through boundary conditions on \(\Gamma_I\) (see Fig. 1).

Let \(\{\mathbf{u}_1, \mathbf{w}_1\}\) and \(\{\mathbf{u}_2, \mathbf{w}_2\}\) be the solid and relative fluid displacements defined separately on \(\Omega^{\mathrm{RD}}\) and \(\Omega^{\mathrm{PML}}\), respectively. The hybrid PML formulation for the solid and fluid displacements, pore pressure, stress history, and pore pressure history reads: find \(\{\mathbf{u}_1, \mathbf{w}_1, p\}\) and \(\{\mathbf{u}_2, \mathbf{w}_2, \pi, \mathbf{S}\}\) satisfying:

\[
\begin{array}{rcll}
\rho \ddot{\mathbf{u}}_1 + \rho_f \ddot{\mathbf{w}}_1 &=& \nabla \cdot \boldsymbol{\sigma}(\mathbf{u}_1, p) & \text{in } \Omega^{\mathrm{RD}} \times T \quad (23a)\\
-\nabla p &=& \rho_f \ddot{\mathbf{u}}_1 + \rho_w \ddot{\mathbf{w}}_1 + \frac{\eta}{\kappa} \dot{\mathbf{w}}_1 & \text{in } \Omega^{\mathrm{RD}} \times T \quad (23b)\\
-\dot{p} &=& \nabla \cdot \{M(\alpha \dot{\mathbf{u}}_1 + \dot{\mathbf{w}}_1)\} & \text{in } \Omega^{\mathrm{RD}} \times T \quad (23c)\\
\rho \mathcal{I} \mathbf{u}_2 + \rho_f \mathcal{I} \mathbf{w}_2 &=& \nabla \cdot \boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S}, \pi) & \text{in } \Omega^{\mathrm{PML}} \times T \quad (23d)\\
\mathcal{D}(\mathcal{I} \mathbf{S}) &=& \frac{1}{2} \left\{ \nabla \dot{\mathbf{u}}_2 \boldsymbol{\Lambda}_p + \boldsymbol{\Lambda}_p (\nabla \dot{\mathbf{u}}_2)^T + \nabla \mathbf{u}_2 \boldsymbol{\Lambda}_e + \boldsymbol{\Lambda}_e (\nabla \mathbf{u}_2)^T \right\} & \text{in } \Omega^{\mathrm{PML}} \times T \quad (23e)\\
-\nabla \dot{\pi} &=& \left( \tilde{\boldsymbol{\Lambda}}_e \frac{\partial}{\partial t} + \tilde{\boldsymbol{\Lambda}}_p \right) \left( \rho_f \dot{\mathbf{u}}_2 + \rho_w \dot{\mathbf{w}}_2 + \frac{\eta}{\kappa} \mathbf{w}_2 \right) & \text{in } \Omega^{\mathrm{PML}} \times T \quad (23f)\\
-\mathcal{I} \pi &=& \nabla \cdot \left\{ M \left( \boldsymbol{\Lambda}_e \frac{\partial}{\partial t} + \boldsymbol{\Lambda}_p \right) (\alpha \mathbf{u}_2 + \mathbf{w}_2) \right\} & \text{in } \Omega^{\mathrm{PML}} \times T \quad (23g)
\end{array}
\]

subject to zero initial values and the Dirichlet and Neumann boundary conditions listed below (see Fig. 1(b)).

\[
\begin{array}{rcll}
\boldsymbol{\sigma}(\mathbf{u}_1, p) \cdot \mathbf{n}_1 &=& \mathbf{g} & \text{in } \Gamma_g \times T \quad (24a)\\
\boldsymbol{\sigma}(\mathbf{u}_1, p) \cdot \mathbf{n}_1 = p\mathbf{I} \cdot \mathbf{n}_1 &=& \mathbf{0} & \text{in } \Gamma^{\mathrm{RD}}_N \times T \quad (24b)\\
\boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S}, \pi) \cdot \mathbf{n}_2 = \dot{\pi} \mathbf{I} \cdot \mathbf{n}_2 &=& \mathbf{0} & \text{in } \Gamma^{\mathrm{PML}}_N \times T \quad (24c)\\
\mathbf{u}_2 &=& \mathbf{0} & \text{in } \Gamma^{\mathrm{PML}}_D \times T \quad (24d)\\
\mathbf{w}_2 &=& \mathbf{0} & \text{in } \Gamma^{\mathrm{PML}}_D \times T \quad (24e)\\
\pi &=& 0 & \text{in } \Gamma^{\mathrm{PML}}_D \times T \quad (24f)
\end{array}
\]

where \(\mathbf{n}_1\) and \(\mathbf{n}_2\) are outward pointing normal vectors to \(\Omega^{\mathrm{RD}}\) and \(\Omega^{\mathrm{PML}}\) (\(\mathbf{n}_1 = -\mathbf{n}_2\) in \(\Gamma_I\)), respectively (see Fig. 1(a)).

Finally, to couple both equations, the continuity of displacements, tractions, and pressures must be imposed on the interface as follows:

\[
\begin{array}{rcll}
\boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S}, \pi) \mathbf{n}_1 + \boldsymbol{\sigma}(\mathbf{u}_1, p) \mathbf{n}_2 &=& \mathbf{0} & \text{in } \Gamma_I \times T \quad (25a)\\
\mathbf{u}_1 &=& \mathbf{u}_2 & \text{in } \Gamma_I \times T \quad (25b)\\
\mathbf{w}_1 &=& \mathbf{w}_2 & \text{in } \Gamma_I \times T \quad (25c)\\
p &=& \dot{\pi} & \text{in } \Gamma_I \times T \quad (25d)
\end{array}
\]

### 4.1. Extension to M-PML

The previous formulations are prone to time instabilities depending on the media properties and the form of the stretching functions. To transform the PML layer in the hybrid problem to the multiaxial case, the functions \(\alpha\) and \(\beta\) (see Eq. (7)) must be redefined as:

\[
\begin{array}{l}
\alpha^{*}(x) = \alpha^{*}(y) = 1\\
\beta^{*}(x,y) = \beta(x) + p^{(y/x)}\beta(y)\\
\beta^{*}(x,y) = \beta(y) + p^{(x/y)}\beta(x)
\end{array} \quad (26)
\]

where \(p^{(y/x)}\) and \(p^{(x/y)}\) are constant parameters that allow the fine-tuning of the M-PML layer. This modified definition of the scaling and attenuation functions does not introduce changes to the previous formulation.

### 4.2. Variational formulation

The interface conditions (25a) and (25b) are fulfilled by using a single continuous function space for \(\mathbf{u}_1\) and \(\mathbf{u}_2\). Similarly, (25c) is fulfilled by using a single continuous function space for \(\mathbf{w}_1\) and \(\mathbf{w}_2\). A Lagrange multiplier is used for the condition (25d).

Thus, using the following function spaces

\[
\begin{array}{l}
V = \{\mathbf{u} \in [H^1(\Omega)]^2 \text{ s.t. } \mathbf{u} = \mathbf{0} \text{ on } \Gamma_D^{\mathrm{PML}} \times T\}\\
Q = \{p \in L^2(\Omega^{\mathrm{RD}})\}\\
Q_0 = \{p \in L^2(\Omega^{\mathrm{PML}}) \text{ s.t. } p = 0 \text{ on } \Gamma_D^{\mathrm{PML}} \times T\}\\
\mathcal{T} = \{\mathbf{S} \in [L^2(\Omega^{\mathrm{PML}})]^{2\times 2}\}\\
L = \{l \in L^2(\Gamma_I)\}
\end{array} \quad (27)
\]

the weak form of the system of PDEs given in (23) reads: find \((\mathbf{u},\mathbf{w},p,\pi,\mathbf{S},\lambda_p)\) for all \((\tilde{\mathbf{u}},\tilde{\mathbf{w}},\tilde{p},\tilde{\pi},\tilde{\mathbf{S}},\tilde{\lambda}_p) \in V \times V \times Q \times Q_0 \times \mathcal{T} \times L\) solution to:

\[
\begin{array}{l}
\displaystyle \int_{\Omega^{\mathrm{RD}}} (\rho \ddot{\mathbf{u}} + \rho_f \ddot{\mathbf{w}}) \cdot \tilde{\mathbf{u}} \, d\Omega + \int_{\Omega^{\mathrm{RD}}} \boldsymbol{\sigma}(\mathbf{u},p) : \nabla \tilde{\mathbf{u}} \, d\Omega = \int_{\Gamma_g} \mathbf{g} \cdot \tilde{\mathbf{u}} \, d\Gamma \\
\displaystyle \quad \int_{\Omega^{\mathrm{PML}}} (\rho \mathcal{I}\mathbf{u} + \rho_f \mathcal{I}\mathbf{w}) \cdot \tilde{\mathbf{u}} \, d\Omega + \int_{\Omega^{\mathrm{PML}}} \boldsymbol{\sigma}^{\mathrm{PML}}(\mathbf{S},\pi) : \nabla \tilde{\mathbf{u}} \, d\Omega = 0 \\
\displaystyle \quad \int_{\Omega^{\mathrm{RD}}} \left( \rho_f \ddot{\mathbf{u}} + \rho_w \ddot{\mathbf{w}} + \frac{\eta}{\kappa} \dot{\mathbf{w}} \right) \cdot \tilde{\mathbf{w}} \, d\Omega - \int_{\Omega^{\mathrm{RD}}} p \nabla \cdot \tilde{\mathbf{w}} \, d\Omega + \int_{\Gamma_g} p \mathbf{n} \cdot \tilde{\mathbf{w}} \, d\Gamma = 0 \\
\displaystyle \quad \int_{\Omega^{\mathrm{PML}}} \left( \tilde{\boldsymbol{\Lambda}}_e \frac{\partial}{\partial t} + \tilde{\boldsymbol{\Lambda}}_p \right) \left( \rho_f \dot{\mathbf{u}} + \rho_w \dot{\mathbf{w}} + \frac{\eta}{\kappa} \mathbf{w} \right) \cdot \tilde{\mathbf{w}} \, d\Omega - \int_{\Omega^{\mathrm{PML}}} \nabla \dot{\pi} \cdot \tilde{\mathbf{w}} \, d\Omega = 0 \\
\displaystyle \quad \int_{\Omega^{\mathrm{RD}}} \dot{p} \tilde{p} \, d\Omega + \int_{\Omega^{\mathrm{RD}}} \nabla \cdot \{ M(\alpha \dot{\mathbf{u}} + \dot{\mathbf{w}}) \} \tilde{p} \, d\Omega = 0 \\
\displaystyle \quad \int_{\Omega^{\mathrm{PML}}} \mathcal{I} \pi \, \tilde{\pi} \, d\Omega + \int_{\Omega^{\mathrm{PML}}} \nabla \cdot \left\{ M \left( \boldsymbol{\Lambda}_e \frac{\partial}{\partial t} + \boldsymbol{\Lambda}_p \right) (\alpha \mathbf{u} + \mathbf{w}) \right\} \tilde{\pi} \, d\Omega = 0 \\
\displaystyle \quad \int_{\Omega^{\mathrm{PML}}} \mathcal{D}(\mathcal{I}\mathbf{S}) : \tilde{\mathbf{S}} \, d\Omega - \frac{1}{2} \int_{\Omega^{\mathrm{PML}}} \left\{ \nabla \dot{\mathbf{u}} \boldsymbol{\Lambda}_p + \boldsymbol{\Lambda}_p (\nabla \dot{\mathbf{u}})^T + \nabla \mathbf{u} \boldsymbol{\Lambda}_e + \boldsymbol{\Lambda}_e (\nabla \mathbf{u})^T \right\} : \tilde{\mathbf{S}} \, d\Omega = 0 \\
\displaystyle \quad \int_{\Gamma_I} \tilde{\lambda}_p (p - \dot{\pi}) \, d\Gamma + \int_{\Gamma_I} (\tilde{p} - \tilde{\pi}) \lambda_p \, d\Gamma = 0
\end{array} \quad (28)
\]

where \(\lambda_p\) is a Lagrange multiplier used to impose the coupling condition (25d).

---

## 5. Numerical experiments

Four experiments were developed to evaluate the performance and accuracy of the proposed hybrid PML formulation. The first experiment (referred to as Experiment 1) considers a homogeneous poroelastic half-space where the water table is located at the free surface. In the second (Experiment 2), three horizontally layered media over a half-space and water table below the second layer was considered (Fig. 2(b)). The third (Experiment 3) adds a more realistic stratification, including an outcropping (see Fig. 2(c)). Finally, the fourth experiment (Experiment 4) is a case consisting of a homogeneous media where the three types of poroelastic waves predicted by the Biot theory propagate with similar energy. The geometry in this case is the same as in Experiment 2, but assuming a homogeneous medium. The material parameters used in the four cases are listed in Table 1. Set 1 corresponds approximately to a soft rock, while sets 2 to 5 represent standard soil parameters, from loose to dense sands. Set 6 is the same as Set 3 but with a reduced fluid viscosity and an increased shear modulus of the solid phase. These parameter modifications were tuned to generate three distinct wave propagation velocities (i.e., clearly differentiated) while also ensuring an equivalent energy level (for improved visualization).

To obtain a reference solution, we solved Biot's equations given in (1) on an extended domain with dimensions large enough to avoid wave reflections at the exterior boundaries of \(\Omega^{\mathrm{RD}}\). A zero displacement Dirichlet condition was imposed to both \(\mathbf{u}\) and \(\mathbf{w}\) on the exterior boundary \(\Gamma_D\). A simulation time of 2.0 s (half of the runtime used for the PML experiments) was enough for comparison purposes in the extended domain simulations. For layered medium, the regular domain was embedded in the extended domain, and layers were extended to the exterior boundary (see Fig. 2(d)). Additionally, a first-order paraxial boundary condition [32] was implemented to solve Biot's equations also for comparison purposes. The PML domain shown in Fig. 2 was removed for paraxial simulations and replaced by this ABC at this boundary.

An alternative approach for solving the proposed PML method would be to solve Eq. (22) directly. In this approach, all variables are solved on the combined domain \(\Omega = \Omega^{\mathrm{RD}} \cup \Omega^{\mathrm{PML}}\), which would result in increased computational costs during the solving stage. However, in contrast to the hybrid form, the direct version is simpler and easier to implement since no special treatment is required for \(\Omega^{\mathrm{RD}}\), and no coupling conditions are necessary at \(\Gamma_I\). Therefore, for the purpose of comparison and to highlight the advantages of the hybrid method, the direct approach was also solved (this solution will be referred as direct PML in the following sections).

**Table 1**
Material parameters, characteristic frequency of the medium, and wave propagation velocities of the porous media considered in the experiments. Set 1 was taken from Tables 1 and 2 in [6], whereas the authors defined sets from 2 to 5 to obtain soil-like wave velocities. Set 6 was tuned to obtain a high contrast between \(f_r\) and \(f_c\) to ensure the propagation of the three types of poroelastic waves with similar energies.

| Parameter | Set 1 | Set 2 | Set 3 | Set 4 | Set 5 | Set 6 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| \(\rho_s\) (kg/m³) | 2650 | 2600 | 2600 | 2600 | 2600 | 2600 |
| \(\rho_f\) (kg/m³) | 900 | 1.29 | 1000 | 1000 | 1.29 | 1000 |
| \(K_s\) (N/m²) | \(12 \times 10^9\) | \(2.3 \times 10^8\) | \(2.3 \times 10^8\) | \(2.5 \times 10^8\) | \(4.6 \times 10^8\) | \(2.3 \times 10^8\) |
| \(K_f\) (N/m²) | \(2 \times 10^9\) | \(1.4 \times 10^5\) | \(2 \times 10^9\) | \(2 \times 10^9\) | \(1.4 \times 10^5\) | \(2 \times 10^{19}\) |
| \(K_b\) (N/m²) | \(10 \times 10^9\) | \(1.5 \times 10^8\) | \(1.5 \times 10^8\) | \(1.7 \times 10^8\) | \(4.2 \times 10^8\) | \(1.5 \times 10^8\) |
| \(\mu_b\) (N/m²) | \(5 \times 10^9\) | \(1.33 \times 10^8\) | \(1.33 \times 10^8\) | \(2.84 \times 10^8\) | \(4.44 \times 10^8\) | \(1.73 \times 10^8\) |
| \(\tau\) | 1.2 | 1 | 1 | 1 | 1 | 1 |
| \(\kappa\) (m²) | \(1 \times 10^{-12}\) | \(8 \times 10^{-9}\) | \(8 \times 10^{-9}\) | \(1 \times 10^{-10}\) | \(1 \times 10^{-11}\) | \(8 \times 10^{-9}\) |
| \(\phi\) | 0.3 | 0.2 | 0.2 | 0.2 | 0.2 | 0.3 |
| \(\eta\) (Pa s) | \(1 \times 10^{-3}\) | \(2 \times 10^{-5}\) | \(1 \times 10^{-3}\) | \(1 \times 10^{-3}\) | \(2 \times 10^{-5}\) | \(2 \times 10^{-4}\) |
| \(f_r\) (Hz) | 15 | 15 | 15 | 15 | 15 | 30 |
| \(f_c\) (Hz) | 44 | 2106 | 24 | 31 | 849 | 3500.8 |
| \(c_s\) (m/s) | 960 | 300 | 300 | 400 | 500 | 288 |
| \(c_{1p}\) (m/s) | 2366 | 474 | 705 | 746 | 367 | 55402 |
| \(c_{2p}\) (m/s) | 775 | 329 | 425 | 513 | 330 | 552 |

\(f_r\): central frequency of the Ricker pulse, \(f_c\): characteristic frequency of the medium, \(c_{1p}\): fast primary wave velocity, \(c_{2p}\): slow primary wave velocity, \(c_s\): shear wave velocity.

A vertical load defined by a Ricker wavelet was applied on a \(0.3\,\text{m}\) width stripe of the free surface in the three experiments (see Fig. 2(a)). The expression defining the source is as follows:

\[
\mathbf{g}(t) = A\begin{bmatrix} 0 \\ S(t) \end{bmatrix} \quad \text{with} \quad S(t) = \frac{(0.25u^2 - 0.5)e^{-0.25u^2} - 13e^{13.5}}{0.5 + 13e^{13.5}} \quad \text{and} \quad 0 \leq t \leq \frac{6\sqrt{6}}{\omega_r} \quad (29)
\]

where \(A\) denotes the pulse amplitude and \(u = \omega_r t - 3\sqrt{6}\). In the previous expression, \(\omega_r = 2\pi f_r\) is the characteristic central angular frequency of the pulse. In all the experiments, \(A\) was set to \(10^4\,\text{N/m}\), and \(f_r\) was set to \(15\,\text{Hz}\), except for the Experiment 4, where \(f_r\) was set to \(30\,\text{Hz}\). In Experiments 1 to 3, the frequency of the source was adjusted to obtain a frequency spectrum similar to those obtained from geophysical seismic surveys, which makes it suitable for the simulations.

Four out of the six media used in the simulations, characterized by the physical parameters in Sets 1, 2, 4, and 5, demonstrate dispersive behavior. This is confirmed by the characteristic frequency of the medium \((f_c)\), which satisfies the inequality \(f_c = (\eta \phi) / (2\pi \rho_f \tau \kappa) > f_r\) (see Table 1). Consequently, the slowest primary wave does not propagate in these media [5,6]. On the other hand, the media related to the physical parameters in Set 3 and 6 (see Table 1) are non-dispersive, i.e., three propagating waves are present. However, their shear wave velocity is slower than its slowest volumetric wave, and therefore, the discretization was chosen considering only the fastest primary wave and the shear wave velocities for all media. The domain was discretized using triangular cells, and the element size \((\Delta x)\) was adjusted to have a minimum of 12 elements per shortest wavelength, with the biggest possible element being chosen. In all simulations, elements in the vicinity of the surface load were refined to a size of \(\Delta x = 0.15\,\text{m}\).

A third-order version of the scaling and attenuation profiles given in (8) (with \(m = 3\)) was used for the simulations. The width of the PML, denoted by \(L_{\mathrm{PML}}\), was chosen to be ten times the element size once \(\Delta x\) was fixed (see Table 1). The constant \(\beta_0\) was calculated using the expression (9) with \(R = 10^{-4}\), while \(\alpha_0\) was fixed at 5. For all the experiments considering M-PML stretching functions, \(p^{(y/x)}\) and \(p^{(x/y)}\) were used as 0.01. A summary of the discretization and PML parameters is presented in Table 2.

The Newmark-\(\beta\) method with \(\beta = 1/4\) and \(\gamma = 1/2\) (i.e., without numerical damping) was used for time discretization [33]. The time-step \(\Delta t\) was calculated using the Courant-Friedrichs-Lewy criteria:

\[
\Delta t < \mathrm{CFL} \frac{\Delta x}{c_{1p}} \quad (30)
\]

where \(c_{1p}\) represents the velocity of the fast primary wave and \(\mathrm{CFL} = 0.75\) is the Courant-Friedrichs-Lewy number.

**Table 2**
Element sizes and PML parameters used in the three experiments with homogeneous and layered media.

| Experiment | General parameters | | | PML parameters | | |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| | \(\Delta x_{\text{global}}\) (m) | \(\Delta x_{\text{source}}\) (m) | \(\Delta t\) (s) | \(R\) | \(L_{\mathrm{PML}}\) (m) | \(\alpha_0\) |
| 1 | 7.8 | 0.15 | \(10^{-3}\) | \(10^{-4}\) | 78 | 5 |
| 2 and 3 | 1.4 | 0.15 | \(10^{-3}\) | \(10^{-4}\) | 14 | 5 |
| 4 | 0.8 | 0.15 | \(10^{-3}\) | \(10^{-4}\) | 8 | 5 |

\(\Delta x_{\text{global}}\): global element size; \(\Delta x_{\text{local}}\): element size refinement near to external source
\(\Delta t\): time step; \(R\): reflection coefficient; \(L_{\mathrm{PML}}\): width of the PML layer.

### 5.1. Metrics for performance evaluation

To evaluate the performance of the tested methods, the poroelastic energy on \(\Omega^{\mathrm{RD}}\) was estimated and compared against the reference solution using the following expression:

\[
\begin{aligned}
E(t_k) = & \frac{1}{2} \int_{\Omega^{\mathrm{RD}}} \rho \dot{\mathbf{u}}(\mathbf{x}, t_k) \cdot \dot{\mathbf{u}}(\mathbf{x}, t_k) \, d\Omega
+ \frac{1}{2} \int_{\Omega^{\mathrm{RD}}} \mathbf{C} \mathbf{e}(\mathbf{u}(\mathbf{x}, t_k), t_k) : \mathbf{e}(\mathbf{u}(\mathbf{x}, t_k), t_k) \, d\Omega \\
& + \frac{1}{2} \int_{\Omega^{\mathrm{RD}}} \rho_w \dot{\mathbf{w}}(\mathbf{x}, t_k) \cdot \dot{\mathbf{w}}(\mathbf{x}, t_k) \, d\Omega
+ \frac{1}{2} \int_{\Omega^{\mathrm{RD}}} \frac{1}{M} p(\mathbf{x}, t_k) \cdot p(\mathbf{x}, t_k) \, d\Omega \\
& + \int_{\Omega^{\mathrm{RD}}} \rho_f \dot{\mathbf{u}}(\mathbf{x}, t_k) \cdot \dot{\mathbf{w}}(\mathbf{x}, t_k) \, d\Omega
\end{aligned} \quad (31)
\]

Finally, we obtained traces of \(\mathbf{u}\), \(\mathbf{w}\), and \(p\) at different locations \(\mathbf{x}_i\) in \(\Omega^{\mathrm{RD}}\) (see Fig. 2). Normalized error metric are defined as:

\[
\begin{array}{l}
e_{\mathbf{u}}(\mathbf{x}_i, t_k) = \dfrac{\|\mathbf{u}_{\text{ref}}(\mathbf{x}_i, t_k) - \mathbf{u}(\mathbf{x}_i, t_k)\|_2}{\max_{t_k} \|\mathbf{u}_{\text{ref}}(\mathbf{x}_i, t_k)\|_2} \quad (32a)\\
e_{\mathbf{w}}(\mathbf{x}_i, t_k) = \dfrac{\|\mathbf{w}_{\text{ref}}(\mathbf{x}_i, t_k) - \mathbf{w}(\mathbf{x}_i, t_k)\|_2}{\max_{t_k} \|\mathbf{w}_{\text{ref}}(\mathbf{x}_i, t_k)\|_2} \quad (32b)\\
e_{p}(\mathbf{x}_i, t_k) = \dfrac{|p_{\text{ref}}(\mathbf{x}_i, t_k) - p(\mathbf{x}_i, t_k)|}{\max_{t_k} |p_{\text{ref}}(\mathbf{x}_i, t_k)|} \quad (32c)
\end{array}
\]

where \(|\cdot|\) is the absolute value and \(\|\mathbf{f}\|_2 = \sqrt{\sum_i f_i^2}\) the vector 2-norm. The sub-index \(()_{\text{ref}}\) denotes the reference solution obtained in the extended domain simulations.

### 5.2. Implementation

All experiments were solved using the open-source computing platform FEniCS [34,35]. To implement the hybrid PML problem, Multiphenics [36] was used as a complementary tool. The finite element meshes were generated using the Frontal-Delaunay algorithm in Gmsh [37]. Discretization of the displacements (\(\mathbf{u}\) and \(\mathbf{w}\)) and pressures (\(p\) and \(\pi\)) was carried out using continuous Lagrange polynomials of second and first order, respectively. The stress history \(\mathbf{S}\) was discretized using discontinuous Lagrange polynomials of first order.

The Multifrontal Massively Parallel sparse direct Solver (MUMPS) and the Generalized Minimal Residual Method (GMRES) solvers were used to solve the linear systems obtained after assembling the discrete weak forms of the problems [38-40]. FEniCS is built with PETSc as linear algebra backend [38-40] and supports both solvers by default. The iterative solver was used only for the extended domain simulations with relative and absolute tolerances of \(10^{-7}\) and \(10^{-9}\), respectively [38]. To accelerate the convergence, the linear system of the extended domain simulation was right-preconditioned using the Parallel ILU preconditioner HYPRE-Euclid [41]. The direct solver was used with default parameters [34].

---

## 6. Results

In the upcoming sections, we will present and analyze energy graphs, traces, error traces, and snapshots of the propagating waves for all three experiments. Through these analyses, we aim to provide a comprehensive understanding of the simulations and their outcomes.

**Table 3**
Number of degrees of freedom (DOFs) solved for extended, paraxial, direct PML, and hybrid PML problems.

| Experiment | Extended | Paraxial | Direct PML | Hybrid PML |
| :--- | :--- | :--- | :--- | :--- |
| 1 | 8,550,901 | 425,518 | 1,054,006 | 607,424 |
| 2 | 13,043,765 | 462,513 | 1,137,430 | 651,758 |
| 3 | 12,759,950 | 480,229 | 1,181,065 | 676,810 |
| 4 | 28,176,668 | 1,392,081 | 3,182,758 | 1,715,750 |

### 6.1. Experiment 1: homogeneous half-space

Fig. 3(a) shows the poroelastic energy estimated using (31) for extended, paraxial, and PML simulations in a homogeneous half-space domain. The results demonstrate good agreement between the PML and reference solutions in the first 2 s of simulation time, although small differences are observed. The energy decay rate of the PML simulations is greater than that of the paraxial case, which consistently decays but at a lower rate. The energy obtained from the hybrid PML formulation and the direct PML implementation showed no observable differences, indicating that both methods provide equivalent solutions (see Fig. 3(a)). Additionally, the number of degrees of freedom (DOFs) in the hybrid case is approximately 1.7 times less than the direct PML approach (as shown in Table 3) because the tensor \(\mathbf{S}\) does not need to be solved in \(\Omega^{\mathrm{RD}}\). Consequently, the hybrid formulation is significantly less computationally expensive than the direct PML implementation while maintaining the same properties.

Fig. 3(b) compares the energies obtained from the hybrid PML and M-PML formulations. During the first second of the simulation, the results are similar, and only small differences are observed. As the simulation progresses, the energy obtained with M-PML stretching functions is slightly larger than that obtained with uniaxial functions, indicating slightly worse performance. However, during the last second of simulation, the energy of the uniaxial case stops decaying showing even a slightly increase, while the energy of the M-PML case continues decaying. The reduced performance of M-PML (compared to PML) during the first seconds of simulation is because the stretching functions are not perfectly matched in \(\Gamma_I\), generating spurious reflections. In fact, a M-PML absorbing layer can be interpreted as a sponge rather than a PML [42]. This is because the coupling of two damping directions causes the loss of the perfectly matched layer characteristic of Berenger's technique [7]. Thus, the theoretical reflection coefficient for an infinite M-PML is not longer zero prior the discretization.

Traces of the solutions and errors calculated using (32) for the two locations shown in Fig. 2(a) are presented in Figs. 4 and 5. There is a good match between the hybrid PML and the extended domain simulations at both locations, in contrast to the paraxial case where reflections are observed. Looking at the errors, the superior performance of the hybrid PML method is evident, showing an improvement of at least three orders of magnitude compared to the paraxial case. In the same figure, the error obtained by using M-PML stretching functions in the hybrid simulation is depicted. As mentioned previously, results obtained using M-PML compared to PML are slightly worse between 0.5 and 1.5 s approximately because the stretching functions are not perfectly matched in this case. However, the results are still better than those obtained with the paraxial boundary conditions and more stable in time compared to the hybrid PML simulation (see Fig. 3). Finally, screenshots of the solutions at different times can be found in Fig. 6.

### 6.2. Experiment 2: horizontally-layered domain

Fig. 7(a) shows the poroelastic energies for extended, paraxial, and PML simulations in the horizontally-layered domain. Similarly, the results indicate good agreement between the PML and reference solutions within the first 2 s of simulation, although some small differences are observed. The energy decay rate of the PML simulations is greater than that of the paraxial case, which decays slowly but consistently. In this experiment, it is not evident from the plots that outgoing waves are reaching the boundaries of \(\Omega^{\mathrm{RD}}\), as the plateaus observed in Fig. 3 are absent. This is due to reflections at the interfaces between layers (see Figs. 2(b) and 10).

The energies obtained using the hybrid formulation and the direct PML implementation produced almost the same results, illustrating that both methods provide similar solutions (see Fig. 3(a)) as expected. Additionally, the number of DOFs in the hybrid case is approximately 1.8 times fewer than the direct approach (as shown in Table 3). Therefore, the hybrid formulation is less computationally expensive than the direct implementation while providing equivalent solutions. Regarding Fig. 7(b), in terms of decay rate, during the first 2.5 s of the simulation, the results obtained using PML and M-PML in the hybrid formulation are similar, and slightly greater errors are observed with M-PML between 1.5 and 2.5s. The worst performance of M-PML in this time window is because the stretching functions are not in perfectly matched at \(\Gamma_I\), as mentioned in previous paragraphs, which generates spurious reflections [42]. However, during the last second, the energy of the uniaxial case decays more slowly than the M-PML case, which keeps a constant rate.

Figs. 8 and 9 present traces of the solutions and errors calculated using (32) for two locations shown in Fig. 2(b). The results show a good match between the hybrid PML and extended domain simulations at both locations, in contrast to the paraxial case where reflections are observed. The error analysis confirms the marginally superior performance of the hybrid PML method with respect to the M-PML case, but shows a considerable improvement compared to the paraxial case. The difference between PML and M-PML results is due to imperfect matching of stretching functions in the last case. Nevertheless, the results obtained with M-PML are still superior to those obtained with paraxial boundary conditions and more stable over time than the hybrid PML simulation (see Fig. 7). Finally, Fig. 10 presents screenshots that depict the propagation of waves at different times. These snapshots clearly show the transitions between different layers and the behavior of waves within each medium. For example, the layer with lower permeabilities (as shown in Table 1 and Fig. 2(b)) exhibited smaller relative fluid displacements ( \(\mathbf{w}\) ) due to the small hydraulic conductivity. Additionally, media filled with air had lower pressures. Despite the complex behavior of waves in different media, the PML effectively absorbed and attenuated the waves during the simulation.

### 6.3. Experiment 3: layered domain with outcropping

As observed in the previous experiments, the hybrid PML implementation produced superior results compared to the paraxial case. The energy plots for this experiment also shows an excellent agreement between the PML and reference solutions within the initial 2 s of simulation, with negligible differences, as illustrated in Fig. 11(a). Although the decay of energy was consistent over time in all methods, PML exhibited faster energy decay compared to the paraxial case. The energy plots indicate that the solutions obtained using the hybrid and direct PML implementations are almost indistinguishable, indicating that both methods produce nearly identical results as expected. This observation is consistent with the results of the previous experiments. Additionally, the hybrid case had a reduction in the number of DOFs by almost 1.8 times compared to the direct PML implementation (as shown in Table 3).

Fig. 11(b) compares the energies obtained using the hybrid PML and M-PML formulations. The results show that the two methods provide almost identical solutions up to 2 seconds of the simulation, while some differences can be observed between 2 and 3 s, where PML slightly surpasses the multiaxial solution. Interestingly, after 3 s of the simulation, the energy decay rate using M-PML is better with respect to the uniaxial case, highlighting the improved time-stability of M-PML stretching functions.

Traces of the solutions and their corresponding errors, calculated using (32), for the two locations illustrated in Fig. 2(c), are presented in Figs. 12 and 13, respectively. These figures confirm the excellent agreement between the hybrid PML and the extended domain simulation at both locations, with no observable reflections unlike the paraxial case. The error analysis shows a considerable improvement compared to the paraxial case and slightly superior performance of the hybrid PML method with respect to the M-PML case. However, it is worth noting that the results obtained with M-PML are still better than those obtained with paraxial boundary conditions and exhibit better stability over time than the hybrid PML simulation. The S-wave is not visible in the pressure trace and screenshots because the fluid only propagates volumetric waves, as expected.

Fig. 14 shows snapshots of the solutions at different times. Due to the complex interfaces between materials in the medium, the waves interact in intricate ways generating complex patterns of reflected waves at interfaces as can be seen in the figure. Despite this complexity, the PML effectively absorbed the waves at the boundary of \(\Omega^{\mathrm{RD}}\), and no unwanted reflections were observed.

### 6.4. Experiment 4: homogeneous domain under high frequency loading

Propagating waves were effectively absorbed at the boundaries using the proposed hybrid PML and adapted M-PML formulations (see Fig. 15), showing minimal differences when compared to the reference solution. These findings are consistent with previous results. In contrast to the previous experiments, the three types of waves predicted by the Biot theory are in this case easily distinguishable. For example, in the traces presented in Fig. 16, the fast P-wave arrives first (as expected), followed by the slow P-wave, and finally, the S-wave arrives in third place, because their wave velocities according to Table 1. These waves can also be distinguished in Fig. 17, which presents screenshots of the propagating waves at different time steps.

Finally, the errors of the traces are presented in Fig. A.1 in Appendix A.

---

## 7. Conclusions

We proposed a hybrid formulation of the PML method for the simulation of poroelastic waves in truncated domains. Compared to other methods, the proposed method introduces only three additional scalar unknowns, i.e., the components of the symmetric stress-history tensor, reducing the number of unknowns with respect to current split- and unsplit-field formulations. Due to the special treatment of the interior domain, the hybrid PML formulation significantly reduced the number of DOFs compared to the direct PML implementation and even more dramatically compared to the extended domain simulations. Compared to the paraxial case, the hybrid PML formulation slightly increased the number of DOFs. It is worth noting that the reduction in the number of DOFs greatly depends on the characteristic of the studied case, particularly in the dimensions of \(\Omega^{\mathrm{RD}}\). However, although effective for some applications, paraxial boundary conditions are not ideal for absorbing surface waves. Therefore, they are not recommended in problems where surface waves are predominant as the solved in this investigation.

The proposed hybrid PML formulation is susceptible to the same issues as other PML formulations and may experience instability over time under certain circumstances. However, a significant advantage of the proposed method is that it allows for the redefinition of scaling and attenuation functions using stretching functions with superior absorbing properties, such as M-PML, without modifying the underlying PDEs. Additionally, the time-integration scheme employed in the simulations did not account for numerical damping, making the problem conditions even more demanding compared to other methods [18,19]. Only minor time instabilities were observed, which were successfully resolved using M-PML.

In terms of discretization, the element size was selected as large as possible to reduce the number of DOFs. Additionally, only P2-P2-P1 finite element triplets were utilized to represent the solutions to the problems. Despite the demanding finite element discretization, the approach proposed in this article considerably reduced the computational effort required to solve the problem, unlike in previous studies. Consequently, the hybrid PML formulation has demonstrated robustness for demanding discretization and physical media conditions.

The following steps regarding this investigation are: (1) to extend the hybrid PML formulation to the 3D case to simulate, for instance, more realistic scenarios in complex seismic geophysical or biomedical applications. (2) To consider variable and discontinuous porosities in space which would introduce discontinuities in the relative fluid displacement. (3) To apply the 2D and 3D formulations to solve inverse problems in porous media. And (4), to simulate and understand the propagation of poroelastic waves in the human body, which is particularly interesting in the fields of Magnetic Resonance Imaging and Ultrasound and is closely related to the elastography problem [43,44].

---

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

---

## Data availability

The code to reproduce the results will be freely available after publication at https://github.com/hmella/POROUS-HYBRID-PML.

---

## Acknowledgments

HM and JM acknowledge the financial support given by ANID through the projects Fondecyt Postdoctorado 2022 #3220266 and Fondecyt Regular 2023 #1230864, respectively. ES was partially funded by a grant from the Research Center for Integrated Disaster Risk Management CIGIDEN Project 1522A0005 FONDAP 2022.

---

## Appendix A. Error plot in the traces of experiment 4

See Fig. A.1.

---

## Appendix B. Supplementary data

Supplementary material related to this article can be found online at https://doi.org/10.1016/j.cma.2023.116386.

---

## References

[1] O.C. Zienkiewicz, C.T. Chang, P. Bettess, Drained, undrained, consolidating and dynamic behaviour assumptions in soils, Geotechnique 30 (4) (1980) 385-395.  
[2] O.C. Zienkiewicz, T. Shiomi, Dynamic behaviour of saturated porous media; The generalized Biot formulation and its numerical solution, Int. J. Numer. Anal. Methods Geomech. 8 (1) (1984) 71-96.  
[3] J. Gwinner, E. Stephan, Advanced Boundary Element Methods: Treatment of Boundary Value, Transmission and Contact Problems, Springer, 2018.  
[4] E. Kausel, Local transmitting boundaries, J. Eng. Mech. 114 (6) (1988) 1011-1027.  
[5] M.A. Biot, Theory of propagation of elastic waves in a fluid-Saturated porous solid. I. Low-frequency range, J. Acoust. Soc. Am. 28 (2) (1956) 168-178.  
[6] N.F. Dudley Ward, T. Lähivaara, S. Eveson, A discontinuous Galerkin method for poroelastic wave propagation: The two-dimensional case, J. Comput. Phys. 350 (2017) 690-727.  
[7] J.P. Berenger, A perfectly matched layer for the absorption of electromagnetic waves, J. Comput. Phys. 114 (2) (1994) 185-200.  
[8] Q. Qi, T.L. Geers, Evaluation of the perfectly matched layer for computational acoustics, J. Comput. Phys. 139 (1) (1998) 166-183.  
[9] W. Chew, Q. Liu, Perfectly matched layers for elastodynamics: a new absorbing boundary condition, J. Comput. Acoust. 04 (04) (1996) 341-359.  
[10] Y.Q. Zeng, J.Q. He, Q.H. Liu, The application of the perfectly matched layer in numerical modeling of wave propagation in poroelastic media, Geophysics 66 (4) (2001) 1258-1266.  
[11] T. Wang, X. Tang, Finite-difference modeling of elastic wave propagation: A nonsplitting perfectly matched layer approach, Geophysics 68 (5) (2003) 1749-1755.  
[12] U. Basu, A.K. Chopra, Perfectly matched layers for transient elastodynamics of unbounded domains, Internat. J. Numer. Methods Engrg. 59 (8) (2004) 1039-1074.  
[13] S. Kucukcoban, L.F. Kallivokas, Mixed perfectly matched layers for direct transient analysis in 2D elastic heterogeneous media, Comput. Methods Appl. Mech. Engrg. 200 (1) (2011) 57-76.  
[14] S. Kucukcoban, L.F. Kallivokas, A symmetric hybrid formulation for transient wave simulations in PML-truncated heterogeneous media, Wave Motion 50 (1) (2013) 57-79.  
[15] F.X. Zhou, Q. Ma, B.B. Gao, Efficient unsplit perfectly matched layers for finite-element time-domain modeling of elastodynamics, J. Eng. Mech. 142 (11) (2016) 1-12.  
[16] R. Song, J. Ma, K. Wang, The application of the nonsplitting perfectly matched layer in numerical modeling of wave propagation in poroelastic media, Appl. Geophys. 2 (4) (2005) 216-222.  
[17] R. Martin, D. Komatitsch, A. Ezziani, An unsplit convolutional perfectly matched layer improved at grazing incidence for seismic wave propagation in poroelastic media, Geophysics 73 (4) (2008) T51-T61.  
[18] Y. He, T. Chen, J. Gao, Unsplit perfectly matched layer absorbing boundary conditions for second-order poroelastic wave equations, Wave Motion 89 (2019) 116-130.  
[19] Y. He, T. Chen, J. Gao, Perfectly matched absorbing layer for modelling transient wave propagation in heterogeneous poroelastic media, J. Geophys. Eng. (2019).  
[20] D. Komatitsch, R. Martin, An unsplit convolutional perfectly matched layer improved at grazing incidence for the seismic wave equation, Geophysics (2007).  
[21] F.H. Drossaert, A. Giannopoulos, Complex frequency shifted convolution PML for FDTD modelling of elastic waves, Wave Motion 44 (7-8) (2007) 593-604.  
[22] Y.Q. Zeng, Q.H. Liu, A staggered-grid finite-difference method with perfectly matched layers for poroelastic wave equations, J. Acoust. Soc. Am. 109 (6) (2001) 2571-2580.  
[23] A. Ezziani, Modelisation mathematique et numerique de la propagation d'ondes dans les milieux viscoelastiques et poroelastiques (Ph.D. thesis), Paris 9, 2005.  
[24] A. Fathi, B. Poursartip, L.F. Kallivokas, Time-domain hybrid formulations for wave simulations in three-dimensional PML-truncated heterogeneous media, Internat. J. Numer. Methods Engrg. 101 (3) (2015) 165-198.  
[25] M.A. Biot, Generalized theory of acoustic propagation in porous dissipative media, J. Acoust. Soc. Am. 34 (9A) (1962) 1254-1264.  
[26] C. Morency, J. Tromp, Spectral-element simulations of wave propagation in porous media, Geophys. J. Int. 175 (1) (2008) 301-345.  
[27] M. Kuzuoglu, R. Mittra, Frequency dependence of the constitutive parameters of causal perfectly matched anisotropic absorbers, IEEE Microw. Guid. Wave Lett. 6 (12) (1996) 447-449.  
[28] D. Correia, Jian-Ming Jin, On the development of a higher-order PML, IEEE Trans. Antennas and Propagation 53 (12) (2005) 4157-4163.  
[29] K.C. Meza-Fajardo, A.S. Papageorgiou, A nonconventional, split-field, perfectly matched layer for wave propagation in isotropic and anisotropic elastic media: Stability analysis, Bull. Seismol. Soc. Am. 98 (4) (2008) 1811-1836.  
[30] S. François, H. Goh, L.F. Kallivokas, Non-convolutional second-order complex-frequency-shifted perfectly matched layers for transient elastic wave propagation, Comput. Methods Appl. Mech. Engrg. 377 (2021) 113704.  
[31] F. Collino, C. Tsogka, Application of the perfectly matched absorbing layer model to the linear elastodynamic problem in anisotropic heterogeneous media, Geophysics 66 (1) (2001) 294-307.  
[32] T. Akiyoshi, K. Fuchida, H.L. Fang, Absorbing boundary conditions for dynamic analysis of fluid-saturated porous media, Soil Dyn. Earthq. Eng. 13 (6) (1994) 387-397.  
[33] N.M. Newmark, A Method of computation for structural dynamics, J. Eng. Mech. Div. 85 (3) (1959) 67-94.  
[34] M.S. Alnaes, J. Blechta, J. Hake, A. Johansson, B. Kehlet, A. Logg, C. Richardson, J. Ring, M.E. Rognes, G.N. Wells, The FEniCS project Version 1.5, Arch. Numer. Softw. 3 (100) (2015) 9-23.  
[35] A. Logg, K.A. Mardal, G. Wells (Eds.), Automated Solution of Differential Equations by the Finite Element Method, Springer, 2012.  
[36] F. Ballarin, Multiphenics - Easy prototyping of multiphysics problems in FEniCS, 2016.  
[37] C. Geuzaine, J.-F. Remacle, Gmsh: A 3-D finite element mesh generator with built-in pre- and post-processing facilities, Internat. J. Numer. Methods Engrg. 79 (11) (2009) 1309-1331.  
[38] S. Balay, S. Abhyankar, M.F. Adams, J. Brown, P. Brune, K. Buschelman, L. Dalcin, A. Dener, V. Eijkhout, W.D. Gropp, D. Kaushik, M.G. Knepley, D.A. May, L.C. McInnes, R.T. Mills, T. Munson, K. Rupp, P. Sanan, B.F. Smith, S. Zampini, H. Zhang, H. Zhang, PETSc web page, 2020.  
[39] S. Balay, W.D. Gropp, L.C. McInnes, B.F. Smith, Efficient management of parallelism in object oriented numerical software libraries, in: E. Arge, A.M. Bruaset, H.P. Langtangen (Eds.), Modern Software Tools in Scientific Computing, Birkhauser Press, 1997, pp. 163-202.  
[40] S. Balay, S. Abhyankar, M.F. Adams, J. Brown, P. Brune, K. Buschelman, L. Dalcin, A. Dener, V. Eijkhout, W.D. Gropp, D. Kaushik, M.G. Knepley, D.A. May, L.C. McInnes, R.T. Mills, T. Munson, K. Rupp, P. Sanan, B.F. Smith, S. Zampini, H. Zhang, H. Zhang, PETSc Users Manual, Technical Report ANL-95/11 - Revision 3.13, Argonne National Laboratory, 2020.  
[41] D. Hysom, A. Pothen, A scalable parallel algorithm for incomplete factor preconditioning, SIAM J. Sci. Comput. 22 (6) (2001) 2194-2215.  
[42] R. Martin, D. Komatitsch, S. Gedney, E. Bruthiaux, A high-order time and space formulation of the unsplit perfectly matched layer for the seismic wave equation using auxiliary differential equations (ADE-PML), CMES Comput. Model. Eng. Sci. 56 (1) (2010) 17-42.  
[43] M.M. Doyle, Model-based elastography: a survey of approaches to the inverse elasticity problem, Phys. Med. Biol. 57 (3) (2012) R35.  
[44] A. Kolipaka, K.P. McGee, A. Manduca, A.J. Romano, K.J. Glaser, P.A. Araoz, R.L. Ehman, Magnetic resonance elastography: Inversions in bounded media, Magn. Reson. Med. 62 (6) (2009) 1533-1542.