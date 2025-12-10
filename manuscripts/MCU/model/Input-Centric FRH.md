---
parent:
  - "[[Research]]"
Time: 2025-10-29T16:40:00
---

## Introduction

#### Steady-state, dynamic stability, and homeostasis

Mean-field theory provides self-consistency equations that define the fixed points of a system, which are states where the network's output is in equilibrium with its input. 
$$
E_{ss} = f_E \left( W_{EE}E_{ss} - W_{EI}I_{ss} \right)
$$
 However, the existence of a fixed point does not guarantee its stability, as even minor perturbations can cause the system to diverge. An Inhibition-Stabilized Network (ISN) represents a specific parameter regime that provides rapid, dynamic stability. In an ISN, the excitatory sub-network is inherently unstable, with recurrent connections strong enough to cause runaway activity. 
 $$
 W_{EE}g_E-1 > 0
 $$
This explosive excitation is controlled by powerful and fast inhibitory feedback, which is recruited by any increase in excitatory firing and acts to quickly restore equilibrium. This tense balance establishes a highly robust form of stability on a fast timescale. However, if the network is subjected to a sustained perturbation, such as a continuous external current, it will simply settle at a new stable fixed point that balances the new input. The ISN mechanism has no intrinsic goal or setpoint; it finds the nearest stable equilibrium for a given set of parameters and inputs. In contrast, homeostasis is the distinct, slower process that actively restores a target firing rate.

Current models of homeostasis face numerous challenges, particularly when applied to recurrent neural networks. A primary issue is the tendency for simple homeostatic rules to become unstable in strongly coupled circuits. For example, in an ISN, any attempt by an excitatory neuron to increase its activity also recruits stronger feedback inhibition. This feedback paradoxically decreases the neuron's firing rate, pushing it further from its target set point. Furthermore, many theoretical derivations that achieve stability do so by invoking mechanisms whose biological implementations are not yet clear. In particular, the process of sensing firing rates or determining its set points are often abstracted, leaving a gap in biological plausibility.


#### Empirical constraints for homeostasis models

Experimental findings suggest:
1. FRH is a network-level phenomenon, allowing individual firing rates to drift. 
2. Interneurons are exempt from FRH.
3. The intrinsic excitability of PV cells dynamically modulates the activity setpoint for the excitatory population.
4. Homeostatic compensation is primarily mediated by intrinsic plasticity, rather than synaptic plasticity.

#### FRH is not cell-autonomous 
The first neurons to respond to a homeostatic challenge, temporarily over-shoot their set point by >50%. This implies a driving force tuned to the network, rather then the single-cell. 
![[Thesis_Fig 16A.png|300]]

Firing rates drift both in vivo and in vitro, even during spontaneous activity. 
![[Ruggiero_Fig 1E.png|300]]

Stability analysis repeatedly shows statistical robustness at the network, but not single-cell level. 
![[Ruggiero_Fig 1B.png|300]]

#### pINT firing is not homeostatically regulated

During chronic perturbation with baclofen, pINT firing rates gradually increase. 
![[Thesis_Fig 14B.png|300]]

During chronic perturbation with Dlx-Gq, pINT firing rates are persistently reduced. 
![[Atsmon_Fig 2E.png|300]]

#### PV modulate activity set points

NMDAR inhibition reduces the set point of network activity. This is causally mediated by an increase in PV intrinsic excitability. 

![[Ruggiero_Fig 5 - PV set point.png|600]]

#### FRH involves intrinsic plasticity

FRH in CA1 in vivo following chronic excitation of interneurons involves intrinsic but not synaptic plasticity (examined in slices).
![[Atsmon_Fig 4B.png|300]]

FRH following chronic application of ketamine in eEF2K-KO hippocampal networks in vitro involves intrinsic but not synaptic plasticity.
![[Ruggiero_Fig 4E.png|300]]


## Analytic derivation

The following model provides a computational framework to understand how network-level homeostasis can emerge from cell-autonomous rules that are both stable and biologically plausible. The central hypothesis is that neurons regulate the statistics of their local synaptic inputs rather than their own output firing rate. 

The model considers a two-population network of interconnected excitatory (E) and inhibitory (I) neurons. The core of the model is built on two fundamental equilibrium conditions that operate on different timescales.

#### Fast network dynamics equilibrium
On a fast timescale, the firing rates of the neuronal populations reach a steady state where the synaptic input current received by a population is precisely the amount required to produce its output firing rate. For the excitatory population, this is described by:
$$
E_{ss} = f_E(W_{EE}E_{ss} - W_{EI}I_{ss})
$$
Where $E_{ss}$ and $I_{ss}$ are the steady-state firing rates, $W_{xy}$ represents the synaptic weight from population $y$ to $x$, and $f_E$ is the transfer function of the excitatory population.
The equilibrium can also be expressed using the inverse transfer function, $f_E^{-1}(E_{ss})$, which represents the specific amount of synaptic input required to produce the excitatory population firing rate. 
$$
f_E^{-1}(E_{ss}) = W_{EE}E_{ss} - W_{EI}I_{ss}
$$

#### Slow homeostatic plasticity equilibrium
On a slow timescale, a homeostatic mechanism adjusts each excitatory neuron's intrinsic properties to ensure that the input current it *requires* to sustain its firing rate ($f_E^{-1}(E_{ss})$) matches a fixed, internal setpoint ($h_{set,E}$):
$$
\tau_f \frac{df_E}{dt} \propto (h_{set,E} - f_E^{-1}(E))
$$
When the homeostatic process reaches its equilibrium, the error term is zero and the neuron's required input must equal its setpoint:
$$
f_E^{-1}(E_{ss}) = h_{set,E}
$$

For the entire system to be stable, both the fast network dynamics and the slow homeostatic plasticity must be at equilibrium simultaneously. Combining the two equilibriums reveals a self-consistency equation, whereby the actual synaptic input provided by the network must equal the neuron's internal setpoint.
$$
W_{EE}E_{ss} - W_{EI}I_{ss} = h_{set,E}
$$

This equation demonstrates that the network's firing rate is not an explicit target but an emergent property that arises from the interplay between recurrent connectivity and the local, input-centric homeostatic rule.

#### Network set point and inhibitory intrinsic excitability 

The following derivation demonstrates how the intrinsic properties of  inhibitory neurons can dynamically modulate the steady-state firing rate of the excitatory population ($E_{ss}$). The goal is to find an expression for $E_{ss}$ that depends on the parameters of the inhibitory transfer function ($g_I$ and $\theta_I$).

The transfer function of the inhibitory population, $f_I$, is modeled as a rectified-linear unit (ReLU), defined by:
$$
f_I(h) = g_I \cdot [h - \theta_I]_+
$$
Where $g_I$ is the gain and $\theta_I$ is the firing threshold of the inhibitory population.

The derivation assumes that both the fast firing rate dynamics and the slow homeostatic plasticity have reached a stable equilibrium.

It further assumes the network operates in a regime where the inhibitory population is active ($I_{ss} > 0$), meaning its input is above the firing threshold. This allows the transfer function to be treated as a linear function:
$$
I_{ss} = g_I (W_{IE}E_{ss} - W_{II}I_{ss} - \theta_I)
$$

First, we rearrange the linear equation for the inhibitory population's firing rate to express $I_{ss}$ as a function of $E_{ss}$. This defines how the steady-state inhibitory activity is determined by the excitatory activity.
$$
I_{ss} = g_I W_{IE}E_{ss} - g_I W_{II}I_{ss} - g_I \theta_I
$$
$$
I_{ss} + g_I W_{II}I_{ss} = g_I W_{IE}E_{ss} - g_I \theta_I
$$
$$
I_{ss}(1 + g_I W_{II}) = g_I (W_{IE}E_{ss} - \theta_I)
$$
$$
I_{ss} = \frac{g_I}{1 + g_I W_{II}} (W_{IE}E_{ss} - \theta_I)
$$

Next, we substitute this expression for $I_{ss}$ into the homeostatic equilibrium equation for the excitatory population. This yields a single equation with $E_{ss}$ as the only unknown variable.
$$
W_{EE}E_{ss} - W_{EI}I_{ss} = h_{set,E}
$$
$$
W_{EE}E_{ss} - W_{EI} \left( \frac{g_I}{1 + g_I W_{II}} (W_{IE}E_{ss} - \theta_I) \right) = h_{set,E}
$$

Finally, we perform algebraic manipulation to solve for $E_{ss}$.

$$
W_{EE}E_{ss} - \left( \frac{W_{EI} g_I W_{IE}}{1 + g_I W_{II}} \right) E_{ss} + \frac{W_{EI} g_I \theta_I}{1 + g_I W_{II}} = h_{set,E}
$$

$$
E_{ss} \left( W_{EE} - \frac{W_{EI} W_{IE} g_I}{1 + g_I W_{II}} \right) = h_{set,E} - \frac{W_{EI} g_I \theta_I}{1 + g_I W_{II}}
$$

This leads to the final analytical expression:
$$
E_{ss} = \frac{h_{set,E} - \frac{W_{EI} g_I \theta_I}{1 + g_I W_{II}}}{W_{EE} - \frac{W_{EI} W_{IE} g_I}{1 + g_I W_{II}}}
$$

This equation provides a formal and quantitative link between the intrinsic properties of inhibitory neurons ($g_I$, $\theta_I$) and the stable, emergent firing rate of the excitatory population ($E_{ss}$). It demonstrates that the network's excitatory setpoint is not a fixed parameter but is dynamically determined by the biophysical state of the inhibitory cells. For example, increasing the gain of the inhibitory neurons ($g_I$) strengthens the feedback inhibition within the network, which generally acts to lower $E_{ss}$. Hence, the inhibitory population dynamically calibrates the global, network-level firing rate that the excitatory population must achieve to satisfy its own fixed, local homeostatic requirements. This allows for global network stability and adaptability while maintaining a simple and robust homeostatic rule at the single-cell level.

## Comments

#### Inverse Transfer Function

A central feature of this model is the use of the inverse transfer function, $f^{-1}(E)$, to formalize the homeostatic learning rule. This function is critical because it distinguishes the synaptic input a neuron *requires* to sustain its current firing rate, from the input the neuron actually *receives* ($h_E = W_{EE}E - W_{EI}I$). This distinction provides an analytical bridge between the two different timescales at which the system finds balance: the fast equilibrium of network dynamics and the slow equilibrium of homeostatic plasticity. At the global steady state, both the input provided by the network ($h_E$) and the neuron's internal goal ($h_{set,E}$) must converge on this single value of required input, $f^{-1}(E)$. Therefore, they must be equal to each other. This formalism provides the mathematical justification for the model's core finding - that the network's recurrent input will organize to match the cell's internal setpoint.


## Biological Basis 

The input-centric model requires a physical mechanism that can measure synaptic input, compare it to an internal setpoint, and generate an error signal to drive plasticity. The proposed biological implementation for this computation is the mitochondrion.

The foundation of this system is a stable reference point, the model's homeostatic setpoint ($h_{set,E}$). This role is filled by the basal concentration of mitochondrial calcium ($[Ca²⁺]_{mito}$). Unlike cytosolic calcium, which fluctuates rapidly with synaptic activity, basal $[Ca²⁺]_{mito}$ is highly buffered and stable over long timescales, providing a reliable reference for the cell's baseline metabolic and energetic state.

The system's sensor is the Mitochondrial Calcium Uniporter (MCU). Due to its low affinity for calcium, the MCU functions as a high-pass filter, preferentially responding to the large, localized calcium transients generated by high-frequency burst firing while ignoring isolated spikes. This ensures that the homeostatic system responds to salient, energy-intensive patterns of activity rather than noisy fluctuations. The calcium influx through MCU thus serves as a direct biophysical proxy for the input current the neuron requires to produce its current output, a value corresponding to the model's inverse transfer function, $f_E^{-1}(E)$.

The error signal itself is generated by the comparison between the MCU-mediated calcium influx (the measurement of demand) and the basal mitochondrial calcium level (the setpoint). When a persistent mismatch occurs, downstream effectors within the mitochondrion generate a sustained corrective signal proportional to the calculated error, $(h_{set,E} - f_E^{-1}(E))$. This signal then engages cellular machinery that adjusts ion channel expression, reshaping the neuron's transfer function to close the homeostatic loop. One molecular candidate for this role is the sodium-calcium exchanger (NCLX), whose efflux rate is dependent on the matrix calcium concentration.

This mechanism reframes firing rate homeostasis as a process of metabolic accounting. In this view, the required input, $f_E^{-1}(E)$, represents the real-time energetic *cost* of sustaining a given firing rate. This cost, which integrates the total ionic load from both synaptic and somatic currents, serves as a direct proxy for the ATP consumption required to maintain ionic gradients. The cost is physically realized by the mitochondrial calcium signal, which is compared to a stable, cell-autonomous setpoint ($h_{set,E}$) that represents the neuron's sustainable energetic *budget*. This budget, embodied by the basal mitochondrial calcium concentration, is dictated by the cell's underlying metabolic capacity - its mitochondrial density and health, the efficiency of its respiratory chain, and its access to fuel sources. Homeostatic plasticity, therefore, is the process of adjusting a neuron's intrinsic excitability to ensure its operational energy demand does not persistently exceed its long-term budget, thereby linking network stability directly to the fundamental principles of cellular bioenergetics.


## Future Refinements

#### Inhibitory homeostasis (unresolved)

While pINT firing rates themselves do not appear to be homeostatically regulated, their intrinsic excitability profoundly modulates the network's activity setpoint. This suggests that PV cells engage in a distinct form of homeostasis, adjusting their properties to maintain a stable level of excitatory drive *from* the local network, rather than regulating their own output.

The biological sensor for this excitatory drive is the NMDAR on the PV interneurons. Due to their voltage-dependency and slow kinetics, NMDARs effectively integrate sustained, correlated activity from the excitatory population. Their activation serves as a reliable measure of the overall excitatory network state. The homeostatic learning rule for the PV cell's transfer function, $f_I$, can thus be formalized as a function of its NMDAR-mediated input:
$$
\Delta f_I \propto ({NMDA}_{set} - {NMDA}_{actual}) \approx \alpha_I ([r_E \cdot W_{IE}]_{set} - [r_E \cdot W_{IE}]_{actual})
$$

This creates a paradoxical effect. For example, a sustained increase in excitatory drive onto PV cells will decrease their intrinsic excitability. This is expected to amplify the initial perturbation, and might even interfere with the network's dynamical state (ISN).  

However, similar mechanisms of "cross-homeostasis," where the plasticity of one population ($I$) is tuned to the activity of a different population ($E$), has been shown to facilitate the emergence of ISN properties when applied to synaptic plasticity [[2022-soldado-magraner-paradoxical]]


#### Burstiness

Burstiness of excitatory units predicts the rate of their homeostatic response ($\alpha$). This could represent greater sensitivity of the mitochondria to $f^{-1}_E$. 

#### Predictions

#### Multi-unit modeling 

The input-centric homeostasis focuses on intrinsic excitability and assumes firing rates are preserved only at the population level. Perhaps, adding synaptic plasticity will allow neurons to preserve individual firing rate values (as in V1).

Moreover, the model attempts to bridge single-cell mechanisms with network stability, including the mean and variance of firing rates. Thus, it most likely requires simulations on heterogenous neurons (multi-unit or spiking models).
