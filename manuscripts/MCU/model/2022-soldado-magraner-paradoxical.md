---
year: 2022
journal: "Proceedings of the National Academy of Sciences"
title: "Paradoxical self-sustained dynamics emerge from orchestrated excitatory and inhibitory homeostatic plasticity rules"
authors:
  - "Saray Soldado-Magraner"
  - "Michael J. Seay"
  - "Rodrigo Laje"
  - "Dean V. Buonomano"
keywords:
doi: "https://doi.org/10.1073/pnas.2200621119"
zotero: "zotero://select/library/items/S49HXCC7"
---

## Abstract
Self-sustained neural activity maintained through local recurrent connections is of fundamental importance to cortical function. Converging theoretical and experimental evidence indicates that cortical circuits generating self-sustained dynamics operate in an inhibition-stabilized regime. Theoretical work has established that four sets of weights (WE←E, WE←I, WI←E, and WI←I) must obey specific relationships to produce inhibition-stabilized dynamics, but it is not known how the brain can appropriately set the values of all four weight classes in an unsupervised manner to be in the inhibition-stabilized regime. We prove that standard homeostatic plasticity rules are generally unable to generate inhibition-stabilized dynamics and that their instability is caused by a signature property of inhibition-stabilized networks: the paradoxical effect. In contrast, we show that a family of “cross-homeostatic” rules overcome the paradoxical effect and robustly lead to the emergence of stable dynamics. This work provides a model of how—beginning from a silent network—self-sustained inhibition-stabilized dynamics can emerge from learning rules governing all four synaptic weight classes in an orchestrated manner.

## My Summary
%% begin notes %%

### Findings

* Standard homeostatic plasticity rules, which adjust synapses based on a neuron's own activity error, are mathematically proven to be unstable in the paradoxical regime that defines Inhibition-Stabilized Networks (ISNs).
* This instability arises because an attempt to increase the activity of an under-active inhibitory neuron by strengthening its excitatory inputs causes a paradoxical decrease in its activity.
* A novel "cross-homeostatic" plasticity rule, where excitatory and inhibitory populations adjust their synapses based on each other's activity errors, is robustly stable and can guide a network to an ISN state.
* The simple cross-homeostatic rule only ensures that population-average firing rates converge to their setpoints, not the rates of individual neurons.
* A hybrid "two-term" rule, combining both standard and cross-homeostatic components, achieves stable convergence for both population-level averages and individual neuron firing rates.
* This two-term rule is shown to be effective across different model types, successfully guiding simple rate-based models, multiunit networks, and large-scale spiking neural networks from a silent state to a fully formed, self-sustained ISN.

### Interpretation

* The emergent dynamic properties of cortical circuits, such as the paradoxical effect, impose fundamental constraints on the types of learning rules the brain can use to self-organize.
* The brain likely employs counter-intuitive or "paradoxical" learning rules to form and maintain stable, high-performance recurrent networks.
* Cross-homeostatic plasticity is presented as a biologically plausible and computationally robust mechanism for the unsupervised development of cortical ISNs.


The authors propose that the counter-intuitive "cross-homeostatic" rule can be implemented through a plausible local mechanism within each neuron. They suggest that a neuron separates the functions of sensing and acting by using two different classes of receptors. A postsynaptic neuron, for example an excitatory cell, would use its slow-acting metabotropic receptors (like GABA-B) to integrate and effectively measure the average activity of its presynaptic inhibitory partners. This measured activity is then compared to an internal setpoint for how much inhibitory input that neuron should be receiving. If a mismatch or "error" is detected, the neuron triggers internal plasticity mechanisms that adjust the number of its fast-acting ionotropic receptors (like AMPA and GABA-A), thereby changing its synaptic weights to help correct the network-wide error. This elegant solution transforms the seemingly non-local rule into a local process where each neuron adjusts its excitability based on a continuous reading of its own presynaptic environment.


### Comments

#### Time Constants

In the differential equations that govern the firing rates of the excitatory and inhibitory populations, the time constants ($τ_E$) and ($τ_I$) represent the characteristic time it takes for the firing rate of a population to change in response to a change in its input. A smaller time constant implies a faster response. 
In the provided paper, as is common in such models, the time constant for the inhibitory population (2 ms) is significantly shorter than that of the excitatory population (10 ms).
This difference is not arbitrary but is rooted in the underlying biophysics of the neurons and synapses they are meant to represent:

Synaptic Kinetics: The primary reason for the different time constants lies in the kinetics of the neurotransmitter receptors. Excitatory synapses typically use glutamate, which acts on receptors like AMPA and NMDA. While AMPA receptors are fast, the persistent activity in recurrent networks often involves the slower NMDA receptors. Inhibitory synapses, on the other hand, primarily use GABA, which acts on GABA-A receptors. These receptors mediate a very fast chloride current, causing a rapid change in the postsynaptic neuron's membrane potential.

Membrane Properties: The intrinsic properties of the neurons themselves also contribute. Inhibitory interneurons, particularly parvalbumin-positive (PV) basket cells which are crucial for network stabilization, have faster membrane time constants than excitatory pyramidal neurons. This means they can integrate and respond to synaptic inputs more quickly.

The faster response of the inhibitory population is a critical feature for maintaining stability in a network with strong recurrent excitation. If excitatory neurons excite each other, their activity can grow uncontrollably, leading to epileptiform activity. Fast, strong inhibition is necessary to counterbalance this positive feedback and keep the network activity stable, a regime known as an inhibition-stabilized network. The different time constants in the model are a crucial element in capturing this fundamental aspect of cortical dynamics.


#### Coupled Dynamic Systems

The behavior of the network is described as a system of two coupled dynamic subsystems, each operating on a different timescale. This framework is essential for analyzing how the network can self-organize through the interaction of neural activity and synaptic learning.

The Fast Neural Subsystem (Neural Activity) governs the moment-to-moment firing rates of the excitatory and inhibitory neuron populations. It operates on a fast timescale of milliseconds, determined by the membrane time constants ($\tau$) of the neurons.

The dynamics are described by a set of coupled differential equations:

$$
\tau_E \frac{dE}{dt} = -E + f_E(W_{EE}E - W_{EI}I)
$$

$$
\tau_I \frac{dI}{dt} = -I + f_I(W_{IE}E - W_{II}I)
$$

* $E, I$: The average firing rates of the excitatory and inhibitory populations, respectively.
* $\tau_E, \tau_I$: The time constants for each population, defining how quickly their activity changes.
* $-E, -I$: A decay term, ensuring that without input, activity returns to zero.
* $W_{XY}$: The average synaptic weight from population $Y$ to population $X$. For example, $W_{EI}$ is the weight from the inhibitory population to the excitatory population.
* $f_E(...), f_I(...)$: Activation functions that convert the total synaptic input into an output firing rate.

For any fixed set of synaptic weights ($W$), this fast subsystem will evolve and settle into a stable state, known as a fixed point, where $dE/dt = 0$ and $dI/dt = 0$. This represents a self-consistent state where the network's activity is balanced by its own recurrent inputs.

The Slow Synaptic Plasticity Subsystem (Learning) describes how the synaptic weights ($W$) change over time. It operates on a much slower timescale, representing the process of learning and adaptation. The changes are not instantaneous but are driven by the average neural activity computed over longer periods.

The dynamics are described by a set of learning rules. For the "cross-homeostatic" rule, the equations are:
$$
\Delta W_{EE} \propto \alpha_E E (I_{set} - I)
$$
$$
\Delta W_{EI} \propto -\alpha_E I (I_{set} - I)
$$
$$
\Delta W_{IE} \propto -\alpha_I E (E_{set} - E)
$$
$$
\Delta W_{II} \propto \alpha_I I (E_{set} - E)
$$

* $\Delta W_{XY}$: The change in the synaptic weight.
* $\alpha_E, \alpha_I$: The learning rates for plasticity onto the excitatory and inhibitory populations.
* $E_{set}, I_{set}$: The target firing rates, or "setpoints," for each population.
* $(E_{set} - E), (I_{set} - I)$: The "error" terms that drive the plasticity. The goal of the learning is to reduce these errors to zero.

The two subsystems are coupled because the state of each directly influences the evolution of the other, creating a continuous feedback loop:

1. The current synaptic weights ($W$) from the slow subsystem define the parameters of the fast neural subsystem.
2. The fast neural subsystem quickly converges to a stable firing rate fixed point ($E$, $I$) determined by those weights.
3. These resulting firing rates ($E$$, $I$) are then used by the slow plasticity subsystem to calculate the error and determine how the weights should be updated ($\Delta W$).
4. The updated weights create a new set of parameters for the fast subsystem, which then finds a new fixed point.

This cycle repeats, with the network's activity continuously shaping its own connectivity, and the connectivity in turn shaping the activity. The analytical power of this approach, known as a **quasi-steady-state approximation**, comes from the separation of timescales. One can assume the fast system is always at its stable fixed point and then analyze how that fixed point slowly moves as the weights evolve, allowing for a prediction of whether the entire coupled system will converge to the desired target state.


> [!NOTE]
> Due to spatial and temporal integration, mitoCa2+ is a better realization of the slow homeostatic kinetics, relative to cytoCa2+. 

#### Loss Function

The paper defines a loss function:
$$L(E, I) = \frac{1}{2}(E_{set}-E)^{2}+\frac{1}{2}(I_{set}-I)^{2}$$

And a gradient descent function in relation to synaptic weights:
$$
\Delta \mathbf{W} \propto -\nabla_{\mathbf{W}} L \quad \text{or} \quad \frac{d\mathbf{W}}{dt} = -\alpha \nabla_{\mathbf{W}} L
$$


However, our data consistently shows that intrinsic plasticity is  involved in FRH, even in the absence of synaptic plasticity.

Hence, applying this principle to intrinsic plasticity (parameters like  gain $g_E$, threshold $\Theta_E$, or adaptation currents), we can write:
$$\frac{dg_{E}}{dt} = -\alpha \frac{\partial L}{\partial g_{E}}$$


See [[LLM - Analytical Derivation of E-Only Homeostasis]] for a similar derivation of learning rules if the loss function only contains the excitatory population. In sum, the derivation seems to work, the four learning rules of synaptic weights will be homeostatic with respect to E. However, two main issues:
1. The learning rules seem biologically implausible (eg, the update for $\Delta W_{EE}$ requires knowledge of $W_{II}$).
2. The functional failure is that the synaptic weights (slow dynamics) will reach stability but the neural (fast) dynamics of E and I will not. 

#### Variability Emerges

Under the cross-homeostatic rule, error signals are determined by the mean activity of all presynaptic partners. Because the error signal is averaged, it is "blind" to the individual variations that constitute it. As a result, when this rule is simulated, the system naturally settles into a state with a non-zero standard deviation. The variability is a result of a rule that only targets the mean.


> [! Cross-Homeostatic Rule:]
> This rule adjusts synaptic weights based on the average error of the presynaptic population. The plasticity of a neuron is driven by the needs of its input partners, not its own activity level. In these equations, `i` is the postsynaptic neuron, `j` is the presynaptic neuron, and the summation over `k` represents the average over the relevant presynaptic population.
> $$
> \begin{aligned}
> \Delta W_{EE}^{ij} &= +\alpha_E E_j \sum_{k=1}^{N_I} (I_{set} - I_k)/N_I \\
> \end{aligned}
> $$



The authors propose a solution to this "problem" in the from of a two-term rule, which adds a local component that is not based on an average but on the specific error of each neuron $i$ ($E_{set} - E_i$). This term actively works to reduce the error of every single neuron to zero. By doing so, it forces the standard deviation of the population's firing rates to also approach zero. The lack of variability is a result of a rule that targets each individual neuron.


> [!Two-Term Cross-Homeostatic Rule]
> This rule is a hybrid that combines the cross-homeostatic rule with a standard homeostatic rule. The first term represents the "personal error" of the postsynaptic neuron `i` (a homeostatic influence), while the second term represents the average error of its presynaptic partners `k` (the cross-homeostatic influence). This combination allows both the population average and individual neuron firing rates to converge to their setpoints.
> $$
> \begin{aligned}
> \Delta W_{EE}^{ij} &= +\alpha E_j(E_{set} - E_i) + \alpha E_j \sum_{k=1}^{N_I} (I_{set} - I_k)/N_I \\
> \end{aligned}
> $$

Assuming FRH is implemented at the network, and not single-cell level, the cross-homeostatic rule is not only sufficient, but also superior as it naturally explains how the variability of the FR distribution emerges. 


#### Perturbation

The network, by its very design, is bistable. It has two fundamentally different stable states:
- The "Down-state": A silent state where E=0 and I=0. This is a trivial but very stable fixed point. If the network starts with zero activity and receives no input, its rate equations dictate that it will stay at zero activity forever. 
- The "Up-state": The desired state of self-sustained, persistent activity where E and I have non-zero firing rates.

Since linear stability analysis assumes small changes, it is only a local property. 
A successful learning rule must do more than just make the target fixed point stable; it must create a large and robust basin of attraction around it.

Hence they numerically simulated a large perturbation (or "kick"). 
The kick is an explicit test of this basin. By starting the system in the Down-state (which is very far from the Up-state in the system's state space) and applying a transient push, the researchers are checking if this push is sufficient to land the system somewhere within the Up-state's basin of attraction. If, after the kick, the network reliably converges to the Up-state, it demonstrates that the basin is large enough to "catch" the system. A learning rule that creates only a tiny, shallow basin of attraction would fail this test; the kick would cause the system to overshoot it, and it would fall back to the more dominant Down-state.

his experiment is detailed in Supplementary Figure S6. The goal was to determine if the learning rule could maintain overall network stability while preserving specific, learned changes in synaptic strength, akin to a memory.
The test was performed as follows:
- First, the multi-unit network was trained using the two-term rule until the firing rates of the neurons stabilized at their target setpoints. At this point, the network is in a balanced, stable state.
- At Trial 200, a large, structural perturbation was introduced. The researchers simulated a "Hebbian-like" event by artificially strengthening the recurrent connections among a small subgroup of neurons (the first 10 excitatory neurons). All of their recurrent weights were increased by a constant factor. This can be thought of as imprinting a "memory" or a cell assembly into the network's connectivity matrix.
- This abrupt change in synaptic weights immediately drove the network's activity away from its stable setpoints, as the artificially strengthened neurons began to fire more.
The results showed that the two-term cross-homeostatic learning rule "re-engaged" in response to this perturbation. Over subsequent trials, the rule successfully guided the average firing rates of the network back to the original target setpoints. Remarkably, it achieved this functional recovery while largely preserving the underlying structural change. The weight matrix after re-stabilization showed that the connections among the perturbed group of neurons remained stronger than their neighbors. This demonstrates that the learning rule is capable of maintaining network stability by adapting the surrounding weights to compensate for the perturbation, thereby allowing the network to both remain stable and retain learned information.

#### Asynchronous Irregular

The asynchronous irregular (AI) regime characterizes a specific pattern of neural activity observed in the awake cortex, where individual neurons fire without temporal correlation to the population and with highly variable inter-spike intervals. This state is understood to emerge from a balance between excitatory and inhibitory inputs. In contrast, the inhibition-stabilized network (ISN) describes the underlying circuit dynamics that can produce this balance. An ISN is defined by recurrent excitatory connections that are strong enough to render the excitatory subnetwork inherently unstable, requiring rapid and potent feedback inhibition for moment-to-moment stabilization.

The relationship between these two concepts is that of mechanism to outcome. The dynamic properties of an ISN provide a robust mechanism for generating and maintaining an AI state of activity. In an ISN, neurons are subject to a constant barrage of strong excitatory and inhibitory synaptic inputs. This balanced, high-conductance state causes the membrane potential of each neuron to fluctuate near its firing threshold, leading naturally to the irregular and asynchronous spike timing that defines the AI regime.

While closely linked, the terms are not interchangeable. An ISN refers to the functional architecture and stability properties of the circuit, whereas the AI regime describes the statistical pattern of its output. It is theoretically possible for an AI state to arise in a balanced network that is not an ISN (i.e., one with stable recurrent excitation). Conversely, a circuit with ISN dynamics may be transiently driven into a more synchronous or oscillatory activity pattern by a strong, coherent stimulus, even though its fundamental stabilizing mechanism remains unchanged.




%% end notes %%



%% Import Date: 2025-12-04T17:33:00.111+02:00 %%
