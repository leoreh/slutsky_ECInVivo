# -*- coding: utf-8 -*-
"""
Brian 2 Implementation of Model 2: A Spiking Network Model of 
Input-Current Homeostasis based on the provided strategic report.

This script simulates a recurrent network of Adaptive Exponential 
Integrate-and-Fire (AdEx) neurons. It can be run in two modes:
1.  'I-regulation': The novel hypothesis where neurons adjust their intrinsic
    properties to stabilize their average total synaptic input current.
2.  'F-regulation': The canonical model where neurons stabilize their
    output firing rate.

The simulation includes a perturbation phase (reduced external drive) to
test the homeostatic mechanisms.
"""

from brian2 import *
import numpy as np
import matplotlib.pyplot as plt

def run_homeostasis_simulation(regulation_type='I', seed_value=42):
    """
    Runs the full network simulation for either I-regulation or F-regulation.

    Args:
        regulation_type (str): 'I' for Input-centric regulation or 'F' for 
                               Output-centric (Firing rate) regulation.
        seed_value (int): Seed for reproducibility.
    """
    # Set the seed for reproducibility
    seed(seed_value)
    np.random.seed(seed_value)
    
    # Use numpy code generation for simplicity
    prefs.codegen.target = 'numpy'
    
    # =========================================================================
    # 1. Model Parameters (as per report and standard AdEx literature)
    # =========================================================================
    
    # --- Network Parameters
    N_E = 800          # Number of excitatory neurons
    N_I = 200          # Number of inhibitory neurons
    N = N_E + N_I      # Total number of neurons
    p_connect = 0.1    # Connection probability
    
    # --- Neuron Parameters (AdEx Model)
    C = 281 * pF              # Membrane capacitance
    gL = 30 * nS              # Leak conductance
    EL = -70.6 * mV           # Leak reversal potential (resting potential)
    VT = -50.4 * mV           # Spike threshold
    DeltaT = 2 * mV           # Slope factor
    V_peak = 20 * mV          # Spike cutoff
    Vr = EL                   # Reset potential
    
    # Adaptation parameters
    tau_w = 144 * ms          # Adaptation time constant
    b = 0.0805 * nA           # Spike-triggered adaptation increment
    a_baseline = 4 * nS       # Baseline subthreshold adaptation parameter
    
    # --- Synaptic Parameters
    tau_exc = 5 * ms          # Excitatory synapse time constant
    tau_inh = 10 * ms         # Inhibitory synapse time constant
    E_exc = 0 * mV            # Excitatory reversal potential
    E_inh = -75 * mV          # Inhibitory reversal potential
    
    # Synaptic weights
    w_ee = 0.5 * nS
    w_ei = 0.5 * nS
    w_ie = 2.0 * nS
    w_ii = 2.0 * nS
    w_ext = 1.5 * nS # External input weight

    # --- Homeostasis Parameters (as per report's description)
    # Timescale for averaging Isyn or Fout
    tau_homeo = 100 * ms      
    # Timescale for intrinsic plasticity (slow change in parameter 'a')
    tau_a = 2000 * ms 
    eta = 0.5                 # Dimensionless learning rate for parameter 'a'
    
    # Set-points
    I_set = 25 * pA           # Target input current for I-regulation
    F_set = 5 * Hz            # Target firing rate for F-regulation

    # =========================================================================
    # 2. Model Equations
    # =========================================================================

    eqs_base = '''
    dV/dt = (-gL*(V-EL) + gL*DeltaT*exp((V-VT)/DeltaT) - w + I_syn_exc + I_syn_inh) / C : volt (unless refractory)
    dw/dt = (a*(V-EL) - w) / tau_w : amp

    I_syn_exc = g_exc * (E_exc - V) : amp
    I_syn_inh = g_inh * (E_inh - V) : amp

    dg_exc/dt = -g_exc / tau_exc : siemens
    dg_inh/dt = -g_inh / tau_inh : siemens
    '''

    I_total_syn = 'g_exc*(E_exc-V) + g_inh*(E_inh-V)'

    if regulation_type == 'I':
        eqs_homeo = '''
        davg_Isyn/dt = (({i_total}) - avg_Isyn) / tau_homeo : amp
        eps_I = I_set - avg_Isyn : amp
        da/dt = -eta * a_baseline * (eps_I / I_set) / tau_a : siemens
        '''.format(i_total=I_total_syn)
        full_eqs = eqs_base + eqs_homeo

    elif regulation_type == 'F':
        eqs_homeo = '''
        dF_out/dt = -F_out / tau_homeo : Hz 
        eps_F = F_set - F_out : Hz
        da/dt = -eta * a_baseline * (eps_F / F_set) / tau_a : siemens
        '''
        full_eqs = eqs_base + eqs_homeo

    else:
        raise ValueError("regulation_type must be 'I' or 'F'")

    # =========================================================================
    # 3. Network Construction
    # =========================================================================
    
    # Neuron Groups
    neurons = NeuronGroup(N, full_eqs,
                          threshold='V > V_peak',
                          reset='V = Vr; w += b',
                          refractory=5*ms,
                          method='exponential_euler')
    
    # Split into excitatory and inhibitory populations
    P_E = neurons[:N_E]
    P_I = neurons[N_E:]

    # --- Initialize Neuron State Variables
    neurons.V = EL
    neurons.w = 0 * nA
    neurons.g_exc = 0 * nS
    neurons.g_inh = 0 * nS
    neurons.a = a_baseline

    if regulation_type == 'I':
        neurons.avg_Isyn = I_set # Start near the set-point
    elif regulation_type == 'F':
        neurons.F_out = F_set # Start near the set-point
        # For F-regulation, we need to increment F_out on each spike
        # This is done via a self-synapse.
        syn_F_update = Synapses(neurons, neurons, on_pre='F_out_post += 1.0/tau_homeo', delay=0*ms)
        syn_F_update.connect(j='i')

    # --- External Input (Poisson process)
    # The rate of this input will be changed to simulate perturbation
    ext_input_rate = 10 * Hz
    P_ext = PoissonGroup(N_E, rates=ext_input_rate)
    
    # --- Synapses
    S_ee = Synapses(P_E, P_E, on_pre='g_exc_post += w_ee')
    S_ei = Synapses(P_E, P_I, on_pre='g_exc_post += w_ei')
    S_ie = Synapses(P_I, P_E, on_pre='g_inh_post += w_ie')
    S_ii = Synapses(P_I, P_I, on_pre='g_inh_post += w_ii')
    S_ext = Synapses(P_ext, P_E, on_pre='g_exc_post += w_ext')

    S_ee.connect(p=p_connect)
    S_ei.connect(p=p_connect)
    S_ie.connect(p=p_connect)
    S_ii.connect(p=p_connect)
    S_ext.connect(j='i') # One-to-one connection

    # =========================================================================
    # 4. Monitors
    # =========================================================================
    
    spike_mon_E = SpikeMonitor(P_E)
    spike_mon_I = SpikeMonitor(P_I)
    
    pop_rate_mon_E = PopulationRateMonitor(P_E)
    
    # Monitor a few sample neurons
    state_mon = StateMonitor(P_E, ['V', 'a'], record=[0, 1, 2])
    
    if regulation_type == 'I':
        homeo_state_mon = StateMonitor(P_E, 'avg_Isyn', record=[0, 1, 2])
    elif regulation_type == 'F':
        homeo_state_mon = StateMonitor(P_E, 'F_out', record=[0, 1, 2])
        
    # =========================================================================
    # 5. Simulation Protocol
    # =========================================================================
    
    # Define simulation phases
    T_baseline = 5 * second
    T_deprivation = 20 * second
    T_recovery = 5 * second
    
    print(f"--- Running {regulation_type}-Regulation Simulation ---")
    
    # Phase 1: Baseline activity
    print(f"1. Running baseline for {T_baseline}...")
    run(T_baseline, report='text')
    
    # Phase 2: Sensory Deprivation (reduce external input)
    print(f"2. Simulating deprivation for {T_deprivation}...")
    P_ext.rates = 2 * Hz
    run(T_deprivation, report='text')
    
    # Phase 3: Recovery (restore external input)
    print(f"3. Simulating recovery for {T_recovery}...")
    P_ext.rates = ext_input_rate
    run(T_recovery, report='text')

    print("--- Simulation Complete ---")

    # =========================================================================
    # 6. Plotting
    # =========================================================================
    
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(4, 1, height_ratios=[2, 1, 1, 1])
    
    # --- Plot 1: Raster Plot
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(spike_mon_E.t/second, spike_mon_E.i, '.k', markersize=2, label='Excitatory')
    # ax1.plot(spike_mon_I.t/second, spike_mon_I.i + N_E, '.r', markersize=2, label='Inhibitory')
    ax1.set(ylabel='Neuron Index', title=f'Network Activity ({regulation_type}-Regulation)', xticklabels=[])
    ax1.legend(loc='upper right')
    
    # --- Plot 2: Population Firing Rate
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    # Smooth the rate for better visualization
    ax2.plot(pop_rate_mon_E.t/second, pop_rate_mon_E.smooth_rate(width=50*ms)/Hz, 'k-')
    ax2.set(ylabel='Excitatory Rate (Hz)', xticklabels=[])
    
    # --- Plot 3: Key Homeostatic Variable (Isyn or Fout)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    if regulation_type == 'I':
        ax3.plot(homeo_state_mon.t/second, homeo_state_mon.avg_Isyn.T/pA)
        ax3.axhline(I_set/pA, color='r', linestyle='--', label=f'I_set = {I_set/pA:.1f} pA')
        ax3.set(ylabel='Avg. Input Current (pA)')
    elif regulation_type == 'F':
        ax3.plot(homeo_state_mon.t/second, homeo_state_mon.F_out.T/Hz)
        ax3.axhline(F_set/Hz, color='r', linestyle='--', label=f'F_set = {F_set/Hz:.1f} Hz')
        ax3.set(ylabel='Avg. Firing Rate (Hz)')
    ax3.legend(loc='upper right')
    ax3.set(xticklabels=[])
    
    # --- Plot 4: Intrinsic Parameter 'a'
    ax4 = fig.add_subplot(gs[3], sharex=ax1)
    ax4.plot(state_mon.t/second, state_mon.a.T/nS)
    ax4.axhline(a_baseline/nS, color='r', linestyle='--', label=f'a_baseline = {a_baseline/nS:.1f} nS')
    ax4.set(xlabel='Time (s)', ylabel='Parameter a (nS)')
    ax4.legend(loc='upper right')
    
    # Add vertical lines to delineate phases
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvline(T_baseline/second, color='gray', linestyle=':')
        ax.axvline((T_baseline + T_deprivation)/second, color='gray', linestyle=':')
        ax.grid(True, linestyle='--', alpha=0.5)

    fig.suptitle(f"Model 2: Spiking Network with {regulation_type}-Regulation", fontsize=16, y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()


if __name__ == '__main__':
    # Run the simulation for the I-regulation hypothesis
    run_homeostasis_simulation(regulation_type='I')
    
    # Run the simulation for the F-regulation control model
    run_homeostasis_simulation(regulation_type='F')