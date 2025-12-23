To analyze the temporal dynamics of bursting across a long recording, we must transition from aggregate summary statistics to time-varying signals. This process involves converting discrete event data into continuous functions. For each of your unit-specific metrics, the strategy for temporal extension depends on whether the metric describes an **occurrence** (e.g., burst rate) or a **property** of the event itself (e.g., duration or intensity).

### Temporal Extension Strategies

The most robust way to study these dynamics is to categorize your metrics into two mathematical frameworks: **Density Estimation** and **Interpolated Property Tracking**.

#### 1. Density Estimation (For Count-Based Metrics)

Metrics like **Burst Rate** and **Burstiness (Bspks)** describe how often events occur or how much of the total activity they consume.

- **Goal:** To visualize the "velocity" or "prevalence" of bursting across the recording.
    
- **Method:** You treat each event as a point in time and apply a smoothing kernel to create a continuous rate function. For Burst Rate, you convolve the burst onset times. For Burstiness, you convolve the individual spikes that were identified as "in-burst."
    
- **Why Convolve?** While a moving average filter on a binary vector is mathematically similar, a Gaussian convolution provides a smoother, more biologically plausible estimation of the underlying state without the "boxcar" artifacts of simple binning.
    

#### 2. Interpolated Property Tracking (For Structural Metrics)

Metrics such as **nspks (spikes per burst)**, **brstDur (duration)**, **freq (intra-burst frequency)**, and **ibi (inter-burst interval)** describe the internal architecture of the bursts.

- **Goal:** To determine if the "shape" or "intensity" of the bursts changes during the perturbation.
    
- **Method:** You should not convolve these. Instead, you assign the metric value to the time of the burst's occurrence and use a moving median or a low-pass filter across the sequence of events.
    
- **Implementation Note:** If a unit stops bursting for a long period, the property values are undefined; you must decide whether to leave these gaps empty or interpolate between the last known burst and the next.
    

---

### Mapping Your Metrics to Dynamics

The following table outlines how each of your existing metrics can be extended into a temporal timecourse.

|**Metric**|**Dynamics Category**|**Temporal Interpretation**|
|---|---|---|
|**Burst Rate**|Density Estimation|Onsets per second (Hz).|
|**Bspks**|Density Estimation|Time-varying fraction of spikes.|
|**nspks**|Property Tracking|Mean spikes per burst.|
|**brstDur**|Property Tracking|Mean burst length (s).|
|**Freq**|Property Tracking|Internal firing intensity (Hz).|
|**IBI**|Property Tracking|Time between successive bursts.|

---

### Detailed Implementation Strategy

To implement this efficiently without re-running the heavy detection logic, you should adopt a "Post-Processing" workflow after your initial `brst_maxInt` run.

Step 1: Event Extraction

Extract the timestamps of every burst onset and every spike marked as a "burst spike." For property tracking, create a vector of the metric values (e.g., all dur values) paired with their corresponding onset times.

Step 2: Signal Generation

For Density Estimation, create a time vector $T$ with your desired resolution (e.g., 1s bins). For each unit, place a "1" at the onset times and convolve this vector with a Gaussian kernel. This results in a continuous estimate of the Burst Rate. You repeat this for "Burst Spikes" to get the Bursting Fraction over time.

Step 3: Smoothing and Filtering

For Property Tracking, you can use a sliding window of $N$ bursts (e.g., a 10-burst moving average). This is often better than a temporal window (e.g., 5 minutes) because it ensures the "smoothness" is consistent with the unit's activity level, regardless of how often it fires.