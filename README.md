# Machine Learning for Wireless Communication Project

<h3 align="center">Keshav Goyal(IMT2022560)</h3>
<h3 align="center">R. Harshavardhan(IMT2022515)</h3>

# A Bayesian Approach to Characterize Unknown Interference Power in Wireless Networks

## Issue
In MIMO networks, a large number of user equipments (UEs) are served on the same time-frequency resources.  
Outage performance is dominated by unknown interference arising from scheduling variations in neighboring cells.

## Solution
- Characterize the distribution of total unknown interference power from UEs in neighboring cells after observing a certain number of samples.
- Compute the **sample mean** and **sample variance**, and estimate the **α** and **β** (shape and scale parameters) of the **Inverse-Gamma distribution** to capture the statistical behavior of unknown interference.
- Compute the **interference threshold** for a **target outage probability**, enabling rate adaptation that guarantees reliability.

---

# Results

## Network Simulation Parameters

| Parameter                  | Value            |
|-----------------------------|------------------|
| Bandwidth                   | 20 MHz           |
| Number of BS antennas (N)   | 16               |
| Pathloss exponent (α)       | 3.76             |
| UL transmit power (pᵢ)       | 100 mW           |
| UL noise power              | -94 dBm          |
| Coherence block length (τc) | 200              |
| Pilot sequence length (τp)  | 10               |

## **Result 1**
   - **The CDF of SINR(dB) for RZF and MR combining**<br> ![Result1](https://github.com/keshavv79/MLforWirelessProject/blob/main/resultsForGit/CDFSINR.png)
## **Result 2**
   - **Site Map and UE Distribution Example**<br> ![Result1](https://github.com/keshavv79/MLforWirelessProject/blob/main/resultsForGit/UEBS.png)
## **Result 3**
   - **Epsilon Outage Proability vs Spectral Efficiency Curve**<br> ![Result1](https://github.com/keshavv79/MLforWirelessProject/blob/main/resultsForGit/epsilonOutage.png)

