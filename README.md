# tuninator
Automatic analysis pipeline for FRA data collected through my other TDT/single-unit extraction pipeline.
Contains the automatic tuning analysis for others ease of use, also contains other analyses post processing that were used to make the figures in an associated paper/dataset.

These accessory custom analysis scripts used to analyze the data may be useful for others.
---
#### Table of Contents

1. [Set up](#setup)
2. [Pipeline Overview](#pipeline)
    - [Response Properties](#main_func)
    - [Thresholds](#thresholds)
    - [Bandwidth](#bw)
    - [Precision](#prec)
3. [Usage](#output)

---



<a name="setup" />

## Set up

Make sure machine has anaconda installed
To set up, 
download the github files to your computer either by manually clicking the download button or by cloning this repository.
download the needed packages (main dependency lies with pandas being at least version 2.2.1 to use the new map function).

```
conda activate your_env
conda install pandas,matplotlib,seaborn,scipy,numpy
pip install ipynb
```
Once you have needed packages proceed with opening the files

```
cd 'your/repo/path'
jupyter notebook 
```
naviagte to .ipynb file, and open.

---
<a name="pipeline" />

## Tuning pipeline overview

This analysis pipeline takes a csv in a structure of dBxkHz for each unit. This is not a typical data structure, but is very useful for some visualization tools also included in the accesory analyses folder. There is another program LINK WILL GO HERE, that will produce this matrix from either TDT multi-unit data or spike-sorted single-unit data. 

![polley_probken](https://github.com/electro-phys/tuninator/assets/155123673/943e61b3-d86e-4aa9-b348-65d4dcc6fceb)

---
<a name="main_func" />

### Main Response Properties Function

Once data is in this structure, the main function will calculate the prestimulus average spike count and std deviation. Note you may need to adjust the timing windows for your analysis given the way the auditory stimulation was performed. In the example case, we used a 20 ms prestimulus window and a 100 ms post-stimulus window for our calculations. The spike counts are then binned into 5 ms bins, where the peak firing bin and latency to peak bin are readily saved. The other information pertains to the first bin to go above $`4*\sigma(prestim)`$ and the last bin to go below $`4*\sigma(prestim)`$, or 20 ms after peak bin, whichever comes first[^1]. From this we can get response duration.

---


<a name="thresholds" />

### Automatic Threshold Calculation

First, the dBxKhz matrix for a given unit is smoothed using a 3x3 meidan filter to remove some of the jaggedness for the later smoothing step[^1]. This matrix is then summed across frequency (i.e. columns). This produces a dBx1 array that is essentially an input/output function losing the frequency domain. So if you graph this as a rate-level function (firing over sound intensity), you get a line. Thresholds are calculated based on a Savitzky–Golay filter with a 5-point qudratic filter similar to Guo 2012[^1]. These thresholds were checked against hand-scored data.

![image](https://github.com/electro-phys/tuninator/assets/155123673/58af363d-6aa9-43f0-8611-6f7a333605ad)

---

### Automatic Best Frequency Calculation

This proceeds exactly opposite to the Threshold calculation. Instead of summing across columns, we sum up the rows and keep the column information. This produces a curve with a maximum value. This maximum value corresponds to the Best Frequency (BF) (doesn't differ too much from the characteristic frequency in my experience and it's easier to calculate). This value is importnat for showing the distribution of BFs between groups to ensure they are equally represented as the BF can change the inherent bandwidth. This value will be used later for q-values.


---



<a name="bw" />

### Automatic bandwidth calculation

Bandwidth is calculated at each intensity above threshold for a given unit, which was calculated in the threshold step. Here, we simply go through each intensity row above threshold and grab the first and last sound-evoked cells, which was caclulated using the [resposne level function](#main_func), where the response crossing the set threshold line is considered to be the edge. There can be mulitple crossing so we take the first and the last to keep the bandwidth estimate conservative.  These become the bounds for the bandwidth calculation where $a$ is the upper bound frequency (i.e. 40000 Hz) and $b$ is the lower bound frequency (i.e. 20000 Hz)  which is simply $bw = log{_2}{(\frac{a}{b})}$. This is then the bandwidth for that intensity above threshold in octaves.

![badnwidth](https://github.com/electro-phys/tuninator/assets/155123673/6df71a76-d542-4f34-b9a8-5f57efa74ed7)

---

#### Q-values

The Q-value is a way to account for frequency specific bandwidth discrepancies. Higher frequencies have a narrow initial bandwidth than lower, for example. Make sure that frequency is in kHz for this step just in case (shouldn't matter but I like to be consistent). We follow the standard form at a given intensity $db$, $q_db = \frac{BF}{bw}$.

![tuning_change_by_position](https://github.com/electro-phys/tuninator/assets/155123673/4ec70710-6e9b-4ba1-a66f-befc0521dc7c).


---
#### d-prime 

Calculating d-prime takes everything we have produced so far and, in a way, summarizes it for each units dBxkHz tuning matrix. In short, it takes the firing properties inside a defiend tone-evoked region and compares to a defined non-evoked region. The farther apart these values are, the 'better' the d-prime. Another words the easier it is to tell these areas apart, or the information encoding is more optimal for higher d-primes. $d' = \frac{&mu;_1 - &mu;_2}{\frac{\sigma_1 - \sigma_2}{2}}$ . 

##### 1st method (Crude)

d-prime was first calculated by taking the previously calculated threshold for the given unit and assuming that every cell below that in intensity is non-evoked. Then for every evoked classified above that, it had to be connected to at least one other evoked cell at a higher intensity. However, our evoked criteria of exceeding $4*\sigma(prestim)$ gives a broad tuning map. So the d-prime calculated here, while sound is much too broad compared to prior results.

##### 2nd method (I/O)

d-prime was also calculated using the same smoothing I/O function described above, but for each frequency. This gave a threshold for each frequency such that anything below the threhsold for each column was deemed non-evoked, while everything above was deemed evoked on a per frequency (column) basis. This yielded more reasonable d-prime results.

##### 3rd method (d-prime per intensity)

d-prime was then calculated on a per intensity basis as this was being accessed for bandwidth calculations anyway. It takes a tuning curve like shown above in the bandwidth function for the given intensity, and takes all columns corresponding to being above the threshold value as evoked and all on the outside as non-evoked. This was consistent with the 2nd method and also produces a d-prime for each intensity, which is interesting. 



#### Precision metric

Another metric that was made is the precision metric. Maybe needs a more descriptive name, but it calculates the area under the best frequency tuning curve (i.e. across all intensities). To account for tuning curves with similar area but differing number peaks, where we defined more peaks as being less precise we use 

$$P_c = M * \int_a^b f(x) dx $$

$P_c$ denotes the Precision of a given curve where $M$ is the number of local maxima after smoothing with a Savitzky-Golay filter, and $f(x)$ being the function of the smoothed curve.

<a name="prec" />

---





<a name="output" />



### Usage
Go to the jupyter tutorial notebook for how functions are used 


[^1]: Guo et al. 2012,Robustness of Cortical Topography across Fields, Laminae, Anesthetic States, and Neurophysiological Signal Types [paper link](https://www.jneurosci.org/content/32/27/9159).
Cite the I/O paper here
