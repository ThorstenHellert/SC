SC - Toolkit for Simulated Commissioning
----------------------------------------

We present the Toolkit for Simulated Commissioning (SC), which allows for realistic commissioning simulations of storage ring light sources by taking into account a multitude of error sources as well as diligently treating beam diagnostic limitations. 

Installation
------------

The SC toolkit requires the Matlab based Accelerator Toolbox (AT) which can be found [here](https://github.com/atcollab/at).

1. Install AT

2. Download the latest version of SC, for example by cloning the git repository
```
$ git clone https://github.com/ThorstenHellert/SC.git
```

3. Include the downloaded folder into the Matlab path, for example with
```
>> addpath('sc')
```

4. You're good to go. Please take a look at the [usage example](https://sc.lbl.gov/main.html#sec:example) or at the more complex [ALS-U Accumulator Ring scripts](https://github.com/ThorstenHellert/SC/tree/master/applications/ALSU_AR) to get startet.

Online documentation
--------------------

Please have a look at the [web site](https://sc.lbl.gov) and the [manual](https://sc.lbl.gov/main.html) for more details.
