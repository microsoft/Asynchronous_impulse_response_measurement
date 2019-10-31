# Acoustic impulse response measurement without clock synchronization
This is a collection of Matlab scripts illustrating the clock drift estimation algorithm introduced in 

```
@InProceedings{gamper2017clock,
  author = {Gamper, Hannes},
  title = {Clock drift estimation and compensation for asynchronous impulse response measurements},
  booktitle = {Proc. Workshop on Hands-free Speech Communication and Microphone Arrays (HSCMA)},
  year = {2017},
  month = {March},
  url = {https://www.microsoft.com/en-us/research/publication/clock-drift-estimation-compensation-asynchronous-impulse-response-measurements/},
}
```

The scripts can be used to estimate the (acoustic) impulse response of a device-under-test (DUT) in the presence of clock drift between the measurement setup and the DUT. 

For more information on the algorithm, see https://www.microsoft.com/en-us/research/publication/clock-drift-estimation-compensation-asynchronous-impulse-response-measurements/

# Running the example
Copy the scripts in this project to a folder and run RUNME.m from within that folder to run through the following steps:
* Generate a test signal for impulse response measurements
* Simulate a noisy, reverberant, asynchronous recording of the test signal
* Estimate the clock drift and impulse response

# Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.
