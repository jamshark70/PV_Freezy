# PV_Freezy

James Harkins (based on work by Joshua Parmenter and Sam Pluta)

This package adds a variant on Joshua Parmenter's PV_Freeze plug-in, from [sc3-plugins](https://github.com/supercollider/sc3-plugins).

The new variant, PV_Freezish, introduces two coefficients, one for FFT bins whose magnitude is increasing, and the other for bins whose magnitude is decreasing, so that the sound can be partially frozen and partially moving -- blurred in time.

Repository layout is modeled on Sam Pluta's [BufFFT](https://github.com/spluta/BufFFT) library. That's just to get the cmake stuff right.

PV_Freezish code directly quotes PV_Freeze -- accordingly, PV_Freezish is GPL-licensed and is shared back to the community as a derivative product, with full credit given to the original author.


### Building and Installation

The [SuperCollider source code repository](https://github.com/supercollider/supercollider) must be available in your system.

Also, minimum cmake 2.8.

In the project root directory, run:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSC_PATH="<Path to SC Source>"
cmake --build . --config Release
```

as you normally would when building SC plugins.

Then, copy or symlink the project directory into your user Extensions directory.
