#How to make extractRaw for your flavour of Linux

To compile the executable you will need to have Modo-develop installed on your machine. The two dll files, `ThermoFisher.CommonCore.Data.dll` and `ThermoFisher.CommonCore.RawFileReader.dll`, must be in the same directory as `extractRaw.cs`.  

Then, from within this directory execute:  
`mono-csc extractRaw.cs -r:ThermoFisher.CommonCore.Data.dll -r:ThermoFisher.CommonCore.RawFileReader.dll`  

You can now use extractRaw with Mono to extract spectra from raw files like this:  
`mono extractRaw.exe file_name.raw`

To make a bundle that does not depend on mono being installed on the target machine execute:  
`mkbundle -o extractRaw --simple extractRaw.exe --no-machine-config --no-config` 

This will create `extractRaw` executable that has mono bundled in. It can now be distributed to other machines.

