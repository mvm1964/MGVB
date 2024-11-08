To compile the executable you will need to have Modo-develop installed on your machine. The tw dll file must be in the same directory as the extractRaw.cs.  
Then, from within this directory execute:  
`mono-csc extractRaw.cs -r:ThermoFisher.CommonCore.Data.dll -r:ThermoFisher.CommonCore.RawFileReader.dll`  
You can now use extractRaw with Mono to extract spectra from raw files like this:  
`mono extractRaw.exe file_name.raw`
