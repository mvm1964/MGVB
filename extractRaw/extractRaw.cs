/*
 * Raw file extractor in C# to be used with rawExtract
 * Uses mono and Therno dll libraries
 * compile with:
 * csc extractRaw.cs -r:ThermoFisher.CommonCore.Data.dll -r:ThermoFisher.CommonCore.RawFileReader.dll
 *
 * Author: Metodi V. Metodiev
 * 30.01.2021
 * Minchini House, Ravna, Bulgaria
*/

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.ExceptionServices;

using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

public class ExtractRaw 
{
    public static void Main(string[] args)
    {
        string filename = string.Empty;
        filename = args[0];
        var rawFile = RawFileReaderAdapter.FileFactory(filename);
        Console.WriteLine ("Extracting {0}...", filename);

        if (!rawFile.IsOpen || rawFile.IsError)
        {
           Console.WriteLine("Unable to access the RAW file using the RawFileReader class!");
                    
           return;
        }

        // Check for any errors in the RAW file
        if (rawFile.IsError)
        {
            Console.WriteLine("Error opening ({0}) - {1}", rawFile.FileError, filename);

            return;
        }

        //Console.WriteLine("H\tThe RAW file has data from {0} instruments", rawFile.InstrumentCount);

        rawFile.SelectInstrument(Device.MS, 1);

        int firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum;
        int lastScanNumber = rawFile.RunHeaderEx.LastSpectrum;
        
        // open files to write to
        System.IO.StreamWriter file = new System.IO.StreamWriter(@filename + ".ms2");
        System.IO.StreamWriter file_ms = new System.IO.StreamWriter(@filename + ".ms1");        

        file.WriteLine("H   Instrument model: " + rawFile.GetInstrumentData().Model);
        file.WriteLine("H   Instrument name: " + rawFile.GetInstrumentData().Name);
        file.WriteLine("H   Serial number: " + rawFile.GetInstrumentData().SerialNumber);
        file.WriteLine("H   Software version: " + rawFile.GetInstrumentData().SoftwareVersion);
        file.WriteLine("H   Firmware version: " + rawFile.GetInstrumentData().HardwareVersion);
        file.WriteLine("H   Units: " + rawFile.GetInstrumentData().Units);
        file.WriteLine("H   Mass resolution: {0:F3} ", rawFile.RunHeaderEx.MassResolution);
        file.WriteLine("H   Number of scans: {0}", rawFile.RunHeaderEx.SpectraCount);
        file.WriteLine("H   Scan range: {0} - {1}", firstScanNumber, lastScanNumber);
        file.WriteLine("H   Mass range: {0:F4} - {1:F4}", rawFile.RunHeaderEx.LowMass, rawFile.RunHeaderEx.HighMass);


        for (int scan = firstScanNumber; scan <= lastScanNumber; scan++)
        {
            double time = rawFile.RetentionTimeFromScanNumber(scan);
            var scanFilter = rawFile.GetFilterForScanNumber(scan);
            var scanEvent = rawFile.GetScanEventForScanNumber(scan);
            
            string analyzer = scanFilter.ToString().Substring(0, 4);

            if (scanFilter.MSOrder == MSOrderType.Ms2)
            {
                /* Todo: should check if mass analyser ITMS then proceed as below
                   else, if it is FTMS should read centroid stream instead:
           
                   if (scanFilter.MassAnalyzer == "MassAnalyzerLQIT) {

                   }
                   else if (scanFilter.MassAnalyzer == "MassAnalyzerFTMS") 
                  

                */                
                var reaction = scanEvent.GetReaction(0);

                double precursorMass = reaction.PrecursorMass;
                int precursorCharge = 0;
                double collisionEnergy = reaction.CollisionEnergy;
                double isolationWidth = reaction.IsolationWidth;
                double monoisotopicMass = 0.0;
                //int masterScan = 0;
                var ionizationMode = scanFilter.IonizationMode;
                var order = scanFilter.MSOrder;
                var trailerData = rawFile.GetTrailerExtraInformation(scan);
                for (int i = 0; i < trailerData.Length; i++)
                {
                    //file.WriteLine(trailerData.Labels[i]);
                    if (trailerData.Labels[i] == "Monoisotopic M/Z:")
                    {
                        monoisotopicMass = Convert.ToDouble(trailerData.Values[i]);
                    }
                    /*
                    if ((trailerData.Labels[i] == "Master Scan Number:") || (trailerData.Labels[i] == "Master Scan Number") || (trailerData.Labels[i] == "Master Index:"))
                    {
                        masterScan = Convert.ToInt32(trailerData.Values[i]);
                    }
                    */
                    if (trailerData.Labels[i] == "Charge State:")
                    {
                        precursorCharge = Convert.ToInt32(trailerData.Values[i]);
                    }

                }
                file.WriteLine(
                            "S\tScan number {0} @ time {1:F2}, Ionization mode={2}, MS Order={3}, Precursor mass={4:F4}, Monoisotopic Mass = {5:F4}, Collision energy={6:F2}, Isolation width={7:F2}, Charge={8}",
                            scan, time, ionizationMode, order, precursorMass, monoisotopicMass, collisionEnergy, isolationWidth, precursorCharge);
                file.WriteLine("I\tFilter: {0}", scanFilter.ToString());
                file.WriteLine("I\tMS Analyzer: {0}", analyzer);
                file.WriteLine("I\tRetTime {0:F2}", time);
                //file.WriteLine("Z\t{0}\t{1}", precursorCharge, (monoisotopicMass*precursorCharge) - (precursorCharge - 1)*1.007276);

                /*
                   Update 7.07.2021: some scans do not have proper monoisotopic mass. Should test if monoisotopic mass is 0 
                   and if so, replace it with precursor mass. 
                */

                if (monoisotopicMass == 0.0){
                    if (precursorMass*precursorCharge <= 2100) monoisotopicMass = precursorMass;
                    else if (precursorMass*precursorCharge > 2100 || precursorMass*precursorCharge < 4200) monoisotopicMass = precursorMass - 1.007276/precursorCharge;
                    else monoisotopicMass = precursorMass - 2*1.007276/precursorCharge;
                }
                file.WriteLine("I\tPrecursorMass\t{0}", (monoisotopicMass*precursorCharge) - (precursorCharge - 1)*1.007276);
                file.WriteLine("Z\t{0}\t{1}", precursorCharge, (monoisotopicMass*precursorCharge) - (precursorCharge - 1)*1.007276);

                var scanStatistics = rawFile.GetScanStatsForScanNumber(scan);
                /*
                if (!scanStatistics.IsCentroidScan)
                {
                    Console.WriteLine("centroid data is present");
                    var centroidStream = rawFile.GetCentroidStream(scan, false);

                    Console.WriteLine("Spectrum (centroid/label) {0} - {1} points", scan, centroidStream.Length);

                    for (int i = 0; i < centroidStream.Length; i++)
                    {
                        Console.WriteLine("  {0} - {1:F4}, {2:F0}, {3:F0}", i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]);
                    }

                }*/

                //else
                //if (Convert.ToString(scanFilter.MassAnalyzer) == "MassAnalyzerITMS") 
                if (scanFilter.MassAnalyzer == MassAnalyzerType.MassAnalyzerITMS)
                {
                    // Get the segmented (low res and profile) scan data
                    var segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scan, scanStatistics);

                    //Console.WriteLine("Spectrum (normal data) {0} - {1} points", scan, segmentedScan.Positions.Length);

                    // Print the spectral data (mass, intensity values)
                    
                    for (int i = 0; i < segmentedScan.Positions.Length; i++)
                    {
                        file.WriteLine("{0:F4} {1:F0}", segmentedScan.Positions[i], segmentedScan.Intensities[i]);
                    }
                    
                }
                //Console.WriteLine();

                else if (scanFilter.MassAnalyzer == MassAnalyzerType.MassAnalyzerFTMS)
                {
                    //Console.WriteLine("centroid data is present");
                    var centroidStream = rawFile.GetCentroidStream(scan, false);

                    //Console.WriteLine("Spectrum (centroid/label) {0} - {1} points", scan, centroidStream.Length);
                    
                    for (int i = 0; i < centroidStream.Length; i++)
                    {
                        //Console.WriteLine("  {0} - {1:F4}, {2:F0}, {3:F0}", i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]);

                        // the extra column might be causing problems with parseMS. Will remove it and try
                        //file.WriteLine("{0:F4} {1:F0} {2:F0}", centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]);
                        file.WriteLine("{0:F4} {1:F0}", centroidStream.Masses[i], centroidStream.Intensities[i]);
                    }

                }                

   
            }
            
            else if (scanFilter.MSOrder == MSOrderType.Ms)
            {
                var scanDependents = rawFile.GetScanDependents(scan, 5);
                file_ms.WriteLine("I  Scan number {0}", scan);
                file_ms.WriteLine("I  MS Analyzer: {0}", analyzer);
                if (scan != lastScanNumber) file_ms.WriteLine("I  Number dependent scans={0}", scanDependents.ScanDependentDetailArray.Length);
                file_ms.WriteLine("I  RetTime {0:F2}", time);
                //var scanStatistics = rawFile.GetScanStatsForScanNumber(scan);         
                //if (scanStatistics.IsCentroidScan)
                {
                    //Console.WriteLine("centroid data is present");
                    var centroidStream = rawFile.GetCentroidStream(scan, false);

                    //Console.WriteLine("Spectrum (centroid/label) {0} - {1} points", scan, centroidStream.Length);

                    for (int i = 0; i < centroidStream.Length; i++)
                    {
                        file_ms.WriteLine("{0:F4} {1:F0} {2:F0}", centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]);
                    }

                }


            }

        }
        
        rawFile.Dispose();
        file.Close();
        file_ms.Close();
    }
}

