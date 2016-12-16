using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Levinson_Durbin_Algorithm;
namespace CleanVocoder
{
    class pitch
    {
        #region Pitch Parameters
        //Frame Length
        static int FLenPitch = 60;
        //Prediction Order
        static int PPitch = 4;
        // frame overlap
        static int ShiftPitch = 20;
        //freq range
        //fs/range
        static int pdlow = 6;
        static int pdhigh = 32;
        #endregion
        
        static public double[] computePitch(string testfile)
        {
            #region reading speech from file
            ////Read the Speech file
            double[] totalspeech;

            
            //Reading the input file
            using (BinaryReader b = new BinaryReader(File.Open(testfile + ".raw", FileMode.Open)))
            {

                //b.BaseStream.Length returns length in bytes
                Int16[] ispeech = new Int16[(long)b.BaseStream.Length / 2];
                for (long i = 0; i < ispeech.Length; i++)
                {
                    ispeech[i] = b.ReadInt16();
                }
                totalspeech = new double[ispeech.Length];
                for (int i = 0; i < totalspeech.Length; i++)
                    totalspeech[i] = (double)ispeech[i];


            }
            #endregion

            #region decimation filter coefficients
            //1KHz cutof, Fs=16KHz 
            //Numerator
            double[] B = {0.00016930232490151135,
                        -0.001238761236660607,
                        0.0040991633524885794,
                        -0.0077941999516641528,
                        0.0086415761242930236,
                        -0.0038767913077029925,
                        -0.0038767913077029925,
                        0.0086415761242930236,
                        -0.0077941999516641528,
                        0.0040991633524885794,
                        -0.001238761236660607,
                        0.00016930232490151135};
            //Denominator
            double[] A = {1,
                        -10.124736342158684,
                        47.097676411086375,
                        -132.82310687048624,
                        252.2606428662661,
                        -338.71142486654685,
                        328.03938767001358,
                        -229.1378574974072,
                        113.12211463822148,
                        -37.591292847078094,
                        7.5679264430264492,
                        -0.69932902632565963};
            #endregion

            #region Variable Declarations
            //for storing Inverse filter memory
            double[] memoryInv = new double[PPitch + 1];
            //LPC Coefficients
            double[] LPCP = new double[(int)totalspeech.Length * (PPitch + 1) / FLenPitch];
            double[] templpcp = new double[PPitch + 1];
            double[] tempresd = new double[FLenPitch];
            double[] Residue = new double[totalspeech.Length / 8];
            //double[] Residue = new double[(totalspeech.Length / (8 * Shift)) * FrameLength];
            double[] Correlation = new double[(totalspeech.Length / (8 * ShiftPitch)) * FLenPitch];
            Int16[] FinalPitch = new Int16[totalspeech.Length / 8];
            double[] Pitches = new double[(totalspeech.Length / (8 * ShiftPitch))];

            double[] PeakValues = new double[totalspeech.Length / 8];
            double[] PeakLocs = new double[totalspeech.Length / 8];
            int Decimate = 8;
            #endregion

            #region filtering and downsampling

            double[] PaddedSpeech = new double[totalspeech.Length + B.Length];
            double[] FilteredSpeech = new double[totalspeech.Length + B.Length];
            for (int i = 0; i < PaddedSpeech.Length; i++)
                if (i < totalspeech.Length)
                    PaddedSpeech[i] = totalspeech[i];
                else
                    PaddedSpeech[i] = 0;

            //Calculate decimation filter output                
            IIRDFII(PaddedSpeech, B, A, ref FilteredSpeech);
            //discard n samples from output
            for (int i = B.Length; i < FilteredSpeech.Length; i++)
                totalspeech[i - B.Length] = FilteredSpeech[i];

            //Decimation
            double[] DownSampledSpeech = new double[totalspeech.Length / Decimate];
            for (int i = 0; i < totalspeech.Length - Decimate; i += Decimate)
                DownSampledSpeech[i / Decimate] = totalspeech[i];

            #endregion

            //for each frame            
            int fcenter = FLenPitch / 2;
            int l = 0;
            double P0 = 0, P1 = 0, P2 = 0, P3 = 0;
            while (fcenter + FLenPitch / 2 < DownSampledSpeech.Length)
            {
                #region windowing and preemphasis

                double[] buff = new double[FLenPitch];
                double[] buffW = new double[FLenPitch];
                double prev = 0;
                double preemp = 0.9;
                for (int j = (fcenter - FLenPitch / 2); j < fcenter + FLenPitch / 2; j++)
                {
                    buffW[j - (fcenter - FLenPitch / 2)] = (0.54 - 0.46 * Math.Cos(2 * Math.PI * (j - (fcenter - FLenPitch / 2)) / (FLenPitch - 1))) * (DownSampledSpeech[j] - preemp * prev);
                    prev = DownSampledSpeech[j];
                    buff[j - (fcenter - FLenPitch / 2)] = DownSampledSpeech[j];
                }

                #endregion

                #region Calculate LPC parameters
                double E = 0;
                templpcp = LPCAnalysis.LDAlgo(buffW, PPitch, ref E);
                #endregion

                #region get the residual signal
                //tempresd = FIRDFI(buff, templpc);
                tempresd = FIRDFI(buff, templpcp, ref memoryInv);
                //apply window to residual
                for (int j = 0; j < tempresd.Length; j++)
                    tempresd[j] = tempresd[j] * (0.54 - 0.46 * Math.Cos(2 * Math.PI * j / (FLenPitch - 1)));
                #endregion

                #region calcute autocorrelation
                double[] AutoCorr = new double[tempresd.Length];
                //Compute Correlation Coefficients
                for (int j = 0; j < tempresd.Length; j++)
                {
                    for (int k = 0; k < tempresd.Length - j; k++)
                    {
                        AutoCorr[j] += tempresd[k] * tempresd[j + k];
                    }
                }
                double[] NAutoCorr = new double[tempresd.Length];
                //normalize the correlation
                for (int j = 0; j < tempresd.Length; j++)
                    NAutoCorr[j] = AutoCorr[j] / AutoCorr[0];


                #endregion

                #region find the peak
                int peakloc = 0;
                double peakval = 0;
                for (int j = pdlow; j <= pdhigh; j++)
                {
                    if (Math.Abs(NAutoCorr[j]) > peakval)
                    {
                        peakval = NAutoCorr[j];
                        peakloc = j;
                    }
                }
                #endregion
                PeakValues[fcenter] = peakval;
                PeakLocs[fcenter] = peakloc;

                if (peakval == 0 | NAutoCorr[peakloc] < NAutoCorr[peakloc - 1])
                    P0 = 0;
                else
                {
                    #region Interpolation
                    double a = (NAutoCorr[peakloc - 1] + NAutoCorr[peakloc + 1]) / 2 - NAutoCorr[peakloc];
                    double b = (NAutoCorr[peakloc + 1] - NAutoCorr[peakloc - 1]) / 4;
                    double ap = NAutoCorr[peakloc] - b * b / a;
                    double al = peakloc - b / a;
                    peakval = ap;   //in the fortran this is divided by ABUF(1)
                    #endregion

                    #region apply variable threshold
                    double threshold = 1;
                    //18 is 9ms for 2khz
                    if (peakloc < 19)
                        threshold = -1 * (peakloc - pdlow) / 13 + 2;
                    else
                        threshold = -1 * (peakloc - 19) / 13 + 1;
                    //scale the peak values
                    peakval /= threshold;
                    #endregion

                    if (peakval > 0.35)
                        P0 = al;
                    else
                        if (P1 != 0 & peakval > 0.3)
                            P0 = al;
                        else
                            P0 = 0;
                }


                #region Decisions

                if (Math.Abs(P1 - P3) < 0.375 * P3)
                    P2 = (P1 + P2) / 2;

                if (P3 == 0 & P2 != 0 & Math.Abs(P0 - P1) <= 0.2 * P1)
                    P2 = 2 * (P1 - P0);

                if (P1 != 0)
                {
                    P3 = P2;
                    P2 = P1;
                    P1 = P0;
                }
                else
                {
                    if (Math.Abs(P2 - P3) > 0.375 * P3)
                        P2 = 0;
                    P3 = P2;
                    P2 = P1;
                    P1 = P0;
                }
                #endregion

                if (l > 2)
                    if (P3 != 0 & Math.Abs((Int16)(2000 / P3)) < 500)
                        Pitches[l - 3] = P3 * 8;
                    else
                        Pitches[l - 3] = 0;
                l++;
                fcenter += ShiftPitch;
            }
            return Pitches;
        }
        
        static public void IIRDFII(double[] inputspeech, double[] B, double[] A, ref double[] FilteredSpeech)
        {

            int L = inputspeech.Length;
            int N = B.Length;
            double temp = 0;
            double[] memory = new double[N];
            for (int i = 0; i < L; i++)
            {
                temp = 0;
                memory[0] = inputspeech[i];
                for (int j = 1; j < N; j++)
                    memory[0] -= A[j] * memory[j];
                for (int j = 0; j < N; j++)
                    temp += (B[j] * memory[j]);
                for (int j = N - 1; j > 0; j--)
                    memory[j] = memory[j - 1];
                FilteredSpeech[i] = temp;

            }

            return;
        }

        static public double[] FIRDFI(double[] inputspeech, double[] B, ref double[] memoryInv)
        {

            int L = inputspeech.Length;
            int N = B.Length;
            double temp = 0;
            double[] filtered = new double[L];
            //double[] memoryInv = new double[N];
            //because the lpc coefficients are with opposite sign
            for (int k = 1; k < B.Length; k++)
                B[k] = -B[k];

            for (int i = 0; i < L; i++)
            {

                temp = 0;
                memoryInv[0] = inputspeech[i];
                for (int j = 0; j < N; j++)
                    temp += B[j] * memoryInv[j];
                for (int j = N - 1; j > 0; j--)
                    memoryInv[j] = memoryInv[j - 1];
                filtered[i] = temp;
            }

            return filtered;
        }
        

            

    }
}
