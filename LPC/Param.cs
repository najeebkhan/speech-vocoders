#define Rec

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Levinson_Durbin_Algorithm;
using Cepstrum;

namespace CleanVocoder
{
   
    class Param
    {

        public static double[] ParamVector(string testfile, ref double[] Sigma, int FrameLength, int Shift, int P, int L )
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

                       
            #region Variable Declarations            
            //LPC Coefficients
            double[] LPC = new double[(int)(totalspeech.Length * (P + 1) / Shift)];
            Sigma= new double[(int)(totalspeech.Length/Shift)];                      
            #endregion
                       
           //for each frame            
           int fcenter = FrameLength / 2;
           int l = 0;
            while (fcenter + FrameLength / 2 < totalspeech.Length)
            {
                #region Framing
                double[] buff = new double[FrameLength];                                
                for (int j = (fcenter - FrameLength / 2); j < fcenter + FrameLength / 2; j++)
                    buff[j - (fcenter - FrameLength / 2)] = (0.54 - 0.46 * Math.Cos(2 * Math.PI * (j - (fcenter - FrameLength / 2)) / (FrameLength - 1))) * totalspeech[j];
                #endregion

                #if (Rec)               
                //method 1
                #region Calculate LPC from Cepstrum
                double E = 0;

                //Compute Real Cepstrum
                double[] cep = Cep.Cepstrum(buff, FrameLength, complex: false);
                // Liftered Cepstrum
                double[] liftcep = Cep.LifterSym(cep, L);
                double[] lpcc = Cep.C2LPC(liftcep, P, ref E);
                //save lpc
                for (int j = 0; j <= P; j++)
                    LPC[l * (P + 1) + j] = lpcc[j];
                //save sigma
                #endregion

                
                
                #else
                

                //Compute Real Cepstrum
                double[] cep = Cep.Cepstrum(buff, FrameLength, complex: false);
                // Liftered Cepstrum
                double[] liftcep = Cep.LifterSym(cep, L);

                //step 1: remove the ifft at last stage of cepstrum

                #region variables
                int N = buff.Length;
                Complex[] x = new Complex[N + 1];
                
                #endregion

                #region prepare data for fft
                for (int i = 1; i <= N; i++)
                {
                    x[i] = new Complex();
                    x[i].Re = liftcep[i - 1];
                }
                #endregion

                #region Compute FFT
                Cep.ComputeFFT(N, 0, ref x);
                #endregion
                
                //step 2: remove the log of cepstrum
                                
                for (int j = 0; j < N; j++)
                {
                    x[j + 1].Re = Math.Exp(x[j + 1].Re);
                    //x[j + 1].Re = Math.Exp(Math.Log(Math.Sqrt(Math.Pow(x[j + 1].Re, 2) + Math.Pow(x[j + 1].Img, 2))));
                    x[j + 1].Img = 0;                    
                }
                
                //step 3: remove the fft in first step of cepstrum to get the v(n) of vocal tract                                
                
                //step 4: take the fft agnitude to get V(w)

                // 3 and 4 cancel out

                //step 5: V(w)=1/A(w) ---> A(w)=1/V(w)

                for (int j = 0; j < N; j++)
                {
                    x[j + 1].Re = 1.0 / x[j + 1].Re;
                    x[j + 1].Img = 0;
                }

                //step 6: get the lpc coefficients
                Cep.ComputeFFT(N, 1, ref x);

                double tmp = x[1].Re;
                for (int i = 1; i < x.GetLength(0); i++)
                    x[i].Re /= tmp;
                x[1].Re = 1;
                double E = 1 / tmp;
                
                
                //save lpc
                for (int j = 0; j <= P; j++)
                    LPC[l * (P+1) + j] = x[j+1].Re;
                
                
                #endif
               
                //if(Double.IsNaN(E)==false)                     
                //    Sigma[l] = Math.Sqrt(E / buff.Length);               
                if (Double.IsNaN(E) == false)
                    Sigma[l] = E;        
                    l++;
                fcenter += Shift;
            }

            return LPC;

        }


    }
}
