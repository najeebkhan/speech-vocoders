using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Levinson_Durbin_Algorithm;
namespace CleanVocoder
{
    class Program
    {
        #region Parameters
        //Frame Length
        static int N = 256; //240 15ms
        static int flen = 160; //160 10ms
        //Prediction Order
        static int P = 255;

        //post emphasis
        static double PE = 0.0;

        #endregion

        #region Variables
        static Random rand = new Random();
        static double[] LPCParams;
        static double[] Sigmas;
        static double[] Pitches;
        static int pc = 0, fc = 0, nb = 0, Nblock = 0, vlast = 0, vL = 1, vR = 0, nn = 0;
        static double pitchL = 0, gainL = 0;
        static double pitchR = flen, gainR = 0;
        static double pitch = flen, sigma, wlc, gain = 1, drv = 0, drvn = 0, yprev;
        static double[] LPCL = new double[P + 1];
        static double[] LPCR = new double[P + 1];
        static double[] LPC = new double[P + 1];
        static double[] yout = new double[flen];
        static double[] memory = new double[P + 1];
        static double[] synth;
        #endregion


        static void Main(string[] args)
        {           
            //get the pitch lpc and sigma
            string testfile = "haa";
            Pitches = CleanVocoder.pitch.computePitch(testfile);
            LPCParams = CleanVocoder.Param.ParamVector(testfile, ref Sigmas, N, flen, P:P, L: 100);

            Nblock = LPCParams.Length / (P + 1);
            synth = new double[Nblock * flen]; 
            process();

            #region normalize
            double max = 0;
            for (long i = 0; i < synth.Length; i++)
                if (Math.Abs(synth[i]) > max)
                    max = Math.Abs(synth[i]);
            for (long i = 0; i < synth.Length; i++)
                synth[i] /= max;
            #endregion

            #region write files
            using (BinaryWriter b = new BinaryWriter(File.Open(testfile + "Synth.dbl", FileMode.Create)))
            {
                //convert int16
                Int16[] outfile = new Int16[synth.Length];
                for (int i = 0; i < outfile.Length; i++)
                    outfile[i] = (Int16)(synth[i]);
                // write shorts to file
                foreach (double i in synth)
                    b.Write(i);
            }

            using (BinaryWriter b = new BinaryWriter(File.Open(testfile + "LPCParams.dbl", FileMode.Create)))
            {
                
                // write shorts to file
                foreach (double i in LPCParams)
                    b.Write(i);
            }
            using (BinaryWriter b = new BinaryWriter(File.Open(testfile + "sigma.dbl", FileMode.Create)))
            {

                // write shorts to file
                foreach (double i in Sigmas)
                    b.Write(i);
            }

            #region debug
            
            string[] stringVal = new string[synth.Length];

            for (int k = 0; k < synth.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(synth[k]);
            }

            File.WriteAllLines("Synth", stringVal);
            
            stringVal = new string[LPCParams.Length];

            for (int k = 0; k < LPCParams.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(LPCParams[k]);
            }

            File.WriteAllLines("LPCParams", stringVal);
            /*
            stringVal = new string[Pitches.Length];

            for (int k = 0; k < Pitches.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(Pitches[k]);
            }

            File.WriteAllLines("Pitches", stringVal);
            */

            stringVal = new string[Sigmas.Length];

            for (int k = 0; k < Sigmas.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(Sigmas[k]);
            }

            File.WriteAllLines("Sigmas", stringVal);
            
            #endregion
            
            #endregion
            
            Console.WriteLine("done");
            Console.Read();
            
            #region Debugging
            /*
            double[] x={ 76,   -38,  -132,    10,   110,    27,   -19,   -42,   127,    91};
            double e=0;
            double[] pram = LPCAnalysis.LDAlgo(x, 3, ref e);
            foreach(double i in pram)
                Console.WriteLine(i);
            Console.WriteLine(e);
            Console.Read();
            */
            #endregion
        }


        static void process()
        {
            while (nb < Nblock)
            {
                if (pc < pitch)
                    synthesize();
                else
                {
                    pc = 0;
                    if (fc >= flen)
                        boundary();
                    Interpolation();                    
                }
            }

        }

        static void boundary()
        {
            fc = fc - flen;
            for (int j = 0; j <= P; j++)
                LPCL[j] = LPCR[j];
            pitchL = pitchR;
            gainL = gainR;
            vlast = vL;
            vL = vR;
            //going from v to uv
            if (vlast == 1 | vL == 0)
                for (int j = 0; j <= P; j++)
                    memory[j] = 0;
            //if (vlast != vL)
            //    for (int j = 0; j <= P; j++)
            //        memory[j] = 0;

            ///////////////////////////////////////////////////////           
            //Read parameters
            for (int j = 0; j <= P; j++)
                LPCR[j] = LPCParams[j + nb * (P + 1)];

            pitchR = Pitches[nb];
            sigma = Sigmas[nb];
            //////////////////////////////////////////////////////
            nb++;
            if (pitchR <= 0)
            {
                gainR = sigma * Math.Sqrt(3 / (double)N) * 10;
                vR = 0;
            }
            else
            {
                gainR = sigma / Math.Sqrt(N);
                vR = 1;
            }

        }

        static void Interpolation()
        {
            if (vR == 1 && vL == 1)
            {
                wlc = (double)(fc) / (double)(flen);                
                //wlc = 0;
                pitch = (pitchR - pitchL) * wlc + pitchL;
            }
            else
            {
                pitch = pitchL;
                if (vL == 0)
                    pitch = flen - fc; // +1;
                wlc = 0;
            }
            drv = 1;
            drvn = -1 / (pitch - 1);

            for (int j = 0; j <= P; j++)
                LPC[j] = (LPCR[j] - LPCL[j]) * wlc + LPCL[j];
            gain = (gainR - gainL) * wlc + gainL;
            if (vL == 1)
                gain = gain * Math.Sqrt(pitch);

        }

        static void synthesize()
        {
            if (vL == 0)
                drv = Gausian(0, 0.1, rand);
            double y = AllPoleFilter(gain * drv, LPC, ref memory);

            yout[nn] = y + PE * yprev;
            yprev = yout[nn];
            drv = drvn;
            fc++; nn++; pc++;
            if (nn >= flen)
            {
                nn = 0;
                for (int t = 0; t < yout.Length; t++)
                    synth[t + nb * flen] = yout[t];
            }
        }

        static double Gausian(double mean, double sigma, Random rand)
        {

            double u1 = rand.NextDouble(); //these are uniform(0,1) random doubles
            double u2 = rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            return (mean + sigma * randStdNormal); //random normal(mean,stdDev^2)

        }

        static public double AllPoleFilter(double input, double[] A, ref double[] memory)
        {
            double FilteredSpeech = new double();
            int N = A.Length;

            FilteredSpeech = input;
            for (int j = 1; j < N; j++)
                FilteredSpeech = FilteredSpeech - A[j] * memory[j];
            memory[0] = FilteredSpeech;
            for (int j = N - 1; j > 0; j--)
                memory[j] = memory[j - 1];

            return FilteredSpeech;
        }        

        
    }
}
