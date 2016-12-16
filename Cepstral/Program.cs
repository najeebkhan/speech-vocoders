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
        static int M = 70;
        

        //post emphasis
        static double PE = 0.0;

        
        #endregion

        #region Variables
        static Random rand = new Random();
        static double[,] coeff;       
        static double[] Pitches;
        static int pc = 0, fc = 0, nb = 0, Nblock = 0, vL = 1, nn = 0;       
        static double pitch = flen, drv = 0, drvn = 0, yprev;        
        static double[] yout = new double[flen];
        static double[,] memory = new double[M,2];
        static double[][] mem = new double[M][];
        static double[] synth;
        static double[,] coeffL=new double[M,2];
        
        #endregion


        static void Main(string[] args)
        {           
            //get the pitch lpc and sigma
            string testfile = "haa";
            Pitches = CleanVocoder.pitch.computePitch(testfile);            
            coeff = CleanVocoder.Param.ParamVector(testfile, N, flen,M);

            Nblock = Pitches.Length;
            synth = new double[Nblock * flen];
            
            mem[0] = new double[1];
            for (int m = 1; m < M;m++ )
                mem[m] = new double[2 * m];

            mem[0][0] = 0;
            for (int m = 1; m < M; m++)
                for (int n = 0; n < mem[m].Length; n++)
                    mem[m][n] = 0;

            
            process();

            #region normalize
            //double max = 0;
            //for (long i = 0; i < synth.Length; i++)
            //    if (Math.Abs(synth[i]) > max)
            //        max = Math.Abs(synth[i]);
            //for (long i = 0; i < synth.Length; i++)
            //    synth[i] /= max;
            #endregion

            #region write files
            using (BinaryWriter b = new BinaryWriter(File.Open(testfile + "Synth.raw", FileMode.Create)))
            {
                //convert int16
                Int16[] outfile = new Int16[synth.Length];
                for (int i = 0; i < outfile.Length; i++)
                    outfile[i] = (Int16)(synth[i]*2000);
                // write shorts to file
                foreach (Int16 i in outfile)
                    b.Write(i);
            }                    

            #region debug
           
            string[] stringVal = new string[synth.Length];

            for (int k = 0; k < synth.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(synth[k]);
            }

            File.WriteAllLines("Synth", stringVal);
            /*
                       
            stringVal = new string[Pitches.Length];

            for (int k = 0; k < Pitches.Length; k++)
            {
                stringVal[k] = System.Convert.ToString(Pitches[k]);
            }

            File.WriteAllLines("Pitches", stringVal);           

            */
            #endregion            
            #endregion
            
            Console.WriteLine("done");
            Console.Read();           
            
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
            int LastVL = vL;
            ///////////////////////////////////////////////////////           
            //Read parameters
            for (int j = 0; j < M; j++)
            {
                coeffL[j, 0] = coeff[j + nb * M, 0];
                coeffL[j, 1] = coeff[j + nb * M, 1];
            }
            pitch = Pitches[nb];
            //////////////////////////////////////////////////////
            vL = 1;
            if(pitch==0)
            {
                vL = 0;
                pitch = flen - fc;                
            }           
            
            //from v to uv reset memory
            if (LastVL == 1 | vL == 0)
            {
                mem[0] = new double[1];
                for (int m = 1; m < M; m++)
                    mem[m] = new double[2 * m];
            }

            nb++;
        }


        static void Interpolation()
        {           
            drv = 1;
            drvn = -1 / (pitch - 1);          

        }


        static void synthesize()
        {
            if (vL == 0)
                drv = Gausian(0, 0.1, rand);
            double y = EDFilter(drv, coeffL);
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
        
        static double EDFilter(double inp, double[,] coef)
        {
            double y = 0, u;

            //stage 0
            //cascading L stages
            int L = 10;
            double factor = 1 / L;
            double a = coef[0, 0] * factor;
            double b = coef[0, 1] * Math.Pow(factor, 2);
            for (int k = 1; k <= L;k++ )
            {
                y = inp * (1 + a + b) / (1 - a + b);
                inp = y;  //cascaded next input = current output
            }


           
            //stage 1
            //cascading L stages
            //L = 1;
            //factor = 1 / L;
            //a = coef[1, 0] * factor;
            //b = coef[1, 1] * Math.Pow(factor, 2);
            //for (int k = 1; k <= L; k++)
            //{
            //    u = inp + a * mem[1][0] - b * mem[1][1];
            //    y = u + a * mem[1][0] + b * mem[1][1];
            //    inp = y;  //cascaded next input = current output
            //    //shift the delayline
            //    mem[1][1] = mem[1][0];
            //    mem[1][0] = u;
            //}


            ////stage 2
            ////cascading L stages
            //L = 20;
            //factor = 1 / L;
            //a = coef[2, 0] * factor;
            //b = coef[2, 1] * Math.Pow(factor, 2);
            //for (int k = 1; k <= L; k++)
            //{
            //    u = inp + a * mem[2][1] - b * mem[2][3];
            //    y = u + a * mem[2][1] + b * mem[2][3];
            //    inp = y;  //cascaded next input = current output
            //    //shift the delayline
            //    mem[2][3] = mem[2][2];
            //    mem[2][2] = mem[2][1];
            //    mem[2][1] = mem[2][0];
            //    mem[2][0] = u;
            //}

            //for stage 3 on
            for (int m = 1; m < M; m++) 
            {
                u = inp + coef[m, 0] * mem[m][m-1] - coef[m, 1] * mem[m][2 * m-1];
                y = u + coef[m, 0] * mem[m][m-1] + coef[m, 1] * mem[m][2 * m-1];
                inp = y;  //cascaded next input = current output
                //shift the delayline
                for (int i = 2 * m-2; i >= 0; i--)
                    mem[m][i+1] = mem[m][i];
                mem[m][0] = u;

            }
            return y;
        }

        
    }
}
