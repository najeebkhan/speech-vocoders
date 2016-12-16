using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Levinson_Durbin_Algorithm;
using System.IO;
namespace Cepstrum
{
    class Cep
    {
        public static void ComputeFFT(long npt, int inv, ref Complex[] x)
        {

            int sign;
            long m, l, le, le1, k, ip, i, j;
            double ur, ui, wr, wi, tr, ti, temp;

            # region in place bit reversal
            j = 1;
            for (i = 1; i < npt; ++i)
            {
                if (i < j)
                {
                    tr = x[j].Re; ti = x[j].Img;
                    x[j].Re = x[i].Re;
                    x[j].Img = x[i].Img;
                    x[i].Re = tr;
                    x[i].Img = ti;
                    k = npt / 2;
                    while (k < j)
                    {
                        j = j - k;
                        k = k / 2;
                    }
                }
                else
                {
                    k = npt / 2;
                    while (k < j)
                    {
                        j = j - k;
                        k = k / 2;
                    }
                }
                j = j + k;
            }
            #endregion

            #region calculate the number of stages and direction
            m = (long)Math.Log(npt, 2);
            if (inv == 1)
                sign = 1;
            else
                sign = -1;
            #endregion
            
            #region perform the fft computation for each stage

            for (l = 1; l <= m; l++)
            {
                le = (long)Math.Pow(2, l);
                le1 = le / 2;
                ur = 1.0; ui = 0;
                wr = Math.Cos(Math.PI / le1);
                wi = sign * Math.Sin(Math.PI / le1);
                for (j = 1; j <= le1; ++j)
                {
                    i = j;
                    while (i <= npt)
                    {
                        ip = i + le1;
                        tr = x[ip].Re * ur - x[ip].Img * ui;
                        ti = x[ip].Img * ur + x[ip].Re * ui;
                        x[ip].Re = x[i].Re - tr;
                        x[ip].Img = x[i].Img - ti;
                        x[i].Re += tr;
                        x[i].Img += ti;
                        i = i + le;
                    }
                    temp = ur * wr - ui * wi;
                    ui = ui * wr + ur * wi;
                    ur = temp;

                }

            }
            #endregion

            #region normalize for inverse fft

            if (inv == 1)
            {
                for (i = 1; i <= npt; ++i)
                {
                    x[i].Re = x[i].Re / npt;
                    x[i].Img = x[i].Img / npt;

                }
            }
            #endregion



        }

        public static double[] Cepstrum(double[] buffW,int C, bool complex)
        {
            #region variables
            int N = buffW.Length;
            Complex[] x = new Complex[N + 1];
            double[] Cepstrum = new double[C];
            double[] PSD = new double[N];
            #endregion

            #region prepare data for fft
            for (int i = 1; i <= N; i++)
            {
                x[i] = new Complex();
                x[i].Re = buffW[i - 1];
            }
            #endregion

            #region Compute FFT
            Cep.ComputeFFT(N, 0, ref x);
            #endregion

            #region Compute Log PSD
            double[] ang = PhaseUnWrap(x);
                  
            for (int j = 0; j < N; j++)
            {
                if (complex)
                {
                    x[j + 1].Re = Math.Log(Math.Sqrt(Math.Pow(x[j + 1].Re, 2) + Math.Pow(x[j + 1].Img, 2)));
                    x[j + 1].Img = ang[j];
                }
                else
                {
                    x[j+1].Re = Math.Log(Math.Sqrt(Math.Pow(x[j + 1].Re, 2) + Math.Pow(x[j + 1].Img, 2)), Math.E);
                    x[j+1].Img = 0;
                }
            }
            #endregion
            
            #region Compute Inverse FFT
            Cep.ComputeFFT(N, 1, ref x);
            #endregion
            
            #region save Cepstrum
            for (int j = 0; j < C; j++)
            {
                Cepstrum[j] = x[j + 1].Re;
            }
            #endregion

            return Cepstrum;
        }

        public static double[] C2LPC(double[] C, int P, ref double E)
        {            
            double[] LPC=new double[P+1];
            LPC[0] = 1;        
            for (int n = 1; n <= P; n++ )
            {
                double sum = 0;
                for(int k=1; k<=n-1; k++)                                
                    sum = sum + (n - k) * C[n - k] * LPC[k];                
                LPC[n] = -C[n] - sum / n;  
            }
            E = Math.Exp(C[0]);
            return LPC;

        }
        
        public static double[] LPC2C(double[] LPC, double r0, int n)
        {
            LPC[0] = 1;            
            int P = LPC.Length - 1;
            double[] C = new double[n];

            for (int m = 1; m <= P; m++)
            {
                double sum = 0;
                for (int k = 1; k <= m - 1; k++)
                    sum = sum - (m - k) * C[m - k] * LPC[k];
                C[m] = -LPC[m] + sum/m;
            }
            for (int m = P + 1; m < n; m++)
            {
                double sum = 0;
                for (int k = 1; k <= P; k++)
                    sum = sum - (m - k) * C[m - k] * LPC[k] / m;
                C[m] = sum;
            }
            C[0] = Math.Log(r0);
            return C;
        }
        
        public static double[] LifterSym(double[] C, int L)
        {
            int N = C.Length;
            double[] LifteredCep = new double[N];
            for (int i = 0; i < N/2; i++)
            {
                if (i < L)
                {
                    LifteredCep[i] = C[i];
                    if (i < L - 1)
                        LifteredCep[N - 1 - i] = C[N - 1 - i];                   

                }
                else
                    LifteredCep[i] = 0;
            }
            return LifteredCep;
        }


        public static double[] PhaseUnWrap(Complex[] x)
        {
            double e=0.01;
            int N = x.Length-1;
            double[] arg = new double[N];
            double[] ARG = new double[N];
            int[] r = new int[N];

            //compute the principal angle
            for (int k = 0; k < ARG.Length; k++ )
                ARG[k] = Math.Atan(x[k+1].Img / x[k+1].Re);

            r[0] = 0;
            for (int k = 1; k < r.Length/2; k++)
            {
                if (ARG[k] - ARG[k - 1] > 2 * Math.PI - e)
                    r[k] = r[k - 1] - 1;
                else
                    if (ARG[k] - ARG[k - 1] < -(2 * Math.PI - e))
                        r[k] = r[k - 1] - 1;
                    else
                        r[k] = r[k - 1];                
            }


            // Compute the unwraped angle
            for (int k = 0; k < arg.Length; k++)
            {
                if (k < arg.Length / 2)
                    arg[k] = ARG[k] + 2 * Math.PI * r[k];
                else
                    if (k > arg.Length / 2)
                        arg[k] = -arg[arg.Length - k];
                    else
                        arg[k] = 0;
            }
            
            return ARG;
        }
    }

    class Complex
    {
        public double Re = 0, Img = 0;
        public Complex()
        {
            Re = 0;
            Img = 0;
        }


    }

}
