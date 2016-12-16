using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Levinson_Durbin_Algorithm
{
    public class LPCAnalysis
    {
        /// <summary>
        /// This method uses the Levinsen Durbin's Algorithm to 
        /// calculate the LPC parameters of the Input speech samples
        /// </summary>
        /// <param name="X">Input speech</param>
        /// <param name="P">Order of Predictor</param>
        /// <returns>Returns the LPC parameters</returns>
        /// 
        static public double[] LDAlgo(double[] X, int P, ref double E)
        {
            int N = X.Length;
            //Correlation Coefficients
            double[] R = new double[P + 1];
            //LPC Coefficients
            double[] A = new double[P + 1];
            double[] tmp = new double[P + 1];


            //Compute Correlation Coefficients R(i) for l=0...p-1
            for (int i = 0; i <= P; i++)
            {
                for (int j = 0; j < N - i; j++)
                {
                    R[i] += X[j] * X[i + j];
                }
            }

            A[0] = 1;
            E = R[0];

            // Levinson Durbin Algorithm Recursion
            for (int i = 1; i <= P; i++)
            {
                double K = 0;
                // Calculate K(j)
                for (int j = 1; j < i; j++)
                {
                    K += A[j] * R[i - j];
                }
                K = (R[i] - K) / E;

                A[i] = K;

                for (int j = 1; j < i; j++)
                {
                    tmp[j] = A[j] - K * A[i - j];
                }
                for (int j = 1; j < i; j++)
                {
                    A[j] = tmp[j];
                }

                E = (1 - (K * K)) * E;
            }


            return A;

        }

        static public double[] CepLDAlgo(double[] X, int P, ref double E )
        {
            int N = X.Length;
            //Correlation Coefficients
            double[] R = new double[P+1];
            //LPC Coefficients
            double[] A = new double[P+1];
            double[] tmp = new double[P+1];
            

            //Compute Correlation Coefficients R(i) for l=0...p-1
            for (int i = 0; i <= P; i++)
            {
                for (int j = 0; j < N - i; j++)
                {
                    R[i] += X[j] * X[i + j];
                }
            }

            A[0] = 1;
            E = R[0];
            
            // Levinson Durbin Algorithm Recursion
            for (int i = 1; i <= P; i++)
            {
                double K = 0;
                // Calculate K(j)
                for (int j = 1; j < i; j++)
                {
                    K += A[j] * R[i - j];
                }
                K = (R[i] - K) / E;

                A[i] = K;

                for (int j = 1; j < i; j++)
                {
                    tmp[j] = A[j] - K * A[i - j];
                }
                for (int j = 1; j < i; j++)
                {
                    A[j] = tmp[j];
                }

                E = (1 - (K * K)) * E;
            }

            for (int i = 1; i < A.Length; i++)
                A[i] = -A[i];

            E = R[0];
            return A;

        }
    }
}
