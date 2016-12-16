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

        public static double[,] ParamVector(string testfile, int FrameLength, int Shift, int M)
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
            int x = 1 + FrameLength / 2;
            //filter Coefficients           
            double[,] coeff = new double[(int)(totalspeech.Length * x / Shift), 2];                   
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

                double[] cep=Cep.Cepstrum(buff,false );
                for (int n = 1; n < M; n++)
                {
                    cep[n] *= 2;
                }
               
                for (int n = 0; n < M; n++ )
                {
                    coeff[n + l * M, 0] = cep[n]/2;
                    coeff[n + l * M, 1] = cep[n] * cep[n]/12;
                }
               
                l++;
                fcenter += Shift;

            }

            return coeff;

        }


    }
}
