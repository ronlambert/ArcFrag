using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ArcFragTools
{
    public class ArcFrag : IArcFrag
    {
        public static string sVersion = "1.0";
        public static double mdum, md;
        public static double[] ma = new double[2];
        public static int[] mma = new int[56];
        public static int typ, nflag, gtpflag, i, j;
        public static int[] tnmat = new int[3];
        public static string inpfile = "";
        public static string sMSD_output = "MSDout.csv";
        public static double vrel;
        public static string sMSDcsv = "Mass_(kg),Length_(m),Width_(m),Thickness_(m)";
        public static double[,] rho = new double[3, 11];
        public static double[,] trho = new double[3, 11];
        public static double[] fnum = new double[101];
        public static double en = 0.0;
        public static double[] masses = new double[101];
        public static double[] size = new double[101];
        public static double[,] fmat = new double[3, 11];
        public static double[,] tfmat = new double[3, 11];
        public static int ntot, velflag, rotflag;
        public static int imax, sizflag, outflag, ncount;
        public static int imin = 19;
        public static int[] nmat = new int[3];
        public static double eninex;
        public static double minsize;
        public static int maxnum, totnum, n;
        public static double eacc;
        public static double alpha, beta, mean;
        public static double seed;
        public static int inext;
        public static int inextp;
        public static int onlyin;
        public static double press;
        public static double vol;
        public static double e1con;
        public static double e1, e2, com;
        public static int[,] npart = new int[imin + 1, 21];
        public static double tappdiv;
        public static double[] bsize = new double[101];
        public static double t1;
        public static double[] t = new double[101];
        public static double[] length = new double[11];
        public static double[] width = new double[11];
        public static double[] th = new double[11];
        public static int[] indind = new int[11];
        public static double[] eheat = new double[3];
        public static double[] espread = new double[3];
        public static double[] eadd = new double[3];
        public static double f, fheat, fspread;
        public static double[] vpeak = new double[101];
        public static double[] vchar = new double[101];
        public static double ke, exp;
        public static double[,] vsp = new double[41, 101];
        public static double[,] dfnum = new double[41, 101];
        public static double[] p = new double[101];
        public static int[] rvtot2 = new int[101];
        public static int[,] rvmat = new int[31, 101];
        public static double[,] vind = new double[11, 4];
        public static double[,] postvel = new double[3, 4];
        public static double[,] vp = new double[2252, 4];
        public static int[] rvnum = new int[101];
        public static int rvnbin;
        public static int iflag;
        public static int[] nrho = new int[11];
        public static double[,,] rvsize = new double[3, 101, 4];
        public static double[] tapp = new double[11];
        public static double[] pp = new double[502];
        public static double x1, x2;
        public static double[,] rind = new double[11, 5];
        public static int[] rindind = new int[11];
        public static int iff = 0;

        public bool MSD()
        {
            param();
            rndstrt();
            getmass();
            expl();
            zero();
            return true;
        }
        public bool MSD(string sMSD_input, string sMSD_output_specified)
        {
            inpfile = sMSD_input;
            sMSD_output = sMSD_output_specified;
            bool bTest = input();
            if (bTest == false)
            {
                Console.WriteLine("\nMissing or opened MSD.in file.");
                return false;
            }
            return MSD();
        }
        public bool MSD(double dInert_mass_kg, int iNumMaterials,
            double[] dMaterialFraction, double[] dMaterialDensity_kg_per_m3,
            double dBreakupEnergy_J, string sMSD_output_specified)
        {
            sMSD_output = sMSD_output_specified;
            ma[1] = dInert_mass_kg;
            nmat[1] = iNumMaterials;
            eninex = dBreakupEnergy_J;
            fmat[1, 1] = dMaterialFraction[0];
            rho[1, 1] = dMaterialDensity_kg_per_m3[0];
            for (int i = 0; i < iNumMaterials; i++)
            {
                fmat[1, i + 1] = dMaterialFraction[i];
                rho[1, i + 1] = dMaterialDensity_kg_per_m3[i];
                if (fmat[1, i + 1] == 0.0 || fmat[1, i + 1] == double.NaN
                  || rho[1, i + 1] == 0.0 || rho[1, i + 1] == double.NaN)
                {
                    return false;
                }
            }
            return MSD();
        }

        // Read in  data from the input file in either old MSD format or new comma deliminted format
        private static bool input()
        {
            System.IO.StreamReader inputReader;
            try
            {
                inputReader = System.IO.File.OpenText(inpfile);
            }
            catch
            {
                return false;
            }
            string sLine = inputReader.ReadLine();                  // Read first line
            if (sLine.Contains("MASS OF TARGET (KG), NUM"))         // Old MSD format has first line as header discarded
            {
                char[] separators = new char[] { ' ' };
                sLine = inputReader.ReadLine();                     // Read mass, number of materials, fractions of naterials
                string[] values = sLine.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                ma[1] = Convert.ToDouble(values[0].Trim());
                nmat[1] = Convert.ToInt32(values[1].Trim());
                for (int i = 2; i < nmat[1] + 2; i++)
                {
                    fmat[1, i - 1] = Convert.ToDouble(values[i].Trim());
                }
                sLine = inputReader.ReadLine();                     // Read densities of naterials
                values = sLine.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 1; i < nmat[1] + 1; i++)
                {
                    rho[1, i] = Convert.ToDouble(values[i - 1].Trim());
                }
                sLine = inputReader.ReadLine();                     // Read and discard second header line
                sLine = inputReader.ReadLine();                     // Read energy of breakup
                values = sLine.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                eninex = Convert.ToDouble(values[0].Trim());
            }
            else                                                    // New comma delimited input format
            {
                string[] values = sLine.Split(',');                 // Read inert mass to be broken up
                ma[1] = Convert.ToDouble(values[0].Trim());
                sLine = inputReader.ReadLine();                     // Read number of materials
                values = sLine.Split(',');
                nmat[1] = Convert.ToInt32(values[0].Trim());
                for (int i = 0; i < nmat[1]; i++)
                {
                    sLine = inputReader.ReadLine();                 // Read material fractions and densities
                    values = sLine.Split(',');
                    fmat[1, i + 1] = Convert.ToDouble(values[0].Trim());
                    rho[1, i + 1] = Convert.ToDouble(values[1].Trim());
                }
                sLine = inputReader.ReadLine();                     // Read energy of breakup
                values = sLine.Split(',');
                eninex = Convert.ToDouble(values[0].Trim());
            }
            inputReader.Close();
            return true;
        }

        // Set many of the parameters used
        private static void param()
        {
            nflag = 0;
            ntot = 0;
            typ = 2;                        // Target type, 1=booster, 2=pbv, 3=satellite
            velflag = 0;
            sizflag = 1;
            rotflag = 0;
            minsize = 0.00010;
            maxnum = 300000;
            totnum = 600000;
            onlyin = 0;
            minsize = 0.00010;
            maxnum = 300000;
            totnum = 600000;
            n = 100;
            seed = 8.42149;                 // Seed for the random number generator
            eacc = 0.05;                    // The accuracy used when conserving 
            vrel = 100000.0;                //   spreading energy (MAY BE DELETED)
            alpha = 1.5;                    // Distribution parameters for the beta distribution 
            beta = 4.5;                     // used to derive tapp from t
            mean = alpha / (alpha + beta);  // mean value of the beta distribution. 
        }


        // Initialize the random number generator using the seed that is provided
        // The routine calls the random number generator a certain number of times
        // depending on the seed number.Once this initialization has occurred future
        // calls can be made directly to the uniform random number generator, 'random'. 
        private static void rndstrt()
        {
            int idum = -Math.Abs((int)(1000.0 * seed));
            double x = ran3(idum);
            x = x + x;
        }

        // Portable randon number uniform random number generator from Press and Flannery,"Numerical Recepies"                                                *****
        private static double ran3(int idum)
        {
            int i, ii, k, mk;
            int mbig = 1000000000;
            int mseed = 161803398;
            int mz = 0;
            int mj = 1;
            double fac = 1.0 / Convert.ToDouble(mbig);

            if (idum < 0 || iff == 0)
            {
                iff = 1;
                mj = mseed - Math.Abs(idum);
                mj = mj % mbig;
                mma[55] = mj;
                mk = 1;
                for (i = 1; i < 55; i++)
                {
                    ii = (21 * i) % 55;
                    mma[ii] = mk;
                    mk = mj - mk;
                    if (mk < mz)
                    {
                        mk = mk + mbig;
                    }
                    mj = mma[ii];
                }
                for (k = 1; k < 5; k++)
                {
                    for (i = 1; i < 56; i++)
                    {
                        mma[i] = mma[i] - mma[1 + ((i + 30) % 55)];
                        if (mma[i] < mz)
                        {
                            mma[i] = mma[i] + mbig;
                        }
                    }
                }
                inext = 0;
                inextp = 31;
                idum = 1;
            }
            inext++;
            if (inext == 56)
            {
                inext = 1;
            }
            inextp++;
            if (inextp == 56)
            {
                inextp = 1;
            }
            mj = mma[inext] - mma[inextp];
            if (mj < mz)
            {
                mj += mbig;
            }
            mma[inext] = mj; // Convert.ToDouble(mj);
            double ran3 = Convert.ToDouble(mj) * fac;
            return ran3;
        }

        // Initialize the mass array
        private static void getmass()
        {
            masses[1] = 10000.0;
            masses[2] = 5000.0;
            masses[3] = 2000.0;
            masses[4] = 1000.0;
            masses[5] = 500.0;
            masses[6] = 200.0;
            masses[7] = 100.0;
            masses[8] = 50.0;
            masses[9] = 20.0;
            masses[10] = 10.0;
            masses[11] = 5.0;
            masses[12] = 2.0;
            masses[13] = 1.0;
            masses[14] = 0.5;
            masses[15] = 0.2;
            masses[16] = 0.1;
            masses[17] = 0.05;
            masses[18] = 0.02;
            masses[19] = 0.01;
            masses[20] = 0.005;
            masses[21] = 0.002;
            masses[22] = 0.001;
            masses[23] = 0.0005;
            masses[24] = 0.0002;
            masses[25] = 0.0001;
            masses[26] = 0.00005;
            masses[27] = 0.00002;
            masses[28] = 0.00001;
            masses[29] = 0.000005;
            masses[30] = 0.000002;
            masses[31] = 0.000001;
        }

        // Coordinate the processing of explosion
        private static void expl()
        {
            transfer();
            getexe(ma[1]);
            pbv(1);
        }

        // Move the information stored in the variables 'nmat', 'rho', and 'fmat' into the new 
        // variable names and is called by the explosion routines.  This is done in place of routine 'sdmmat' 
        private static void transfer()
        {
            tnmat[1] = nmat[1];
            for (int i = 1; i < nmat[1] + 1; i++)
            {
                trho[1, i] = rho[1, i];
                tfmat[1, i] = fmat[1, i];
            }
        }

        // Calculate the energy released in a pressure induced explosion, given
        // the initial volume and presure of the gas just before the explosion
        private static void ptoen()
        {
            double gamma = 1.4;
            eninex = (press * vol) / (gamma - 1.0);
        }

        // determine the exponent used in the exponential law. 
        // It is an empiricle function of the energy available per unit mass of the target.
        private static void getexe(double maa)
        {
            e1con = 0.0000102;
            com = 0.0;
            e2 = 0.0;
            if (eadd[1] >= 0.0)
            {
                e1 = e1con * Math.Sqrt(eninex / maa);
            }
            else
            {
                e1 = Math.Abs(eadd[1]);
                eadd[1] = 0.0;
            }
        }

        // Coordinate the fragmentation modeling for a Post-Boost Vehicle or a similar structure
        private static void pbv(int iin)
        {
            mdum = ma[iin];
            massdist(typ, mdum);
            getimax();
            mksizes(iin);
            vspread(en);
            ntot = ggntot();
            sMSDcsv = "Mass_(kg),Length_(m),Width_(m),Thickness_(m)\n";
            getvel(iin, ntot);
        }

        // Produce a mass distribution with a total mass within 'cutoff' of the given mass, 'md'
        // The flag, 'flag', indicates which type of mass distribution is to be used
        // flag=0 indicates that a power law based distribution should be used
        // flag=1 indicated that an exponential law distribution should be used
        private static void massdist(int flag, double maa)
        {
            double cutoff, cmin, cmax, cmid, massmin;
            double massmid = 0.0;
            double massmax = 0.0;
            md = maa;
            cutoff = 0.1 * md;
            cmin = 0.0;
            cmax = 3.0;
            massmin = md;
            massmax = emtot(cmax);

        // Iterate to get coefficient value cmid of distribution whose total mass
        // is within cutoff (0.1 * md) of md
        label_10:
            cmid = (cmax + cmin) / 2.0;
            massmid = emtot(cmid);
            if (Math.Abs(massmid) > cutoff
                && Math.Abs(cmax - cmid) > 0.000001
                && Math.Abs(cmid - cmin) > 0.000001)
            {
                if (massmin * massmid > 0.0)
                {
                    cmin = cmid;
                    massmin = massmid;
                }
                else
                {
                    cmax = cmid;
                    massmax = massmid;
                }
                goto label_10;
            }
            mnmbal(maa);
        }


        // Compute total mass in a mass distribution using an exponential law
        // and a scaling constant 'c'
        //  inputs:
        //    c - scaling constant for cumulative fragments equation
        //     com - always zero?
        //     e1 - exponent constant for cumulative fragments equation
        //     md - vehicle total mass minus large fragments and rv's
        //     masses(100) - mass(i) is mass in kg.Of fragments in bin i
        //  
        //  outputs:
        //      fnum(100) - mass distribution: fnum(i) is number of frags of  mass masses(i)
        //      imin  - index of minimum mass
        //      imax  - index of maximum mass
        //      dmass - differenct between desired mass, md, and total in distibution
        //  
        //   Intermediate variables:
        //      nn(100) - cumulative number of fragments less than masses(i)
        private static double emtot(double c)
        {
            double dmass, mmax, mtotal;
            double[] nn = new double[101];
            double[] bb = new double[101];

            int num;

            // Compute maximum fragment mass(mmax)
            md = md * 1000.0;
            mmax = Math.Pow(-1.0 / e1 * Math.Log(2.0 / (c * md)), 2.0);
            mmax = mmax / 1000.0;
            md = md / 1000.0;
            num = 1;

            // Compute index of maximum fragment mass (imax)
            while (masses[num] > mmax)
            {
                num++;
            }
            imax = num;

            // Compute cumulative number of fragments whose mass is less than bb(I), nn(I)
            bb[num] = mmax;
            if (bb[num] > com / 1000.0)
            {
                nn[imax] = c * md * 1000.0 * Math.Exp(-e1 * Math.Sqrt(bb[imax] * 1000.0));
            }
            else
            {
                nn[imax] = c * Math.Exp(Math.Sqrt(com) * (e2 - e1)) * md * Math.Exp(-e2 * Math.Sqrt(bb[imax] * 1000.0));
                md = md / 1000.0;
            }
            for (int i = imax + 1; i < imin + 2; i++)
            {
                bb[i] = (masses[i - 1] + masses[i]) / 2.0;
            }
            for (int i = imax; i < imin + 2; i++)
            {
                if (bb[i] > com / 1000.0)
                {
                    nn[i] = c * md * 1000.0 * Math.Exp(-e1 * Math.Sqrt(bb[i] * 1000.0));
                }
                else
                {
                    md *= 1000.0;
                    nn[i] = c * Math.Exp(Math.Sqrt(com) * (e2 - e1)) * md * Math.Exp(-e2 * Math.Sqrt(bb[i] * 1000.0));
                    md /= 1000.0;
                }
            }

            // Set number of fragments in bins of mass greater than max mass equal to zero
            mtotal = 0.0;
            for (int i = 1; i < imax + 1; i++)
            {
                fnum[i] = 0.0;
            }

            // Compute number of fragments in each bin
            for (i = imax; i < imin + 1; i++)
            {
                fnum[i] = Convert.ToDouble(Convert.ToInt32(nn[i + 1] - nn[i]));
                mtotal += masses[i] * fnum[i];
            }
            dmass = md - mtotal;
            return dmass;
        }

        // Balance the mass and the momentum in the given mass distribution
        // Ensure that the total mass in the distribution is not greater than
        // that all fragment mass bins have even numbers of fragments
        // This makes it possible to produce paired sets of spread velocities
        // whose net momentums are zero
        private static void mnmbal(double mdd)
        {
            int iter, scale, flagg2, num;
            double per, dmm, mtot;

            per = 0.003;
            iter = 0;

        label_5:
            mtot = 0.0;
            for (int i = imax; i < imin + 1; i++)
            {
                mtot += masses[i] * fnum[i];
            }
            dmm = mdd - mtot;

            // Avoid wiping out heaviest fragment
            if (Math.Abs(dmm) < fnum[imax] * masses[imax] || iter > 10)
            {
                goto label_205;
            }
            iter++;
            for (int i = imax; i < imin + 1; i++)
            {
                fnum[i] = Convert.ToInt32(fnum[i] * mdd / mtot);
            }
            goto label_5;

        label_205:
            for (int i = imax; i < imin + 1; i++)
            {
                if (Convert.ToInt32(fnum[i]) % 2 == 1)
                {
                    if (dmm < 0.0)
                    {
                        fnum[i]--;
                        dmm += masses[i];
                    }
                    else
                    {
                        if (masses[i] < dmm)
                        {
                            fnum[i]++;
                            dmm -= masses[i];
                        }
                        else
                        {
                            fnum[i]--;
                            fnum[i + 1] += 2.0;
                            if ((i + 1) % 3 == 0)
                            {
                                fnum[i + 2]++;
                            }
                        }
                    }
                }
                else
                {
                    if (2.0 * masses[i] < dmm)
                    {
                        fnum[i] += 2.0;
                        dmm -= 2.0 * masses[i];
                    }
                }
            }

            //  Ensure mass in the distribution is not greater than total pre-breakup mass
            if (dmm < 0.0)
            {
                i = 1;
            label_220:

                if (2.0 * masses[i] > Math.Abs(dmm)
                  && 2.0 * masses[i + 1] < Math.Abs(dmm))
                {
                    flagg2 = i;
                    if (fnum[i] > 0.5)
                    {
                        fnum[i] -= 2.0;
                        dmm += 2.0 * masses[i];
                    }
                }
                else
                {
                    i++;
                    goto label_220;
                }

                if (dmm < 0.0)
                {
                    if (Math.Abs(dmm) < 2.0 * masses[imax])
                    {
                        i = 1;

                    label_230:
                        if (fnum[flagg2 - i] > 0.5)
                        {
                            fnum[flagg2 - i] -= 2.0;
                            dmm += 2.0 * masses[flagg2 - i];
                        }
                        else
                        {
                            i++;
                            goto label_230;
                        }
                    }
                    else
                    {
                    label_240:
                        i++;
                        if (fnum[i] > 0.5)
                        {
                            if (fnum[i] * masses[i] > Math.Abs(dmm))
                            {
                                scale = (int)(Math.Abs(dmm) / masses[i]);
                                if (scale % 2 == 1)
                                {
                                    scale++;
                                }
                                else
                                {
                                    scale += 2;
                                }
                                fnum[i] -= Convert.ToDouble(scale);
                                dmm += Convert.ToDouble(scale) * masses[i];
                            }
                            else
                            {
                                dmm += fnum[i] * masses[i];
                                fnum[i] = 0.0;
                                goto label_240;
                            }
                            goto label_205;
                        }
                        else
                        {
                            i++;
                            goto label_240;
                        }
                    }
                }
            }

        // Make the final mass conservation corrections
        label_300:
            if (dmm > per * mdd && per * mdd > 2.0 * masses[imin])
            {
                num = imax;
                while (2.0 * masses[num] > dmm)
                {
                    num++;
                }
                dmm -= 2.0 * masses[num];
                fnum[num] += 2.0;
                goto label_300;
            }
            while (fnum[imax] < 0.1)
            {
                imax++;
            }
        }

        // Redetermin 'iamx' after RV's or extra large fragments are addd to distribution
        private static void getimax()
        {
            int i = 0;
        label_10:
            i++;
            if (i < imin)
            {
                if (fnum[i] < 0.1)
                {
                    goto label_10;
                }
                else
                {
                    imax = i;
                }
            }
            else
            {
                imax = imin;
            }
        }


        // Coordinate materials distributions and mass/size/shape relationships for fragments
        private static void mksizes(int iin)
        {
            ncount = 1;
            matdis(iin);
            getsize(iin);
        }

        // Determine the material distribution in each mass bin
        // Designate the number of fragments in each bin to be a material type 
        private static void matdis(int iin)
        {
            double mt = 0.0;
            double mtotal, exsum;
            double[] mmass = new double[21];
            int[] mnum = new int[21];
            int[,] exmat = new int[101, 21];
            int count, npar, max;
            for (int i = 1; i < imin + 1; i++)
            {
                for (int j = 1; j < 20 + 1; j++)
                {
                    npart[i, j] = 0;
                }
            }

            for (int i = imax; i < imin + 1; i++)
            {
                mt += masses[i] * fnum[i];
            }

            for (int i = 1; i < tnmat[iin] + 1; i++)
            {
                mmass[i] = tfmat[iin, i] * mt;
            }

            for (int i = imax; i < imin + 1; i++)
            {
                count = 0;
                mtotal = 0.0;
                if (fnum[i] < 0.5)
                {
                    continue;
                }
                npar = Convert.ToInt32(fnum[i]);
                for (int j = 1; j < tnmat[iin] + 1; j++)
                {
                    if (nflag != 1)
                    {
                        if (masses[i] <= mmass[j])
                        {
                            count++;
                            mnum[count] = j;
                            mtotal += mmass[j];
                        }
                    }
                    else
                    {
                        if (Convert.ToInt32(fnum[i]) > 0)
                        {
                            count++;
                            mnum[count] = j;
                            mtotal += mmass[j];
                        }
                    }

                }
                exsum = 0.0;
                for (int k = 1; k < tnmat[iin] + 1; k++)
                {
                    exsum += Convert.ToDouble(exmat[i, k]);
                }

                for (j = 1; j < count + 1; j++)
                {
                    npart[i, mnum[j]] = (int)((mmass[mnum[j]] / mtotal)
                                      * Convert.ToInt32((fnum[i]) - exsum));
                    mmass[mnum[j]] -= masses[i] * Convert.ToDouble(npart[i, mnum[j]]);
                    npar -= npart[i, mnum[j]];
                }

                for (int j = 1; j < tnmat[iin] + 1; j++)
                {
                    npart[i, j] += exmat[i, j];
                    npar -= exmat[i, j];
                }

                if (npar > 0)
                {
                    for (int k = 1; k < npar + 1; k++)
                    {
                        max = 1;
                        for (int j = 1; j < tnmat[iin] + 1; j++)
                        {
                            if (mmass[max] < mmass[j])
                            {
                                max = j;
                            }
                        }
                        npart[i, max] += 1;
                        mmass[max] -= masses[i];
                    }
                }
            }
        }

        // Calculate the characteristic size for each mass bin
        // which represents the average apparent size of the most likley 
        // fragment averaged over the two largest of its dimensions
        private static void getsize(int iin)
        {
            t1 = 0.001;
            double woverl, bound, ts;
            double llength, wwidth, tarea;
            double[] rhoave = new double[101];
            double[] mind = new double[101];
            double[] area = new double[101];
            int c;
            tappdiv = 0.1;
            for (int i = 1; i < imin + 1; i++)
            {
                rhoave[i] = 0;
                if (Convert.ToInt32(fnum[i]) > 0)
                {
                    for (int j = 1; j < tnmat[iin] + 1; j++)
                    {
                        rhoave[i] += Convert.ToDouble(npart[i, j]) * trho[iin, j];
                    }
                    rhoave[i] /= fnum[i];
                }
                else
                {
                    if (i == 1)
                    {
                        rhoave[i] = trho[iin, 1];
                    }
                    else
                    {
                        rhoave[i] = rhoave[i - 1];
                    }
                }
                mind[i] = 0.5 * Math.Log10(masses[i] / (25.0 * rhoave[i] * Math.Pow(t1, 3.0))) + 1.0;
            }

            c = 1;
        label_7:
            if ((int)(mind[c]) != (int)(mind[c + 1]))
            {
                if ((mind[c] - (int)(mind[c])) >
                    ((int)(mind[c + 1]) + 1 - mind[c + 1]))
                {
                    mind[c] = Convert.ToDouble((int)(mind[c])) - 0.33;
                    mind[c + 1] = Convert.ToDouble((int)(mind[c + 1])) + 0.33;
                }
                else
                {
                    mind[c + 1] = Convert.ToDouble((int)(mind[c + 1])) + 0.33;
                    mind[c] = Convert.ToDouble((int)(mind[c])) - 0.33;
                }
                c += 2;
            }
            else
            {
                mind[c] = Convert.ToDouble((int)(mind[c]));
                c++;
            }
            if (c < imin)
            {
                goto label_7;
            }
            else
            {
                mind[c] = Convert.ToDouble((int)(mind[c]));
            }
            for (int i = 1; i < imin + 1; i++)
            {
                if (mind[i] > 3.0)
                {
                    mind[i] = 3.0;
                }
                if ((int)(mind[i]) == 0)
                {
                    t[i] = t1;
                }
                else
                {
                    t[i] = t1 * Math.Pow(100.0, (mind[i] - 1.0) / 3.0);
                }
                area[i] = masses[i] / (rhoave[i] * t[i]);
                woverl = 0.62;
                llength = Math.Sqrt(area[i] / woverl);
                wwidth = area[i] / llength;
                size[i] = (llength + wwidth) / 2.0;

                // Calculate the smallest size in a mass bin
                if (i != imin)
                {
                    bound = (masses[i] + masses[i + 1]) / 2.0;
                }
                else
                {
                    bound = (masses[i] * 1.5 / 2.0);
                }
                if (Convert.ToInt32((mind[i])) == 0)
                {
                    ts = t1;
                }
                else
                {
                    ts = t[i];
                }
                tarea = bound / (rhoave[i] * ts);
                woverl = 0.62;
                llength = Math.Sqrt(tarea / woverl);
                wwidth = tarea / llength;
                bsize[i] = (llength + wwidth) / 2.0;
            }
        }

        private static void output(int iin, int ntot, int outflag, int ncount)
        {
            for (int i = 1; i < ncount + 1; i++)
            {
                sMSDcsv += masses[indind[i]].ToString("0.000000") + ","
                      + length[i].ToString("0.000000") + ","
                      + width[i].ToString("0.000000") + ","
                      + th[i].ToString("0.000000") + "\n";

            }
            System.IO.File.WriteAllText(sMSD_output, sMSDcsv);
        }

        private static void endist()
        {
            eheat[1] = (fheat / f) * eninex;
            espread[1] = (fspread / f) * eninex;
        }

        private static void getf()
        {
            f = fheat + fspread;
        }

        private static void charvel(int iin, double ked)
        {
            double sme, m1, vs, min, mid, max, ratio, vrelf, ve, cvtol, maxrat;
            float dkemax, dkemid, dkemin, ket;
            int maxind;

            for (int i = imax; i < imin + 1; i++)
            {
                if (typ == 2)
                {
                    sme = espread[iin] / fspread;
                }
                else
                {
                    sme = 2.0 * (eheat[1] + espread[1] + eadd[1]);
                    vrel = Math.Sqrt(2.0 * 10.0 * sme / (0.01 * ma[1]));
                }
                double Exp1 = Math.Pow(sme, (1.0 / 3.0));
                double Log1 = Math.Log10(size[i] * 8.010488E+08 / Exp1);
                double Exp2 = Math.Pow(Log1, 2.0);
                double Exp3 = Math.Pow(10.0, -0.12478 - 0.067599 * Exp2);
                vpeak[i] = vrel * Exp3;
                // vpeak(i) = vrel * 10.**(-0.12478 - 0.067599 * ( log10( (size(i) * 8.010488e+8) / sme**(1./3.)) )**2. )

                vpeak[i] *= 2.0 / 3.0;
            }

            // Calculate characteristic velocity
            exp = 4.0;
            ke = ked;
            m1 = masses[imax];
            vrelf = 1.0;
            cvtol = 0.01;
            vs = vrelf * vrel;
            min = 0.0;
            max = vs;

            // Adjust the KE distribution to within 5% of the actual KE if necessary
            dkemax = (float)getdke(max, vs);
            dkemin = (float)getdke(min, vs);
        label_21:
            mid = (max + min) / 2.0;
            dkemid = (float)getdke(mid, vs);
            if (Math.Abs(dkemid / ke) > cvtol)
            {
                if (dkemin * dkemid > 0.0)
                {
                    dkemin = dkemid;
                    min = mid;
                }
                else
                {
                    dkemax = dkemid;
                    max = mid;
                }
                goto label_21;
            }
            vchar[imax] = mid;
            ket = (float)0.0;
            for (int i = imax + 1; i < imin + 1; i++)
            {
                vchar[i] = 0.0;
            }

            // check to verify that vpeak is not too large relative to vchar
            maxrat = 0.0;
            for (int i = imax; i < imin + 1; i++)
            {
                ratio = vpeak[i] / vchar[i];
                if (ratio > maxrat)
                {
                    maxrat = ratio;
                    maxind = i;
                }
            }
            if (maxrat > Math.Sqrt(3.0))
            {
                for (int i = imax; i < imin + 1; i++)
                {
                    if (vpeak[i] / vchar[i] > Math.Sqrt(3.0))
                    {
                        vpeak[i] = vchar[i] * (Math.Sqrt(3.0) - 1.0);
                    }
                }
            }
            for (int i = imax; i < imin + 1; i++)
            {
                if (vpeak[i] / Math.Sqrt(3.0) * vchar[i] < 0.08)
                {
                    vpeak[i] = 0.08 * Math.Sqrt(3.0) * vchar[i];
                }
            }
        }


        // Calculates the total kinetic energy associated with a characteristic 
        // velocity distribution given the velocity of the largest fragment, v1, 
        // and the associated fragment distribution, fnum[].
        // The output, 'dke', is the difference between the actual kinetic energy, 
        // 'ke', and the kinetic energy in the distribution.
        private static double getdke(double v1, double vs)
        {
            double ket, vt, dke;
            ket = 0.5 * fnum[imax] * masses[imax] * Math.Pow(v1, 2.0);
            for (int i = imax + 1; i < imin + 1; i++)
            {
                vt = Math.Pow(masses[imax] / masses[i], 1.0 / exp) * v1;
                if (vt > vs)
                {
                    vt = vs;
                }
                ket += 0.5 * fnum[i] * masses[i] * Math.Pow(vt, 2.0);
            }
            dke = ke - ket;
            return dke;
        }

        // dds in additional energy, if present, to both the heat/light energy and spread energy
        private static void geten(int iin, double en)
        {
            if (eadd[iin] >= 0.0)
            {
                if (fspread > 0.0)
                {
                    en = espread[iin] + eadd[iin] * (fspread / f);
                }
                else
                {
                    if (fheat > 0.0)
                    {
                        en = espread[iin];
                    }
                    else
                    {
                        en = 0.333 * eadd[iin];
                        eheat[iin] = 0.667 * eadd[iin];
                    }
                }
            }
            if (fheat > 0.0)
            {
                eheat[iin] += eadd[iin] * (fheat / f);
            }
        }

        // Distributes the spread velocities in each mass bin using a beta function.
        // This distribution conserves mass, energy and momentum by pairing fragments.
        // All possible spread velocities are stored in array 'vsp' and the number of
        // fragments at each spread velocity are stored in the array 'dfnum'.
        private static void vspread(double ked)
        {
            double ketot, Charea, BalCoef, d, nvpeak, nvchar, x, fac, dket, dke;
            double[] keref = new double[101];
            double[] kedist = new double[101];
            double[] pr = new double[1000];
            int ind;

            ke = ked;

            // Use beta distributions in each size bin
            ketot = 0.0;
            Console.WriteLine("\n       Generating fragments....\n");
            Console.WriteLine("\n       LYFG will continue running after this.\n");
            for (int i = imax; i < imin + 1; i++)
            {
                keref[i] = 0.5 * fnum[i] * masses[i] * Math.Pow(vchar[i], 2.0);
                if (fnum[i] > -2.0)
                {
                    Charea = Math.Pow(size[i], 2.0) / 1.0287 * 3.28084 * 3.28084;
                    BalCoef = masses[i] * 2.20462262 / 0.75 / Charea;
                }
                if (fnum[i] < 0.001)
                {
                    continue;
                }
                if (vpeak[i] < vchar[i])
                {
                    d = (3.0 * Math.Sin(vpeak[i] / vchar[i] * 3.14159265) + Math.Sqrt(3.0)) * vchar[i];
                }
                else
                {
                    d = (-(Math.Sqrt(3.0) - vpeak[i] / vchar[i]) * Math.Sin((vpeak[i] / vchar[i] - 1.0)
                          * 3.14159265 / Math.Sqrt(3.0)) + Math.Sqrt(3.0)) * vchar[i];
                }
                nvpeak = vpeak[i] / d;
                nvchar = vchar[i] / d;
                for (int j = 1; j < n + 1; j++)
                {
                    vsp[i, j] = d / (2.0 * Convert.ToDouble(n)) * (2.0 * Convert.ToDouble(j) - 1.0);
                    dfnum[i, j] = 0.0;
                }
                getbeta(nvpeak, nvchar);

                // Generate the distributions from the probabilities in p[i]
                fac = 2.0 * (1.0 + Convert.ToDouble((int)(fnum[i] / 10000.0)));
                for (int j = 1; j < (int)(fnum[i] / fac) + 1; j++)
                {
                    x = ran3(1);
                    ind = 0;
                label_1010:
                    if (p[ind] >= x)
                    {
                        dfnum[i, ind] += fac;
                    }
                    else
                    {
                        ind++;
                        goto label_1010;
                    }
                }

                // Add on remainder fragments when the fac reduction is used
                if (fac > 2.0)
                {
                    x = ran3(1);
                    ind = 1;
                label_1011:
                    if (x >= p[ind] && x < p[ind + 1])
                    {
                        dfnum[i, ind] += fnum[i] - fac * Convert.ToDouble((int)(fnum[i]
                            / fac));
                    }
                    else
                    {
                        ind++;
                        goto label_1011;
                    }
                }

                dket = 0.0;
                for (int j = 1; j < n + 1; j++)
                {
                    if (dfnum[i, j] > 0.1)
                    {
                        ketot += 0.5 * dfnum[i, j] * masses[i] * Math.Pow(vsp[i, j], 2.0);
                        dket += 0.5 * dfnum[i, j] * masses[i] * Math.Pow(vsp[i, j], 2.0);
                    }
                }
                kedist[i] = dket;
            }
            dke = ke - ketot;
            if (Math.Abs(dke / ke) > eacc)
            {
                ebal(); // ebal(keref, dke, kedist);  // This should never be called
            }
        }

        // Calculate the appropriate beta distribution and cumulative probability distribution
        private static void getbeta(double vpeak, double vchar)
        {
            double d, v, sum, vsum, w, aa, a, b, diff, amin, amax, bb;
            d = 1.0;
            amax = 30.0;
            amin = 1.00001;
        label_100:
            a = (amax + amin) / 2.0;
            bb = ((a - 2.0) * d - (a - 2.0)); // * Convert.ToDouble(vpeak) / Convert.ToDouble(vpeak));
            if (bb > 100.0)
            {
                bb = 100.0;
            }
            sum = 0.0;
            w = 100.0;
            aa = 1.0;
            for (int i = 1; i < Convert.ToInt32(w) + 1; i++)
            {
                v = (Convert.ToDouble(i) - 0.5) * d / w;
                b = aa * Math.Pow(v, (a - 1.0)) * Math.Pow((1.0 - v), (bb - 1.0));
                sum += 1.0 / w + b;
            }
            aa = 1.0 / sum;
            sum = 0.0;
            vsum = 0.0;
            for (int i = 1; i < Convert.ToInt32(w) + 1; i++)
            {
                v = (Convert.ToDouble(i) - 0.5) * d / w;
                b = aa * Math.Pow(v, (a - 1.0)) * Math.Pow((1.0 - v), (bb - 1.0));
                sum += 1.0 / w + b;
                vsum += 1.0 / w * d * Math.Pow(v, 2.0);
            }
            diff = vsum - Convert.ToDouble(Math.Pow(vchar, 2.0));
            if (Math.Abs(diff) > 0.001 * Math.Pow(Convert.ToDouble(vchar), 2.0)
            && (amax - amin) / amin > 1.0E-08)
            {
                if (vpeak < vchar)
                {
                    if (Convert.ToDouble(vchar) < d * Math.Sqrt(1.0 / 3.0))
                    {
                        if (diff > 0.0)
                        {
                            amin = a;
                        }
                        else
                        {
                            amax = a;
                        }
                    }
                    else
                    {
                        if (diff > 0.0)
                        {
                            amax = a;
                        }
                        else
                        {
                            amin = a;
                        }
                    }
                }
                else
                {
                    if (diff > 0.0)
                    {
                        amax = a;
                    }
                    else
                    {
                        amin = a;
                    }
                }
                goto label_100;
            }
            p[1] = 0.0;
            sum = 0.0;
            for (int i = 1; i < Convert.ToInt32(w) + 1; i++)
            {
                v = (Convert.ToDouble(i) - 0.5) * d / w;
                b = aa * Math.Pow(v, (a - 1.0)) * Math.Pow((1.0 - v), (bb - 1.0));
                sum += 1.0 / w + b;
                if (sum > 1.0E-30)
                {
                    p[i] = sum;
                }
                else
                {
                    p[i] = 0.0;
                }
            }
        }

        private static void ebal() // If this is called, there is an error
        {
            Console.WriteLine("\n Error in balancing kinetic energy.\n");
            Console.WriteLine("\n Press any key to exit.\n");
            Console.ReadKey();
            Environment.Exit(0);
        }

        // produces a set of velocity vectors in ECI for all fragments meeting certain criteria.
        // These include a minimum fragment size cutoff, a maximum total number of fragments,
        // and a maximum number of fragments in a mass bin.
        private static void getvel(int iin, int ntot)
        {

            double rys, mu, var, temp, x, ftapp, dens, area; //s, 
            double woverl, rrho, wave, sigma, kk, s;
            double omega = 0.0;
            double sprev = 0.0;
            double[] d = new double[4];
            int nbin, mb, nn;//, totnum ;
            int inbin = 1;
            int ndf = 0;
            int indr = 1;
            int ptotal = 0;
            int nrvs, fflag, rvpro, ind;
            int[] np = new int[20];
            rys = 0.001 * 2.5e+8;
            rvnbin = 1;
            nrvs = 1;
            nbin = 0;
            ncount = 0;
            rvpro = 0;
            mu = 0.0;
            var = 2.0;
            if (onlyin > 0)
            {
                mb = onlyin;
            }
            else
            {
                mb = imax;
            }
            fflag = 1;
            for (int L = 1; L < Convert.ToInt32(Convert.ToDouble(ntot / 2)) + 1; L++)
            {
                if (nbin == fnum[mb] || fflag == 1)
                {
                    if (fflag == 0)
                    {
                        rvpro += rvtot2[mb];
                        mb += 1;
                        while (fnum[mb] < 0.5)
                        {
                            mb += 1;
                        }
                    }
                    fflag = 0;
                    nbin = 0;
                    inbin = 1;
                    ndf = 0;
                    if (nflag == 1 && iin == 1)
                    {
                        ptotal = 2;
                        np[1] = 2;
                    }
                    else
                    {
                        ptotal = Convert.ToInt32(fnum[mb]) - rvtot2[mb];
                        for (int i = 1; i < tnmat[iin] + 1; i++)
                        {
                            np[i] = npart[mb, i];
                        }
                        for (int i = 1; i < rvtot2[mb] + 1; i++)
                        {
                            np[rvmat[mb, i]] = np[rvmat[mb, i]] - 1;
                        }
                    }

                    // RV setup
                    rvnbin = 1;
                }

            label_30:
                if (ndf == Convert.ToInt32(dfnum[mb, inbin]))
                {
                    inbin++;
                    ndf = 0;
                    goto label_30;
                }
                ndf += 2;

                // Get velocity vectors
                if (velflag > 0)
                {
                    x = ran3(L);
                    ind = (int)(2252 * x) + 1;
                    for (int i = 1; i < 4; i++)
                    {
                        vind[ncount + 1, i] = postvel[iin, i] + vsp[mb, inbin] * vp[ind, i];
                        vind[ncount + 2, i] = postvel[iin, i] - vsp[mb, inbin] * vp[ind, i];
                        indind[ncount + 1] = mb;
                        indind[ncount + 2] = mb;
                    }
                }
                // Get fragment dimensions
                for (int k = 1; k < 3; k++)
                {
                    nbin++;
                    if (rvtot2[mb] > 0 && rvnbin < rvtot2[mb])
                    {
                        i = 1;
                    label_60:
                        if (nbin != rvnum[i])
                        {
                            i++;
                        }
                        if (i <= rvtot2[mb])
                        {
                            goto label_60;
                        }
                        else
                        {
                            nrho[ncount + k] = rvmat[mb, i];
                            if (sizflag == 1)
                            {
                                if (rvsize[iin, nrvs - rvnbin + 1, i] > 0.0001)
                                {
                                    length[ncount + k] = rvsize[iin, nrvs - rvnbin + i, 1];
                                    width[ncount + k] = rvsize[iin, nrvs - rvnbin + i, 2];
                                    th[ncount + k] = rvsize[iin, nrvs - rvnbin + i, 3];
                                    tapp[ncount + k] = rvsize[iin, nrvs - rvnbin + i, 3];
                                    nrvs++;
                                    rvnbin++;
                                    goto label_50;
                                }
                                else
                                {
                                    th[ncount + k] = t[mb];
                                    tapp[ncount + k] = th[ncount + k];
                                }
                            }
                            nrvs++;
                            rvnbin++;
                            goto label_41;
                        }
                    }

                    // Get sizes
                    x = ran3(1);
                    ind = (int)(ptotal * x) + 1;
                    j = 0;
                    nn = 0;
                label_43:
                    j++;
                    nn += np[j];
                    if (ind <= nn)
                    {
                        nrho[ncount + k] = j;
                        ptotal--;
                        np[j]--;
                    }
                    else
                    {
                        goto label_43;
                    }
                label_41:
                    if (sizflag == 1)
                    {
                        if (bsize[mb] > mean * 25.0 * t1)
                        {
                            if (tapp[ncount + k] < 0.000001)
                            {
                                th[ncount + k] = t[mb];

                                // Test code to use two districutions for tapp
                                if (L != 1)
                                {
                                    s = sprev;
                                }
                                else
                                {
                                    s = 1;
                                }
                                if (s > tappdiv)
                                {
                                    iflag = 1;
                                }
                                else
                                {
                                    iflag = 0;
                                }
                                sprev = (length[ncount + k] + width[ncount + k] + tapp[ncount + k]) / 3.0;
                                // End of test code

                                ftapp = gettapp(iflag);
                                tapp[ncount + k] = ftapp * th[ncount + k];
                            }
                            dens = trho[iin, nrho[ncount + k]];
                            area = masses[mb] / (dens * th[ncount + k]);
                            gauss(0.62, 0.24);
                            if (x1 > 1.0)
                            {
                                woverl = 1.0;
                            }
                            else
                            {
                                woverl = x1;
                            }
                            if (woverl < 0.0)
                            {
                                woverl = 0.05;
                            }
                            length[ncount + k] = Math.Sqrt(area / woverl);
                            width[ncount + k] = area / length[ncount + k];
                        }
                        else
                        {
                            if (tapp[ncount + k] < 0.000001)
                            {
                                th[ncount + k] = 2.0 * Math.Pow((3.0 * masses[mb])
                                    / (4.0 * Math.PI * trho[iin, nrho[ncount + k]]), 3.0);
                            }
                            for (int j = 1; j < 4; j++)
                            {
                                gauss(mu, var);
                                if (x1 < 0.0)
                                {
                                    d[j] = th[ncount + k] / (Math.Abs(x1) + 1.0);
                                }
                                else
                                {
                                    d[j] = (x1 + 1.0) * th[ncount + k];
                                }
                            }
                            for (int j = 1; j < 3; j++)
                            {
                                for (int i = 2; i < 4; i++)
                                {
                                    if (d[i] > d[j])
                                    {
                                        temp = d[i];
                                        d[i] = d[j];
                                        d[j] = temp;
                                    }
                                }
                            }
                            length[ncount + k] = d[1];
                            width[ncount + k] = d[2];
                            if (tapp[ncount + k] < 0.000001)
                            {
                                tapp[ncount + k] = d[3];
                            }
                        }
                    }
                label_50:
                    if (velflag == 0)
                    {
                        indind[ncount + k] = mb;
                    }

                    // Get rotation information
                    if (rotflag == 1)
                    {
                        if (k == 1)
                        {
                            rrho = trho[iin, nrho[ncount + k]];
                            wave = Math.Sqrt(rys / Math.Pow((rrho * size[mb]), 2.0));
                            sigma = wave / 2.575;
                            x = ran3(1);
                            indr = (int)(2252.0 * x) + 1;
                        label_70:
                            gauss(wave, sigma);
                            if (x1 < 0.0)
                            {
                                if (x2 < 0.0)
                                {
                                    goto label_70;
                                }
                                else
                                {
                                    omega = x2;
                                }
                            }
                            else
                            {
                                omega = x1;
                            }
                        }
                        if (k == 1)
                        {
                            kk = 1.0;
                        }
                        else
                        {
                            kk = -1.0;
                        }
                        for (int i = 1; i < 4; i++)
                        {
                            rind[ncount + k, i] = kk * vp[indr, i];
                        }
                        rind[ncount + k, 4] = omega;
                        rindind[ncount + k] = mb;
                    }

                }
                ncount += 2;
                if (ncount == 10 && (nflag != 1 || iin == 2))
                {
                    outflag = 3;
                    output(iin, ntot, outflag, ncount);
                    ncount = 0;
                    for (int i = 1; i < 11; i++)
                    {
                        tapp[i] = 0.0;
                    }
                }
            }
            if (ncount > 0 && (nflag != 1 || iin == 2))
            {
                outflag = 3;
                output(iin, ntot, outflag, ncount);
                ncount = 0;
            }
            if (nbin != 0)
            {
                for (int i = 1; i < rvtot2[mb] + 1; i++)
                {
                    if (rvnum[i] <= nbin)
                    {
                        rvpro++;
                    }
                }
            }
        }

        // Generate the factor tapp/t from a beta distribution derived from experimental data
        private static double gettapp(int iflag)
        {
            double max, min, aa, bb, sum, v, x;
            double ftapp = 0;
            int ind;
            if (iflag != 1)
            {
                max = 25.0;
            }
            else
            {
                max = 100.0;
            }
            min = 1.0;
            aa = 11.66916;
            if (gtpflag == 0)
            {
                sum = 0.0;
                for (int i = 1; i < 502; i++)
                {
                    v = Convert.ToDouble((i - 1) / 500.0);
                    if (i == 1 || i == 501)
                    {
                        bb = 0.0;
                        goto label_20;
                    }
                    bb = aa * Math.Pow(v, (alpha - 1)) * Math.Pow((1.0 - v), (beta - 1.0));
                label_20:
                    sum += (1.0 / 500.0) * bb;
                    pp[i] = sum;
                }
                gtpflag = 1;
            }
            x = ran3(1);
            ind = (int)(344.8526 * Math.Pow(((x + 0.91069) / 1.866896), (1.0 / 0.224392)));
            if (ind < 1)
            {
                ind = 1;
            }
            if (ind > 500)
            {
                ind = 500;
            }
            if (x >= pp[ind])
            {
            label_50:
                if (x < pp[ind + 1])
                {
                    ftapp = (max - min) * ((Convert.ToDouble(ind - 1) + 0.5) / 500.0) + min;
                    return ftapp;
                }
                else
                {
                    ind++;
                    goto label_50;
                }
            }
            else
            {
                ind--;
            label_60:
                if (x >= pp[ind])
                {
                    ftapp = (max - min) * ((Convert.ToDouble(ind - 1) + 0.5) / 500.0) + min;
                    return ftapp;
                }
                else
                {
                    ind--;
                    goto label_60;
                }
            }
        }

        // Generate a gaussian distribution
        // 'ave' is the mean value of the resulting distribution
        // 'sd' is the standard deviation (sigma)
        private static void gauss(double ave, double sd)
        {
            double u1, u2, x;
            x = ran3(1);
            u1 = Math.Sqrt(-2.0 * Math.Log(x));
            x = ran3(1);
            u2 = x * 6.2831853072;
            x1 = u1 * Math.Cos(u2);
            x2 = u1 * Math.Sin(u2);
            x1 = ave + sd * x1;
            x2 = ave + sd * x2;
            //return x1;
        }

        // Calculate the number of individually treated fragments that will be
        // processed by either getvel or getrot.  The flage tflag' tells if the
        // number for getvel or getrot is to be determined.
        private static int ggntot()
        {
            int gntot;
            if (velflag + sizflag + rotflag == 0)
            {
                gntot = 0;
            }
            else
            {
                if (onlyin > 0)
                {
                    if (Convert.ToInt32(fnum[onlyin]) > totnum)
                    {
                        gntot = totnum;
                    }
                    else
                    {
                        gntot = Convert.ToInt32(fnum[onlyin]);
                    }
                }
                else
                {
                    gntot = 0;
                    for (int i = imax; i < imin + 1; i++)
                    {
                        if (Convert.ToInt32(fnum[i] - rvtot2[i]) > maxnum)
                        {
                            goto label_20;
                        }
                        if (gntot + Convert.ToInt32(fnum[i]) > totnum)
                        {
                            gntot = totnum;
                            goto label_20;
                        }
                        if (bsize[i] < minsize)
                        {
                            if (i != 1)
                            {
                                if (bsize[i - 1] > minsize)
                                {
                                    gntot += Convert.ToInt32(fnum[i]);
                                }
                                else
                                {
                                    gntot = 0;
                                }
                            }
                            else
                            {
                                gntot += Convert.ToInt32(fnum[i]);
                            }
                            goto label_20;
                        }
                        gntot += Convert.ToInt32(fnum[i]);
                    label_20:;
                    }
                }
            }
            if (gntot != 2 * Convert.ToInt32(Convert.ToDouble(gntot / 2.0))
                && nflag != 1)
            {
                gntot--;
            }
            return gntot;
        }

        private static void zero()
        {
            for (int i = 0; i < 101; i++)
            {
                fnum[i] = 0;
                masses[i] = 0.0;
                size[i] = 0;
                bsize[i] = 0.0;
                t[i] = 0.0;
                bsize[i] = 0.0;
                vpeak[i] = 0.0;
                vchar[i] = 0.0;
                vsp[1, i] = 0.0;
                dfnum[1, i] = 0.0;
                p[i] = 0.0;
                rvtot2[i] = 0;
                rvmat[1, i] = 0;
                rvnum[i] = 0;

                if (i < 56)
                {
                    mma[i] = 0;
                }

                if (i < 11)
                {
                    fmat[1, i] = 0.0;
                    tfmat[1, i] = 0.0;
                    rho[1, i] = 0.0;
                    trho[1, i] = 0.0;
                    length[i] = 0.0;
                    width[i] = 0.0;
                    th[i] = 0.0;
                    indind[i] = 0;
                    nrho[i] = 0;
                    tapp[i] = 0;
                }

                if (i < 21)
                {
                    npart[1, i] = 0;
                }
            }

            for (int i = 0; i < 501; i++)
            {
                pp[i] = 0.0;
            }
            //        ma[1] = 0.0;
            nmat[1] = 0;
            //        eninex = 0.0;
            mdum = 0.0;
            md = 0.0;
            tnmat[1] = 0;
            vrel = 0.0;
            en = 0.0;
            ntot = 0;
            totnum = 0;
            n = 0;
            alpha = 0.0;
            beta = 0.0;
            mean = 0.0;
            inext = 0;
            inextp = 0;
            onlyin = 0;
            e1 = 0.0;
            e2 = 0.0;
            com = 0.0;
            tappdiv = 0.0;
            t1 = 0.0;
            x1 = 0.0;
            x2 = 0.0;
            iff = 0;

        }

        // LFrag/HAZX parameters
        // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        // Constant definitions
        private const int MAX_FRAG = 100;                // Maximum number of fragments groups (bins)
        private const int MAX_RECORDS = 9975;            // Maximum possible number of records
        private const int MAX_TOTAL = 1000000;           // Maximum total number of fragments
        private const double KG_to_LB = 2.204622622;       // Conversion constants
        private const double M_to_FT = 3.280839895;
        private const double FT_to_M = 1.0 / M_to_FT;
        private const double M2_to_FT2 = 10.76391042;
        private const double NperSQM_to_LBperSQFT = 0.020885472;
        private const double N_to_LB = 0.224808943;
        private const double J_to_FTLB = 0.7375621;
        private const double GRAMS_TNT_to_J = 4184.0;
        private const double LB_to_GRAMS = 453.592;
        private const double LB_TNT_to_J = LB_to_GRAMS * GRAMS_TNT_to_J;
        private const double J_to_LB_TNT = 1.0 / LB_TNT_to_J;
        private const double M_PI = Math.PI;
        private const double M_PI_2 = Math.PI / 2.0;
        private const double M_PI_4 = Math.PI / 4.0;
        private const double Maxwell_Mean_to_3sigma = 2.372111;
        private const double GRAVITY_ft_per_s2 = 32.17404;

        public class FRAGMENT                           // Description of each mass bin
        {
            public double Vel3Sig_fps, LDRatio;         // 3-sigma induced velocity and L/D ratio
            public double Vel97, Vel_fps;               // 97th percentile and mean indiv spread vel
            public double Area_ft2, AreaMax;            // Reference area and maximum projected area
            public double Mass_kg, Weight_lb;           // Mass and weight
            public double Beta_psf, CD;                 // Ballistic coefficient and drag coef
            public double Length, Width, Thickness;     // Fragment overall dimensions, ft
            public double Distance;                     // Distance of fragment from center of explosion, ft
            public Shape shape;                         // Shape of fragment
            public string Description = "";             // Description of fragment (optional)
            public bool bIsIntact = false;              // Flag indicating intact fragment
            public bool bIncludeVrel = false;           // Flag indicating that relative velocity components are in intact components file
            public bool bIncludeBeta = false;           // Flag indicating that beta value (N/m2) is included in intact components file
            public double Number;                       // Number of fragments
            public int RecordNumber;                    // Record number to which the fragment belongs
        }

        // HAZX related flags
        private bool bValidOpen = true;
        private bool bFileSaved = false;
        private bool bHistogramGenerated = false;

        // Fragment generation

        // Parameters used for HAZX and FRG2 files
        private double dFragment_Minimum_Distance_ft = 0.0;
        private double dFragment_Maximum_Distance_ft = 0.0;
        private bool bCompute_HAZX_Fragments = false;
        private double dSpread_Velocity_Efficiency_pct = 0.0;
        public string s_Fragment_Bin_Counts_TP14 = "";
        public double dMax_Energy_to_Mass_Ratio_J_per_kg = 500.0E+06;
        public bool bEnergy_mass_ratio_exceeded = false;

        // Parameters used for FRG2 files only
        private bool bInclude_Unreacted_Liquids = true;
        private bool bMostly_Subsonic_Fallback = false;
        private bool bInclude_Intact_Components = false;
        private string sIntact_Components_File = null;
        private bool bIntact_Components_Only = false;
        private bool bSort_Intact_Components_Into_Groups = false;
        private double dMinimum_Beta_psf = 0.0;
        private double dSpread_Velocity_Energy_J = 0.0;

        public enum Shape { box, plate, thin_plate };   // Fragment shape
        private const double DVefficiency = 0.01;
        private string[] values;                         // Individual paramters in line of data
        private ArrayList FragmentList;                         // List of fragment objects
        private double SumAR3;
        private double dEnergy_J;
        private double dTotalFragmentMass_lbm = 0;                          // Running total fragment mass, lbm
        private double dTotalIntactMass_lbm = 0.0;
        private int iTotalNumberOfFragments = 0;
        private double[] MSD_Mass = { 0.002, 0.005, 0.01, 0.02,  0.05,   0.1,   0.2,    0.5,    1.0,    2.0,
                                       5.0,  10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0 };
        private double[] dRmax = new double[20];
        private double[] dVel = new double[20];
        private int[] iCount = new int[20];
        private double[] V0 = { 0.0, 200.0,  400.0,  600.0,  800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0,
                                   3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0, 5600.0, 5800.0, 6000.0 };
        private double[,] Rmax = {
            { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0 },
            { 421.8, 477.7, 537.3, 620.5, 685.6, 750.4, 834.1, 893.9, 950.1,1017.0,1061.5, 1100.4, 1143.5, 1170.4, 1193.4, 1217.7, 1232.4, 1244.3, 1256.7 },
            { 644.1, 755.0, 881.9,1076.6,1245.3,1431.2,1703.2,1926.9,2165.3,2491.8,2745.4, 2996.4, 3318.8, 3548.1, 3760.4, 4010.4, 4177.1, 4322.4, 4482.8 },
            { 768.3, 911.7,1078.9,1344.2,1580.6,1851.2,2264.2,2621.5,3021.6,3606.7,4094.8, 4613.6, 5335.4, 5897.3, 6459.5, 7191.1, 7716.6, 8207.4, 8793.4 },
            { 851.0,1015.8,1210.6,1522.8,1806.3,2135.4,2648.6,3103.0,3624.2,4411.1,5091.5, 5840.9, 6934.5, 7829.4, 8773.9,10075.1,11069.2,12061.8,13325.0 },
            { 909.9,1090.3,1304.6,1651.2,1968.0,2339.8,2925.8,3452.0,4063.5,5003.5,5832.3, 6763.1, 8162.9, 9340.7,10630.2,12465.2,13944.4,15470.5,17511.5 },
            { 953.1,1144.6,1373.6,1744.8,2086.9,2489.6,3129.3,3708.8,4386.1,5442.5,6381.4, 7456.8, 9089.4,10502.3,12065.4,14365.5,16263.6,18284.2,21073.4 },
            { 985.9,1186.4,1425.9,1816.6,2176.7,2604.2,3284.3,3905.1,4634.8,5778.6,6807.5, 7988.2, 9814.5,11403.1,13198.3,15871.6,18134.4,20588.0,24094.7 },
            { 1011.9,1219.2,1467.3,1873.2,2248.5,2694.0,3408.0,4059.4,4831.1,6044.5,7143.5, 8414.6,10388.4,12129.8,14100.6,17096.6,19651.2,22477.8,26588.5 },
            { 1033.5,1246.4,1501.8,1919.9,2307.7,2768.4,3509.6,4188.3,4992.5,6266.2,7419.8, 8765.9,10864.4,12729.6,14860.8,18114.4,20934.1,24064.6,28717.7 },
            { 1053.0,1270.7,1532.6,1961.8,2361.4,2836.3,3601.4,4305.0,5138.5,6465.9,7673.5, 9082.4,11298.8,13270.5,15546.9,19040.8,22096.3,25528.8,30660.6 },
            { 1069.2,1290.9,1558.5,1996.9,2406.7,2899.5,3684.6,4412.7,5274.7,6649.0,7906.7, 9375.6,11697.8,13777.1,16177.3,19903.2,23168.5,26879.8,32473.5 },
            { 1086.2,1311.6,1585.0,2033.9,2451.6,2957.3,3763.6,4511.5,5401.1,6819.1,8123.0, 9650.8,12067.7,14247.0,16769.1,20703.6,24181.7,28135.8,34178.6 },
            { 1101.8,1331.5,1610.0,2067.7,2494.5,3011.8,3838.5,4605.0,5519.4,6981.5,8325.1, 9908.0,12416.3,14686.2,17326.9,21451.9,25129.8,29330.1,34178.6 },
            { 1116.8,1350.6,1633.8,2100.3,2535.5,3062.9,3908.7,4693.1,5631.0,7134.5,8516.4,10150.2,12748.7,15099.6,17832.5,21451.9,26022.8,29330.1,34178.6 },
            { 1130.7,1368.3,1655.8,2130.6,2573.5,3111.1,3974.2,4777.3,5735.6,7277.7,8698.6,10377.0,13059.9,15492.7,18343.4,22840.1,26022.8,29330.1,34178.6 },
            { 1143.6,1384.6,1676.5,2158.7,2609.8,3156.5,4035.9,4855.3,5833.6,7411.3,8868.3,10589.2,13349.6,15860.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1155.8,1399.9,1696.1,2185.0,2643.3,3199.3,4093.6,4928.6,5926.4,7536.2,9027.0,10790.9,13620.5,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1167.3,1414.8,1714.6,2210.1,2675.4,3240.0,4147.8,4997.4,6014.1,7635.8,9176.2,10980.3,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1178.1,1428.7,1732.2,2233.7,2704.9,3277.4,4199.1,5062.4,6069.2,7764.5,9316.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1188.4,1441.2,1748.8,2255.7,2733.2,3313.6,4247.4,5123.9,6174.4,7871.2,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1198.4,1453.8,1764.3,2277.2,2760.1,3347.8,4294.9,5182.1,6248.8,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1207.7,1465.6,1779.3,2297.9,2785.8,3380.2,4339.4,5237.0,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1216.9,1476.9,1794.0,2317.1,2809.9,3411.5,4381.6,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1225.4,1487.6,1806.8,2335.2,2833.4,3441.2,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1233.7,1498.1,1820.2,2353.9,2856.0,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1241.1,1508.0,1833.2,2371.0,2877.7,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1249.2,1517.8,1844.8,2387.7,2877.7,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1256.7,1526.9,1856.6,2387.7,2877.7,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1264.0,1536.0,1856.6,2387.7,2877.7,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 },
            { 1270.5,1536.0,1856.6,2387.7,2877.7,3469.7,4422.0,5290.3,6319.0,7972.8,9449.7,11159.0,13874.9,16204.5,18801.2,22840.1,26022.8,29330.1,34178.6 } };

        private double[,] Rmax_TP14 = {
            { 0,0,0,0,0,0,0,0,0,0 },
            { 1089.1,1037.1,977.7,909.5,836.5,756.7,675.7,597.3,524,441.3 },
            { 2920.4,2602.4,2292.7,1990.6,1711.4,1450.5,1218.5,1020.1,852.8,682.1 },
            { 4452.6,3814.8,3245.5,2725.7,2277.1,1879.8,1542.2,1266.2,1040.3,817 },
            { 5605.2,4698.5,3921.7,3237.5,2664.9,2170.5,1760.2,1430.5,1165.3,906.8 },
            { 6467.7,5352.5,4415.4,3609.9,2944.5,2379.8,1916.3,1548.1,1254.8,970.6 },
            { 7114.4,5836.3,4676.7,3882.9,3149.9,2532.9,2030.9,1634.3,1320.2,1017.7 },
            { 7610.7,6209.4,5059.7,4092.6,3306.1,2650,2117.5,1700.3,1369.9,1053.5 },
            { 8007.7,6502.9,5282.3,4257.2,3431,2742.3,2187,1752.1,1409.4,1081.8 },
            { 8334.3,6744.1,5465.2,4257.2,3533.6,2818.2,2244.5,1794.8,1442.2,1105.1 },
            { 8628.3,6969.5,5629.7,4393.7,3626.2,2887.7,2295.8,1833.4,1471.5,1126.1 },
            { 8903.5,7173,5782.6,4518.2,3711,2952,2343.3,1869,1498.5,1145.6 },
            { 9159.1,7362,5925.9,4738.4,3790.3,3011.1,2387.3,1902.8,1523.7,1163.7 },
            { 9398.1,7539.3,6060.1,4837.7,3866.3,3066.9,2428.4,1934,1547.4,1180.9 },
            { 9623.2,7709.2,6186.3,4931.4,3937.1,3119.5,2468.1,1964,1570.1,1197.4 },
            { 9833.9,7868.3,6304.7,5020.6,4003.3,3168.3,2505.2,1991.6,1591.2,1212.5 },
            { 10033.7,8016.3,6415.2,5104.2,4065.5,3215,2539.7,2017.7,1610.9,1226.6 },
            { 10179.6,8155.1,6518.5,5182.2,4123.8,3258.3,2572.2,2041.9,1629.6,1226.6 },
            { 10397.4,8285.6,6617.3,5255.4,4178.7,3300.3,2602.8,2064.8,1647.3,1252.4 },
            { 10563.3,8408.4,6711.1,5324.9,4230.4,3339.1,2631.9,2086.7,1664.2,1264.4 },
            { 10720.4,8524.5,6799.1,5390.6,4278.7,3375.7,2659.1,2106.9,1679.6,1275.6 },
            { 10869.2,8636.8,6883.1,5452.8,4326.7,3410.8,2685,2127.1,1694.5,1286.3 },
            { 11010.8,8744,6962.8,5512.2,4371.7,3444.4,2710.5,2145.8,1708.8,1296.3 },
            { 11145.6,8845.9,7038.9,5568.2,4414.4,3475.4,2733.2,2164,1722.3,1306.2 },
            { 11275.7,8934.6,7111.6,5622.3,4455.2,3506.2,2756.2,2180.4,1735.2,1315.5 },
            { 11401.5,9036.5,7181,5674.7,4494.3,3535.1,2777.6,2197.2,1748,1324.7 },
            { 11522.1,9125.7,7247.6,5725.2,4532.2,3562.8,2798.5,2213.3,1759.9,1332.7 },
            { 11637.9,9211.6,7311.6,5773.5,4568.4,3590.1,2819.7,2228.5,1771.1,1341.4 },
            { 11749.2,9294.3,7373.1,5820.6,4603.1,3616,2839.2,2243.2,1782.4,1349.9 },
            { 11856.4,9373.6,7432.8,5865.2,4637.1,3640.4,2857.3,2257.1,1792.8,1357.4 },
            { 11959.8,9450.5,7491.2,5898.5,4669,3664.2,2875.7,2270.8,1803.2,1364.6 } };

        // TP-14 Bin data
        private double[] dCount_TP14 = new double[10];
        private double[] dWeight_TP14 = { 35.7, 14.9, 6.34, 2.66, 1.13, 0.473, 0.199, 0.0852, 0.0379, 0.0142 };
        private double[] dVel_TP14 = new double[10];
        private double[] dRmax_TP14 = new double[10];
        private bool bComputeTP14bins = true;

        // Static yield parameters
        private double dWeight = 0.0;
        private string sStaticFuelLoc = "";
        private double dStaticFuelDia = 0.0;
        private double dStaticFuelTankLength = 0.0;
        private double dStaticFuelFillHeight = 0.0;
        private double dStaticFuelPressure = 0.0;
        private double dStaticFuelWt = 0.0;
        private double dStaticFuelDensity = 0.0;
        private string sStaticOxidizerLoc = "";
        private double dStaticOxidizerDia = 0.0;
        private double dStaticOxidizerTankLength = 0.0;
        private double dStaticOxidizerFillHeight = 0.0;
        private double dStaticOxidizerPressure = 0.0;
        private double dStaticOxidizerWt = 0.0;
        private double dStaticOxidizerDensity = 0.0;
        private double dStaticMixedHt = 0.0;
        private double dTotalWeightLeaked = 0.0;
        private double dStaticTankHeight = 0.0;
        private double dStaticTotalProp = 0.0;
        private double dStaticHoleDiameter = 0.0;
        private double dStaticIgnitionTime = 0.0;
        private double dStaticMixedDen = 0.0;
        private double Cd = 0.95;
        private double dLowerUllageHeight = 0.0;
        private double dLowerTankFillHeight = 0.0;
        private double dUpperTankRadius = 0.0;
        private double dLowerTankRadius = 0.0;
        private double dPressureUpperTankTop = 0.0;
        private double dPressureLowerTankTop = 0.0;
        private double dUpperTankFillHeight = 0.0;
        private double dLowerTankLength = 0.0;
        private double dLowerTankWeight = 0.0;
        private double dMaximumPropellantLeaked = 9.9E15;
        private double dStaticUpperDensity = 0.0;
        private double dStaticLowerDensity = 0.0;
        private double dLCH4_OFratio = 4.00;
        private double dRP1_OFratio = 3.40;  // See ref. [1]
        private double dLH2_OFratio = 7.94;
        private double dMaxOFratio = 0.0;
        private double dHypergolStaticYieldFactor = 0.05;  //  5% based on DESR
        private double dHypergolDynamicYieldFactor = 0.10; // 10% based on DESR
        // Ref: DESR 6055.09_AFMAN91-201, 28 May 2020, Table V5.E4.T5 "Energetic Liquid Equivalent Weights"
        // Defense Explosives Safety Regulation
        // Subject: "Department of the Air Force Guidance Memorandum to Defense Explosives Safety
        //           Regulation 6055.09 Air Force Manual 91-201, Explosives Safety Standards"
        private bool bValidStaticRun = false;
        private bool bMessageGivenOnDefaultCd = false;
        private bool bMessageGivenOnAreaRatio = false;
        private bool bLeakAtDownComer = false;
        private double[,] dYieldvsOF_LCH4 = { { 0.0, 0.0},       // Yield as a function of LOX/LCH4 O/F ratio
                                            { 1.0, 0.7 },
                                            { 4.0, 1.6 },
                                            {14.0, 0.7 },
                                            {22.705, 0.01 },
                                            {100.0, 0.01 } };
        private double[,] dYieldvsOF_RP1 =  { {0.0, 0.0000},       // Yield as a function of LOX/RP-1 O/F ratio
                                            {0.5, 1.0577},
                                            {1.0, 1.5158},
                                            {1.5, 1.7796},
                                            {2.0, 1.9491},
                                            {2.5, 2.0687},
                                            {3.0, 2.1572},
                                            {3.5, 2.1777},
                                            {3.6, 2.1293},
                                            {3.8, 2.0603},
                                            {4.0, 1.9540},
                                            {4.5, 1.7696},
                                            {5.0, 1.6154},
                                            {5.5, 1.4847},
                                            {6.0, 1.3725},
                                            {6.5, 1.2751},
                                            {7.0, 1.1897},
                                            {7.5, 1.1144},
                                            {8.0, 1.0474},
                                            {8.5, 0.9873},
                                            {9.0, 0.9333},
                                            {9.5, 0.8843},
                                            {10.0,0.8398},
                                            {10.5,0.7991},
                                            {11.0,0.7618},
                                            {11.5,0.7275},
                                            {12.0,0.6957},
                                            {13.0,0.6391},
                                            {14.0,0.5899},
                                            {15.0,0.5469},
                                            {20.0,0.3928},
                                            {30.0,0.2529},
                                            {30.01,0.01 } };
        private double[,] dYieldvsOF_LH2 = { { 0.00, 0.0000},        // Yield as a function of LOX/LH2 O/F ratio
                                            {0.5, 0.6468},
                                            {1.0, 1.3311},
                                            {1.5, 1.7417},
                                            {2.0, 2.0155},
                                            {2.5, 2.2110},
                                            {3.0, 2.3403},
                                            {3.5, 2.4169},
                                            {4.0, 2.4731},
                                            {4.5, 2.5138},
                                            {5.0, 2.5421},
                                            {5.5, 2.5610},
                                            {6.0, 2.5716},
                                            {6.5, 2.5752},
                                            {7.0, 2.5734},
                                            {7.5, 2.5658},
                                            {8.0, 2.5345},
                                            {8.5, 2.4151},
                                            {9.0, 2.3049},
                                            {9.5, 2.2041},
                                            {10.0,2.1117},
                                            {10.5,2.0266},
                                            {11.0,1.9480},
                                            {11.5,1.8753},
                                            {12.0,1.8077},
                                            {13.0,1.6861},
                                            {14.0,1.5796},
                                            {15.0,1.4855},
                                            {17.0,1.3268},
                                            {20.0,1.1426},
                                            {25.0,0.9259},
                                            {30.0,0.7767},
                                            {35.0,0.6676},
                                            {40.0,0.5783},
                                            {50.0,0.4501},
                                            {60.0,0.3639},
                                            {60.01,0.01 } };


        // Generate Fragments by running MSD, then either creating fragments for HAZX or LFrag
        public string Generate_Fragments(string sInput_File_Name, string sMSD_Output_Specified,
                                         string sOutput_File_Name, bool bSingle_Record,
                                         double dUser_specified_max_energy_to_mass_ratio_MJ_per_kg)
        {
            dMax_Energy_to_Mass_Ratio_J_per_kg =
                dUser_specified_max_energy_to_mass_ratio_MJ_per_kg * 1.0E+06;
            return Generate_Fragments(sInput_File_Name, sMSD_Output_Specified,
                                      sOutput_File_Name, bSingle_Record);
        }
        // Generate Fragments by running MSD, then either creating fragments for HAZX or LFrag
        // Default max energy to mass ratio is 5 MJ/kg unless called by overloaded method
        public string Generate_Fragments(string sInput_File_Name, string sMSD_Output_Specified,
                                         string sOutput_File_Name, bool bSingle_Record)
        {
            System.IO.StreamReader inputReader;
            try
            {
                inputReader = System.IO.File.OpenText(sInput_File_Name);
            }
            catch
            {
                return " Error: unable to opent input file.";
            }

            string sLine = inputReader.ReadLine();              // Read inert mass to be broken up, kg
            string[] values = sLine.Split(',');
            try
            {
                ma[1] = Convert.ToDouble(values[0].Trim());
                sLine = inputReader.ReadLine();                     // Read number of materials
                values = sLine.Split(',');
                nmat[1] = Convert.ToInt32(values[0].Trim());
                for (int i = 0; i < nmat[1]; i++)
                {
                    sLine = inputReader.ReadLine();                 // Read material fractions and densities (kg/m3)
                    values = sLine.Split(',');
                    fmat[1, i + 1] = Convert.ToDouble(values[0].Trim());
                    rho[1, i + 1] = Convert.ToDouble(values[1].Trim());
                }
                sLine = inputReader.ReadLine();                     // Read energy of breakup, J
                values = sLine.Split(',');
                eninex = Convert.ToDouble(values[0].Trim());

                // Limit energy-to-mass ratio to no more than 5E-6 J/kg
                double dEnergy_mass_ratio = eninex / ma[1];
                if (dEnergy_mass_ratio > dMax_Energy_to_Mass_Ratio_J_per_kg)
                {
                    eninex = ma[1] * dMax_Energy_to_Mass_Ratio_J_per_kg;
                    bEnergy_mass_ratio_exceeded = true;
                }

                sLine = inputReader.ReadLine();                     // Read fragment minimum distance, m
                values = sLine.Split(',');
                dFragment_Minimum_Distance_ft = Convert.ToDouble(values[0].Trim());
                sLine = inputReader.ReadLine();                     // Read fragment maximum distance, m
                values = sLine.Split(',');
                dFragment_Maximum_Distance_ft = Convert.ToDouble(values[0].Trim());

                sLine = inputReader.ReadLine();                     // Read flag to compute TP14 fragments (not FRG2)
                values = sLine.Split(',');
                bCompute_HAZX_Fragments = Convert.ToBoolean(values[0].Trim().ToLower());
            }
            catch
            {
                inputReader.Close();
                string sError = Write_Error_Message();
                try
                {
                    System.IO.File.WriteAllText(sOutput_File_Name, sError);
                }
                catch
                {
                    System.IO.File.WriteAllText(sOutput_File_Name, "Error: MSD parameters not read correctly");
                }
                return sError;
            }

            if (!bCompute_HAZX_Fragments)                      // Skip remaining inputs if fragments are for HAZX (not FRG2)
            {
                try
                {
                    sLine = inputReader.ReadLine();                 // Read spread velocity energy, J
                    values = sLine.Split(',');
                    dSpread_Velocity_Energy_J = Convert.ToDouble(values[0].Trim());

                    sLine = inputReader.ReadLine();                 // Read spread efficiency, %
                    values = sLine.Split(',');
                    dSpread_Velocity_Efficiency_pct = Convert.ToDouble(values[0].Trim());

                    sLine = inputReader.ReadLine();                 // Read flag to include effects of unreacted liquids
                    values = sLine.Split(',');
                    bInclude_Unreacted_Liquids = Convert.ToBoolean(values[0].Trim().ToLower());

                    sLine = inputReader.ReadLine();                 // Read flag to treat fallback as mostly subsonic
                    values = sLine.Split(',');
                    bMostly_Subsonic_Fallback = Convert.ToBoolean(values[0].Trim().ToLower());

                    sLine = inputReader.ReadLine();                 // Read flag to include intact components
                    values = sLine.Split(',');
                    bInclude_Intact_Components = Convert.ToBoolean(values[0].Trim().ToLower());

                    sLine = inputReader.ReadLine();                 // Read intact components input file
                    values = sLine.Split(',');
                    sIntact_Components_File = values[0].Trim();

                    sLine = inputReader.ReadLine();                 // Read flag to use intact components only
                    values = sLine.Split(',');
                    bIntact_Components_Only = Convert.ToBoolean(values[0].Trim().ToLower());

                    sLine = inputReader.ReadLine();                 // Read flag to sort intact components into groups
                    values = sLine.Split(',');
                    bSort_Intact_Components_Into_Groups = Convert.ToBoolean(values[0].Trim().ToLower());

                    sLine = inputReader.ReadLine();                 // Read minimum beta to output, lb/ft2
                    values = sLine.Split(',');
                    dMinimum_Beta_psf = Convert.ToDouble(values[0].Trim());
                    inputReader.Close();
                }
                catch
                {
                    inputReader.Close();
                    string sError = Write_Error_Message();
                    try
                    {
                        System.IO.File.WriteAllText(sOutput_File_Name, sError);
                    }
                    catch
                    {
                        return sError;
                    }
                    return sError;
                }
            }
            inputReader.Close();

            // Run MSD
            sMSD_output = sMSD_Output_Specified;
            bool bSuccessful_MSD_Run = MSD();
            if (!bSuccessful_MSD_Run)
            {
                return " Unsucssessful MSD run.\n";
            }

            // Compute HAZX or FRG2 fragments
            if (bCompute_HAZX_Fragments)
            {
                return Compute_HAZX_Fragments(sOutput_File_Name, bSingle_Record);
            }
            else
            {
                return Compute_FRG2_Fragments("Fragments.FRG2.csv");
            }
        }

        private string Compute_FRG2_Fragments(string sFRG2_file)
        {





            return "";
        }



        private string Compute_HAZX_Fragments(string sFragmentBinFile, bool bSingle_Record)
        {
            bool bAdjustInertWeightsBasedOnFillHeight = false;
            double dFillRatio = 0.5;
            double dEnergy_J = eninex;
            double dInertMass = ma[1];
            string sFragment_bins = "";
            double dBinYield = eninex * J_to_LB_TNT;
            double dYield_lb_TNT = dBinYield;

            // Initialized arrays
            for (int i = 0; i < 20; i++)
            {
                dVel[i] = 0;
                dRmax[i] = 0;
            }
            for (int i = 0; i < 10; i++)
            {
                dCount_TP14[i] = 0;
                dVel_TP14[i] = 0;
                dRmax_TP14[i] = 0;
            }

            bool bInvalidMSD = false;
            double dWeightMax = 0.0;
            double dWeightMin = 1.0e9;
            FragmentList = new ArrayList();
            System.IO.StreamReader inputReader;
            inputReader = System.IO.File.OpenText(sMSD_output);
            string sLine = inputReader.ReadLine();                     // Skip header
            while ((sLine = inputReader.ReadLine()) != null)
            {
                if (bInvalidMSD) break;
                FRAGMENT fragment = new FRAGMENT();
                values = sLine.Split(',');                      // Parse line of data
                try
                {   // Add fragment inputs to list
                    fragment.Mass_kg = Double.Parse(values[0]);
                    fragment.Weight_lb = fragment.Mass_kg * KG_to_LB;
                    fragment.Length = Double.Parse(values[1]) * M_to_FT;
                    fragment.Width = Double.Parse(values[2]) * M_to_FT;
                    fragment.Thickness = Double.Parse(values[3]) * M_to_FT;
                    fragment.Number = 1.0;
                    FragmentList.Add(fragment);

                    // Record running max and min fragment weight
                    dWeightMax = Math.Max(dWeightMax, fragment.Weight_lb);
                    dWeightMin = Math.Min(dWeightMin, fragment.Weight_lb);
                }
                catch
                {
                    bInvalidMSD = true;
                }
            }
            inputReader.Close();
            if (bInvalidMSD || FragmentList.Count == 0)          // Invalid MSD run
            {
                string sErrorMessage = "\n MSD inputs do not seem to be valid.  Please review them.\n";
                if (dEnergy_J / (dInertMass / KG_to_LB) <= 1e4)
                {
                    sErrorMessage += " Often the problem is a very low energy-to-mass ratio.\n";
                }
                else
                {
                    sErrorMessage += " Often the problem is a very low or high energy-to-mass ratio.\n";
                }
                return sErrorMessage;
            }

            // Initialized arrays
            for (int i = 0; i < 20; i++)
            {
                iCount[i] = 0;
            }

            // Sort fragments by MSD mass and find max V0
            foreach (object o in FragmentList)
            {
                FRAGMENT fragment = (FRAGMENT)o;
                if (fragment.Mass_kg < MSD_Mass[0])
                {
                    fragment.Mass_kg = MSD_Mass[0];   // Round up unlikely very small fragments
                }
                if (fragment.Mass_kg > MSD_Mass[19])
                {
                    fragment.Mass_kg = MSD_Mass[19];  // Round down unlikely very large fragments
                }
                for (int i = 0; i < 20; i++)
                {
                    if (MSD_Mass[i] == fragment.Mass_kg)
                    {
                        iCount[i] += 1;
                        break;
                    }
                }
            }

            // Scale up fragment counts based on dFillRatio
            if (bAdjustInertWeightsBasedOnFillHeight)
            {
                for (int i = 0; i < 20; i++)
                {
                    iCount[i] *= Convert.ToInt32(iCount[i] / dFillRatio);
                }
            }

            // Compute bin counts
            // Per J. Chrostowski email dated 10/20/2024:
            // TP Bin 1 = LFRAG bins 13 thru 20
            // TP Bin 2 = LFRAG bins 11 + 12
            // TP Bin 3 = LFRAG bins 10
            // TP Bin 4 = LFRAG bins 9
            // TP Bin 5 = LFRAG bins 8
            // TP Bin 6 = LFRAG bins 7
            // TP Bin 7 = LFRAG bins 6
            // TP Bin 8 = LFRAG bins 5
            // TP Bin 9 = LFRAG bins 4
            // TP Bin 10 = LFRAG bins 1 - 3
            if (bComputeTP14bins)
            {
                // Compute TP14 bin counts using bins 13 thru 20
                double dBin_1_wt = 0.0;
                for (int i = 12; i < 20; i++)
                {
                    dBin_1_wt += MSD_Mass[i] * KG_to_LB * (double)iCount[i];
                }
                dCount_TP14[0] = dBin_1_wt / dWeight_TP14[0];

                // Compute bin counts for TP14 bin 2 using bins 11 and 12
                double dBin_2_wt = (MSD_Mass[10] * (double)iCount[10] + MSD_Mass[11] * (double)iCount[11]) * KG_to_LB;
                dCount_TP14[1] = dBin_2_wt / dWeight_TP14[1];

                // Compute bin counts for TP14 bins 3 through 9 using bins 10 through 4
                dCount_TP14[2] = MSD_Mass[9] * KG_to_LB * (double)iCount[9] / dWeight_TP14[2];
                dCount_TP14[3] = MSD_Mass[8] * KG_to_LB * (double)iCount[8] / dWeight_TP14[3];
                dCount_TP14[4] = MSD_Mass[7] * KG_to_LB * (double)iCount[7] / dWeight_TP14[4];
                dCount_TP14[5] = MSD_Mass[6] * KG_to_LB * (double)iCount[6] / dWeight_TP14[5];
                dCount_TP14[6] = MSD_Mass[5] * KG_to_LB * (double)iCount[5] / dWeight_TP14[6];
                dCount_TP14[7] = MSD_Mass[4] * KG_to_LB * (double)iCount[4] / dWeight_TP14[7];
                dCount_TP14[8] = MSD_Mass[3] * KG_to_LB * (double)iCount[3] / dWeight_TP14[8];

                // Compute bin counts for TP14 bin 10
                double dBin_10_wt = MSD_Mass[0] * KG_to_LB * (double)iCount[0];
                dBin_10_wt += MSD_Mass[1] * KG_to_LB * (double)iCount[1];
                dBin_10_wt += MSD_Mass[2] * KG_to_LB * (double)iCount[2];
                dCount_TP14[9] = dBin_10_wt / dWeight_TP14[9];


                if (bSingle_Record)
                {
                    sFragment_bins = "Yield_lb_TNT,";
                    sFragment_bins += "Propbability,Bin 1,Bin 2,Bin 3,Bin 4,Bin 5,Bin 6,Bin 7,Bin 8,Bin 9,Bin 10\n";
                    sFragment_bins += Math.Round(dYield_lb_TNT, 3).ToString("0.000E+00");
                    sFragment_bins += ",1,";
                }
                else
                {
                    sFragment_bins = ",";
                }

                // Write line of TP14 fragment bin counts
                if (dBinYield > 0.0)
                {
                    sFragment_bins += Convert.ToInt32(dCount_TP14[0]).ToString();
                    for (int i = 1; i < 10; i++)
                    {
                        sFragment_bins += "," + Convert.ToInt32(dCount_TP14[i]).ToString();
                    }
                }
                else
                {
                    sFragment_bins += ",0,0,0,0,0,0,0,0,0,0";
                }
            }
            else
            {
                // Write line of MSD fragment bin counts
                sFragment_bins += iCount[0].ToString();
                for (int i = 1; i < 20; i++)
                {
                    sFragment_bins += "," + iCount[i].ToString();
                }
            }
            sFragment_bins += "\n";

            // Compute imparted velocities
            double dTankRadius = dFragment_Minimum_Distance_ft;                 // Radius of Tank, ft
            double dTankArea = 2.0 * M_PI * dTankRadius * dTankRadius;          // Area of Tank * 2, ft2
            double Dmin_ft = dTankRadius;
            double Dmax_ft = dFragment_Maximum_Distance_ft;
            double LogWtMax = Math.Log10(dWeightMax);
            double LogWtMin = Math.Log10(dWeightMin);
            double Slope = (Dmax_ft - Dmin_ft) / (LogWtMax - LogWtMin);
            double Offset = Dmax_ft - Slope * LogWtMax;

            // Compute distributed fragment characteristics
            foreach (object o in FragmentList)
            {
                FRAGMENT fragment = (FRAGMENT)o;
                fragment = ComputeFragmentCharacteristics(fragment);
                fragment.Distance = Offset + Slope * Math.Log10(fragment.Weight_lb);  // Distributed fragment distances from center of explosion
            }

            // Compute parameters for imparted velocities
            SumAR3 = 0.0;
            bool bIncludedUnexploded = true;
            foreach (object o in FragmentList)                          // Compute stats for distributed fragments
            {
                FRAGMENT fragment = (FRAGMENT)o;
                SumAR3 += fragment.Number * fragment.Area_ft2           // Compute running sum of nA/R3, ft-1
                    / fragment.Distance / fragment.Distance / fragment.Distance;
            }
            if (bIncludedUnexploded)                                    // Include effects of unexploded propellant
            {
                double PropDist = dTankRadius / 2.0;                       // Propellant modeled distance from center of explosion, ft (estimated as half of tank radius)
                SumAR3 += dTankArea / PropDist / PropDist / PropDist;      // Add running sum
            }
            s_Fragment_Bin_Counts_TP14 = sFragment_bins;
            sFragment_bins += ',';

            // Write bin data for this yield
            if (bSingle_Record && sFragmentBinFile != "")
            {
                try
                {
                    sFragment_bins += Compute_Fragment_Bins(dYield_lb_TNT);
                    System.IO.File.WriteAllText(sFragmentBinFile, sFragment_bins);
                }
                catch
                {
                    sFragment_bins = "Error: fragment bins not written.";
                }
            }
            return sFragment_bins;
        }


        public FRAGMENT ComputeFragmentCharacteristics(FRAGMENT fragment)
        {
            // Characterize fragment shape as a box, plate or thin plate
            double L1 = fragment.Length; // Abreviated dimension names
            double L2 = fragment.Width;
            double L3 = fragment.Thickness;
            if (L3 < 2.0 * 0.0254 && L3 / L1 < 0.20) fragment.shape = Shape.plate;
            else fragment.shape = Shape.box;
            if (fragment.shape == Shape.plate && L2 / L1 < 0.10) fragment.shape = Shape.thin_plate;

            // Compute drag coefficient (depends on the shape)
            if (fragment.shape == Shape.box) fragment.CD = 0.75;
            else if (fragment.shape == Shape.plate) fragment.CD = 0.75;
            else /* shape=thin plate */
            {
                if (L2 / L1 < 0.05) fragment.CD = (2.0 - 10.0 * L2 / L1) / 1.56;
                else fragment.CD = (1.7 - 4.0 * L2 / L1) / 1.56;
            }

            // Estimate lift/drag based on shape
            if (fragment.shape == Shape.box) fragment.LDRatio = 0.01;
            else if (fragment.shape == Shape.plate) fragment.LDRatio = 0.03;
            else /* shape=thin plate */ fragment.LDRatio = 0.05;

            // Compute maximum projected area, reference area and ballistic coefficient
            fragment.AreaMax = L1 * L2;
            fragment.Area_ft2 = (L1 * L2 + L1 * L3 + M_PI_2 * L2 * L3) / M_PI_4 / M_PI;
            fragment.Beta_psf = fragment.Weight_lb / (fragment.CD * fragment.Area_ft2);
            return fragment;
        }


        public string Compute_Fragment_Bins(double dYield_lb_TNT)
        {
            //   string sFragmentBinFile = textBoxSingleFragListName.Text;
            double DVenergy_FTLB = dYield_lb_TNT * LB_TNT_to_J * J_to_FTLB;             // Convert TNT Yield to foot-lb
            foreach (object o in FragmentList)                                          // Calculate imparted velocities for distributed fragments
            {
                FRAGMENT fragment = (FRAGMENT)o;
                double AR3 = fragment.Area_ft2 / fragment.Distance / fragment.Distance / fragment.Distance;
                fragment.Vel3Sig_fps = Math.Sqrt(64.4 * DVenergy_FTLB * DVefficiency * AR3
                    / fragment.Weight_lb / SumAR3);
            }

            // Initialized arrays
            for (int i = 0; i < 20; i++)
            {
                dVel[i] = 0;
            }

            // Sort fragments by MSD mass and find max V0
            foreach (object o in FragmentList)
            {
                FRAGMENT fragment = (FRAGMENT)o;
                if (fragment.Mass_kg < MSD_Mass[0])
                {
                    fragment.Mass_kg = MSD_Mass[0];   // Round up unlikely very small fragments
                }
                if (fragment.Mass_kg > MSD_Mass[19])
                {
                    fragment.Mass_kg = MSD_Mass[19];  // Round down unlikely very large fragments
                }
                for (int i = 0; i < 20; i++)
                {
                    if (MSD_Mass[i] == fragment.Mass_kg)
                    {
                        dVel[i] = Math.Max(dVel[i], fragment.Vel3Sig_fps);
                        break;
                    }
                }
            }

            // Compute average TP14 bin velocities
            // See comments in method Compute_Fragment_Masses for explanation of bin indices
            if (bComputeTP14bins)
            {
                // Compute average velocities for TP14 bins 13 thru 20
                double dBin_1_vel = 0.0;
                double dCount_bins_13_thru_20 = 0.0;
                for (int i = 12; i < 20; i++)
                {
                    dBin_1_vel += dVel[i] * (double)iCount[i];
                    dCount_bins_13_thru_20 += (double)iCount[i];
                }
                if (dCount_bins_13_thru_20 > 0.0)
                {
                    dVel_TP14[0] = dBin_1_vel / dCount_bins_13_thru_20;
                }
                else
                {
                    dVel_TP14[0] = 0.0;
                }

                // Compute average velocity for TP14 bin 2
                double dBin_2_vel = dVel[10] * (double)iCount[10] + dVel[11] * (double)iCount[11];
                double dCount_bins_11_and_12 = (double)iCount[10] + (double)iCount[11];
                if (dCount_bins_11_and_12 > 0.0)
                {
                    dVel_TP14[1] = dBin_2_vel / dCount_bins_11_and_12;
                }
                else
                {
                    dVel_TP14[1] = 0.0;
                }

                // Compute velocities for TP14 bins 3 through 9
                for (int i = 2; i < 9; i++)
                {
                    dVel_TP14[i] = dVel[11 - i];
                }

                // Compute average velocities for TP14 bin 10
                double dBin_10_vel = dVel[0] * (double)iCount[0] + dVel[1] * (double)iCount[1] + dVel[2] * (double)iCount[2];
                double dCount_bins_1_thru_3 = (double)(iCount[0] + iCount[1] + iCount[2]);
                if (dCount_bins_1_thru_3 > 0.0)
                {
                    dVel_TP14[9] = dBin_10_vel / dCount_bins_1_thru_3;

                }
                else
                {
                    dVel_TP14[9] = 0.0;
                }
            }

            // Compute max throw distance for each mass bin
            if (bComputeTP14bins)
            {
                for (int i = 0; i < 10; i++)
                {
                    if (Convert.ToInt32(dCount_TP14[i]) > 0)
                    {
                        dRmax_TP14[i] = Interpolate_Rmax(i, dVel_TP14[i]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < 20; i++)
                {
                    if (iCount[i] > 0)
                    {
                        dRmax[i] = Interpolate_Rmax(i, dVel[i]);
                    }
                }
            }

            // Write bin data for this yield
            string sFragment_bins = "";
            if (bComputeTP14bins)
            {
                for (int i = 0; i < 10; i++)
                {
                    sFragment_bins += "," + Math.Round(dRmax_TP14[i], 0);
                }
                sFragment_bins += "\n,";
                for (int i = 0; i < 10; i++)
                {
                    sFragment_bins += "," + Math.Round(dVel_TP14[i], 0);
                }
            }
            else
            {
                for (int i = 0; i < 20; i++)
                {
                    sFragment_bins += "," + Math.Round(dRmax[i], 0);
                }
                sFragment_bins += "\n,";
                for (int i = 0; i < 20; i++)
                {
                    sFragment_bins += "," + Math.Round(dVel[i], 0);
                }
            }
            sFragment_bins += "\n";

            return sFragment_bins;
        }

        private double Interpolate_Rmax(int iMass, double dVelocity)
        {
            // Interpolates max horizontal throw distance (dRmax) data to the initial throw velocity for a given mass index
            int j;
            double dRmax = 0.0;
            if (dVelocity == 0.0)
            {
                return 0.0;
            }
            for (j = 0; j < 30; j++)
            {
                if (dVelocity <= V0[j])
                {
                    break;
                }
            }
            j--;
            if (bComputeTP14bins)
            {
                dRmax = Rmax_TP14[j, iMass] + (Rmax_TP14[j + 1, iMass] - Rmax_TP14[j, iMass]) * (dVelocity - V0[j]) / (V0[j + 1] - V0[j]);
            }
            else
            {
                dRmax = Rmax[j, iMass] + (Rmax[j + 1, iMass] - Rmax[j, iMass]) * (dVelocity - V0[j]) / (V0[j + 1] - V0[j]);
            }
            return dRmax;
        }

        private string Write_Error_Message()
        {
            string sOut = " Error reading input file.  Correct format:\n\n";
            sOut += "<value>,         Inert mass broken up (kg)\n";
            sOut += "<value>,         Number of materials\n";
            sOut += "<value>,<value>, Material fraction, Density (kg/m3)\n";
            sOut += "<value>,<value>, Material fraction, Density (kg/m3)\n";
            sOut += ".\n";
            sOut += ".\n";
            sOut += ".\n";
            sOut += "<value>,<value>, Material fraction, Density (kg/m3)\n";
            sOut += "<value>,         Energy of breakup (J)\n";
            sOut += "<value>,         Fragment Minimum distance from center (m)\n";
            sOut += "<value>,         Fragment Maximum distance from center (m)\n";
            sOut += "<value>,         Compute HAZX fragments instead of FRG2 (true/false)\n";
            sOut += "<value>,         Spread velocity energy (J)\n";
            sOut += "<value>,         Spread velocity efficiency (%)\n";
            sOut += "<value>,         Include effects of unreacted liquids (true/false)\n";
            sOut += "<value>,         Mostly subsonic fragment fallback (true/false)\n";
            sOut += "<value>,         Include intact components (true/false)\n";
            sOut += "<value>,         Intact components input file name\n";
            sOut += "<value>,         Intacts only without running MSD (true/false)\n";
            sOut += "<value>,         Sort intact components into groups (true/false)\n";
            sOut += "<value>,         Minimum Beta (lb/ft2)\n";
            return sOut;
        }



    }
}
