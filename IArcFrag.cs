using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Runtime.InteropServices;

namespace ArcFragTools
{
    [ComVisible(true)]
    public interface IArcFrag
    {
        [DispId(2001)]
        string Generate_Fragments(string sInput_File_Name, string sMSD_Output_Specified, string sOutput_File_Name, bool bSingle_Record, double dUser_specified_max_energy_to_mass_ratio_MJ_per_kg);
    }
}
