from pyiron_base import TemplateJob, ImportAlarm, Project
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from string import Template

### TODO 

"""
Quantities to be obtained from lower length scale simulation

1. Isotropic elastic constants (mu and nu) (from first principles)
2. Mobilities of screw and edge dislocation (from MD simuations)

"""

class ParaDis(TemplateJob):

    def __init__(self, project, job_name):
        super(ParaDis, self).__init__(project, job_name)
        self._ctrlfilename  = None 
        self._datafilename  = None
        self._ctrlargs = None
        self._dataargs = None
        self._dataparams = {}
        self._ctrlparams = {}
        self._write_default_datafile = 0
        self.SetDefaultParams()

    def validate_ready_to_run(self):
        #        if (self._write_default_datafile):
        #    args = f"paradisgen -type {self._dataparams['type']} -cubel {self._dataparams['cubel']} -nfrsrcs {self._dataparams['nfrsrcs']} -frlen {self._dataparams['frlen'][0]} -maxseg {self._dataparams['maxseglen']} -outfile {self._datafilename} > paradisgen.log"
        #    subprocess.run(args, shell=True, capture_output=True)
            #subprocess.run("cat paradisgen.log", shell=True, capture_output=True)

        if self._datafilename is None:
            raise ValueError("job._datafilename is not set")
        if self._ctrlfilename is None:
            raise ValueError("job._ctrlfilename is not set")
        self.executable = f"paradis -d {self._datafilename} {self._ctrlfilename} > paradis.log"

    def write_input(self):
        self.write_ctrlfile()
        self.write_datafile()

    @staticmethod
    def translate_dict(d):
        """
        Convert a dictionary to a format that can be written to a file

        Args:
            d (dict): Dictionary to be converted
        """
        d_result = {}
        for key, value. in d.items():
            # if value is list or numpy array, convert it to a string
            if isinstance(value, (list, np.ndarray)):
                for ii, item in enumerate(value):
                    d_result[f"{key}_{ii}"] = item
            else:
                d_result[key] = value
        return d_result

    def write_ctrlfile(self):
        with open(os.path.join(self.working_directory, self._ctrlfilename), "w") as f:
            template = Template(f.read())

        # Substitute the placeholders with actual values
        output_str = template.substitute(self.translate_dict(self._ctrlparams))

        with open(os.path.join(self.working_directory, self._ctrlfilename), "w") as f:
            f.write(output_str)

    def write_datafile(self):
        #file_path = os.path.join(self.working_directory,"frank_read_src.data")
        if (not self._write_default_datafile):

            str = """
dataFileVersion =   4
numFileSegments =   1
minCoordinates = [
    -1.750000e+04
    -1.750000e+04
    -1.750000e+04
]
maxCoordinates = [
    1.750000e+04
    1.750000e+04
    1.750000e+04
]
nodeCount =   3
dataDecompType =   1
dataDecompGeometry = [
    1
    1
    1
    ]
#
#  END OF DATA FILE PARAMETERS
#

domainDecomposition =
     -17500.0000
     -17500.0000
     -17500.0000
      17500.0000
      17500.0000
      17500.0000
 nodalData =
#  Primary lines: node_tag, x, y, z, num_arms, constraint
#  Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz
0,0      -4000.0000   500.0000        6000.0000   1   7
        0,1       0.5773503000     0.5773503000    -0.5773503000
                    0.              1.              1.
0,1       17.0000      500.0000        6000.0000  2   0
        0,0     -0.5773503000    -0.5773503000     0.5773503000
                    0.              1.              1.
        0,2      0.5773503000     0.5773503000    -0.5773503000
                    0.              1.              1.
            0,2       4000.0000    500.0000        6000.0000   1   7
        0,1     -0.5773503000    -0.5773503000     0.5773503000
                0.              1.              1.
"""
            self._dataargs = str
            with open(os.path.join(self.working_directory, self._datafilename), "w") as f:
                f.write(self._dataargs)

        else:
           
            #    pass
            os.chdir(self.working_directory)
            args = f"paradisgen -type {self._dataparams['type']} -cubel {self._dataparams['cubel']} -nfrsrcs {self._dataparams['nfrsrcs']} -frlen {self._dataparams['frlen'][0]} -maxseg {self._dataparams['maxseglen']}"
            subprocess.run(args, shell=True, capture_output=True)
            subprocess.run(f"mv paradis.data {self._datafilename}", shell=True, capture_output=True)



    def SetDefaultParams(self):
        
        #Setting data file parameters
        self._dataparams["cubel"] = 4000
        self._dataparams["hexl"] = 50
        self._dataparams["type"] = "frank-read-src"
        self._dataparams["npf"] = 3
        self._dataparams["nfrsrcs"] = 1
        self._dataparams["frlen"] = [400,500] # [min, max] 
        self._dataparams["seed"] = -784574883
        self._dataparams["nchains"] = 20
        self._dataparams["maxseglen"] = 200
        self._dataparams["nloops"] = 10 
        self._dataparams["radius"] = self._dataparams["maxseglen"] / 2.0
        self._dataparams["xsurf"] = [0,0]
        self._dataparams["ysurf"] = [0,0]
        self._dataparams["zsurf"] = [0,0]
        self._dataparams["outfile"] = "frank_read_src.data"
        self._dataparams["cOVERa"] = 1.568 # for beryllium

        #Setting control file parameters 
        self._ctrlparams["dirname"] = "results"
        self._ctrlparams["numdomains"] = [1, 1, 1] #if more than 1 then parallel version of paradis should be compiled
        self._ctrlparams["numcells"] = [4,4,4]

        ## Fast multipole parameters 
        self._ctrlparams["fmEnabled"] = 1 #toggle 1: yes 0: no
        self._ctrlparams["fmMPOrder"] = 2 
        self._ctrlparams["fmTaylorOrder"] = 5 
        self._ctrlparams["fmCorrectionTbl"] = "../../../../ParaDiS/inputs/fm-ctab.Ta.600K.0GPa.m2.t5.dat"
        self._ctrlparams["fmEnableAnisotropy"] = 0 # by default disabled 

        ## Discretization and topological changes parameters  
        self._ctrlparams["remeshRule"] = 2 
        self._ctrlparams["minSeg"] = 100.0 
        self._ctrlparams["maxSeg"] = 300.0
        self._ctrlparams["splitMultiNodeFreq"] = 1 
        self._ctrlparams["collisionMethod"] = 2 

        ## Simulation timestepping controls 
        self._ctrlparams["cycleStart"] = 0 
        self._ctrlparams["timeNow"] = 0.0 
        self._ctrlparams["timeStart"] = 0.0 
        self._ctrlparams["timestepIntegrator"] = "trapezoid" 
        self._ctrlparams["trapezoidMaxIterations"] = 2 
        self._ctrlparams["deltaTT"] = 0.0 
        self._ctrlparams["maxDT"] = 1.0e-07 
        self._ctrlparams["nextDT"] = 1.0e-07 
        self._ctrlparams["dtIncrementFact"] = 1.2
        self._ctrlparams["dtDecrementFact"] = 0.5
        self._ctrlparams["dtExponent"] = 4.0 
        self._ctrlparams["dtVariableAdjustment"] = 1 
        self._ctrlparams["maxstep"] = 500 
        self._ctrlparams["rTol"] = 1.0  
        self._ctrlparams["rmax"] = 100 

        ## Materials and Mobility parameters
        self._ctrlparams["shearModulus"] = 6.488424e+10 #shear modulus units=Pa for Ta temp 300K p = 0GPa 
        self._ctrlparams["pois"] = 3.327533e-01 # Poisson's ratio for Ta temp 300K p=0 GPa  
        self._ctrlparams["YoungModulus"] = 172.95e+09 # Youngs modulus unit=Pa for Ta temp=300K p=0GPa 
        self._ctrlparams["burgMag"] = 2.875401e-10  # Burgers vector unit=m temp 300 K p =0 GPa 
        self._ctrlparams["rc"] = 5 #core radius of dislocation in units of Burgers vector magnitude  
        self._ctrlparams["MobScrew"] = 10.0  # unit = 1/(Pa.s)
        self._ctrlparams["MobEdge"] = 10.0   # unit = 1/(Pa.s)
        self._ctrlparams["mobilityLaw"] = "BCC_0"

        ## Loading parameters 
        self._ctrlparams["loadType"] = 0 
        self._ctrlparams["edotdir"] = [0.0, 0.0, 1.0]
        self._ctrlparams["appliedStress"] = [0.0, 0.0, 5.0e+08, 0.0, 0.0, 0.0]
        self._ctrlparams["TempK"] = 600.0 # in Kelvins 
        self._ctrlparams["pressure"] = 0.0
        self._ctrlparams["eRate"] = 1.0
        self._ctrlparams["indxErate"] = 0 
        self._ctrlparams["useLabFrame"] = 1 
        self._ctrlparams["rotationMatrix"] = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]

        ## I/O parameters 

        ### arm files 
        self._ctrlparams["armfile"] = 0 #by default armfiles are not written. toggle on to enable it. 
        self._ctrlparams["armfilefreq"] = 100  
        self._ctrlparams["armfiledt"] = -1  
        self._ctrlparams["armfiletime"] = 0  
        self._ctrlparams["armfilecounter"] = 0  

        ### flux files 
        self._ctrlparams["fluxfile"] = 1 #by default fluxfiles are written 
        self._ctrlparams["fluxfreq"] = 100  
        self._ctrlparams["fluxdt"] = -1  
        self._ctrlparams["fluxtime"] = 0  
        self._ctrlparams["fluxcounter"] = 0  

        ### gnuplot files 
        self._ctrlparams["gnuplot"] = 1 
        self._ctrlparams["gnuplotfreq"] = 100 
        self._ctrlparams["gnuplotdt"] = -1 
        self._ctrlparams["gnuplottime"] = 0 
        self._ctrlparams["gnuplotcounter"] = 0 

        ### restart files 
        self._ctrlparams["savecn"] = 1 
        self._ctrlparams["savecnfreq"] = 100 
        self._ctrlparams["savecndt"] = -1 
        self._ctrlparams["savecntime"] = 0 
        self._ctrlparams["savecncounter"] = 0 

        ### povray files 
        self._ctrlparams["povray"] = 1 
        self._ctrlparams["povrayfreq"] = 100 
        self._ctrlparams["povraydt"] = -1 
        self._ctrlparams["povraytime"] = 0 
        self._ctrlparams["povraycounter"] = 0 

        ### properties files 
        self._ctrlparams["saveprop"] = 1 
        self._ctrlparams["savepropfreq"] = 2
        self._ctrlparams["savepropdt"] = -1
        self._ctrlparams["saveproptime"] = 0

        ### tecplot files 
        self._ctrlparams["tecplot"] = 0 
        self._ctrlparams["tecplotfreq"] = 100 
        self._ctrlparams["tecplotdt"] = -1 
        self._ctrlparams["tecplottime"] = 0 
        self._ctrlparams["tecplotcounter"] = 0 

        ### polefigure files 
        self._ctrlparams["polefigfile"] = 0 
        self._ctrlparams["polefigfreq"] = 100 
        self._ctrlparams["polefigdt"] = -1 
        self._ctrlparams["polefigtime"] = 0 
        self._ctrlparams["polefigcounter"] = 0 

        ### nodal velocity files 
        self._ctrlparams["velfile"] = 0 
        self._ctrlparams["velfilefreq"] = 100 
        self._ctrlparams["velfiledt"] = -1 
        self._ctrlparams["velfiletime"] = 0 
        self._ctrlparams["velfilecounter"] = 0 

        ### nodal force files 
        self._ctrlparams["writeForce"] = 0 
        self._ctrlparams["writeForceFreq"] = 100 
        self._ctrlparams["writeForceDT"] = -1 
        self._ctrlparams["writeForceTime"] = 0 
        self._ctrlparams["writeForceCounter"] = 0 

        ### postscript files 
        self._ctrlparams["psfile"] = 0 
        self._ctrlparams["psfilefreq"] = 100 
        self._ctrlparams["psfiledt"] = -1 
        self._ctrlparams["psfiletime"] = 0 

    def collect_output (self):
        pass

    def plot_results(self, properties):

        numofplots = len(properties)
        fig, ax = plt.subplots(numofplots)
        plotnum = 0

        for propid in properties:
            if (propid == "density"):

                data = np.loadtxt(self.working_directory + "/results" + "/properties/density")
                plstr = data[:,0]
                density = data[:,2]
                ax[plotnum].set_xlabel(r'$\varepsilon_p$', fontsize=16)
                ax[plotnum].set_ylabel(r'$\rho (m^{-2})$', fontsize=16)
                ax[plotnum].plot(plstr, density, linestyle='-', linewidth='2.5')
                plotnum = plotnum + 1

            elif (propid == "epsdot"):

                data = np.loadtxt(self.working_directory + "/results" + "/properties/epsdot")
                time = data[:,0]
                epsdot = data[:,1]
                ax[plotnum].set_xlabel(r'time', fontsize=16)
                ax[plotnum].set_ylabel(r'$\dot{\epsilon} (s^{-1})$', fontsize=16)
                ax[plotnum].plot(time, epsdot, linestyle='-', linewidth='2.5')
                plotnum = plotnum + 1

            elif (propid == "time_pl_str"):
                data = np.loadtxt(self.working_directory + "/results" + "/properties/time_Plastic_strain")
                time = data[:,0]
                plstr = data[:,1]
                ax[plotnum].set_xlabel(r'time', fontsize=16)
                ax[plotnum].set_ylabel(r'$\varepsilon_p$', fontsize=16)
                ax[plotnum].plot(time, plstr, linestyle='-', linewidth='2.5')
                plotnum = plotnum + 1

        fig.tight_layout()
