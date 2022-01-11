import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as NC
import sys
import glob

DrawEns=False

z=np.loadtxt('data/init/z.txt')

#nvar=17
#variables=['P1_c', 'P1_n', 'P1_p', 'P1_Chl', 'P1_s', 'P2_c', 'P2_n', 'P2_p', 'P2_Chl', 'P3_c', 'P3_n', 'P3_p', 'P3_Chl', 'P4_c', 'P4_n', 'P4_p', 'P4_Chl']

with open('data/init/names.txt') as f:
    lines = f.readlines()

nvar=0
variables=[]
for line in lines:
    for variable in line.split():
        if variable!="":
            nvar+=1
            variables.append(variable)
            
nz=z.size

#ens_size=10
ens_size=len(glob.glob("data/analysis/ens_??_ana.txt"))

secs=(31+29+15)*24*60*60

def readstate(infile, dimensions):
    
    statesize=1
    for dimension in dimensions:
        statesize*=dimension
    
    flatstate=np.zeros([statesize])

    with open(infile) as f:
        lines = f.readlines()

    c=0
    for line in lines:
        for number in line.split():
            if number!="":
                flatstate[c]=float(number)
                c+=1
                
    if not c==statesize:
        print("ERROR with ", infile," size! c=",c)
        sys.exit(1)
        
    return np.reshape(flatstate, dimensions)

chlindexes=[variables.index("P1_Chl"),variables.index("P2_Chl"),variables.index("P3_Chl"),variables.index("P4_Chl")]

def ens2state(ensemble):
    chl=np.sum(ensemble[:,chlindexes,:],axis=1)
    return np.mean(ensemble,axis=0), np.mean(chl,axis=0), np.std(chl, axis=0, ddof=1), chl

if ens_size:
    instate_ens=np.zeros([ens_size,nvar,nz])
    for i in range(ens_size):
        instate_ens[i]=readstate('data/forecast/ens_{:02d}.txt'.format(i+1), [nvar,nz])

    instate, inchl, instd, inchl_ens = ens2state(instate_ens)
    
    outstate_ens=np.zeros([ens_size,nvar,nz])
    for i in range(ens_size):
        outstate_ens[i]=readstate('data/analysis/ens_{:02d}_ana.txt'.format(i+1), [nvar,nz])

    outstate, outchl, outstd, outchl_ens = ens2state(outstate_ens)

else:
    instate=readstate('data/diag/state_init.txt', [nvar,nz])
    inchl=np.sum(instate[chlindexes,:],axis=0)
    
    outstate=readstate('data/analysis/state_ana.txt', [nvar,nz])
    outchl=np.sum(outstate[chlindexes,:],axis=0)
    
sat = np.loadtxt('data/obs/sat.txt')
sat_std = np.loadtxt('data/obs/sat_std.txt')
eofs = np.loadtxt('data/init/eof.txt')
climstd= np.linalg.norm(eofs, axis=0)

chl=np.array([inchl, outchl])

argo =readstate('data/obs/argo.txt', [1,nz]) 
argo_std = readstate('data/obs/argo_std.txt', [1,nz])

def writenc(outfile,statepre, statepost):
    with NC.netcdf_file(outfile,"w") as ncOUT:

        ncOUT.createDimension('time',2)
        ncOUT.createDimension('lev',nz)
        ncOUT.createDimension('lon',1)
        ncOUT.createDimension('lat',1)

        ncvar = ncOUT.createVariable('time','i',('time',))
        ncvar[:] = np.array([secs-1,secs])
        ncvar.units = 'seconds since 2016-01-01 00:00:00'
        ncvar = ncOUT.createVariable('lev','f',('lev',))
        ncvar[:] = z
        ncvar.units = 'm'
        ncvar = ncOUT.createVariable('lon','f',('lon',))
        ncvar[:] = np.array([0.0])
        ncvar.units = 'deg'
        ncvar = ncOUT.createVariable('lat','f',('lat',))
        ncvar[:] = np.array([0.0])
        ncvar.units = 'deg'
        for i,var in enumerate(variables):
            ncvar = ncOUT.createVariable(var,'f',('time',"lev","lat","lon"))
            ncvar[0,:,0,0] = statepre[i,:]
            ncvar[1,:,0,0] = statepost[i,:]

outfile = 'postproc/state_ana.nc'
writenc(outfile,instate, outstate)

def drawstate(fig, ax, variable, climstd, color, name, 
              varens=None, ensstd=None):
    
    lineensdev=None
    linehybriddev=None
    lineens=None
    
    line,=ax.plot(variable,z,color,label=name)
    lineclimdev,=ax.plot(variable+climstd,z,":"+color,label=name+' clim. st.dev.')
    ax.plot(variable-climstd,z,":"+color)
    
    if ens_size:
    
        lineensdev,=ax.plot(variable+ensstd,z,"--"+color,label=name+' ens. st.dev.')
        ax.plot(variable-ensstd,z,"--"+color)
        linehybriddev,=ax.plot(variable + np.sqrt(0.5*(ensstd**2+climstd**2)), z, "-."+color, label=name+' hybrid st.dev.')
        ax.plot(variable - np.sqrt(0.5*(ensstd**2+climstd**2)), z, "-."+color)

        if DrawEns:
            for i in range(ens_size):
                lineens,=ax.plot(varens[i],z,color,label=name+' ens. member', alpha=0.5)
                
    return (line, lineclimdev,
            lineensdev, linehybriddev, lineens)

fig,ax=plt.subplots()

(lineforecast, lineforecastclimdev,
 lineforecastensdev, lineforecasthybriddev, lineforecastens
 ) = drawstate(fig, ax, inchl, climstd, "b", "forecast",
               inchl_ens, instd)
          

(lineanal, lineanalclimdev,
 lineanalensdev, lineanalhybriddev, lineanalens
 )=drawstate(fig, ax, outchl, climstd, "r", "analysis",
             outchl_ens, outstd)
          

ax.plot(sat,[0.0],"^g")
ax.plot([sat-sat_std, sat+sat_std], [0.0,0.0], "g")
ax.plot([sat-sat_std], [0.0], "|g")
ax.plot([sat+sat_std], [0.0], "|g")

lineargo,=ax.plot(argo[0],z,"g",label='argo')
lineargodev,=ax.plot(argo[0]+argo_std[0],z,"--g",label='argo st.dev.')
ax.plot(argo[0]-argo_std[0],z,"--g")

plt.ylim([400.0,-40.0])
limx=plt.xlim()

ax.plot(limx,[0.0,0.0],"k",linewidth=1)
linesat,=ax.plot(sat,[0.0],"^g",label='sat')
ax.plot([sat-sat_std, sat+sat_std], [0.0,0.0], "g")
ax.plot([sat-sat_std], [0.0], "|g")
ax.plot([sat+sat_std], [0.0], "|g")
plt.xlim(limx)

plt.title("Chl")
plt.xlabel("Chl (mg/m^3)")
plt.ylabel("Depth (m)")

if ens_size:
    if DrawEns:
        plt.legend([lineforecast, lineforecastclimdev, lineforecastensdev, lineforecasthybriddev, lineforecastens,
                    lineanal, lineanalclimdev, lineanalensdev, lineanalhybriddev, lineanalens,
                    linesat, lineargo, lineargodev],
                   ['forecast', 'forecast clim. st.dev.', 'forecast ens. st.dev.', 'forecast hybrid st.dev.', "forecast ens. member",
                    'analysis', 'analysis clim. st.dev.', 'analysis ens. st.dev.', 'analysis hybrid st.dev.', "analysis ens. member",
                    'sat', 'argo', 'argo st.dev.'])
    else:
        plt.legend([lineforecast, lineforecastclimdev, lineforecastensdev, lineforecasthybriddev,
                    lineanal, lineanalclimdev, lineanalensdev, lineanalhybriddev,
                    linesat, lineargo, lineargodev],
                   ['forecast', 'forecast clim. st.dev.', 'forecast ens. st.dev.', 'forecast hybrid st.dev.',
                    'analysis', 'analysis clim. st.dev.', 'analysis ens. st.dev.', 'analysis hybrid st.dev.',
                    'sat', 'argo', 'argo st.dev.'])
else:
    plt.legend([lineforecast, lineforecastclimdev,
                lineanal, lineanalclimdev,
                linesat, lineargo, lineargodev],
               ['forecast', 'forecast clim. st.dev.',
                'analysis', 'analysis clim. st.dev.',
                'sat', 'argo', 'argo st.dev.'])

plt.savefig("postproc/chl.png")

plt.show()
