import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as NC
import sys

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

ens_size=10
DrawForecastEns=False

secs=(31+29+15)*24*60*60

def readstate(infile):

    flatstate=np.zeros([nz*nvar])

    with open(infile) as f:
        lines = f.readlines()

    c=0
    for line in lines:
        for number in line.split():
            if number!="":
                flatstate[c]=float(number)
                c+=1
                
    if not c==nz*nvar:
        print("ERROR with ", infile," size! c=",c)
        sys.exit(1)
        
    return np.reshape(flatstate, [nvar,nz])

#instate=readstate('data/forecast/phyto.txt')
instate_ens=np.zeros([ens_size,nvar,nz])
for i in range(ens_size):
    instate_ens[i]=readstate('data/forecast/ens_{:02d}.txt'.format(i+1))
    
chlindexes=[variables.index("P1_Chl"),variables.index("P2_Chl"),variables.index("P3_Chl"),variables.index("P4_Chl")]

def ens2state(ensemble):
    chl=np.sum(ensemble[:,chlindexes,:],axis=1)
    return np.mean(ensemble,axis=0), np.mean(chl,axis=0), np.std(chl, axis=0, ddof=1), chl

instate, inchl, instd, inchl_ens = ens2state(instate_ens)
outstate=readstate('data/analysis/state_ana.txt')
    

obs = np.loadtxt('data/obs/sat.txt')
obsdev = np.loadtxt('data/obs/sat_std.txt')
eofs = np.loadtxt('data/init/eof.txt')
stdev= np.linalg.norm(eofs, axis=0)

chl=np.array([inchl, np.sum(outstate[chlindexes,:],axis=0)])


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


fig,ax=plt.subplots()

if DrawForecastEns:
    for i in range(ens_size):
        lineforecastens,=ax.plot(inchl_ens[i],z,"b",label='forecast ens. member', alpha=0.5)

lineforecast,=ax.plot(chl[0],z,"b",label='forecast')
lineforecastclimdev,=ax.plot(chl[0]+stdev,z,":b",label='forecast clim. st.dev.')
ax.plot(chl[0]-stdev,z,":b")
lineforecastensdev,=ax.plot(chl[0]+instd,z,"--b",label='forecast ens. st.dev.')
ax.plot(chl[0]-instd,z,"--b")
lineforecasthybriddev,=ax.plot(chl[0] + np.sqrt(0.5*(instd**2+stdev**2)), z, "-.b", label='forecast hybrid st.dev.')
ax.plot(chl[0] - np.sqrt(0.5*(instd**2+stdev**2)), z, "-.b")
lineanal,=ax.plot(chl[1],z,"r",label='analysis')
#ax.plot(chl[1]+stdev,z,":r")
#ax.plot(chl[1]-stdev,z,":r")
ax.plot(obs,[0.0],"^g")
ax.plot([obs-obsdev, obs+obsdev], [0.0,0.0], "g")
ax.plot([obs-obsdev], [0.0], "|g")
ax.plot([obs+obsdev], [0.0], "|g")
plt.ylim([400.0,-40.0])
limx=plt.xlim()
ax.plot(limx,[0.0,0.0],"k",linewidth=1)
lineobs,=ax.plot(obs,[0.0],"^g",label='obs')
ax.plot([obs-obsdev, obs+obsdev], [0.0,0.0], "g")
ax.plot([obs-obsdev], [0.0], "|g")
ax.plot([obs+obsdev], [0.0], "|g")
plt.xlim(limx)

plt.title("Chl")
plt.xlabel("Chl (mg/m^3)")
plt.ylabel("Depth (m)")

if DrawForecastEns:
    plt.legend([lineforecast,lineforecastens,       lineforecastclimdev,      lineforecastensdev,      lineforecasthybriddev,    lineobs, lineanal],
               ['forecast', "forecast ens. member", 'forecast clim. st.dev.', 'forecast ens. st.dev.', 'forecast hybrid st.dev.', 'obs','analysis'])
else:
    plt.legend([lineforecast, lineforecastclimdev,    lineforecastensdev,      lineforecasthybriddev,    lineobs, lineanal],
               ['forecast', 'forecast clim. st.dev.', 'forecast ens. st.dev.', 'forecast hybrid st.dev.', 'obs','analysis'])

plt.savefig("postproc/chl.png")

plt.show()
