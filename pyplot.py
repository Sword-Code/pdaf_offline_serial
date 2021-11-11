import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as NC

z=np.loadtxt('data/z.txt')

nz=z.size
nvar=17
variables=['P1_c', 'P1_n', 'P1_p', 'P1_Chl', 'P1_s', 'P2_c', 'P2_n', 'P2_p', 'P2_Chl', 'P3_c', 'P3_n', 'P3_p', 'P3_Chl', 'P4_c', 'P4_n', 'P4_p', 'P4_Chl']
obsdev=0.1

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
        print("ERROR! c=",c)
        
    return np.reshape(flatstate, [nvar,nz])

instate=readstate('data/phyto.txt')
outstate=readstate('state_ana.txt')

obs = np.loadtxt('data/sat.txt')
eofs = np.loadtxt('data/eof.txt')
stdev= np.linalg.norm(eofs, axis=0)

chl=np.array([np.sum(instate[[3,8,12,16],:],axis=0),np.sum(outstate[[3,8,12,16],:],axis=0)])


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

outfile = 'state_ana.nc'
writenc(outfile,instate, outstate)


fig,ax=plt.subplots()
lineforcast=ax.plot(chl[0],z,"b",label='forecast')
lineforcastdev=ax.plot(chl[0]+stdev,z,":b",label='forecast st.dev.')
ax.plot(chl[0]-stdev,z,":b")
lineanal=ax.plot(chl[1],z,"r",label='analysis')
#ax.plot(chl[1]+stdev,z,":r")
#ax.plot(chl[1]-stdev,z,":r")
ax.plot(obs,[0.0],"^g")
ax.plot([obs-obsdev, obs+obsdev], [0.0,0.0], "g")
ax.plot([obs-obsdev], [0.0], "|g")
ax.plot([obs+obsdev], [0.0], "|g")
plt.ylim([400.0,-40.0])
limx=plt.xlim()
ax.plot(limx,[0.0,0.0],"k",linewidth=1)
lineobs=ax.plot(obs,[0.0],"^g",label='obs')
ax.plot([obs-obsdev, obs+obsdev], [0.0,0.0], "g")
ax.plot([obs-obsdev], [0.0], "|g")
ax.plot([obs+obsdev], [0.0], "|g")
plt.xlim(limx)

plt.title("Chl")
plt.xlabel("Chl (g/m^3)")
plt.ylabel("Depth (m)")

ax.legend([lineforcast,lineforcastdev,lineobs,lineanal],['forecast','forecast st.dev.','obs','analysis'])

plt.savefig("chl.png")

plt.show()
