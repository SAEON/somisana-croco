import numpy as np
import netCDF4 as netcdf
import sys

#######################################################
#Transfert a field at rho points to psi points
#######################################################
def rho2psi(var_rho):
    '''
    Root function to transfert croco variable from horizontal
    CROCO rho to psi coordinates

    Input:
      var_rho        var (2D,3D or 4D) on rho-grid

    Output:
      var_psi        var on psi-grid
    '''

    if np.ndim(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)
    return var_psi
###########
def rho2psi_2d(var_rho):
    '''
    Convert 2D variables from CROCO rho to psi grid

    Input:
      var_rho       2D var on rho-grid

    Output:
      var_psi       2D var on psi-grid
    '''

    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])
    return var_psi
###########
def rho2psi_3d(var_rho):
    '''
    Convert 3D variables from CROCO rho to psi grid

    Input:
      var_rho       3D var on rho-grid

    Output:
      var_psi       3D var on psi-grid
    '''
    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])
    return var_psi
###########
def rho2psi_4d(var_rho):
    '''
    Convert 4D variables from CROCO rho to psi grid

    Input:
      var_rho       4D var on rho-grid

    Output:
      var_psi       4D var on psi-grid
    '''
    var_psi = 0.25*(var_rho[:,:,1:,1:]+var_rho[:,:,1:,:-1]+var_rho[:,:,:-1,:-1]+var_rho[:,:,:-1,1:])
    return var_psi

#######################################################
#Transfert a field at rho points to u points
#######################################################
def rho2u(var_rho):
    '''
    Root function to transfert croco variable from horizontal
    CROCO rho to u coordinates

    Input:
      var_rho        var (2D,3D or 4D) on rho-grid

    Output:
      var_u          var on u-grid
    '''

    if np.ndim(var_rho)==1:
        var_u = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:
        var_u = rho2u_2d(var_rho)
    elif np.ndim(var_rho)==3:
        var_u = rho2u_3d(var_rho)
    else:
        var_u = rho2u_4d(var_rho)
    return var_u
###########
def rho2u_2d(var_rho):
    '''
    Convert 2D variables from CROCO rho to u grid

    Input:
      var_rho       2D var on rho-grid

    Output:
      var_u         2D var on u-grid
    '''
    var_u = 0.5*(var_rho[:,1:]+var_rho[:,:-1])
    return var_u
###########
def rho2u_3d(var_rho):
    '''
    Convert 3D variables from CROCO rho to u grid

    Input:
      var_rho       3D var on rho-grid

    Output:
      var_u         3D var on u-grid
    '''

    var_u = 0.5*(var_rho[:,:,1:]+var_rho[:,:,:-1])
    return var_u
###########
def rho2u_4d(var_rho):
    '''
    Convert 4D variables from CROCO rho to u grid

    Input:
      var_rho       4D var on rho-grid

    Output:
      var_u         4D var on u-grid
    '''

    var_u = 0.5*(var_rho[:,:,:,1:]+var_rho[:,:,:,:-1])
    return var_u

#######################################################
#Transfert a field at rho points to v points
#######################################################
def rho2v(var_rho):
    '''
    Root function to transfert croco variable from horizontal
    CROCO rho to v coordinates

    Input:
      var_rho        var (2D,3D or 4D) on rho-grid

    Output:
      var_v          var on v-grid
    '''

    if np.ndim(var_rho)==1:
        var_v = 0.5*(var_rho[1:]+var_rho[:-1])
    elif np.ndim(var_rho)==2:
        var_v = rho2v_2d(var_rho)
    elif np.ndim(var_rho)==3:
        var_v = rho2v_3d(var_rho)
    else:
        var_v = rho2v_4d(var_rho)
    return var_v
###########
def rho2v_2d(var_rho):
    '''
    Convert 2D variables from CROCO rho to v grid

    Input:
      var_rho       2D var on rho-grid

    Output:
      var_v         2D var on v-grid
    '''
    var_v = 0.5*(var_rho[1:,:]+var_rho[:-1,:])
    return var_v
###########
def rho2v_3d(var_rho):
    '''
    Convert 3D variables from CROCO rho to v grid

    Input:
      var_rho       3D var on rho-grid

    Output:
      var_v         3D var on v-grid
    '''
    var_v = 0.5*(var_rho[:,1:,:]+var_rho[:,:-1,:])
    return var_v
###########
def rho2v_4d(var_rho):
    '''
    Convert 4D variables from CROCO rho to v grid

    Input:
      var_rho       4D var on rho-grid

    Output:
      var_v         4D var on v-grid
    '''
    var_v = 0.5*(var_rho[:,:,1:,:]+var_rho[:,:,:-1,:])
    return var_v

#######################################################
#Transfert a field at u points to the rho points
#######################################################
def v2rho(var_v):
    '''
    Root function to transfert croco variable from horizontal
    CROCO v to rho coordinates

    Input:
      var_v        var (2D,3D or 4D) on v-grid

    Output:
      var_rho      var on rho-grid
    '''

    if np.ndim(var_v)<3:
        var_rho = v2rho_2d(var_v)
    elif np.ndim(var_v)==3:
        var_rho = v2rho_3d(var_v)
    else:
        var_rho = v2rho_4d(var_v)
    return var_rho
###########
def v2rho_2d(var_v):
    '''
    Convert 2D variables from CROCO v to rho grid

    Input:
      var_v         2D var on v-grid

    Output:
      var_rho       2D var on rho-grid
    '''

    [L,Mp]=var_v.shape
    Lp=L+1
    Lm=LM-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[1:L,:]=0.5*(var_v[0:Lm:,:]+var_v[1:L,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Lp-1,:]=var_rho[L-1,:]
    return var_rho
###########
def v2rho_3d(var_v):
    '''
    Convert 3D variables from CROCO v to rho grid

    Input:
      var_v         3D var on v-grid

    Output:
      var_rho       3D var on rho-grid
    '''

    [N,L,Mp]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Lp,Mp))
    var_rho[:,1:L,:]=0.5*(var_v[:,0:Lm,:]+var_v[:,1:L,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Lp-1,:]=var_rho[:,L-1,:]
    return var_rho
##########
def v2rho_4d(var_v):
    '''
    Convert 4D variables from CROCO v to rho grid

    Input:
      var_v         4D var on v-grid

    Output:
      var_rho       4D var on rho-grid
    '''
    [T,N,L,Mp]=var_v.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((T,N,Lp,Mp))
    var_rho[:,:,1:L,:]=0.5*(var_v[:,:,0:Lm,:]+var_v[:,:,1:L,:])
    var_rho[:,:,0,:]=var_rho[:,:,1,:]
    var_rho[:,:,Lp-1,:]=var_rho[:,:,L-1,:]
    return var_rho

#######################################################
#Transfert a 2 or 2-D field at u points to the rho points
#######################################################
def u2rho(var_u):
    '''
    Root function to transfert croco variable from horizontal
    CROCO u to rho coordinates

    Input:
      var_u        var (2D,3D or 4D) on u-grid

    Output:
      var_rho      var on rho-grid
    '''
    if np.ndim(var_u)==1:
        [M]=var_u.shape
        Mp=M+1
        Mm=M-1
        var_rho=np.zeros((Mp))
        var_rho[1:M]=0.5*(var_u[0:Mm]+var_u[1:M])
        var_rho[0]=var_rho[1]
        var_rho[Mp-1]=var_rho[M-1]
    elif np.ndim(var_u)==2:
        var_rho = u2rho_2d(var_u)
    elif np.ndim(var_u)==3:
        var_rho = u2rho_3d(var_u)
    else:
        var_rho = u2rho_4d(var_u)
    return var_rho
##########
def u2rho_2d(var_u):
    '''
    Convert 2D variables from CROCO u to rho grid

    Input:
      var_u         2D var on u-grid

    Output:
      var_rho       2D var on rho-grid
    '''
    [Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[:,1:M]=0.5*(var_u[:,0:Mm]+var_u[:,1:M])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Mp-1]=var_rho[:,M-1]
    return var_rho
##########
def u2rho_3d(var_u):
    '''
    Convert 3D variables from CROCO u to rho grid

    Input:
      var_u         3D var on u-grid

    Output:
      var_rho       3D var on rho-grid
    '''

    [N,Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Lp,Mp))
    var_rho[:,:,1:M]=0.5*(var_u[:,:,0:Mm]+var_u[:,:,1:M])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Mp-1]=var_rho[:,:,M-1]
    return var_rho
##########
def u2rho_4d(var_u):
    '''
    Convert 4D variables from CROCO u to rho grid

    Input:
      var_u         4D var on psi-grid

    Output:
      var_rho       4D var on rho-grid
    '''

    [T,N,Lp,M]=var_u.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((T,N,Lp,Mp))
    var_rho[:,:,:,1:M]=0.5*(var_u[:,:,:,0:Mm]+var_u[:,:,:,1:M])
    var_rho[:,:,:,0]=var_rho[:,:,:,1]
    var_rho[:,:,:,Mp-1]=var_rho[:,:,:,M-1]
    return var_rho

#######################################################
#Transfert a field at psi points to rho points
#######################################################
def psi2rho(var_psi):
    '''
    Root function to transfert croco variable from horizontal
    CROCO psi to rho coordinates

    Input:
      var_psi      var (2D,3D or 4D) on psi-grid

    Output:
      var_rho      var on rho-grid
    '''
    if np.ndim(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    elif np.ndim(var_psi)==3:
        var_rho = psi2rho_3d(var_psi)
    else:
        var_rho = psi2rho_4d(var_psi)
    return var_rho
###########
def psi2rho_2d(var_psi):
    '''
    Convert 2D variables from CROCO psi to rho grid

    Input:
      var_psi       2D var on psi-grid

    Output:
      var_rho       2D var on rho-grid
    '''
    [L,M]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1
    var_rho=np.zeros((Lp,Mp))
    var_rho[1:L,1:M]=0.25*(var_psi[0:Lm,0:Mm]+var_psi[0:Lm,1:M]+var_psi[1:L,0:Mm]+var_psi[1:L,1:M])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Mp-1]=var_rho[:,M-1]
    var_rho[0,:]=var_rho[1,:]
    var_rho[Lp-1,:]=var_rho[L-1,:]
    return var_rho
###########
def psi2rho_3d(var_psi):
    '''
    Convert 3D variables from CROCO psi to rho grid

    Input:
      var_psi       3D var on psi-grid

    Output:
      var_rho       3D var on rho-grid
    '''
    [Nz,Lz,Mz]=var_psi.shape
    var_rho=np.zeros((Nz,Lz+1,Mz+1))
    for iz in range(0, Nz, 1):
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])
    return var_rho
##########
def psi2rho_4d(var_psi):
    '''
    Convert 4D variables from CROCO psi to rho grid

    Input:
      var_psi       4D var on psi-grid

    Output:
      var_rho       4D var on rho-grid
    '''
    [Tz,Nz,Lz,Mz]=var_psi.shape
    var_rho=np.zeros((Tz,Nz,Lz+1,Mz+1))
    for it in range(0, Tz, 1):
        var_rho[it,:]=psi2rho_3d(var_psi[it,:])
    return var_rho


#######################################################
#Transfert a 3-D field from verical w points to vertical rho-points
#######################################################
def w2rho(var_w):
    '''
    Transfert 3D field from vertical w points to rho points

    Input:
      var_w       3D var on vertical w-points

    Output:
      var_rho     3d var on vertical rho-points
    '''
    [N,L,M]=var_w.shape
    print( '[N,L,M]',[N,L,M])
    var_rho = np.zeros((N-1,L,M))
    for iz in range(1,N-2):
        var_rho[iz,:,:]  = 0.5625*(var_w[iz+1,:,:] + var_w[iz,:,:]) -0.0625*(var_w[iz+2:,:] + var_w[iz-1,:,:])
    var_rho[0,:,:]  = -0.125*var_w[1,:,:] + 0.75*var_w[1,:,:] +0.375*var_w[0,:,:]
    var_rho[N-2,:,:]  = -0.125*var_w[N-3,:,:] + 0.75*var_w[N-2,:,:] +0.375*var_w[N-1,:,:]
    return var_rho


##############################################################
##############################################################
##############################################################

