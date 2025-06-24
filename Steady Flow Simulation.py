import numpy as np
from tabulate import tabulate

def airfoil(cam,pos,thk,crd,pts):
    ang=np.linspace(0,np.pi,pts)
    x=0.5*crd*(1-np.cos(ang))
    ythk=thk*crd*((1.4845*(x/crd)**0.5)-(0.63*(x/crd))-(1.785*(x/crd)**2)+(1.4215*(x/crd)**3)-(0.5075*(x/crd)**4))
    fwd=x<pos*crd
    ycam=np.where(fwd,cam*x/(pos**2)*(2*pos-x/crd),cam*(crd-x)/((1-pos)**2)*(1+x/crd-2*pos))
    slp=np.where(fwd,2*cam/(pos**2)*(pos-x/crd),-2*cam/((1-pos)**2)*(pos-x/crd))
    ang=np.arctan(slp)
    xup,xlo=x-ythk*np.sin(ang),x+ythk*np.sin(ang)
    yup,ylo=ycam+ythk*np.cos(ang),ycam-ythk*np.cos(ang)
    xpts=np.concatenate([xup,xlo[::-1]])
    ypts=np.concatenate([yup,ylo[::-1]])
    return xpts,ypts

def Singularity_Element(xpts,ypts):
    N=len(xpts)-1
    xc,yc=(xpts[:-1]+xpts[1:])/2,(ypts[:-1]+ypts[1:])/2
    dx,dy=xpts[1:]-xpts[:-1],ypts[1:]-ypts[:-1]
    S=np.sqrt(dx**2+dy**2)
    theta=np.arctan2(dy,dx)
    return xc,yc,dx,dy,S,theta

def Grid_Generation(theta,S):
    nx,ny=-np.sin(theta),np.cos(theta)
    tx,ty=np.cos(theta),np.sin(theta)
    return nx,ny,tx,ty

def Influence_Coefficients(xc,yc,xpts,ypts,theta,S):
    N=len(xc)
    A=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i==j:
                A[i,j]=0.5
            else:
                dx=xc[i]-xpts[j]
                dy=yc[i]-ypts[j]
                r2=dx**2+dy**2
                beta=np.arctan2((yc[i]-ypts[j]),(xc[i]-xpts[j]))-theta[j]
                A[i,j]=(1/(2*np.pi))*np.cos(beta)*np.log(r2)
    return A

def RHSVector(nx,ny,U_inf,alpha):
    return -(U_inf*np.cos(alpha)*nx+U_inf*np.sin(alpha)*ny)

def Equation_Solver(A,RHS):
    print("\n Influence Coefficient Matrix A ")
    print(tabulate(np.round(A,4), tablefmt="fancy_grid"))
    print("\n RHS Vector ")
    rhs_table=[[val] for val in np.round(RHS,4)]
    print(tabulate(rhs_table, headers=["RHS"], tablefmt="fancy_grid"))
    gamma=np.linalg.solve(A,RHS)
    print("\n Solved Gamma Vector ")
    gamma_table=[[val] for val in np.round(gamma,4)]
    print(tabulate(gamma_table, headers=["Gamma"], tablefmt="fancy_grid"))
    return gamma

def final_computations(gamma,tx,ty,U_inf,alpha,xc,yc,S):
    Vt=U_inf*np.cos(alpha)*tx+U_inf*np.sin(alpha)*ty+gamma
    cp=1-(Vt/U_inf)**2
    dL=-cp*S
    dLx=dL*np.zeros_like(tx)
    dLy=dL
    total_Lx=np.sum(dLx)
    total_Ly=np.sum(dLy)
    moment=np.sum(xc*dLy - yc*dLx)
    print(f"\nTotal Lift = {total_Lx:.4f}i - {total_Ly:.4f}j N")
    print(f"Moment about leading edge = {moment:.4f} Nm")

U_inf=float(input("Enter freestream velocity (m/s): "))
alpha_deg=float(input("Enter angle of attack (degrees): "))
alpha=np.radians(alpha_deg)

cam=0.02
pos=0.4
thk=0.12
crd=1.0
pts=8

xpts,ypts=airfoil(cam,pos,thk,crd,pts)
xc,yc,dx,dy,S,theta=Singularity_Element(xpts,ypts)
nx,ny,tx,ty=Grid_Generation(theta,S)
A=Influence_Coefficients(xc,yc,xpts,ypts,theta,S)
RHS=RHSVector(nx,ny,U_inf,alpha)
gamma=Equation_Solver(A,RHS)
final_computations(gamma,tx,ty,U_inf,alpha,xc,yc,S)
