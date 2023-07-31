import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

pos = np.array([]) #atomic positions in the supercell
f2 = open('pos.dat','r')
for line in f2:
	pos = np.append(pos,[float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
f2.close()

alat = 3.563 #lattice constant
c11 = 252.0  #elastic constant
c12 = 158.0  #elastic constant
c44 = 142.0  #elastic constant
mass = (58.6934*58.933195*51.9961)**(1/3.0) #averaged atomic mass
n = 30 #grid points along each direction
N = n**3
Basis = 2 * np.pi/alat * np.array([[-1,1,1],[1,-1,1],[1,1,-1]])

K0 = np.array([]) #coordinates of k-points in reciprocal space
for i in range(n):
	for j in range(n):
		for k in range(n):
			K0 = np.append(K0,[i/float(n)-0.5,j/float(n)-0.5,k/float(n)-0.5])

K0 = np.resize(K0,(int(np.size(K0)/3),3))
K = np.dot(K0, Basis)

n1 = 4
N1 = 4*n1**3-1
R = n1 * pos
tmp = np.array([])
for coords in R:
	if abs(coords) >= abs(n1 - coords):
		tmp = np.append(tmp, coords-n1)
	else:
		tmp = np.append(tmp, coords)

R = alat * tmp
R = np.resize(R, (int(len(R)/3),3))

KF1 = np.zeros((N1,3)) #Kanzaki forces of Ni in real space
KF2 = np.zeros((N1,3)) #Kanzaki forces of Co in real space
KF3 = np.zeros((N1,3)) #Kanzaki forces of Cr in real space

KF1[1] = [-0.00112078,  0.11037717,  0.11008989]
KF1[2] = [ 0.08894672,  0.00495242,  0.07261606]
KF1[0] = [ 0.09692353,  0.10657994,  0.00754831]
KF1[13] = [-0.00871956,  0.10056483, -0.09210733]
KF1[14] = [ 0.10127122,  0.00267614, -0.08393542]
KF1[48] = [ 0.09962564, -0.10456461,  0.00519078]
KF1[49] = [-0.00443375, -0.10912036,  0.09480789]
KF1[194] = [-0.09512364, -0.00214386,  0.08974644]
KF1[192] = [-0.09467292,  0.10630453, -0.01923344]
KF1[61] = [-0.00490631, -0.10251756, -0.10031372]
KF1[206] = [-0.07584981,  0.00161756, -0.07777872]
KF1[240] = [-0.10942536, -0.11710836, -0.00040389]


KF2[1] = [0.00643989,  0.04120942,  0.04806208]
KF2[2] = [0.03625378,  0.00666708,  0.04511056]
KF2[0] = [0.08846311,  0.075158,    0.00929356]
KF2[13] = [-0.00885731,  0.07281708, -0.06532908]
KF2[14] = [0.05568353,  0.00172036, -0.05051525 ]
KF2[48] = [0.07539781, -0.07617125,  0.01548419]
KF2[49] = [0.01198822, -0.069614,    0.05844214]
KF2[194] = [-0.05313919,  0.00702181,  0.04735992]
KF2[192] = [-0.08440217,  0.07616861,  0.00274569]
KF2[61] = [0.00264881, -0.07065664, -0.045429]
KF2[206] = [-0.06850808,  0.00612561, -0.04706806]
KF2[240] = [-0.08048906, -0.07192669, -0.02363325]

KF3[1] = [-0.01011914,  0.08665875,  0.07892783]
KF3[2] = [0.09603664, -0.01736403,  0.07735933]
KF3[0] = [0.06022139,  0.08664325, -0.000562]
KF3[13] = [-0.01359422,  0.06517028, -0.06905433]
KF3[14] = [0.08675697, -0.03178975, -0.10122611]
KF3[48] = [0.06530258, -0.07130444, -0.00746328]
KF3[49] = [-0.00277039, -0.08307367,  0.07683253]
KF3[194] = [-0.08532894,  0.01103203,  0.08438619]
KF3[192] = [-0.05985303,  0.09669425,  0.00522111]
KF3[61] = [0.00315786, -0.06904675, -0.08754217]
KF3[206] = [-0.05595583, -0.00613006, -0.06414894]
KF3[240] = [-0.07680408, -0.07661169, -0.01173422]

F1q = np.zeros((N,3), dtype=complex) #Kanzaki forces of Ni in reciprocal space
F2q = np.zeros((N,3), dtype=complex) #Kanzaki forces of Co in reciprocal space
F3q = np.zeros((N,3), dtype=complex) #Kanzaki forces of Cr in reciprocal space

for q in range(N):
	for i in range(N1):
		F1q[q] = F1q[q] + KF1[i] * np.exp(-1j * np.dot(K[q],R[i]))
		F2q[q] = F2q[q] + KF2[i] * np.exp(-1j * np.dot(K[q],R[i]))
		F3q[q] = F3q[q] + KF3[i] * np.exp(-1j * np.dot(K[q],R[i]))

Gq = np.zeros((N, 3, 3)) #Lattice Green's function in reciprocal space
Gtmp = np.zeros((3,3))
D=np.zeros((3,3)) #Dynamical matrix of the host lattice

for q in range(N):
	D[0,0] = alat*c11*(2-np.cos(alat*K[q][0]/2.0)*(np.cos(alat*K[q][1]/2.0)+np.cos(alat*K[q][2]/2.0))) + alat*(2*c44-c11)*(1-np.cos(alat*K[q][1]/2.0)*np.cos(alat*K[q][2]/2.0))
	D[1,1] = alat*c11*(2-np.cos(alat*K[q][1]/2.0)*(np.cos(alat*K[q][0]/2.0)+np.cos(alat*K[q][2]/2.0))) + alat*(2*c44-c11)*(1-np.cos(alat*K[q][0]/2.0)*np.cos(alat*K[q][2]/2.0))
	D[2,2] = alat*c11*(2-np.cos(alat*K[q][2]/2.0)*(np.cos(alat*K[q][1]/2.0)+np.cos(alat*K[q][0]/2.0))) + alat*(2*c44-c11)*(1-np.cos(alat*K[q][1]/2.0)*np.cos(alat*K[q][0]/2.0))
	D[0,1] = alat*(c12+c44)*np.sin(alat*K[q][0]/2.0)*np.sin(alat*K[q][1]/2.0)
	D[0,2] = alat*(c12+c44)*np.sin(alat*K[q][0]/2.0)*np.sin(alat*K[q][2]/2.0)
	D[1,0] = alat*(c12+c44)*np.sin(alat*K[q][1]/2.0)*np.sin(alat*K[q][0]/2.0)
	D[1,2] = alat*(c12+c44)*np.sin(alat*K[q][1]/2.0)*np.sin(alat*K[q][2]/2.0)
	D[2,0] = alat*(c12+c44)*np.sin(alat*K[q][2]/2.0)*np.sin(alat*K[q][0]/2.0)
	D[2,1] = alat*(c12+c44)*np.sin(alat*K[q][2]/2.0)*np.sin(alat*K[q][1]/2.0)
	if np.linalg.norm(D) > 0:
		Gq[q] = 160.21766 * np.linalg.inv(D) 
	else:
		Gq[q] = np.zeros((3,3))

Q11 = 0  #self interactions
Q22 = 0
Q33 = 0
Q12 = 0
Q13 = 0
Q23 = 0


#Strain induced interaction in reciprocal space
Vsq11 = np.zeros((N,1))
Vsq22 = np.zeros((N,1))  
Vsq33 = np.zeros((N,1))
Vsq12 = np.zeros((N,1))
Vsq13 = np.zeros((N,1))
Vsq23 = np.zeros((N,1))

for q in range(N):
	Vsq11[q] = - np.real(np.dot(np.dot(F1q[q], Gq[q]), np.conj(F1q[q])))
	Vsq22[q] = - np.real(np.dot(np.dot(F2q[q], Gq[q]), np.conj(F2q[q])))
	Vsq33[q] = - np.real(np.dot(np.dot(F3q[q], Gq[q]), np.conj(F3q[q])))
	Vsq12[q] = - np.real(np.dot(np.dot(F1q[q], Gq[q]), np.conj(F2q[q])))
	Vsq13[q] = - np.real(np.dot(np.dot(F1q[q], Gq[q]), np.conj(F3q[q])))
	Vsq23[q] = - np.real(np.dot(np.dot(F2q[q], Gq[q]), np.conj(F3q[q])))
	Q11 = Q11 + np.real(1/float(N) * np.dot(np.dot(F1q[q], Gq[q]), np.conj(F1q[q])))
	Q22 = Q22 + np.real(1/float(N) * np.dot(np.dot(F2q[q], Gq[q]), np.conj(F2q[q])))
	Q33 = Q33 + np.real(1/float(N) * np.dot(np.dot(F3q[q], Gq[q]), np.conj(F3q[q])))
	Q12 = Q12 + np.real(1/float(N) * np.dot(np.dot(F1q[q], Gq[q]), np.conj(F2q[q])))
	Q13 = Q13 + np.real(1/float(N) * np.dot(np.dot(F1q[q], Gq[q]), np.conj(F3q[q])))
	Q23 = Q23 + np.real(1/float(N) * np.dot(np.dot(F2q[q], Gq[q]), np.conj(F3q[q])))

for q in range(N):
	Vsq11[q] = Vsq11[q] + Q11
	Vsq22[q] = Vsq22[q] + Q22
	Vsq33[q] = Vsq33[q] + Q33
	Vsq12[q] = Vsq12[q] + Q12
	Vsq13[q] = Vsq13[q] + Q13
	Vsq23[q] = Vsq23[q] + Q23


#Strain induced interactions in reciprocal space
f = open('Vsq.dat', 'w')
f.write('#kx ky kz Ni-Ni Co-Co Cr-Cr Ni-Co Ni-Cr Co-Cr (meV)\n')
for q in range(N):
	f.write('%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (K0[q,0],K0[q,1],K0[q,2],Vsq11[q]*1e3,Vsq22[q]*1e3,Vsq33[q]*1e3,Vsq12[q]*1e3,Vsq13[q]*1e3,Vsq23[q]*1e3))

f.close()

nr = 4

R0 = np.zeros((nr,3)) #coordinates of neighboring shells

R0[0] = np.array([-0.5,0.5,0.0])
R0[1] = np.array([-1.0,0.0,0.0])
R0[2] = np.array([-1.0,0.5,0.5])
R0[3] = np.array([-1.0,0.0,1.0])

R0 = alat * R0

#Strain induced interaction in real space
Vsr11 = np.zeros((nr,1))
Vsr22 = np.zeros((nr,1))
Vsr33 = np.zeros((nr,1))
Vsr12 = np.zeros((nr,1))
Vsr13 = np.zeros((nr,1))
Vsr23 = np.zeros((nr,1))

for s in range(nr):
	for q in range(N):
		Vsr11[s] = Vsr11[s] + 1/float(N) * Vsq11[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))
		Vsr22[s] = Vsr22[s] + 1/float(N) * Vsq22[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))
		Vsr33[s] = Vsr33[s] + 1/float(N) * Vsq33[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))
		Vsr12[s] = Vsr12[s] + 1/float(N) * Vsq12[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))
		Vsr13[s] = Vsr13[s] + 1/float(N) * Vsq13[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))
		Vsr23[s] = Vsr23[s] + 1/float(N) * Vsq23[q] * np.real(np.exp(1j * np.dot(K[q],R0[s])))


f = open('Vsr.dat','w')
f.write("# no. Rx Ry Rz Ni-Ni Co-Co Cr-Cr Ni-Co Ni-Cr Co-Cr (meV)\n")
for i in range(nr):
	f.write("%s %s %s %s %s %s %s %s %s %s\n" % (i,R0[i,0]/alat,R0[i,1]/alat,R0[i,2]/alat,Vsr11[i,0]*1e3,Vsr22[i,0]*1e3,Vsr33[i,0]*1e3,Vsr12[i,0]*1e3,Vsr13[i,0]*1e3,Vsr23[i,0]*1e3))
f.close()
