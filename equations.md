
### Jacobian calculation
```latex
dir1 = (x)*pos1 +c 
dir2 = (x)*pos2 +c
dir3 = (x)*pos3 +c
v_x0 = ax0+by0+cz0+d
v_y0 = ex0+fy0+gz0+h
v_z0 = ix0+jy0+kz0+l
v_x1 = ax1+by1+cz1+d
v_y1 = ex1+fy1+gz1+h
v_z1 = ix1+jy1+kz1+l
v_x2 = ax2+by2+cz2+d
v_y2 = ex2+fy2+gz2+h
v_z2 = ix2+jy2+kz2+l

v_x0 = ax0+by0+cz0+d
v_x1 = ax1+by1+cz1+d
v_x2 = ax2+by2+d

[row of A(x)][col of P(0)]=
A=f1 .. fn

J=dV/dP=A.T
f(x+h)=f(x)+f'(x)*h
V_i=V_0+J*(P_i-P_0)
V_i=J*P_i+(V_0-J*P_0)
shift P_0=>0
V_i=J*P_i+V_0
V=[A t]|P|
       |1|
V=[V1 V2...Vn]

V_i-V=J*(P_i-P)
is the same as
dy=dy/dp*dp

Idea:choose P in center of P_i's, shift all P_i to P_i-P
P=(px,py,pz)
V_0x=a*P_0x+b*P_0y+c*P_0z+E=J*P_0+(V-J*P) => E=V-J*P

P= px0 px1 px2
   py0  #   # 
   pz0  #   # 
V= vx0 vx1 vx2
   vy0  #   # 
   vz0  #   # 
A= [F1 .. ]

P.T*A.T= V.T
AP=V
```

### Vector Approximation
```latex
V(x)=sum(V(x+h_i)-J(x)*h_i)/n                                     //jacobian info
V(x)=sum(w_iV(x+h_i))                                             //interpolation (jacobian is not known at x)
V(x)=sum(w_i(V(x+h_i)-J(x)*h_i))                                  //merged (jacobian is known at x)
V(x)*m=sum(w_iV(x+h_i)-J(x)*h_i) + sum(V(x+h_i)+J(x+h_i)*h_i)     //jacobian is known at x+h_i and at x
V(x)*m=sum(w_iV(x+h_i)) + sum(V(x+h_i)+J(x+h_i)*h_i)              //jacobian is known at x+h_i

// jacobian
v_i=sum(v_j)/n-sum(J_i*(p_j-p_i))/n
J_i*sum(p_j-p_i)=sum(v_j)-v_i*n

// merged
J_i*(sum(w_ij*p_j)-p_i) + sum(w_iu*J_u*(p_u-p_i)) = sum(w_ij*v_j) + sum(w_iu*v_u) - v_i*m       (jacobian at x and x+h_i)
J_i*(sum(w_ij*p_j)-p_i) = sum(w_ij*v_j) - v_i                                         (jacobian at x)
sum(w_iu*J_u*(p_u-p_i)) = sum(w_iu*v_u) - v_i*m                                       (jacobian at x+h_i)
0 = sum(w_ij*v_j) - v_i                                                               (no jacobian info nearby)   
```

different ideas:
   Nguyen's: 
      use jacobian information to add more points to the field 
      then interpolate with upgraded field
   Add Equations(current and not very successful):
      sum the interpolation and jacobian equations 
   Merge Equations(current and not very successful):
      use interpolation weight in jacobian equations() 
   Augment Matrix:
      add rows to the matrix which represent the jacobian equations
      problems: equations can contradict each other
   

problem: 2 adjacent vertices
case 1 has dir: interpolate
case 1 has jacobian: jacobian linear equation
case 2 have jacobian: ^
case 1 has dir and jacobian: 