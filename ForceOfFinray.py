import numpy as np 
import math
import time
import matplotlib.pyplot as plt
import sympy as sym
from scipy.integrate import quad, dblquad
from scipy.optimize import fsolve


def cal_F(u, v, k, f,Am):
	""" 
	Calculate base element of hopf oscillator
	Input:
	Output:
	"""
	return k*(Am*Am-u**2-v**2)*u - 2*np.pi*f*v,  k*(Am*Am-u**2-v**2)*v + 2*np.pi*f*u

def cal_P_head(post_u, post_v, epsilon, psi):
	"""
	"""
	return epsilon * (post_v*np.cos(psi) - post_u*np.sin(psi))

def cal_P_tail(pre_u, pre_v, epsilon, psi):
	"""
	"""
	return epsilon * (pre_u*np.sin(psi) + pre_v*np.cos(psi))

def cal_P_body(pre_u, pre_v, post_u, post_v, epsilon, psi):
	"""
	"""
	return epsilon * (pre_u*np.sin(psi) + pre_v*np.cos(psi) - post_u*np.sin(psi) + post_v*np.cos(psi))
def quadratic_waveform(cnt):
		fx=1*cnt*1/15+0
		return fx
def sin_waveform(fx):
		return fx
def eliptic_waveform(cnt):
		if cnt <=7:
			fx=25*cnt*1/15+5
		elif cnt>7:
			fx=25*(15-cnt)*1/15+5
		return fx
endtime = 10
step = 0.01

def Cal_velocity(force_thuster,v_r):
	m=1
	po=1000 #kg/m^3
	Cv=0.1
	R=75 #mm ban kinh của than robot
	Area=2*np.pi*np.power( R/1000 , 2 )
	v_water=0;  #nuoc tinh
	F_drag_up=1/2*po*Cv*Area*np.power(v_r-v_water,2)
	eqs=force_thuster - F_drag_up - m*v_r/step;
	return eqs

def fn_cal(w,h,theta,post_theta,A,post_A,u,post_u,v,post_v):
	
	D=20/1000
	L=D*15
	p = 1000
	Cn = 2.8
	hmax = 0.1
	f=1
	Am=1
	k=10
	v0=5


	u_dot = k*(Am*Am-u**2-v**2)*u - 2*np.pi*f*v
	post_u_dot = k*(Am*Am-post_u**2-post_v**2)*post_u - 2*np.pi*f*post_v

	u_theta_dot = A*u_dot
	u_post_theta_dot =post_A*post_u_dot

	delta= np.sqrt(  2*D*D*w*(w-1)*( 1-np.cos(post_theta - theta) )+ h*h*np.sin(post_theta - theta)*np.sin(post_theta - theta)+D*D  )
	n0= (1*D)*h/D*np.sin(post_theta-theta)
	vq=  (1/delta)*D*h*(-h/D*v0*np.sin(post_theta-theta)+ (np.cos(theta)*(1-w)+w*np.cos(post_theta))*(np.cos(theta)*u_theta_dot+( np.cos(post_theta)*u_post_theta_dot-np.cos(theta)*u_theta_dot )*w )+(np.sin(theta)*(w-1) -w*np.sin(post_theta) )*(-np.sin(theta)*u_theta_dot+(np.sin(theta)*u_theta_dot -np.sin(post_theta)*u_post_theta_dot)*w ) )
	fn=(-1/2)*p*Cn*vq*vq*n0*np.sign(vq)
	return fn
def intergrate_fn(theta,post_theta,A,post_A,u,post_u,v,post_v):
	hmax=0.1
	return dblquad(lambda h,w:fn_cal(w,h,theta,post_theta,A,post_A,u,post_u,v,post_v),0, 1, 0, hmax)[0]



def CPG_Calulate():
	array_u = np.zeros([1,16])
	array_v = np.zeros([1,16])
	array_theta=np.zeros([1,16])	
	array_time= [0] 
	array_psi_r=[0]
	# initialize parameter
	time=0
	Ax=1
	k = 10
	f = 1
	epsilon = 0.8
	psi = -np.pi/3

	# initial state
	array_u[0][0] = 0
	array_v[0][0] = 0.001

	for idx in range(1,int(endtime/step)):

		time=time+0.01
		state_u = []
		state_v = []
		state_theta=[]
		
		array_theta = np.append(array_theta, np.zeros([1,16]), axis=0)
		array_u = np.append(array_u, np.zeros([1,16]), axis=0)
		array_v = np.append(array_v, np.zeros([1,16]), axis=0)

		for i in range(16):
			#Waveform
			# A=quadratic_waveform(i)
			A=sin_waveform(1)
			# compute the base element
			F_u, F_v = cal_F(array_u[idx-1][i], array_v[idx-1][i], k, f,Ax)
			# compute new state of ith CPG at time idx*step with newton approximation
			new_u = F_u * step + array_u[idx-1][i]
			if i==0:
				new_v = (F_v + cal_P_head(array_u[idx-1,1], array_v[idx-1,1], epsilon, psi)) * step + array_v[idx-1][i]
			elif i==15:
				new_v = (F_v + cal_P_tail(array_u[idx-1,14], array_v[idx-1,14], epsilon, psi)) * step + array_v[idx-1][i]
			else: 
				new_v = (F_v + cal_P_body(array_u[idx-1,i-1], array_v[idx-1,i-1], array_u[idx-1,i+1], array_v[idx-1,i+1], epsilon, psi)) * step + array_v[idx-1][i]
			new_theta=A*new_u
			# create new state vector
			state_u.append(new_u)
			state_v.append(new_v)
			state_theta.append(new_theta)
		# add new state vector to original placeholder
		array_u[idx] = array_u[idx] + state_u
		array_v[idx] = array_v[idx] + state_v
		array_theta[idx]=array_theta[idx]+state_theta

	return array_theta,array_u,array_v

if __name__ == '__main__':
	
	print(" Start CPG !!!")
	theta,u,v = CPG_Calulate()
	print(" CPG calculated !!!")

	ti=0
	arr_time=[]
	arr_force=[]
	arr_vel=[]
	# int(endtime/step)
	for idx in range(1,int(endtime/step)):
			ti=ti+step			
			sum_F=0

			for i in range(15):
				Amp =sin_waveform(1)
				Amp_post= sin_waveform(1)
				# Amp =quadratic_waveform(i)
				# Amp_post= quadratic_waveform(i+1)
				fx=intergrate_fn(theta[idx,i], theta[idx,i+1], Amp , Amp_post, u[idx,i] , u[idx,i+1] ,v[idx,i] , v[idx,i+1],)
				
				sum_F=sum_F+fx
			Vel = fsolve(lambda v_r:Cal_velocity(sum_F,v_r),1.0)
			print(sum_F)
			arr_vel.append(Vel)
			arr_time.append(ti)
			arr_force.append(sum_F)
	print(" Force calculated !!!")





	fig, ax = plt.subplots()
	legend_list = []
	for i in range(16):		
		ax.plot(np.linspace(0, endtime, int(endtime/step)),theta[:,i])
		legend_list.append(i+1)
	ax.legend(legend_list)
	ax.set_title("CPG ")
	ax.set_ylabel("Output of oscillator")
	ax.set_xlabel("time(s)")
	y_axis = ax.axes.get_yaxis()
	# y_axis.set_visible(False)
	# y_axis.set_ticklabels([])
	# ax.set_xlim(0,10)
	# ax.set_ylim(-2,34)
	ax.grid('on')

	fig1, ax1 = plt.subplots()
	ax1.plot(arr_time,arr_force);
	ax1.set_title("Force ")
	ax1.set_ylabel("F(N)")
	ax1.set_xlabel("time(s)")

	fig2, ax2 = plt.subplots()
	ax2.plot(arr_time,arr_vel);
	ax2.set_title(" Velocity ")
	ax2.set_ylabel("V(m/s)")
	ax2.set_xlabel("time(s)")


	plt.show()