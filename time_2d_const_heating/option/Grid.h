integer,parameter::nt=64
real(4),parameter::tlen=1000.0
real(4),parameter::dt=tlen/real(nt)

real(4),parameter::axt(nt)=(/(real(it-nt/2)*dt,it=1,nt)/)
