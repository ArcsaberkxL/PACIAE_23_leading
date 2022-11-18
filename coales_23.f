	subroutine coales(ijk,neve,nnout,nap,nat,nzp,nzt)    
c       A simple coalescence model writen by Sa Ben-Hao on 04/06/2004
c       Its input messages are in 'pyjets'   ! 220822
c       Its storing array is 'pyjets' 
c	Its output message is in 'sa1_h' (in 'pyjets' either)
c	ijk: the run number
c	neve: total number of runs
c       nnout: a internal printing per nnout runs
c	nap and nzp: atomic and charge number of projectile
c       nat and nzt: atomic and charge number of target
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (KSZJ=80000)
      parameter (mplis=80000)
	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c	Those variables in above four statements are only used here and   
c	in subroutine 'decayh','findb' and 'thephi'. 
c	PYDAT1,PYDAT2,PYDAT3 and PYJETS are the subroutines in PYTHIA
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
        common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 280822
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) ! 220822
        common/sbh/nbh,nonh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  napp,natt,nzpp,nztt,pio
	dimension rc(3),b(3),pstr(3),rstr(3),numb(3),jk(20)
	dimension p0(4),pf(20,4),ppt(50,2),peo(5)
	dimension ppr(kszj,4),peoh(5),pnn(kszj,5),rr(3)
        dimension spglu(4) ! 080512, e (p) sum of unsplitted gluons

c271022 Moved from 'coal' to here       
	ithroq=0
	ithrob=0
	ich=0
	do i=1,4
	throe(i)=0.
	enddo
c271022

	rrp=1.16
	nn=0   ! nn: the particle number of this event.
        do i1=1,kszj
          do j1=1,5
          kn(i1,j1)=0.
          pn(i1,j1)=0.
          rn(i1,j1)=0.
	  enddo
	enddo

	nout=nnout
	imc=adj1(13)
        ibc=adj1(14)
	iphas=adj1(21)

c220822 remove junctions
        jb=0
2010    do i1=jb+1,n  ! i1 loop      
        kf=k(i1,2)
        kfab=iabs(kf)
        if(kfab.ne.88)then
        jb=jb+1
        goto 2020
        endif
c       move particle list 'pyjets' one step downward since i1+1 to n
c       do j=i1+1,n
c       do jj=1,5
c       k(j-1,jj)=k(j,jj)
c       p(j-1,jj)=p(j,jj)
c       v(j-1,jj)=v(j,jj)
c       enddo
c       enddo
        call updad_pyj(n,i1+1,1)   ! 090922 'updad_pyj' in sfm_23.f
        n=n-1
        goto 2010
2020    enddo   ! i1 loop        
c220822        

c	Conservation of net baryon.
	netba=0
	do i1=1,nbh ! note: 'sbh' is hadron list before hadronization 060119
	  kf=kbh(i1,2)
	  kfab=iabs(kf)
	  if(kf.gt.0.and.(kf.gt.1000 .and. kf.lt.10000))netba=netba+1
	  if(kf.lt.0.and.(kfab.gt.1000 .and. kfab.lt.10000))netba=netba-1
	enddo
c060813 120214
	if(nap.gt.1 .and. nat.gt.1)then
	netba=nap+nat-netba   ! AA
	elseif(nap.eq.1.and.nat.gt.1.and.ipden.eq.0)then
	netba=nzp+nat-netba   ! p-A (pbar-A)
	elseif(nap.gt.1.and.nat.eq.1)then 
	netba=nap+nzt-netba   ! A-p (A-pbar)
	elseif(nap.eq.1.and.nat.eq.1)then
	netba=nzp+nzt-netba   ! pp (p-pbar,pbar-p)	
        else
 	netba=nat-netba   ! lepton-A
	endif
c060813 120214
c060119 if(ijk.eq.1)then
c       napt=nap+nat
c       nzpt=nzp+nzt
c       sbaryi=float(napt)
c       schgei=3.*float(nzpt)
c060119 endif

888   continue

c220122
c       write(22,*)'in coales n=',n
        call pyedit(2) 
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       move gluons from 'pyjest' to 'sa36'
        call remo_glu
c       write(22,*)'af. remo_glu n,nglu=',n,nglu
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       call prt_sa36(nglu,cc)
c       break-up gluon (with E_g>2E_u in 'sa36') -> qqbar string 
c        (filling in 'pyjets')
        call break_glu
c       write(22,*)'af. break_glu n=',n
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       so far, the parton list ('pyjets') is composed of q and qba only 

	adj12=adj1(12)
	adj16=adj1(16)
	adj17=adj1(17)
c200222 adj17=max(4.0,adj17) ! 070612, yan

c280822 energetic q (qbar) de-excitation according to Field-Feynman model
        n0=n
        igens=0
c       goto 900   ! 280822
700	do 800 i1=1,n0   ! 080512 280822   
        kf0=k(i1,2)
	ee=p(i1,4)
        iflav=1
        if(kf0.lt.0)iflav=-1
c       iflav = 1 : if source parton is quark
c             =-1 : if source parton is antiquark

        if(ee.gt.adj17)then
        call ffm(i1,kf0,igen,iflav)   ! 280822
        igens=igens+1   ! # of 'call ffm'
c       print*,'n0,adj17,i1,igens=',n0,adj17,i1,igens   ! 280822
        endif
c       igen : # of generations per source q (qba)
800     enddo   ! 280822 continue->enddo
c280822 n=n3
900	continue
c       energetic q (qbar) de-excitation, finished.
c       write(22,*)'af. ffm igens,n=',igens,n
c       call pylist(1)
c       call prt_sbh(nbh,cc)
c       print*,'af. ffm 2'

c220122
        numb(1)=0   ! 220822
        numb(2)=0
        numb(3)=0
c       numb(1),(2), and (3): the order # of last g,qba & q in 'pyjets'
c220822	1: refers to g (no gluon at all now), 2: qba, 3: q

c       make the partons in order of qba and q   ! 220822
c        i.e. move q to the end   ! 220822 
        jh=n
        jl=0
2030    continue
        do j=jl+1,jh
        kf=k(j,2)
        kfab=iabs(kf)
        if(kfab.lt.7 .and. kf.gt.0)then   ! q, consider d,u,s,c,b,t only
        n=n+1      
        numb(3)=numb(3)+1  
        do i4=1,4
        k(n,i4)=k(j,i4)
        p(n,i4)=p(j,i4)
        v(n,i4)=v(j,i4)
        enddo
c       move particle list 'pyjets' one step downward since j+1 to n
c       do i2=j+1,n
c       do jj=1,5
c       k(i2-1,jj)=k(i2,jj)
c       p(i2-1,jj)=p(i2,jj)
c       v(i2-1,jj)=v(i2,jj)
c       enddo
c       enddo
        call updad_pyj(n,j+1,1)
        n=n-1
        jh=jh-1
        jl=j-1
        goto 2030
        endif
        enddo
        numb(1)=0
        numb(2)=n-numb(3)
        numb(3)=n
        n1=numb(1)
        n2=numb(2)
        n3=n
c220822
c       write(22,*)'af. ordering numb1,2,3=',numb(1),numb(2),numb(3)       
c       write(22,*)'n1,n2,n3,n=',n1,n2,n3,n      
c       call pylist(1)
c       call prt_sbh(nbh,cc)        
c220822        

c       Order the qba according to energy from the maximal to minimal.
        call eord(n1+1,n2)  
c       Order the q according to energy from the maximal to minimal.
        call eord(n2+1,n3)
c       write(22,*)'af. eord n1,n2,n3,n=',n1,n2,n3,n
c       call pylist(1)
c       call prt_sbh(nbh,cc)

c220822
c       share the energy in 'ppsa' among partons
        sn=dfloat(n)
        spx=ppsa(1)/sn
        spy=ppsa(2)/sn
        spz=ppsa(3)/sn
        see=ppsa(4)/sn
        do i1=1,n
        p(i1,1)=p(i1,1)+spx
        p(i1,2)=p(i1,2)+spy
        p(i1,3)=p(i1,3)+spz
        p(i1,4)=p(i1,4)+see
        enddo
        do i1=1,4
        ppsa(i1)=0.
        enddo

c	Parton coalescence 
	if(ijk.eq.1)call tabhb 
c       ijk is the event number
c       Read the table of hadron (meson: pseudoscalar-spin 0 & vector-spin 1 
c        only, baryon: octet-spin 1/2 & decuplet-spin 3/2 only)

        iqba=n2
	call coal(n3,iqba,ijk,rrp,iphas,netba)
c	n3: total number of partons (qba and q) 
c	iqba: total # of qba (qba is ordered before q) 
c	ithroq : the total # of quarks thrown away
c	ithrob : the total # of antiquarks thrown away
c	throe : total 4-momentum of the partons thrown away
c	ich : total charge of the partons thrown away

	iqba=ithrob   ! not active originally
	n3=ithroq+ithrob   ! not active originally
	if(iphas.eq.1 .and. n3.ge.2)then
	  call coal(n3,iqba,ijk,rrp,0,0)   ! 110905
	endif
 
c150922 ichth=ich   ! 092600

c	'sbh' to 'sa1_h'.
c010518	if(nbh.ge.1)then
c	do l=1,nbh
c	l1=nn+l
c	do m=1,5
c		kn(l1,m)=kbh(l,m)
c		pn(l1,m)=pbh(l,m)
c		rn(l1,m)=vbh(l,m)
c	enddo
c	enddo
c	nn=nn+nbh
c010518	endif

c	'sa1_h' to 'pyjets'.
	n=nn
	do j1=1,nn
        do j2=1,5
          k(j1,j2)=kn(j1,j2)
          p(j1,j2)=pn(j1,j2)
          v(j1,j2)=rn(j1,j2)
        enddo
      enddo

c     Decay of unstable hadrons
        if(adj12.ne.0)call decayh(rrp)    ! 060119 

	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pes(pei,il,ih,peo)
c	sum up momentum and energy.  
c	pei : two dimension array of input momentum and energy
c	il and ih : lower and higher limits of summation
c	peo : one dimension array of output momentum,energy & sqrt(s)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (KSZJ=80000)
      parameter (mplis=80000)
      dimension pei(kszj,4),peo(5)
	do i=1,5
	  peo(i)=0.
	enddo
	do i=il,ih
	  do j=1,4
	    peo(j)=peo(j)+pei(i,j)
	  enddo
	enddo
	peo(5)=peo(4)*peo(4)
	do i=1,3
	  peo(5)=peo(5)-peo(i)*peo(i)
	enddo
	peo5=peo(5)
	if(peo5.lt.0.)then
	  peo5=abs(peo5)
	endif

100	format('            px           py          pz         e      
     c	 sqrt(s)')
200	format(4x,5(1x,f9.3))
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine tabhb
c	Table of primary meson (pseudoscalar (spin 0) and vector (spin 1) 
c        only) and baryon (octet (spin 1/2) and decuplet (spin 3/2) only)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
c       imc (ibc): dimension of meson (baryon) table considered
c       kqh(i,j): flavor code of j-th constituent q (or qba) in i-th quark 
c        configuration of meson
c       kfh(i,1): flavor code of pseudoscalar meson with i-th quark
c        configuration of meson
c       kfh(i,2): flavor code of vector meson with i-th quark
c        configuration of meson
c       proh(i,1): probability of pseudoscalar meson with i-th quark
c        configuration of meson
c       proh(i,2): probability of vector meson with i-th quark
c        configuration of meson
c       amash(i,1): mass of pseudoscalar meson with i-th quark
c        configuration of meson
c       amash(i,2): mass of vector meson with i-th quark
c        configuration of meson
c       kqb(i,j): flavor code of j-th constituent q (or qba) in i-th quark
c        configuration of baryon
c       kfb(i,1): flavor code of octet baryon with i-th quark
c        configuration of baryon
c       kfb(i,2): flavor code of decuplet baryon with i-th quark
c        configuration of baryon
c       prob(i,1): probability of octet baryon with i-th quark
c        configuration of baryon
c       prob(i,2): probability of decuplet baryon with i-th quark
c        configuration of baryon
c       amasb(i,1): mass of octet baryon with i-th quark
c        configuration of baryon
c       amasb(i,2): mass of decuplet baryon with i-th quark
c        configuration of baryon
	data (kqh(i,1),i=1,80)/2*2,2*1,2*3,3*2,3*1,2*3,4,1,4,2,4,3,4,
     c	 2,5,1,5,5,54*0/	
	data (kqh(i,2),i=1,80)/-1,-3,-2,-3,-2,-1,-2,-2,-2,-1,-1,-1,-3,
     c	     -3,-1,-4,-2,-4,-3,-4,-4,-5,-2,-5,-1,-5,54*0/
	data (kfh(i,1),i=1,80)/211,321,-211,311,-321,-311,111,221,331,
     c	     111,221,331,221,331,411,-411,421,-421,431,-431,441,521,
     c	     -521,511,-511,553,54*0/
	data (kfh(i,2),i=1,80)/213,323,-213,313,-323,-313,113,223,0,
     c	     113,223,0,333,0,413,-413,423,-423,433,-433,443,2*0,513,
     c	     -513,0,54*0/
	data (proh(i,1),i=1,80)/6*1.,0.5,2*0.25,0.5,2*0.25,2*0.5,12*1,
     c	     54*0./
	data (proh(i,2),i=1,80)/6*1.,2*0.5,0.,2*0.5,0.,1.,0.,7*1,
     c	     2*0,2*1,0,54*0./

	data (kqb(i,1),i=1,80)/5*2,3*1,3,4*2,1,2,1,3,2,62*0/
	data (kqb(i,2),i=1,80)/3*2,1,3,2*1,2*3,3*1,2,1,3*3,1,62*0/
	data (kqb(i,3),i=1,80)/2,1,3,1,3,1,5*3,6*4,5,62*0/
	data (kfb(i,1),i=1,80)/0,2212,3222,2112,3322,0,3112,3312,0,
     c	     3122,3212,4122,4222,4112,4232,4132,4332,5122,62*0/
	data (kfb(i,2),i=1,80)/2224,2214,3224,2114,3324,1114,3114,3314,
     c	     3334,0,3214,4212,4222,0,4232,4132,2*0,62*0/
	data (prob(i,1),i=1,80)/0.,4*1,0.,2*1.,0.,2*0.5,7*1.,62*0./
	data (prob(i,2),i=1,80)/9*1.,0.,1.,2*1.,0.,2*1.,2*0.,62*0./

	do i1=1,imc
	  kf1=kfh(i1,1)
	  kf2=kfh(i1,2)
	  amash(i1,1)=pymass(kf1)
	  amash(i1,2)=pymass(kf2)
	enddo
	do i1=1,ibc
	  kf1=kfb(i1,1)
	  kf2=kfb(i1,2)
	  amasb(i1,1)=pymass(kf1)
	  amasb(i1,2)=pymass(kf2)
	enddo

	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine eord(ni,nc)
c   Order particle set (ni to nc) according to energy
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (kszj=80000,mplis=80000)   ! 280822
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 280822
	dimension rr(4),p1(4)
	do 100 i1=ni,nc
          ii=i1
	  if(i1.eq.nc)goto 100
c280822   j=ii
	  alar=p(ii,4)
	  do 200 i2=i1+1,nc
c280822 communication between i1 and i2 of which the energy is lagest
c        among i1+1, i1+2, ..., nc 
	    ee=p(i2,4)
	    if(ee.gt.alar)then   ! 280822 .ge. -> .gt.
		j=i2 
		alar=ee
	    endif
200	  enddo   ! continue-> enddo 280822
c280822 now, j: order # of particle with largest energy
	  kii2=k(ii,2)   
	  do jj=1,4
	  p1(jj)=p(ii,jj)
	  rr(jj)=v(ii,jj)
	  enddo
	  pii5=p(ii,5)

	  k(ii,2)=k(j,2) 
	  do jj=1,4
	  p(ii,jj)=p(j,jj)
	  v(ii,jj)=v(j,jj)
	  enddo
	  p(ii,5)=p(j,5)

	  k(j,2)=kii2
	  do jj=1,4
	    p(j,jj)=p1(jj)
	    v(j,jj)=rr(jj)
	  enddo
	  p(j,5)=pii5
100	enddo

	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coal(n1,nqb,ijk,rrp,iphas,netba)
c	Parton coalescence (hadronization)
c	n1 : total # of partons (q & qba only)
c	nqb : total # of qba (qba is ordered before q)
c	ijk : the run number
c	iphas=1: with phase space adjudgment
c	     =0: without phase space adjudgment
c	netba: number of baryons to be first generated keeping
c	       baryon conservation 
c	ithroq : total # of quarks thrown away
c	ithrob : total # of antiquarks thrown away
c	throe : total four momentum of the partons thrown away   ! 090922
c	ich : total charge of the partons thrown away
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000,mplis=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa6_c/ithroq,ithrob,ich,non6_c,throe(4)
      common/sa24/adj1(40),nnstop,non24,zstop
	common/sa37/nth,npadth,kth(kszj,5),pth(kszj,5),vth(kszj,5)   ! 150922
	dimension pc(4),rc(3),iar(3),rcp(3)
	dimension psu(3),peo(5),pnn(kszj,5)
	dimension numb(3)   ! 110905

	nth=0
      do j=1,kszj
      do i=1,5
      kth(j,i)=0
      pth(j,i)=0.
      vth(j,i)=0.
      enddo
      enddo

c	Coalescence
	ibarp=0
	ibarm=0
	imes=0

c	Generate first 'netba' baryons (if 'netba'>0) or -'netba' antibaryons 
c	 (if 'netba'<0) keeping baryon conservation 
	if(netba.gt.0)call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,1)
	if(netba.lt.0)then
302     continue
	  do ii1=1,nqb
	    call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,1)
	    if(isu.eq.1 .and. ibarm.lt.-netba)goto 303
          if(isu.eq.1 .and. ibarm.eq.-netba)goto 304
          if(isu.eq.0) stop 6666
	  enddo
303     goto 302
304     continue
	endif
	
c     Compose first the antiquark, one by one, into hadron 150922
c     A antiquark can be a component of antibaryon or meson.
	jj=0
300	continue
	  jj=jj+1
	  if(nqb.eq.0)goto 100   ! composing baryon
	  if(nqb.lt.3 .and. n1.eq.nqb)goto 301 
c	No q, cannot compose meson, throw away nqb.

	  if(nqb.eq.1)ii1=1
	  if(nqb.gt.1)ii1=nqb*pyr(1)+1

        kf=k(ii1,2)
        relpr=adj1(31)
        relpr=relpr/(1.+relpr)
c       probability of antibaryon
        if(kf.eq.-3)then
          relpr=adj1(31)*adj1(33)
          relpr=relpr/(1.+relpr)
c         probability of strange antibaryon
        endif

	  rand=pyr(1)
        if(rand.le.relpr .and. nqb.ge.3)then   !
          call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)

          if(isu.eq.1 .and. nba.eq.1)then   ! one antibaryon produced
	      goto 300
	    endif   

          if(isu.eq.0)then   ! ii1 can not produce antibryon
            call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
		  if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
			goto 300   
		  endif

            if(isu.eq.0)then   !!! ii1 can not produce meson either, throw away
	      nth=nth+1  
              ithrob=ithrob+1
              kk=k(ii1,2)
              ich=ich+pychge(kk)
              kth(nth,2)=kk
              pth(nth,5)=p(ii1,5)
              do i2=1,4
		    throe(i2)=throe(i2)+p(ii1,i2)
	    	    pth(nth,i2)=p(ii1,i2)
		    vth(nth,i2)=v(ii1,i2)
              enddo

c     Move parton list one steps downward since ii1+1
              call updad_pyj(n1,ii1+1,1)   
              n1=n1-1
              nqb=nqb-1
              goto 300
            endif   
          endif   ! ii1 can not produce antibryon
        elseif(rand.gt.relpr .or. nqb.lt.3)then   !
          call mespro(n1,nqb,ii1,iphas,nme,imes,rrp,isu)
          if(isu.eq.1 .and. nme.eq.1)then   ! one hadron produced
		goto 300   
	    endif
          if(isu.eq.0)then   !! ii1 can not produce meson
		call an_barpro(n1,nqb,ii1,iphas,nba,ibarm,netba,rrp,isu,2)
		if(isu.eq.1 .and. nba.eq.1)then   ! one antibryon produced
		  goto 300   
		endif
            if(isu.eq.0)then !!! can not produce antibryon either, throw away
		  nth=nth+1   
		  ithrob=ithrob+1
		  kk=k(ii1,2)
		  ich=ich+pychge(kk)
        	  kth(nth,2)=kk
		  pth(nth,5)=p(ii1,5)
		  do i2=1,4
		    throe(i2)=throe(i2)+p(ii1,i2)
		    pth(nth,i2)=p(ii1,i2)
		    vth(nth,i2)=v(ii1,i2)
		  enddo

c     Move parton list one steps downward since ii1+1
		  call updad_pyj(n1,ii1+1,1)   
		  n1=n1-1
		  nqb=nqb-1
		  goto 300
            endif   !!!
          endif   !!
	  else   !
        endif  !

c	Throw away those qba remained.
301	continue
	if(nqb.eq.0)goto 100
      do i1=1,nqb
	nth=nth+1   ! 110905
        ithrob=ithrob+1
        kk=k(i1,2)
        ich=ich+pychge(kk)
        kth(nth,2)=kk
        v(nth,5)=v(i1,5)
        do i2=1,4
          throe(i2)=throe(i2)+p(i1,i2)
          pth(nth,i2)=p(i1,i2)
          vth(nth,i2)=v(i1,i2)
        enddo
      enddo

c     Move parton list nqb steps downward since nqb+1
	call updad_pyj(n1,nqb+1,nqb)   
	n1=n1-nqb
	nqb=0
100	continue
   
600	format(20(1x,i3))
400   if(n1-nqb.eq.0)goto 505   ! get out of 'coal'
      if(n1-nqb.lt.3)goto 500   ! throw away those quarks

c     Proceed for baryon production
	call barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,0)
	if(n1-nqb.eq.0)goto 505
500   continue   

c     Throw away remained q   150922
503   continue
	if(n1-nqb.eq.0)goto 505
	ithroq=ithroq+n1-nqb
	do i1=nqb+1,n1
	  nth=nth+1   ! 110905 171122 Lei corrects it from nth=nth+i1
	  j=nth
	  kk=k(i1,2)
	  ich=ich+pychge(kk)
	  kth(nth,2)=k(i1,2)
	  vth(nth,5)=v(i1,5)
	  do i2=1,4
	    throe(i2)=throe(i2)+p(i1,i2)
	    pth(nth,i2)=p(i1,i2)
	    vth(nth,i2)=v(i1,i2)
	  enddo
	enddo
	n1=0   ! 110905
505	continue

c150922	Reconstruct parton list for calling 'coal' again, in case of iphas=1   
	if(iphas.eq.1)then
c110905
c     Make the parton list ('sa37') in order of qba and q
c150922
        numb(1)=0   
        numb(2)=0
        numb(3)=0
c       numb(1),(2), and (3): the order # of last g,qba & q in 'pyjets'
c        1: refers to g (no gluon at all now), 2: qba, 3: q
        jh=nth
        jl=0
2030    continue
        do j=jl+1,jh
        kf=kth(j,2)
        kfab=iabs(kf)
        if(kfab.lt.7 .and. kf.gt.0)then   ! q, consider d,u,s,c,b,t only
        nth=nth+1      
        numb(3)=numb(3)+1  
        do i4=1,5
        kth(nth,i4)=kth(j,i4)
        pth(nth,i4)=pth(j,i4)
        vth(nth,i4)=vth(j,i4)
        enddo
c       move particle list 'sa37' one step downward since j+1 to nth
        call updad_pyj(nth,j+1,1)
        nth=nth-1
        jh=jh-1
        jl=j-1
        goto 2030
        endif
        enddo
        numb(1)=0
        numb(2)=nth-numb(3)
        numb(3)=nth
        n1=numb(1)
        n2=numb(2)
        n3=nth

	  nqb=ithrob
	  n1=nqb+ithroq
	  do i1=1,n1
	    do j1=1,5
            k(i1,j1)=kth(i1,j1)
            p(i1,j1)=pth(i1,j1)
            v(i1,j1)=vth(i1,j1)    
            enddo
          enddo  
        do i1=n1+1,kszj
          do j1=1,5
            k(i1,j1)=0
            p(i1,j1)=0.
            v(i1,j1)=0.
          enddo
        enddo
	endif   ! 150922

	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine findm(kf1,kf2,cm,kfii,amasi,isucc,iflav)
c   Find out the primary meson from mesonic table according to kf1 & kf2
c	cm : invariant mass of kf1 & kf2
c	kfii : flavor code of the primary meson
c	amasi : mass of the primary meson
c	isucc = 1 : success
c             = 0 : fail
c	iflav = 1 : kf1>0,do not need to permute kf1 & kf2
c	      = -1 : kf1<0,need to permute kf1 & kf2 (never used now, 241022) 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
	dimension ikf(2)

	if(iflav.eq.1)then   ! 1
	  if1=kf1
	  if2=kf2
	  do 500 i4=1,imc
          kfi=kqh(i4,1)
          kfj=kqh(i4,2)
	    amas1=amash(i4,1)
	    amas1=abs(cm-amas1)
	    amas2=amash(i4,2)
	    amas2=abs(cm-amas2)
          if(kfi.eq.if1 .and. kfj.eq.if2)then   ! 2 success
	    if(proh(i4,1).eq.0 .and. proh(i4,2).eq.0.)goto 500
	    if(proh(i4,1).eq.0 .and. proh(i4,2).ne.0.)goto 506   ! vector
	    if((proh(i4,1).ne.0.and.proh(i4,2).ne.0.).and.amas2.le.amas1)
     c    goto 506   ! vector

c         Proceed for pseudoscalar
          kfii=kfh(i4,1)
          amasi=amash(i4,1)
          proi=proh(i4,1)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 500
          goto 504   ! success

506       kfii=kfh(i4,2)   ! vector
          amasi=amash(i4,2)
          proi=proh(i4,2)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 500
	    goto 504   ! success

        endif   ! 2
500	  continue
        isucc=0   ! fail
	  return
	endif   ! 1

	ikf(1)=kf1
	ikf(2)=kf2
c	Two body permutation = arrangement (2,2)
	do 501 i1=1,2
	  if1=ikf(i1)
	  do 502 i2=1,2
	    if(i2.eq.i1)goto 502
		if2=ikf(i2) 
	    do 503 i4=1,imc
            kfi=kqh(i4,1)
            kfj=kqh(i4,2)
            amas1=amash(i4,1)
            amas1=abs(cm-amas1)
            amas2=amash(i4,2)
            amas2=abs(cm-amas2)
          if(kfi.eq.if1 .and. kfj.eq.if2)then   ! success
          if(proh(i4,1).eq.0 .and. proh(i4,2).eq.0.)goto 503
          if(proh(i4,1).eq.0 .and. proh(i4,2).ne.0.)goto 505   ! vector
          if((proh(i4,1).ne.0.and.proh(i4,2).ne.0.).and.amas2.le.amas1)
     c    goto 505   ! vector

c         Proceed for pseudoscalar
          kfii=kfh(i4,1)
          amasi=amash(i4,1)
          proi=proh(i4,1)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 503
          goto 504   ! success

505       kfii=kfh(i4,2)   ! vector
          amasi=amash(i4,2)
          proi=proh(i4,2)
          ran1=pyr(1)
          if(ran1.gt.proi)goto 503
	    goto 504

        endif
503       continue
502	  continue
501	continue
      isucc=0   ! fail
	return
504   isucc=1   ! success
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine findb(kf0,kf1,kf2,cm,kfii,amasi,isucc,iflav)
c     Find out the primary baryon (antibaryon) from baryonic table  
c	according to kf0,kf1,and kf2,these flavor codes are all > 0
c	cm: invariant mass of kf0,kf1 & kf2
c	kfii : flavor code of the primary baryon
c	amasi : mass of the primary baryon
c	isucc = 1 : success
c       isucc = 0 : fail
c	iflav = 1 : if composing parton is quark
c             =-1 : if composing parton is antiquark
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
	common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
	dimension ikf(3)

	ikf(1)=kf0
	ikf(2)=kf1
	ikf(3)=kf2

c	Three body permutation = arrangement (3,3)
	do 104 i1=1,3
	  if1=ikf(i1)
	  do 105 i2=1,3
	    if(i2.eq.i1)goto 105
	      if2=ikf(i2)
		do 106 i3=1,3
		  if(i3.eq.i2)goto 106
		  if(i3.eq.i1)goto 106
		  if3=ikf(i3)
		    do 107 i4=1,ibc

	  kfi=kqb(i4,1)
	  kfj=kqb(i4,2)
        kfk=kqb(i4,3)
        amas1=amasb(i4,1)
        amas1=abs(cm-amas1)
        amas2=amasb(i4,2)
        amas2=abs(cm-amas2)
        if(kfi.eq.if1 .and. kfj.eq.if2 .and. kfk.eq.if3)then ! success
	  if(prob(i4,1).eq.0.and.prob(i4,2).eq.0)goto 107 ! fail and try again
        if(prob(i4,1).eq.0.and.prob(i4,2).ne.0.)goto 108   ! 3/2
        if((prob(i4,1).ne.0.and.prob(i4,2).ne.0.).and.amas2.le.amas1)
     c  goto 108   ! 3/2 
c	Goto 108, for spin 3/2 decuplet.
c	Proceed for spin 1/2 octet.
        kfii=kfb(i4,1)
        amasi=amasb(i4,1)
	  if(iflav.eq.-1)then
	    kfii=-kfb(i4,1)
	  endif
        proi=prob(i4,1)
        ran1=pyr(1)
        if(ran1.gt.proi)goto 107 ! fail and try again
	  goto 109   ! success

108     kfii=kfb(i4,2)   ! spin 3/2 decuplet
        amasi=amasb(i4,2)
	  if(iflav.eq.-1)then
	    kfii=-kfb(i4,2)
	  endif
        proi=prob(i4,2)
        ran1=pyr(1)
        if(ran1.gt.proi)goto 107 ! fail and try again
	  goto 109   ! success
	endif

107	      continue
106	    continue
105	  continue
104	continue
	isucc=0   ! fail
	return
109	isucc=1   ! success

	return   
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine barpro(n1,nqb,iphas,nba,ibarp,netba,rrp,iway)
c     To compose baryon
c     n1 : total # of partons (q & qba)
c     nqb : total # of qba (qba is ordered before q)  
c     iphas=1: with phase space adjudgment
c          =0: without phase space adjudgment
c	ibarp: statistic number of baryon
c	netba: number of baryons keeping baryon conservation 
c	iway: a switch
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000,mplis=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension rcp(3)
	dpmax=adj1(27)
        drmax=adj1(28)	
	nba=0
400	continue

	if(n1-nqb.lt.3)return   ! number of quarks less than three
      do 404 i1=nqb+1,n1-2
	 kf1=k(i1,2)   ! 150922
        do 405 i2=i1+1,n1-1
           kf2=k(i2,2)
          do 406 i3=i2+1,n1
             kf3=k(i3,2)

	  sume=p(i1,4)+p(i2,4)+p(i3,4)
	  sump1=p(i1,1)+p(i2,1)+p(i3,1)
	  sump2=p(i1,2)+p(i2,2)+p(i3,2)
          sump3=p(i1,3)+p(i2,3)+p(i3,3)
	  cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	  if(cm.gt.1.e20)cm=1.e20
	  if(cm.le.0.)goto 406   ! fail 071204
          cm=sqrt(cm)

c     Find out the primary baryon from hadron table due to kf1,kf2 & kf3
        call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)

        if(isucc.eq.0)goto 406   ! fail, and keep on cycle, try again.

c	Proceed for success   
c	Phase space adjudgment
        if(iphas.eq.1)then
	    call phas(i1,i2,i3,isucc,3)
          if(isucc.eq.0)goto 406   ! fail
        endif

c	  Proceed for success
        ibarp=ibarp+1
	  nba=nba+1

c       Give proper variables to the primary baryon
        nnol=nn
        nn=nn+1
	kn(nn,1)=1
        kn(nn,2)=kfii
        kn(nn,3)=0
        pn(nn,5)=amasi
	pn(nn,1)=sump1
	pn(nn,2)=sump2
	pn(nn,3)=sump3
       	pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	pnnmm=amasi*amasi+pnnm
	if(pnnmm.gt.1.e20)pnnmm=1.e20
        if(pnnmm.le.0.)pnnmm=1.e-20
	  pnnn=sqrt(pnnmm)
	  pn(nn,4)=pnnn
	  dele=sume-pnnn

c     Produced hadron is arranged among constituent partons randomly.
	  pyrx=pyr(1)
	  pyry=pyr(1)
	  pyrz=pyr(1)
	  rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)+pyrz*v(i3,1)
	  rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)+pyrz*v(i3,2)
	  rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)+pyrz*v(i3,3)

c     Move parton list one step downward from i3+1 to n1
411	  call updad_pyj(n1,i3+1,1)
	  n1=n1-1

c     Move parton list one step downward from i2+1 to n1
	  call updad_pyj(n1,i2+1,1)
	  n1=n1-1

c     Move parton list one step downward from i1+1 to n1
	  call updad_pyj(n1,i1+1,1)
	  n1=n1-1

c     Share the surplus energy. 
	  if(n1+nn.gt.0)then
	    dele=dele/float(n1+nn)
	    if(n1.gt.0)then
		do i4=1,n1
		  p(i4,4)=p(i4,4)+dele
		  if(dele.lt.0.)then
		    if(p(i4,4).lt.0.)p(i4,4)=p(i4,4)-dele
		    pabs=abs(p(i4,3))
		    if(pabs.ge.p(i4,4))p(i4,4)=p(i4,4)-dele
		  endif
		enddo
	    endif
	    if(nn.gt.0)then
		do i4=1,nn
		  pn(i4,4)=pn(i4,4)+dele
		  if(dele.lt.0.)then
		  if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
		  pabs=abs(pn(i4,3))
		  if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
		  endif
		enddo
	    endif
	  endif
	  if(iway.eq.1)then 
	    if(nba.lt.netba)goto 400 ! recycle all the partons remained, do again.
	    if(nba.eq.netba)return
	  endif
        if(iway.eq.2 .and. nba.eq.1)return
	  if(iway.eq.0)goto 400   ! 121204
c   iway=1: when the baryon number equal net baryon, return. Used in creat net baryon
c   iway=2: it can return when there generate one baryon. 
c   iway=0: check all the probability of parton constituent baryon, then return.
406       continue   ! fail
405     continue
404   continue
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine an_barpro(n1,nqb,i1,iphas,nba,ibarm,netba,rrp,isu,
     c	iway)
c   To compose an anti-baryon
c     n1 : total # of partons (q & qba)
c     nqb : total # of qba (qba is ordered before q)
c     i1: antiquark wanted to compose antibaryon   
c     iphas=1: with phase space adjudgment
c          =0: without phase space adjudgment
c     ibarm: statistic number of anti-baryon
c     -netba: number of anti-baryons keeping baryon conservation 
c     isu: =1 success
c          =0 fail
c     iway: a switch
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000,mplis=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension rcp(3)
	dpmax=adj1(27)
	drmax=adj1(28)	
	isu=1
	nba=0

	if(nqb.lt.3)goto 100   ! number of qbar less than three,return
	kf1=k(i1,2)
      do 405 i2=1,nqb
	  if(i2.eq.i1)goto 405
        kf2=k(i2,2)
        do 406 i3=1,nqb
	    if(i3.eq.i1)goto 406
	    if(i3.eq.i2)goto 406
	    kf3=k(i3,2)

          sume=p(i1,4)+p(i2,4)+p(i3,4)
          sump1=p(i1,1)+p(i2,1)+p(i3,1)
          sump2=p(i1,2)+p(i2,2)+p(i3,2)
          sump3=p(i1,3)+p(i2,3)+p(i3,3)
	  cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	  if(cm.gt.1.e20)cm=1.e20
	  if(cm.le.0.)goto 406   ! fail 071204
        cm=sqrt(cm)

c     Find out the primary antibaryon from hadron table according to kf1,kf2 
c	 & kf3
        call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)

        if(isucc.eq.0)goto 406   ! fail

c	Proceed for success   
c	Phase space adjudgment
        if(iphas.eq.1)then
	    call phas(i1,i2,i3,isucc,3)
          if(isucc.eq.0)goto 406   ! fail
        endif

c	Proceed for success
	  goto 400
406     continue   ! fail

405   continue
	goto 100
400	continue
      ibarm=ibarm+1
	nba=nba+1

c     Give proper variables to the primary antibaryon
        nnol=nn
        nn=nn+1
        kn(nn,1)=1
        kn(nn,2)=kfii
        kn(nn,3)=0
        pn(nn,5)=amasi
        pn(nn,1)=sump1
        pn(nn,2)=sump2
        pn(nn,3)=sump3
	pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	pnnmm=amasi*amasi+pnnm
	if(pnnmm.gt.1.e20)pnnmm=1.e20
	if(pnnmm.le.0.)pnnmm=1.e-20
	pnnn=sqrt(pnnmm)
	pn(nn,4)=pnnn
	dele=sume-pnnn

c     Arrange produced particle on the surface of sphere with radius
c     rrp and centered at the center of mass

c	Produced hadron is arranged among contituent partons randomly
	pyrx=pyr(1)
	pyry=pyr(1)
	pyrz=pyr(1)
        rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)+pyrz*v(i3,1)
        rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)+pyrz*v(i3,2)
        rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)+pyrz*v(i3,3)

c     Move parton list one step downward from i3+1 to n1
411	call updad_pyj(n1,i3+1,1)
	if(i1.gt.i3)i1=i1-1
	if(i2.gt.i3)i2=i2-1
	nqb=nqb-1
	n1=n1-1

c     Move parton list one step downward from i2+1 to n1
	call updad_pyj(n1,i2+1,1)
	if(i1.gt.i2)i1=i1-1
	nqb=nqb-1
	n1=n1-1

c     Move parton list one step downward from i1+1 to n1
	call updad_pyj(n1,i1+1,1)
	nqb=nqb-1
	n1=n1-1

c	Share the surplus energy.
	if(n1+nn.gt.0)then
	  dele=dele/float(n1+nn)
	  if(n1.gt.0)then
	    do i4=1,n1
		p(i4,4)=p(i4,4)+dele
		if(dele.lt.0.)then
		  if(p(i4,4).lt.0.)p(i4,4)=p(i4,4)-dele
		  pabs=abs(p(i4,3))
		  if(pabs.ge.p(i4,4))p(i4,4)=p(i4,4)-dele
		endif
	    enddo
	  endif
	  if(nn.gt.0)then
          do i4=1,nn
		pn(i4,4)=pn(i4,4)+dele
		if(dele.lt.0.)then
		  if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
		  pabs=abs(pn(i4,3))
		  if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
		endif
	    enddo
	  endif
	endif
c   iway=1: creat an antibaryon and return
c   iway=2: creat an antibaryon and a baryon, then return
	if(iway.eq.1 .and. nba.eq.1)return
      if(iway.eq.2 .and. nba.eq.1)then
c     An antibaryon generation must be followed immediately a baryon
c     generation keeping baryon conservation
	  call barpro(n1,nqb,iphas,nbaa,ibarp,netba,rrp,2)
	  return
	endif

100	isu=0
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mespro(n1,nqb,i1,iphas,nme,imes,rrp,isu)
c       Compose a meson
c       n1 : total # of partons (q & qba)
c       nqb : total # of qba (qba is ordered before q)
c	i1: antiquark wanted to compose a meson
c       iphas=1: with phase space adjudgment
c            =0: without phase space adjudgment
c	imes: statistic number of meson
c       isu: =1 success
c            =0 fail
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000,mplis=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
      common/sa1_h/nn,non1_h,kn(kszj,5),pn(kszj,5),rn(kszj,5)
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa24/adj1(40),nnstop,non24,zstop
      dimension rcp(3)
	dpmax=adj1(27)
        drmax=adj1(28)
	isu=1
	nme=0

c	write(9,*)'enter mespro'

	if(n1.eq.nqb)goto 100   ! 300105, no quark, return
	kf1=k(i1,2)
	do 102 i2=nqb+1,n1
	kf2=k(i2,2)
	sume=p(i1,4)+p(i2,4)
	sump1=p(i1,1)+p(i2,1)
	sump2=p(i1,2)+p(i2,2)
        sump3=p(i1,3)+p(i2,3)
	cm=sume*sume-sump1*sump1-sump2*sump2-sump3*sump3
	if(cm.gt.1.e20)cm=1.e20
	if(cm.le.0.)goto 102   ! fail 071204
        cm=sqrt(cm)

c     Find out primary meson from hadronic table according to kf2 & kf1
        call findm(kf2,kf1,cm,kfii,amasi,isucc,1)

        if(isucc.eq.0)goto 102   ! fail
c	Proceed for success
c	Phase space adjudgment
	  if(iphas.eq.1)then
	    call phas(i1,i2,0,isucc,2)
          if(isucc.eq.0)goto 102   ! fail
        endif

c	Proceed for success
	  imes=imes+1
	  nme=nme+1

c     Give proper variables to the primary meson
        nnol=nn
        nn=nn+1
        kn(nn,1)=1
        kn(nn,2)=kfii
        kn(nn,3)=0
        pn(nn,5)=amasi
	pn(nn,1)=sump1
	pn(nn,2)=sump2
	pn(nn,3)=sump3
	pnnm=sump1*sump1+sump2*sump2+sump3*sump3
	pnnmm=amasi*amasi+pnnm
        if(pnnmm.gt.1.e20)pnnmm=1.e20
        if(pnnmm.le.0.)pnnmm=1.e-20
	pnnn=sqrt(pnnmm)
	pn(nn,4)=pnnn
c090922
c       do j1=1,5
c       pnnnj1=pn(nn,j1)
c       aabbss=abs(pnnnj1)
c       pnnnj1=pnnnj1/aabbss   ! symbol, i. e. (+ or -)1
c       dpmax=dpmax/1000.0
c       if(aabbss.ge.dpmax)then
c       aabbss=dpmax
c       pn(nn,j1)=aabbss*pnnnj1
c       endif
c       enddo
c090922        
	  dele=sume-pnnn

c	Produced hadron is seded in between contituent partons randomly
	  pyrx=pyr(1)
	  pyry=pyr(1)
	  rn(nn,1)=pyrx*v(i1,1)+pyry*v(i2,1)
	  rn(nn,2)=pyrx*v(i1,2)+pyry*v(i2,2)
	  rn(nn,3)=pyrx*v(i1,3)+pyry*v(i2,3)

c     Move parton list one step downward since i2+1
111	  call updad_pyj(n1,i2+1,1)
        n1=n1-1

c     Move parton list one step downward since i1+1
        call updad_pyj(n1,i1+1,1)
        nqb=nqb-1
        n1=n1-1
	  if(n1+nn.gt.0)then
	    dele=dele/float(n1+nn)
	    if(n1.gt.0)then
		do i4=1,n1
		  p(i4,4)=p(i4,4)+dele
		  if(dele.lt.0.)then
		    if(p(i4,4).lt.0.)p(i4,4)=p(i4,4)-dele
		    pabs=abs(p(i4,3))
		    if(pabs.ge.p(i4,4))p(i4,4)=p(i4,4)-dele
		  endif
		enddo
	    endif
	    if(nn.gt.0)then
            do i4=1,nn
		  pn(i4,4)=pn(i4,4)+dele
		  if(dele.lt.0.)then
		    if(pn(i4,4).lt.0.)pn(i4,4)=pn(i4,4)-dele
		    pabs=abs(pn(i4,3))
		    if(pabs.ge.pn(i4,4))pn(i4,4)=pn(i4,4)-dele
		  endif
            enddo
	    endif
	  endif
	  if(nme.eq.1)return
102	enddo   ! 090922 continue-> enddo
100	isu=0   ! 300105
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conse_c(pp,ps,npl,np,nstep)
c     Adjust four momentum conservation by iteration,no more than
c	 5000 iterations
c       pp : four momentum of particle
c       ps : the four momentum should be conserved to
c	npl : order # of the first particle
c       np : order # of last particle
c	nstep : interval of the step
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 140604 060813
        dimension pp(kszj,5),ps(4),ff(kszj),pxyz(3),arp(3)

        ps4=ps(4)
        do i=1,3
        pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=npl,np,nstep
        es=es+pp(i,4)
        enddo
        fr=es/ps4

        do i=npl,np,nstep
        amas=pp(i,5)
        amas2=amas*amas
        ppm=pp(i,4)
        ppf=ppm/fr
c090700
	den2=ppm*ppm-amas2
	if(den2.le.0.)then
	den2=1.e-15
	endif
c090700
        ff(i)=sqrt(abs(ppf*ppf-amas2)/den2)
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=abs(1.-pxyz(i)/ps(i))
        enddo
        if(abs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)-ps(i)
        pxyz(i)=pxyz(i)/(float(np-npl)/float(nstep)+1)
        enddo
        do i=npl,np,nstep
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
	pp5=pp(i,5)
	pp52=pp5*pp5
        pp(i,4)=sqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
        enddo
        jj=jj+1
        if(jj.lt.5000)goto 100
200     return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccrcccccccccccccccc
	subroutine ffm(ii,kf0,igen,iflav)   ! 280822
c	qqbar pair generation according to Field-Feynman model
c280822 i. e. energetic q (qbar) de-exciatation
c	ii : the order # of source quark (or antiquark)
c	kf0 : flavor code of source quark (or antiquark)
c	igen : # of generations per source q (qba)  ! 280822 
c	iflav = 1 : if source parton is quark (kf0>0)
c             =-1 : if source parton is antiquark (kf0<0)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
      common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp   ! 280822
      common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
      common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension p0(4),p1(4),p1c(4),p00(4),rc(3),rr(3),pnn(kszj,5),
     c	 peo(5),pdec(20,5)   ! 090922

        adj16=adj1(16)   ! # of allowed generations 280822
	adj17=adj1(17)   ! threshold energy
	n0=n

	do i1=1,3
	rc(i1)=v(ii,i1)   ! three coordinate of source q (qba)
	enddo
	do i2=1,4
	p0(i2)=p(ii,i2)   ! four momentum of source q (qba)
        p00(i2)=p0(i2)
	enddo
	kf00=kf0   ! kf code of source q (qba)
	e0=p0(4) ! E
        p00(3)=p0(3) 

	if(e0.lt.0.)return   ! stop generation 280822 

c	qqba creation from vacuum
	igen=0   ! count # of generation
100	continue

c	sample the flavor for generated  q (qba)
	eg=e0
        call break_f(eg,kf,amq)
        amasi=2*amq    ! 220222 mass of created qqbar pair (a object) 280822
	kf1=kf
	kf2=-kf
	if(iflav.eq.1)then
	kf1=-kf
	kf2=kf
	endif
	
104	continue   

c       Sample transverse momentum of created qqba pair by thermodynamic 
c        distribution (exponential distribution) ! 081022        
        call thptdi(pt,px,py)
        p1(1)=px    ! of created qqba pair
        p1(2)=py    ! of created qqba pair
        ppt=pt*pt   ! pt square of created qqba pair

c	Sample z (energy fraction of created qqba pair taking from source
c        q (qba)) by Field-Feymman fragmentation function 081022
        call funcz(z1)

	e1=z1*e0 ! energy of first generation qqba pair
	p1(4)=e1
	p13=-(amasi*amasi+ppt)+e1*e1  ! sqruare p1(3) of created qqba pair
        if(p13.lt.1.d-20)p13=1.d-20
        p1(3)=sqrt(p13)            

c	Fill the generated q & qba into parton list ('pyjets') after n

c       Give four position to generated q & qba        
c        generated q and qba are arranged around sourve parton within 0.5 fm 
c        randumly in each one of the three coordinates and has same fourth 
c        coordinate as sourve parton
	do i=1,3
	rr(i)=pyr(1)*0.5
	v(n+1,i)=rc(i)+rr(i)
        if(pyr(1).gt.0.5)v(n+1,i)=rc(i)-rr(i)
        rr(i)=pyr(1)*0.5
	v(n+2,i)=rc(i)+rr(i)
        if(pyr(1).gt.0.5)v(n+2,i)=rc(i)-rr(i)
	enddo
	v(n+1,4)=v(ii,4) ! ii: order # of sourse parton
        v(n+2,4)=v(ii,4) ! n+1 (n+2): order # of new generated parton

        k(n+1,1)=2   ! 'A' 250922  
        k(n+2,1)=1   ! 'V' 250922 

c	Give four momentum to generated q (qba)   ! 090922
	amq=0.5*amasi
        decsuc=1
        call decmom(p1,pdec,amq,amq,decsuc)  ! for 'decay method'
c090922 p1: four momentum of decaying particle
c       pdec: four momentum of decayed particles   ! 090922
c       amq: mass of decayed particle
c090922 'decmom': in parini_23.f
        if(decsuc.eq.1)then
        do i1=1,4
        p(n+1,i1)=pdec(1,i1)
        enddo
        do i1=1,4
        p(n+2,i1)=pdec(2,i1)
        enddo      
        if(p00(3).lt.0.0)then  ! for the negative direction p_z
        p(n+1,3)=-p(n+1,3)
        p(n+2,3)=-p(n+2,3)
        endif
        goto 300  
        elseif(decsuc.eq.0)then
        goto 400   ! random three momentum method
        else        
        endif
400     do i1=1,3   ! random three momentum method
	pii=pyr(1)*p1(i1)
	p(n+1,i1)=pii
	p(n+2,i1)=p1(i1)-pii
	enddo

        pn11=p(n+1,1)   ! pnn(1,1) 280822
        pn12=p(n+1,2)   ! pnn(1,2) 280822
        pn13=p(n+1,3)   ! pnn(1,3) 280822
        agsq=amq*amq+pn11*pn11+pn12*pn12+pn13*pn13   ! 280822
        if(agsq.lt.1.d-20)agsq=1.d-20   ! 280822
	p(n+1,4)=sqrt(agsq)   ! 280822
        pn21=p(n+2,1)   ! pnn(2,1) 280822
        pn22=p(n+2,2)   ! pnn(2,2) 280822
        pn23=p(n+2,3)   ! pnn(2,3) 280822
        agsq=amq*amq+pn21*pn21+pn22*pn22+pn23*pn23   ! 280822
        if(agsq.lt.1.d-20)agsq=1.d-20   ! 280822
        p(n+2,4)=sqrt(agsq)   ! 280822        

300     ppsa(4)=ppsa(4)+(p0(4)-p(n+1,4)-p(n+2,4))   ! 280822

c	Give other properties to generated q and qba 
	p(n+1,5)=amq
        p(n+2,5)=amq
	k(n+1,2)=kf1
	k(n+2,2)=kf2

c	Give proper variables to the remnant parton
	e1c=e0-e1   ! conservation
	do i3=1,3   ! three momentum conservation
	p1c(i3)=p0(i3)-p1(i3)   
	enddo
        ee=amq*amq+p1c(1)*p1c(1)+p1c(2)*p1c(2)+p1c(3)*p1c(3)
        if(ee.lt.1.d-20)ee=1.d-20
        ee=sqrt(ee)
        ppsa(4)=ppsa(4)+(e1c-ee)
c       p0 refers to original q (qbar), p1 refers to generated qqbar pair, 
c       p1c refers to remnant (original parion after generating qqbr pair)
	p1c(4)=e1c

	kf0=kf00
	if(kf0.gt.0)iflav=1
	if(kf0.lt.0)iflav=-1
	e0=e1c
	do i3=1,4
	p0(i3)=p1c(i3)
	enddo
	igen=igen+1

c       write(9,*)'iii,ii,eorig=',iii,ii,p(ii,4)   ! 090922
c       write(9,*)'px,py,z1=',px,py,z1   ! 090922
c       write(9,*)'p(n+1)=',(p(n+1,i1),i1=1,4)   ! 090922
c       write(9,*)'p(n+2)=',(p(n+2,i1),i1=1,4)   ! 090922
c       write(9,*)'p0(1-4)=',(p0(i1),i1=1,4)   ! 090922

	n=n+2
	if(igen.ge.adj16)goto 106   ! stop generation 280822
	if(e0.lt.adj17)goto 106   ! stop generation

	goto 100

106	continue
	if(igen.eq.0)return
c	Update four momentum of the remnant of ii-th source q (qba) 280822
	do i3=1,4   ! 280822
	p(ii,i3)=p0(i3)
	enddo

c280822 do i1=n0+1,n
c       j1=i1-n0
c       do i2=1,4
c       pnn(j1,i2)=p(i1,i2)
c       enddo
c       pnn(j1,5)=p(i1,5)
c280822 enddo
c280822 ii1=n-n0   ! new generated q and qba from ii-th source q (qba)

c	Incluse remnant of ii-th source q (qba)
c280822 ii1=ii1+1   
c       do i2=1,4
c       pnn(ii1,i2)=p(ii,i2)
c       enddo
c280822 pnn(ii1,5)=p(ii,5)

	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine funcz(zz)
c       Sample zz from Field-Feymman fragmentation function using selecting 
c        sample method	  
c	Distribution function: f(z)dz=[1-a+3*a*(1-z)**2]dz, 0<z<1.
c	The largest value of function: fmax=f(0)=1-a+3*a=1+2*a.
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      common/sa24/adj1(40),nnstop,non24,zstop
	a=adj1(6)   ! 0.77  
	b=3.*a   
	fmax=1-a+b   
100	ran1=pyr(1)
	ran2=pyr(1)
	fm=ran1*fmax
c101204	fz=1.-0.77+3*0.77*(1.-ran2)**2.
c101204	if(fm.le.fz)goto 100
	fz=1.-a+b*(1.-ran2)**2.   ! 101204
      if(fm.gt.fz)goto 100   ! 101204
	zz=ran2  	
	return 
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine phas(i1,i2,i3,isucc,j)
c	The phase space judgement.
c	j=2 for meson
c	j=3 for baryon
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000,mplis=80000)   ! 150922
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 150922
      common/sa24/adj1(40),nnstop,non24,zstop
	dimension ri1(3),ri2(3),ri3(3),pi1(3),pi2(3),pi3(3)
	delc=adj1(22)

	if(j.eq.2)goto 100   ! for meson
c	 proceed for baryon 
        ri1(1)=v(i1,1)
        ri1(2)=v(i1,2)
        ri1(3)=v(i1,3)
        ri2(1)=v(i2,1)
        ri2(2)=v(i2,2)
        ri2(3)=v(i2,3)
        ri3(1)=v(i3,1)
        ri3(2)=v(i3,2)
        ri3(3)=v(i3,3)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        ri131=ri1(1)-ri3(1)
        ri132=ri1(2)-ri3(2)
        ri133=ri1(3)-ri3(3)
        ri231=ri2(1)-ri3(1)
        ri232=ri2(2)-ri3(2)
        ri233=ri2(3)-ri3(3)
	delr12=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
	delr13=sqrt(ri131*ri131+ri132*ri132+ri133*ri133)
	delr23=sqrt(ri231*ri231+ri232*ri232+ri233*ri233)
        pi1(1)=p(i1,1)
        pi1(2)=p(i1,2)
        pi1(3)=p(i1,3)
        pi2(1)=p(i2,1)
        pi2(2)=p(i2,2)
        pi2(3)=p(i2,3)
        pi3(1)=p(i3,1)
        pi3(2)=p(i3,2)
        pi3(3)=p(i3,3)
        pi121=pi1(1)-pi2(1)
        pi122=pi1(2)-pi2(2)
        pi123=pi1(3)-pi2(3)
        pi131=pi1(1)-pi3(1)
        pi132=pi1(2)-pi3(2)
        pi133=pi1(3)-pi3(3)
        pi231=pi2(1)-pi3(1)
        pi232=pi2(2)-pi3(2)
        pi233=pi2(3)-pi3(3)
	delp12=sqrt(pi121*pi121+pi122*pi122+pi123*pi123)
	delp13=sqrt(pi131*pi131+pi132*pi132+pi133*pi133)
	delp23=sqrt(pi231*pi231+pi232*pi232+pi233*pi233)
	del12=delr12*delp12
	del13=delr13*delp13
	del23=delr23*delp23

	if(del12.le.delc.and.del13.le.delc.and.del23.le.delc)then
	  isucc=1
	else
	  isucc=0
	endif
	return

100	continue   ! for meson
        ri1(1)=v(i1,1)
        ri1(2)=v(i1,2)
        ri1(3)=v(i1,3)
        ri2(1)=v(i2,1)
        ri2(2)=v(i2,2)
        ri2(3)=v(i2,3)
        ri121=ri1(1)-ri2(1)
        ri122=ri1(2)-ri2(2)
        ri123=ri1(3)-ri2(3)
        delr=sqrt(ri121*ri121+ri122*ri122+ri123*ri123)
        pi1(1)=p(i1,1)
        pi1(2)=p(i1,2)
        pi1(3)=p(i1,3)
        pi2(1)=p(i2,1)
        pi2(2)=p(i2,2)
        pi2(3)=p(i2,3)
        pi121=pi1(1)-pi2(1)
        pi122=pi1(2)-pi2(2)
        pi123=pi1(3)-pi2(3)
        delp=sqrt(pi121*pi121+pi122*pi122+pi123*pi123)
        delrp=delr*delp
	if(delrp.le.delc)then
	  isucc=1
	else
	  isucc=0
	endif
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_f(eg,kf,amq)
c       sample flavor (mass) of generated qqbar pair  
c       eg: energy of original q or qbar 
c       kf (amq): flavor code (mass) of generated quark
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa38/csp_31,csp_32,csp_41,csp_42,csp_43,csp_51,csp_52,
     c   csp_53,csp_54,csp_61,csp_62,csp_63,csp_64,csp_65   ! 161022
        amd=pymass(1)   ! constituent mass in GeV
        amu=pymass(2)   ! amu=amd
        ams=pymass(3)
        amc=pymass(4)
        amb=pymass(5)
        amt=pymass(6)
        amuu=2*amu
        amdd=2*amd
        amss=2*ams
        amcc=2*amc  
        ambb=2*amb
        amtt=2*amt
        aa=pyr(1)
c	  if(eg.lt.amuu)goto 200   ! throw away amuu
c161022 if(eg.ge.amdd .and. eg.lt.amss)then   ! d,u
	  if(eg.lt.amss)then   ! d,u (with same flavor generation probability)
           if(aa.le.0.5)then
           kf=1   ! d
	   amq=amd
           else  
           kf=2   ! u
	   amq=amu
	   endif
          goto 200        
	  endif

	  if(eg.ge.amss .and. eg.lt.amcc)then   ! d,u,s
	   if(aa.le.csp_31)then
	   kf=1   ! d    
           amq=amd
           elseif(aa.gt.csp_31 .and. aa.le.csp_32)then
           kf=2   ! u  
           amq=amu       
	   else
	   kf=3   ! s
	   amq=ams
	   endif
          goto 200
          endif

        if(eg.ge.amcc .and. eg.lt.ambb)then ! d,u,s,c
         if(aa.le.csp_41)then          
         kf=1   ! d    
         amq=amd
         elseif(aa.gt.csp_41 .and. aa.le.csp_42)then
         kf=2   ! u  
         amq=amu       
         elseif(aa.gt.csp_42 .and. aa.le.csp_43)then
         kf=3
         amq=ams
         else
         kf=4
         amq=amc  
         endif
        goto 200
        endif

        if(eg.ge.ambb .and. eg.lt.amtt)then ! d,u,s,c,b
         if(aa.le.csp_51)then         
         kf=1
         amq=amd
         elseif(aa.gt.csp_51 .and. aa.le.csp_52)then
         kf=2
         amq=amu
         elseif(aa.gt.csp_52 .and. aa.le.csp_53)then
         kf=3
         amq=ams
         elseif(aa.gt.csp_53 .and. aa.le.csp_54)then
         kf=4
         amq=amc
         else
         kf=5
         amq=amb   
         endif
        goto 200
        endif

        if(eg.ge.amtt)then ! d,u,s,c,b,t
         if(aa.le.csp_61)then
         kf=1
         amq=amd
         elseif(aa.gt.csp_61 .and. aa.le.csp_62)then
         kf=2
         amq=amu
         elseif(aa.gt.csp_62 .and. aa.le.csp_63)then
         kf=3
         amq=ams
         elseif(aa.gt.csp_63 .and. aa.le.csp_64)then
         kf=4
         amq=amc
         elseif(aa.gt.csp_64 .and. aa.le.csp_65)then
         kf=5
         amq=amb
         else
         kf=6
         amq=amt
         endif
        endif
200     continue
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tdgaus(v,pmax,np,pp)
c.... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
c.... set pmax < 0 if pmax should be infinity.
c.... np : the total # of particles wanted to sample their transverse momenta
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 220312 240412 131212
      dimension pp(50,2)
      do 30 i=1,np
        p2 = 0
        if(v.le.1.e-8)return
        if(pmax.lt.1.E-9)return
        if(pmax.lt.0)then
          a = 1.
          goto 10
        endif
        aa=-pmax/v
        if(aa.lt.-70)then
          a=1.
          goto 10
        endif
        a = 1. - exp(aa)
10      p2 = -v*log(max(1.e-20,1. - a*pyr(1)))
        if(p2.LT.0.)goto 10
        ps=sqrt(p2)
        fi=2.*3.1415926*pyr(1)
c220312 randomly sample [px,py] on circle of sphere with radius ps
c220312	pp(i,1)=ps*cos(fi)
c220312	pp(i,2)=ps*sin(fi)
c220312 randomly sample [px,py] on circle of ellipsoid with half major axis
c220312 of ps*(1+smadel) and half minor axis of ps*(1-smadel)
        pp(i,1)=ps*cos(fi)*(1+smadel)   ! 220312
        pp(i,2)=ps*sin(fi)*(1-smadel)   ! 220312
c220312 note: ps is not in the dimension list
30    continue
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo_glu   ! 160822
c       moves gluons from 'pyjets' to 'sa36'   
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	nglu=0
	do i1=1,kszj
        do j1=1,5
        kglu(i1,j1)=0
        pglu(i1,j1)=0.
        vglu(i1,j1)=0.
        enddo
        enddo
        jb=0
201     do i1=jb+1,n   
        kf=k(i1,2)
        kfab=iabs(kf)
        eng=p(i1,4)
        if(kfab.ne.21)then   ! stay   
        jb=jb+1
        goto 202
        endif

        nglu=nglu+1
        do i2=1,5
        kglu(nglu,i2)=k(i1,i2)
        pglu(nglu,i2)=p(i1,i2)
        vglu(nglu,i2)=v(i1,i2)
        enddo
        if(i1.eq.n)then
        n=n-1
        goto 203
        endif
c	move particle list one step downward from i1+1 to n
        do j=i1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 201   ! this statement is needless 090122
202     enddo   ! do loop 
203     continue
	return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break_glu   ! 160822
c       break-up gluon (in 'sa36') and give flavor (mass) and four momentum 
c        (position) to broken objects, which are assumed to be a string and 
c        filling in 'pyjets'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 220110
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
c       write(22,*)'enter break_glu n=',n
        amu=pymass(2)   ! 271022
        amuu=2*amu   ! 271022
c       throw away g with energy<amuu
        do i1=1,nglu
        eg=pglu(i1,4)
        if(eg.lt.amuu)then
        do i2=1,4
        ppsa(i2)=ppsa(i2)+pglu(i1,i2)
        enddo
c       move particle list ('sa36') one step downward from i1+1 to nglu
        do j=i1+1,nglu
        do jj=1,5
        kglu(j-1,jj)=kglu(j,jj)
        pglu(j-1,jj)=pglu(j,jj)
        vglu(j-1,jj)=vglu(j,jj)
        enddo
        enddo
        nglu=nglu-1
        endif
        enddo     

c       g (in 'sa36') -> qq_{bar} (as a string filling in 'pyjets')
100     do i1=1,nglu   ! do loop over gluons
c       write(22,*)'break_glu 1 nglu,n=',nglu,n
	eg=pglu(i1,4)
        call break_f(eg,kf,amq)
        kf1=kf                                         
        kf2=-kf
        am1=amq
        am2=amq

c       write(22,*)'break_glu 2 nglu,n=',nglu,n         
        k(n+1,1)=2   ! 'A'
        k(n+2,1)=1   ! 'V'
        k(n+1,2)=kf1
        k(n+2,2)=kf2
        k(n+1,3)=0
        k(n+2,3)=0
        k(n+1,4)=0
        k(n+2,4)=0
        k(n+1,5)=0
        k(n+2,5)=0
c       p(n+1,5)=am1
c       p(n+2,5)=am2

c       give four momentum to the breaked quarks
	call bream_glu(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord_glu(i1)
        if(i1.eq.nglu)then
        n=n+2
        goto 200
        endif
        n=n+2
c       goto 100
        enddo   ! do loop over gluons ended
200     return
        end


        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream_glu(ii,kf1,kf2)   ! 160822
c       give four momentum to the broken quarks
c       ii: line number of initial gluon in 'sa36'
c       kf1,kf2: flavor codes of broken quarks
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   
        am1=pymass(kf1)
        am2=pymass(kf2)
        pp(1,5)=am1   ! mass of first broken quark
        pp(2,5)=am2   ! mass of second broken quark
c       pp : four momenta & mass of broken quarks, local variable 
        do i1=1,4
        ps(i1)=pglu(ii,i1)   
        enddo
c       ps : four momentum of initial gluon, local variable 
	goto 400   ! for 'decay method'
401	do i1=1,3   ! for 'random three momentum method'
        pi(i1)=pyr(1)*pglu(ii,i1)
        pp(1,i1)=pi(i1)
        pp(2,i1)=ps(i1)-pi(i1)
        enddo
c250503
	pp11=pp(1,1)
	pp12=pp(1,2)
	pp13=pp(1,3)
c021005
        pp14=am1*am1+pp11*pp11+pp12*pp12+pp13*pp13
        if(pp14.le.0.)pp14=1.e-20
        pp(1,4)=dsqrt(pp14)
c021005
	pp21=pp(2,1)
	pp22=pp(2,2)
	pp23=pp(2,3)
c021005
        pp24=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
        if(pp24.le.0.)pp24=1.e-20
        pp(2,4)=dsqrt(pp24)
c021005
	goto 300   ! activate it for 'random three momentum method'
c250503
c260503
400	continue   ! for 'decay method'
	decsuc=1
	call decmom(ps,pp,am1,am2,decsuc)
	if(decsuc.eq.0)goto 401   ! return to random three momentum method
300	continue
c       adjust four momentum conservation by iteration,no more than
c        4000 iterations
c	call conser(2,pp,ps)   
c160822
c       do i1=1,2
c       do i2=1,4
c       ppi1i2=pp(i1,i2)
c       pp12=abs(ppi1i2)
c       if(pp12.gt.1.d6)pp12=1.d6
c       pp(i1,i2)=sign(pp12,ppi1i2)
c       enddo
c       enddo
c160822
        do i1=1,4
        p(n+1,i1)=pp(1,i1)
        enddo
        p(n+1,5)=am1
        do i1=1,4
        p(n+2,i1)=pp(2,i1)
        enddo
        p(n+2,5)=am2
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord_glu(ii)   ! 160822
c       give four position to broken quarks
c       first broken quark takes the four position of gluon
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three positions and has same
c        fourth position as gluon
c       ii: order # of gluon in 'sa36'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5) 
        dimension rr(3)

        do i1=1,3
        v(n+1,i1)=vglu(ii,i1)
        rr(i1)=pyr(1)*0.5   ! 261002
        v(n+2,i1)=vglu(ii,i1)+rr(i1)
        if(pyr(1).gt.0.5d0)v(n+2,i1)=vglu(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=vglu(ii,4)
        v(n+2,4)=vglu(ii,4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa36(nn,cc)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa36/nglu,nongu,kglu(kszj,5),pglu(kszj,5),vglu(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(mstu(11),*)i,kglu(i,2),(pglu(i,j),j=1,4)
c	write(9,*)i,kglu(i,2),(pglu(i,j),j=1,4)
c       enddo
        call psum(pglu,1,nglu,peo)
        ich1=0.
        do i1=1,nn
        kf=kglu(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sa36 nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   !
c	write(9,*)peo,ich1/3   !
        return
        end



c*********************************************************************
        subroutine thptdi(pt,px,py)   ! 150822 
c       Generates transverse momentum according to thermodynamic distribution
c        (exponential distribution)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)        
      SAVE /PYDAT1/
 
c       Generate p_T and azimuthal angle, gives p_x and p_y.
        pt=-parj(21)*log(max(1d-10,pyr(1)))
        phi=paru(2)*pyr(1)
c       randomly sample [px,py] on circle of sphere with radius pt
        px=pt*cos(phi)
        py=pt*sin(phi)   
        return
        end 


 
