	subroutine parini(time_ini,parp21,parp22,win,psno,ijk)   ! 081010
c210921  generate partonic initial state in relativistic  
c	 pA,Ap,AA,lp, & lA collision based on 'pythia'   ! 140414
c	it was composed by Ben-Hao Sa on 04/12/2003
c	the intermediate working arraies are in common statement 'sa2'
c110123 arraies in 'saf' is consistent with intermediate working arraies 
c       'saf' to 'pyjets' after call 'scat'   ! 220110 
c110123 output message is in 'pyjets' (partons) and 'sbh' (hadrons) for
c        case of mstpt1=0, but is in 'sbh' for mstpt1=1 and low energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter(kszj=80000)
        parameter(nsize=280000)
	double precision bst(4),bzp,bzt,bbb(3),bb(3)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c	those variables in above common blocks are defined in 'jetset'
        COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
	COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)   ! 221203
c	those variables in above common block are defined in 'pythia'
	COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
c080104
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c	 ,vnlep   ! 260314
        common/sa24/adj1(40),nnstop,non24,zstop   ! 140414
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3,
     c   parj21,parj4,adiv,gpmax,nnc   !   070417 010518
	common/sa30/vneump,vneumt,mstptj   ! 241110 100821 230722
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
c080104
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/count/isinel(600)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c	nap,nat,nzp,nzt,pio
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
	common/sa23/kpar,knn,kpp,knp,kep   ! 200601 060813
        common/sa33/smadel,ecce,secce,parecc,iparres   ! 270312 240412 131212
        common/schuds/schun,schudn,schudsn,sfra   ! She and Lei
c	iii : number of current event
c       csen: e+p total x section in fm^2
c	neve : total number of events 
c	bp : impact parameter
c       'sbe': store initial parton confiquration (with diquark) 
c       'saf': store parton configuration after parton rescattering 
c              (w/o diquark) 
c       c17(i,1-3) :  three position of i-th nucleon (origin is set at the center of
c       target nucleus)
c       tp(i) : time of i-th nucleon counted since collision of two nuclei
c       p17(i,1-4) : four momentum of i-th nucleon 
c       ishp(i)=1 if i-th particle inside the simulated volume
c              =0 if i-th particle outside the simulated volume
c	cspsn : total cross section of J/Psi (Psi') + n
c	cspsm : total cross section of J/Psi (Psi') + meson
c	iabsb = 0 : without J/Psi (Psi') + baryon
c	      = 1 : with J/Psi (Psi') + baryon
c	iabsm = 0 : without J/Psi (Psi') + meson
c	      = 1 : with J/Psi (Psi') + meson
c       ecsen: largest collision distance between lepton and p ! 060813
c060813 120214 note: incident lepton collides with nucleon in nucleus once 
c	 only, because of very low total x-section. that collision is the one with 
c	 lowest minimum approaching distance.
c       sig (fm^2): cross section of pion + pion to kaon + kaon
c       edipi: largest interaction distance between two pions.
c       epin: largest interaction distance between pion and nucleon.
c       ekn: largest interaction distance between kaon and nucleon.
c       ecsnn: largest interaction distance between two nucleons.
c       t0 : average proper formation time at rest.
c       ddt : time accuracy used in parton initiation 
c	time accuracy used in parton cascade is dddt 
c       rou0 : normal nucleon density.
c       rao : enlarged factor for the radius of simulated volume.
c       nap and nzp (nat and nzt) are the mass and charge numbers of 
c        projectile (target) nucleus
c       r0p=rnp     : the standard radius of projectile
c       r0t=rnt     : the standard radius of target
c       nctl: number of collision pairs in current collision list
c       nctl0: number of collision pairs in initial collision list
c180121 nctlm: maxmimum number of collision pairs 
c       noinel(1): statistics of nn elas. colli.;
c	noinel(592): statistics of nn colli. calling pythia
c230121 noinel(593): statistics of nn colli. not calling pythia

	dimension peo(4),pi(4),pj(4),xi(4),xj(4)
	dimension inoin(kszj)
        dimension lc(nsize,5),tc(nsize),tw(nsize)

	kpar=0
        knn=0
        kpp=0
        knp=0
        kep=0   ! 060813
c	kep: to statistics of the # of times calling pythia in 
c	 case of lepton is projectile

c270312 initiation of x,y,xy,x^2,y^2 and sump (statistics of the number of
c	 nucleons in overlap region)   ! 131212
        sumx=0.
        sumy=0.
        sumxy=0.   ! 131212
        sumx2=0.
        sumy2=0.
        sump=0.
        adj130=adj1(30)   ! 121222 Lei
c270312


c	initiates pp (pA,Ap,AA,lp & lA) collision system
c241110
c       creat the initial particle (nucleon) list

c230311 in position phase space
c191110
c       A+B (nucleus-nucleus)   ! 230311
        if(ipden.eq.1 .and. itden.eq.1)then   !! 230311
c       distribute projectile nucleons by Woods-Saxon   ! 060921              
	napt=nap
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0p
	am=suppm   ! upper bound in sampling the radius of projectile nucleon
	ac=suppc   ! maximum radius for projectile
        ratps=vneump/nap   ! ratio of projectile participant nucleons to total 
c	if(iii.eq.10)write(9,*)'nap,vneump,ratps=',nap,vneump,ratps
	do i1=1,nap
        if( ABS(adj130 - 1.).lt.1D-15 )then   ! 121222 Lei i.e. adj130=1, avoids unexpected bugs.
        rann=pyr(1)
        if(rann.lt.ratps)then
c       sample position of projectile nucleon in overlap region of colliding 
c         nuclei
        call arrove(i1,1,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c       alp,r0,am,ac)   ! 270312 131212 101014
        else
c	sample position of projectile nucleon according to Woods-Saxon
c	 distribution
        call woodsax_samp(i1,1,alp,r0,am,ac,1)   ! 230311
c230311 last argument here is 'iway', iway=1: particle i1 must be outside the 
C230311  overlap region of colliding nuclei, iway=0: no more requirement 
        endif
        elseif( ABS(adj130 - 0.).lt.1D-15 )then   ! 121222 Lei i.e. adj130=0
        call woodsax_samp(i1,1,alp,r0,am,ac,0)   ! 230311 060921
        endif
	enddo
c230311
c	distribute target nucleons by Woods-Saxon   ! 060921
	napt=nat
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0t
	am=suptm   ! upper bound in sampling the radius of target
	ac=suptc   ! maximum radius for target
        ratps=vneumt/nat   ! ratio of target participant nucleons to total 
	do i1=1,nat
	i2=i1+nap
        if( ABS(adj130 - 1.).lt.1D-15 )then   ! 121222 Lei i.e. adj130=1
        rann=pyr(1)
        if(rann.lt.ratps)then
c       sample position of target nucleon in overlap region of colliding nuclei
        call arrove(i2,0,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c       alp,r0,am,ac)   ! 270312 131212 101014
        else
c	sample position of target nucleon according to Woods-Saxon
c	 distribution
        call woodsax_samp(i2,0,alp,r0,am,ac,1)
        endif
        elseif( ABS(adj130 - 0.).lt.1D-15 )then   ! 121222 Lei i.e. adj130=0
        call woodsax_samp(i2,0,alp,r0,am,ac,0)   ! 060921
        endif
	enddo
c191110
        do i=1,nap
c050322 c17(i,1)=c17(i,1)+bp
        c17(i,1)=c17(i,1)+0.5*bp ! 050322 move x-component of origin to 0.5*bp
c1512222
c       Calculates eccentricity correctly for both adj(30)=0 and 1.   ! 151222 Lei
        x = c17(i,1)
        y = c17(i,2)
        z = c17(i,3)
c       Relative distance between the projectile nucleon i and the target center (-bp/2., 0, 0)
        rel_dist = SQRT( (x+bp/2.)**2 + y**2 + z**2 )
c       The projectile nucleon i is inside the target, i.e. inside the overlap region.
        if( rel_dist .le. r0t )then
            sumx  = sumx  + x
            sumy  = sumy  + y
            sumxy = sumxy + x*y
            sumx2 = sumx2 + x**2
            sumy2 = sumy2 + y**2
            sump  = sump  + 1.
        end if
c151222
        enddo
c050322
        do i=nap+1,nap+nat
        c17(i,1)=c17(i,1)-0.5*bp
c151222
c       Calculates eccentricity correctly for both adj(30)=0 and 1.   ! 151222 Lei
        x = c17(i,1)
        y = c17(i,2)
        z = c17(i,3)
c       Relative distance between the target nucleon i and the projectile center (+bp/2., 0, 0)
        rel_dist = SQRT( (x-bp/2.)**2 + y**2 + z**2 )
c       The target nucleon i is inside the projectile, i.e. inside the overlap region.
        if( rel_dist .le. r0p )then
            sumx  = sumx  + x
            sumy  = sumy  + y
            sumxy = sumxy + x*y
            sumx2 = sumx2 + x**2
            sumy2 = sumy2 + y**2
            sump  = sump  + 1.
        end if
c151222
        enddo
c050322        
c191110

c       p+A or lepton+A   ! 060813 120214
        elseif((ipden.eq.0.or.ipden.gt.1) .and. itden.eq.1)then !060813 120214
c100821 distribute projectile proton
        do i=1,3
        c17(1,i)=0.
        if(i.eq.1)c17(1,i)=c17(1,i)+0.5*bp   ! 050322 bp->0.5*bp
        enddo
c	distribute target nucleons by Woods-Saxon   ! 180921
c	distribute nat-vneumt target nucleons by Woods-Saxon   ! 100821
c100821 vneumt: # of target participant nucleons
c180921 ineumt=int(vneumt)   ! 100821
	napt=nat   ! -ineumt 180921
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0t
	am=suptm   ! upper bound in sampling the radius of target
	ac=suptc   ! maximum radius for target
	do i1=1,napt   ! 100821 nat->nat-ineumt=napt
	i2=i1+nap
	call woodsax_samp(i2,0,alp,r0,am,ac,0)
	enddo
c240513
c050322
        do i=nap+1,nap+nat
        c17(i,1)=c17(i,1)-0.5*bp
        enddo
c050322

c       A+p
        elseif(ipden.eq.1 .and. itden.eq.0)then   !!
c180921 distribute projectile nucleons by Woods-Saxon
c distribute nap-vneump projectile (spectator) nucleons by Woods-Saxon ! 100821
c100821 vneump: # of projectile participant nucleons 
c180921 ineump=int(vneump)    
	napt=nap   ! 180921 -ineump
	if(napt.lt.27)then
        alpt=0.47
        elseif(napt.gt.27.and.napt.lt.108)then
        alpt=0.488
	else
	alpt=0.54
	endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
        elseif(napt.eq.56)then
        alpt=0.49
        elseif(napt.eq.64)then
        alpt=0.49     
        elseif(napt.eq.108)then
        alpt=0.495
        elseif(napt.eq.184)then
        alpt=0.53
        elseif(napt.eq.197)then
        alpt=0.54
        elseif(napt.eq.207)then
        alpt=0.545
        elseif(napt.eq.238)then
        alpt=0.55
 	endif
 	alp=alpt
	r0=r0p
	am=suppm   ! upper bound in sampling the radius of projectile nucleon
	ac=suppc   ! maximum radius for projectile
	do i1=1,napt
	call woodsax_samp(i1,1,alp,r0,am,ac,0)
	enddo
c191110 100821
        do i=1,napt
        c17(i,1)=c17(i,1)+0.5*bp   ! 050322 bp->0.5*bp
        enddo
c191110 100821
c100821 move x-component of origin to 0.5*bp        
        do i=1,3
        c17(nap+1,i)=0.
        if(i.eq.1)c17(nap+1,i)=-0.5*bp   ! 0.->-0.5*bp
        enddo
c240513
c	p+p or lepton+p   ! 070417
	elseif((ipden.eq.0 .and. itden.eq.0) .or.
     c   (itden.eq.0 .and. ipden.ge.11))then   !! 070417
	do i=1,3
	c17(1,i)=0.
	c17(2,i)=0.
	enddo        
        endif   !!
c230311
	r0pt=r0p+r0t   ! 191110
c270312
        if(sump.ne.0.)then   !!!
        asumx=sumx/sump
        sigmx2=sumx2/sump-asumx*asumx
        asumy=sumy/sump
        sigmy2=sumy2/sump-asumy*asumy
	asumxy=sumxy/sump   ! 131212
        sigmxy=asumxy-asumx*asumy   ! 131212
        sigmsu=sigmy2+sigmx2   ! change from sigmxy to sigmsu 131212
        sigmde=sigmy2-sigmx2   ! 131212
        argu=sigmde*sigmde+4*sigmxy*sigmxy   ! 131212        
c131212
c       participant eccentricity of participant nucleons
        if(argu.gt.0. .and. sigmsu.gt.0.)
     c   ecce=sqrt(argu)/sigmsu !131212
c       calculate \epsilon{2}=\sqrt(<\epsilon_{part}^2>)
cc      ecce=ecce*ecce
c       note, \epsilon{2} should be \sqrt(aecceo), aecceo is a output
c        in paciae_21c.f
c       calculate transverse overlap area
        argu1=sigmx2*sigmy2-sigmxy*sigmxy
        if(argu1.gt.0.)secce=3.1416*sqrt(argu1) ! overlop area 250113
c131212
c       assuming ecce=geometric eccentricity of ellipsoid (\sqrt{(1-b^2/a^2)})
c        with half major axis b=pt*(1+smadel) and half minor axis
c        a=pt*(1-smadel), the resulted smadel=-ecce*ecce/4 (if neglecting
c        the samll term of ecce*ecce*(-2*smadel+smadel*smadel)
	ecc2=ecce*ecce   ! 250113
	smadel_a=parecc*ecc2/4. ! approximated deformation parameter 250113
c250113
        delta1=(2.-ecc2+2.*(1.-ecc2)**0.5)/ecc2
        delta2=(2.-ecc2-2.*(1.-ecc2)**0.5)/ecc2
        if(delta1.le.1.)then
        smadel=parecc*delta1  ! exact deformation parameter
        elseif(delta2.le.1.)then  
        smadel=parecc*delta2  ! exact deformation parameter
        else   
        endif   
c       write(9,*)'ecce,smadel_a,smadel=',ecce,smadel_a,smadel
c250113
c       here a sign change is introduced because of asymmetry of initial
c        spatial space is oppsed to the final momentum space
c       write(9,*)'vneump,vneumt,sump,ecce,smadel=',
c     c  vneump,vneumt,sump,ecce,smadel
        endif   !!!
c270312
c021018	note: psno=0 (bp=0) for pp,lp and lA
c	if(iii.eq.10)then
c	write(9,*)'after woodnat,bp=',bp
c	do i=1,nat+nap
c	write(9,*)c17(i,1),c17(i,2),c17(i,3)
c	enddo
c	endif
c191110
c	the beam direction is identified as the z axis
c	the origin in position space is set on the center of
c	target nucleus and the origin of time is set at the moment of 
c	first nn colission assumed to be 1.e-5
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	continue

c230311 in momentum phase space
	if(ifram.eq.1)then
	ep1=0.5*win   ! energy of projetile particle (if it is proton)
	et1=ep1   ! energy of target particle (if proton)
	ep2=0.5*win   ! energy of projetile particle (if neutron)
	et2=ep2   ! energy of target particle (if neutron)
	pm2=pmas(pycomp(2212),1)**2   ! square mass of proton	
	pp1=dsqrt(ep1*ep1-pm2)   ! momentum of projetile particle (if proton)
	pt1=-dsqrt(et1*et1-pm2)  ! momentum of target particle (if proton)
	pm2=pmas(pycomp(2112),1)**2   ! square mass of nucleon	
	pp2=dsqrt(ep2*ep2-pm2)   ! momentum of projetile particle (if neutron)
	pt2=-dsqrt(et2*et2-pm2)  ! momentum of target particle (if neutron)
c260314 set four momentum and mass for incident lepton 
	if(ipden.ge.11.and.ipden.le.16)then   ! in cms
	pincl(1)=0.
	pincl(2)=0.
	pincl(4)=0.5d0*win
	pincl(5)=pmas(pycomp(ipden),1)
	pincl3=pincl(4)*pincl(4)-pincl(5)*pincl(5)
	pincl3=dmax1(pincl3,1.d-20)
	pincl(3)=dsqrt(pincl3)
	pinch(1)=0.
	pinch(2)=0.
	pinch(4)=0.5d0*win
	pinch(5)=pmas(pycomp(2212),1)
	pinch3=pinch(4)*pinch(4)-pinch(5)*pinch(5)
	pinch3=dmax1(pinch3,1.d-20)
	pinch(3)=dsqrt(pinch3)
	endif
c260314
	endif	
	if(ifram.eq.0)then
	pp1=win   ! momentum of projetile particle (if proton)
	pt1=1.e-20   ! momentum of target particle (if proton)
	pp2=win   ! momentum of projetile particle (if neutron)
	pt2=1.e-20   ! momentum of target particle (if neutron)
	pm2=pmas(pycomp(2212),1)**2   ! square mass of proton	
	ep1=dsqrt(pp1*pp1+pm2)   ! energy of projetile particle (if proton)
	et1=dsqrt(pt1*pt1+pm2)   ! energy of target particle (if proton)
	pm2=pmas(pycomp(2112),1)**2   ! square mass of neutron
	ep2=dsqrt(pp2*pp2+pm2)   ! energy of projetile particle (if neutron)
	et2=dsqrt(pt2*pt2+pm2)   ! energy of target particle (if neutron)
c260314 set four momentum and mass for incident lepton 
	if(ipden.ge.11.and.ipden.le.16)then   ! in lab
	pincl(1)=0.
	pincl(2)=0.
	pincl(3)=win
	pincl(5)=pmas(pycomp(ipden),1)
	pincl4=pincl(3)*pincl(3)+pincl(5)*pincl(5)
	pincl4=dmax1(pincl4,1.d-20)
	pincl(4)=dsqrt(pincl4)
	pinch(1)=0.
	pinch(2)=0.
	pinch(3)=0.
	pinch(5)=pmas(pycomp(2212),1)
	pinch(4)=pinch(5)
	endif
c260314
	endif
c260314
c	if(ipden.ge.11.and.ipden.le.16)then
c	write(9,*)'pincl=',(pincl(i),i=1,5)
c	write(9,*)'pinch=',(pinch(i),i=1,5)
c	endif
c260314	
100	inzp=iabs(nzp)
        inzt=iabs(nzt)
	do i=1,nap
	p17(i,1)=0.   ! four momenta of projectile particle i
	p17(i,2)=0.
	if(i.le.inzp)then  
	p17(i,3)=pp1   ! projectile particle is proton
	p17(i,4)=ep1
	else
	p17(i,3)=pp2   ! projectile particle is neutron
	p17(i,4)=ep2
	endif
	enddo
	napt=nap+nat
	do i=nap+1,napt
	p17(i,1)=0.   ! four momenta of target particle i
	p17(i,2)=0.
	if(i.le.nap+inzt)then   	
	p17(i,3)=pt1   ! target particle is proton
	p17(i,4)=et1
	else
	p17(i,3)=pt2   ! target particle is neutron
	p17(i,4)=et2
	endif	
	enddo

	do i=1,napt
	tp(i)=0.
	enddo

c	calculate the velocity of the CM of collision system in LAB or 
c        in nucleon-nucleon CM system
	bst(1)=p17(1,1)*nap+p17(nap+1,1)*nat
	bst(2)=p17(1,2)*nap+p17(nap+1,2)*nat
	bst(3)=p17(1,3)*nap+p17(nap+1,3)*nat
	bst(4)=p17(1,4)*nap+p17(nap+1,4)*nat
	bst(1)=-bst(1)/bst(4)
	bst(2)=-bst(2)/bst(4)
	bst(3)=-bst(3)/bst(4)

	n=0
        nbe=0   ! 080104
        naf=0   ! 080104
        nsa=0
	idi=0
	idio=0
        do i1=1,kszj
        do j1=1,5
	k(i1,j1)=0
	p(i1,j1)=0.
	v(i1,j1)=0.
        kbe(i1,j1)=0
        pbe(i1,j1)=0.
        vbe(i1,j1)=0.
        kaf(i1,j1)=0
        paf(i1,j1)=0.
        vaf(i1,j1)=0.
        ksa(i1,j1)=0
        psa(i1,j1)=0.
        vsa(i1,j1)=0.
        enddo
        ndiq(i1)=0
        npt(i1)=0
        ifcom(i1)=0   ! 220110
        ishp(i1)=0
        tau(i1)=0.
        enddo

        nctl=0
        do i=1,nsize
        do j=1,5
        lc(i,j)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo

        do i=1,100
        numb(i)=0
        numbs(i)=0
        enddo

c      '1 -> |nzp|' are projectile protons or lepton, '|nzp|+1 -> nap' 
c        are projectile neutrons; 'nap+1 -> nap+nzt' are targer protons, 
c        the rest are target nuctrons in 'pyjets' after nuclear initiation above
	n=napt
	do i=1,n
	k(i,1)=1
	k(i,2)=2112
	p(i,5)=pmas(pycomp(2112),1)
	if((i.le.iabs(nzp).and.ipden.lt.2).or.(i.gt.nap .and. i.le.nap+
     c	 nzt))then   ! 060813 120214
	k(i,2)=2212
	p(i,5)=pmas(pycomp(2212),1)
c060813 120214
        elseif(i.le.nap.and.(ipden.eq.11.and.nzp.eq.-1))then
        k(i,2)=11
        p(i,5)=pmas(pycomp(11),1)
        elseif(i.le.nap.and.(ipden.eq.11.and.nzp.eq.1))then
        k(i,2)=-11
        p(i,5)=pmas(pycomp(11),1)
        elseif(i.le.nap.and.(ipden.eq.12.and.nzp.eq.-1))then
        k(i,2)=12
        p(i,5)=pmas(pycomp(12),1)
        elseif(i.le.nap.and.(ipden.eq.12.and.nzp.eq.1))then
        k(i,2)=-12
        p(i,5)=pmas(pycomp(12),1)
        elseif(i.le.nap.and.(ipden.eq.13.and.nzp.eq.-1))then
        k(i,2)=13
        p(i,5)=pmas(pycomp(13),1)
        elseif(i.le.nap.and.(ipden.eq.13.and.nzp.eq.1))then
        k(i,2)=-13
        p(i,5)=pmas(pycomp(13),1)
	elseif(i.le.nap.and.(ipden.eq.14.and.nzp.eq.-1))then
        k(i,2)=14
        p(i,5)=pmas(pycomp(14),1)
        elseif(i.le.nap.and.(ipden.eq.14.and.nzp.eq.1))then
        k(i,2)=-14
        p(i,5)=pmas(pycomp(14),1)
        elseif(i.le.nap.and.(ipden.eq.15.and.nzp.eq.-1))then
        k(i,2)=15
        p(i,5)=pmas(pycomp(15),1)
        elseif(i.le.nap.and.(ipden.eq.15.and.nzp.eq.1))then
        k(i,2)=-15
        p(i,5)=pmas(pycomp(15),1)
        elseif(i.le.nap.and.(ipden.eq.16.and.nzp.eq.-1))then
        k(i,2)=16
        p(i,5)=pmas(pycomp(16),1)
        elseif(i.le.nap.and.(ipden.eq.16.and.nzp.eq.1))then
        k(i,2)=-16
        p(i,5)=pmas(pycomp(16),1)
        else
c060813 120214
	endif
	do j=1,3
	p(i,j)=p17(i,j)
	v(i,j)=c17(i,j)
	enddo
	p(i,4)=p17(i,4)
	v(i,4)=tp(i)
	enddo
500	continue    ! 031103
c	v, vbh and vsa arraies are the position four vector
c	note: for v etc., we do not take care of their fifth component
c	 for array k, we take care of only first three components
c	write(9,*)'after initializing iii=',iii   !s
c	call psum(p,1,n,peo)   !!
c	write(9,*)'psum, after initializing nucleus-nucleus collision'   !!
c	write(9,*)peo   !!

c	boost PYJETS into cms of initial nucleus-nucleus collision system 
c	 from lab or initial nucleon-nucleon cms system.
c	call pyrobo(1,n,0.0,0.0,bst(1),bst(2),bst(3))
c	Lorentz contract
	bzp3=0.
	bzp4=0.
	bzt3=0.
	bzt4=0.
	do i=1,nap
	bzp3=bzp3+p(i,3)
	bzp4=bzp4+p(i,4)
	enddo
	do i=nap+1,napt
	bzt3=bzt3+p(i,3)
	bzt4=bzt4+p(i,4)
	enddo
	bzp=bzp3/bzp4
	bzt=bzt3/bzt4
	gamp=1.d0/dsqrt(dmax1(1.d-20,(1.0d0-bzp*bzp)))
c060813 120214 no Lowrantz contraction for incident lepton
        if(ipden.ge.2)gamp=1.   ! 060813 120214
	gamt=1.d0/dsqrt(dmax1(1.d-20,(1.0d0-bzt*bzt)))
c	try no lorentz contract for target
c	gamt=1.
	do i=1,nap
	c17(i,3)=c17(i,3)/gamp
	v(i,3)=v(i,3)/gamp	
	enddo
	do i=nap+1,napt
	c17(i,3)=c17(i,3)/gamt
	v(i,3)=v(i,3)/gamt
	enddo
c	filter out those kind of particles wanted to study and make 
c	 the order of proton, neutron, ... (cf. 'filt')
	call filt
c060813 120214
c       since lepton was moved to last position after calling filt, one has to
c        remove it to the fist position
        if(ipden.ge.2)call ltof(n)  
c060813 120214
c161021 'pyjets' to 'sa2' 
	nsa=n
        do i=1,n
        do j=1,5
	ksa(i,j)=k(i,j)
        psa(i,j)=p(i,j)
        vsa(i,j)=v(i,j)
        enddo
	ishp(i)=1
	enddo
        do m=1,kfmax
        numb(m)=numbs(m)
        enddo
c	note: particle list is composed of the arraies in common block 
c	 'sa2', the array 'ishp' in common block 'wz', the array 'tau' in 
c        common block 'sa4', and the array 'numb' in common block 'sa5'
	time=time_ini   ! 081010
	irecon=0

        if(nosc.eq.2 .OR. nosc.eq.3) call oscar(win,0)   !Lei20230214 full_events_history of OSC1999A

c       calculate the position for the center of mass of the
c	non-freeze-out system. The distance of a particle, when checking
c	is it freezing out or not, is measured with respect to this center
	call copl(time)

c       creat the initial collision list, note: be sure that the initial  
c	collision list must not be empty
	call ctlcre(lc,tc,tw)

c	find out colli. pair with least colli. time
	call find(icp,tcp,lc,tc,tw,0)
	if(icp.eq.0)stop 'initial collision list is empty'   !
	time=tcp

c070417 perform classical Newton motion in Lab. system for all particles 
	call his(time,lc,tc,tw,istop)
	do ij=1,nsa
	vsa(ij,4)=0.
	enddo

c070417 move origin of time to collision time of first nucleon-nucleon collision
	do ij=1,nctl
	tc(ij)=tc(ij)-time+1.e-5
	enddo
	time=time_ini   ! 081010
	call copl(time) 
400	continue

c       administrate a nucleus-nucleus collision   ! 180520
        call scat(time,lc,tc,tw,win,parp21,parp22,psno,ijk,ipau,irecon,
     c   gamt)   ! 021207
        if(ijk.eq.1)return   
        time_ini=time   ! 081010
c	write(9,*)'af scat, iii,time,time_ini=',iii,time,time_ini   !
c281121 write(9,*)'af. scat iii,nctl=',iii,nctl0,nctlm,nctl   ! 100821

800	continue
c       'saf' to 'pyjets'
c180520	if(adj1(40).ne.5)call tran_saf   ! 140414 
        if(mstptj.eq.0)call tran_saf   ! 140414 180520 230722
c241110
c       if(iii.eq.5)then
c       write(22,*)'af scat'
c       call pylist(1)
c       endif
c241110
        n00=n   ! 220110
c220110 n00: 'largest line number' in 'pyjets'
c220110 partons above n00 appear after inelastic collision  
c       'sa2' to 'sbh'
        nbh=0
        if(nsa.ge.1)then
        nbh=nsa
        do i1=1,nsa
        do i2=1,5
        kbh(i1,i2)=ksa(i1,i2)
        pbh(i1,i2)=psa(i1,i2)
        vbh(i1,i2)=vsa(i1,i2)
        enddo
        enddo
        endif
	do i1=nbh+1,kszj
        do i2=1,5
        kbh(i1,i2)=0
        pbh(i1,i2)=0.
        vbh(i1,i2)=0.
        enddo
        enddo
c 	P(N,5)=SQRT(MAX(-P(N,1)**2-P(N,2)**2-P(N,3)**2+P(N,4)**2,0.0))
c 	P(N-1,5)=SQRT(MAX(-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
c     &	+P(N-1,4)**2,0.0))
c	call pyboro(1,n,0.0,0.0,-bst(1),-bst(2),-bst(3))
c	boost PYJETS back to lab or nucleon-nucleon cms system.

	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sysini(win)   ! 060813
c       give the initial values to quantities needed in calculation
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP  
        PARAMETER (KSZJ=80000)
        COMMON/PYCIDAT1/KFACOT(100),DISDET(100),ISINELT(600)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/count/isinel(600)
        COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
        common/sa6/kfmaxi,nwhole
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c   ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   ! 221203 250204

        anat=nat
        anap=nap
        param(1)=para1_1   ! 250204 200504
c       write(9,*)'para1_1=',param(1)   ! 250204 200504
c       rou0=PARAM(11)
c       considering the nucleus as a sphere with radii rnt for target
c        and rnp for projectile.
c       rnt=(3.*anat/(4.*3.1415926*rou0))**(0.33333)
c       rnp=(3.*anap/(4.*3.1415926*rou0))**(0.33333)
        rp00=1.12   ! 1.05 to 1.12 070613
        rt00=1.12   ! 1.05 to 1.12 070613
c070613 if(nap.gt.16)rp00=1.16*(1-1.16*anap**(-0.666666)) !rp00=1.122 (nat=208)
c070613 if(nat.gt.16)rt00=1.16*(1-1.16*anat**(-0.666666)) ! rt00=1.12 (nat=197)
        if(itden.eq.0)rnt=rt00*anat**(0.33333)   ! 310805
        if(itden.eq.1)rnt=rt00*anat**(0.33333)   ! +0.54  160511
        if(nat.eq.2 .and. nzt.eq.1)rnt=4.0   ! 2.60 2.095  1.54 2603141
c060813 120214 if(ipden.eq.2)rnt=0.5   ! lepton
        if(ipden.eq.0)rnp=rp00*anap**(0.33333)   ! 310805
        if(ipden.eq.1)rnp=rp00*anap**(0.33333)   ! +0.54  160511
        if(ipden.ge.2)rnp=0.5   ! lepton   ! 060813 120214
        if(nap.eq.2 .and. nzp.eq.1)rnp=4.0   ! 2.60 2.095  1.54
        rou0=3./4./3.1416*anat/(rnt*rnt*rnt)   ! 310805
        r0p=rnp
        r0t=rnt
C       set initial values to some quantities
c       in the program the x-sections are given in a unit of fm^2   ! 060813
        csnn=PARAM(1)*0.1
        cspin=PARAM(2)*0.1
        cskn=PARAM(3)*0.1
        cspipi=PARAM(4)*0.1
        cspsn=PARAM(13)*0.1
        cspsm=PARAM(14)*0.1
        csspn=PARAM(15)*0.1
        csspm=PARAM(16)*0.1
c060813 120214
        if(ipden.ge.2)then
        if(ifram.eq.0)then
        ept=sqrt(win*win+0.938*0.938)
        rots=sqrt((ept+0.938)*(ept+0.938)-win*win)
        endif
        if(ifram.eq.1)rots=win
        call crosep(rots,csen)   ! temporary using e^-p total x-section  
c       if(nzp.lt.0)call crosep(rots,csen)   ! e^-p total x-section
c       if(nzp.ge.0)call crosepp(rots,csen)   ! e^+p total x-section  
        csen=csen*0.1
        endif   
c060813 120214
c       largest collision distance between two colliding particles.
        edipi=dsqrt(cspipi/3.1416)
        epin=dsqrt(cspin/3.1416)
        ekn=dsqrt(cskn/3.1416)
        ecsnn=dsqrt(csnn/3.1416)
        ecspsn=dsqrt(cspsn/3.1416)
        ecspsm=dsqrt(cspsm/3.1416)
        ecsspn=dsqrt(csspn/3.1416)
        ecsspm=dsqrt(csspm/3.1416)
        ecsen=sqrt(csen/3.1416)   ! 060813

        anp=nap**.3333
        ant=nat**.3333
        do ia=1,2
        if(ia.eq.1)napt=nap
        if(ia.eq.2)napt=nat
        if(napt.lt.27)then
        alpt=0.47
        elseif(napt.ge.27.and.napt.le.108)then
        alpt=0.488
        else
        alpt=0.54
        endif
        if(napt.eq.27)then
        alpt=0.478
        elseif(napt.eq.28)then
        alpt=0.48
        elseif(napt.eq.32)then
        alpt=0.49
         elseif(napt.eq.56)then
        alpt=0.49
          elseif(napt.eq.64)then
        alpt=0.49     
           elseif(napt.eq.108)then
        alpt=0.495
           elseif(napt.eq.184)then
        alpt=0.53
           elseif(napt.eq.197)then
        alpt=0.54
           elseif(napt.eq.207)then
        alpt=0.545
            elseif(napt.eq.238)then
        alpt=0.55
        endif
        if(ia.eq.1)alp=alpt
        if(ia.eq.2)alt=alpt
        enddo
        suppc=rp00*anp+2.*alp
        suptc=rt00*ant+2.*alt
        suppm=1.d0/(1+dexp(0.d0-r0p/alp))
        suptm=1.d0/(1+dexp(0.d0-r0t/alt))

        sig=PARAM(5)*0.1
        rcsit=PARAM(6)
        t0=PARAM(7)   ! 230805
c230805t0=0.   ! 221102   
c221102 proper formation time of particle from 'pythia'
        dep=PARAM(9)
        ddt=PARAM(8)
        rao=PARAM(10)
        kfmax=KFMAXT
        do i=1,100
        kfaco(i)=KFACOT(i)
        enddo
        do j=1,600
        isinel(j)=ISINELT(j)
        enddo
        do i=1,100
        do j=1,100
        disbe(i,j)=0.
        enddo
        enddo
        do j=1,kfmax
c       something might be missing here ?
        disbe(1,j)=DISDET(j)
        disbe(2,j)=DISDET(j)
        disbe(3,j)=DISDET(j)
        disbe(4,j)=DISDET(j)
c       disbe(26,j)=DISDET(j)
c       disbe(27,j)=DISDET(j)
c       disbe(28,j)=DISDET(j)
c       disbe(29,j)=DISDET(j)
        enddo
400     do i=1,99
        do j=i+1,100
        disbe(j,i)=disbe(i,j)
        enddo
        enddo
c        kfmaxi=kfmax
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine arrove(ii,jj,sumx,sumy,sumxy,sumx2,sumy2,sump,
     c	 alp,r0,am,ac)   ! 101014)   
c	 ! 191110 270312 131212
c	arrange randomly particle ii in overlap region of colliding nuclei 
c	jj=0 and 1 for target and projectile, respectively  
	PARAMETER (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc 
	common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	b=bp
        iii=0
54      iii=iii+1
        if(iii.eq.10000)then
        write(9,*)'difficult to arrange produced nucleon in'
        write(9,*)'subroutine arrove,infinitive loop may occur'
	goto 55   ! set larget number of try is equal to 10000
        endif
c       sample a point in the unit sphere 
c101014        x=1.-2.*pyr(1)
c        y=1.-2.*pyr(1)
c        z=1.-2.*pyr(1)
c        rr=x*x+y*y+z*z
c101014        if(rr.gt.1) goto 54
	if(jj.eq.0)then   ! ii in target (-b/2.)
c101014	x=x*r0t
c	y=y*r0t
c101014	z=z*r0t
c101014
c       sample a point according to woodsax distribution
        call woodsax_samp(ii,jj,alp,r0,am,ac,0)
        x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)
c101014
c       relative to projectile center, they are b-x, y, and z, respectively 
c	adjudge does (x-b,y,z) is in the sphere of projectile
        r1=sqrt((b-x)*(b-x)+y*y+z*z)
        if(r1.gt.r0p)goto 54
c101014        c17(ii,1)=x
c        c17(ii,2)=y
c101014        c17(ii,3)=z
c270312
c151222 They are done in parini now for both adj(30)=0 and 1.   ! 151222 Lei
c151222        sumx=sumx+x
c151222        sumy=sumy+y
c151222        sumxy=sumxy+x*y   ! 131212
c151222        sumx2=sumx2+x*x
c151222        sumy2=sumy2+y*y
c151222        sump=sump+1.
c270312
	endif
	if(jj.eq.1)then   ! ii in projectile (+b/2.)
c101014	x=x*r0p
c	y=y*r0p
c101014	z=z*r0p
c101014
c       sample a point according to woodsax distribution
        call woodsax_samp(ii,jj,alp,r0,am,ac,0)
        x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)
c101014
c       relative to target center, they are x+b, y, and z, respectively 
c	adjudge does (x+b,y,z) is in the sphere of target
        r1=sqrt((x+b)*(x+b)+y*y+z*z)
        if(r1.gt.r0t)goto 54
c101014	c17(ii,1)=x 
c	 c17(ii,2)=y
c101014        c17(ii,3)=z
c270312
c151222 They are done in parini now for both adj(30)=0 and 1.   ! 151222 Lei
c151222	xb=x+b   ! 101014 chen
c151222        sumx=sumx+xb   ! 101014
c151222        sumy=sumy+y
c151222        sumxy=sumxy+xb*y   ! 131212 101014
c151222        sumx2=sumx2+xb*xb   ! 101014
c151222        sumy2=sumy2+y*y
c151222        sump=sump+1.
c270312
	endif
55	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine woodsax_samp(ii,jj,alp,r0,am,ac,iway)   ! 191110 230311
c       sample position of nucleon ii in nucleus according to
c        Woods-Saxon distribution 
c       jj=0 and 1 for target and projectile, respectively
c	alp: diffusion length
c	r0: radius of nucleus
c	am: upper bound in sampling the radius
c	ac: maximum radius 
c230311 iway=1: ii must be outside overlap region of colliding nuclei
c230311 iway=0: no more requirement
	PARAMETER (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc  
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	b=bp 
c       selecting sample method
c	Woods-Sax distribution
	iii=0
100	iii=iii+1
        if(iii.eq.100000)then
        write(9,*)'difficult to arrange produced nucleon in'
        write(9,*)'subroutine woodsax,infinitive loop may happen'
	goto 200   ! set larget number of try is equal to 100000
        endif
	a1=pyr(1)
	xf=ac*(a1)**(1./3.)
	b1=pyr(1)
	deno2=1.+exp((xf-r0)/alp)
c       if(deno2.eq.0.)deno2=1.e-10
	yf=1./deno2
c	Gaussan distribution
cc	yf=exp(-xf*xf/2./r0)
	if(b1.gt.yf/am) goto 100
	call samp(xf,ii)
c       subroutine 'samp' is sampling the direction according to isotropic
c        distribution
        x=c17(ii,1)
        y=c17(ii,2)
        z=c17(ii,3)
        if(iway.eq.0)goto 200   ! 230311
c	ii must be outside overlap region of colliding nuclei
        if(jj.eq.0)then   ! ii in target (-b/2.)
c       relative to projectile center, above x, y, and z are b-x, y, and z, 
c	 respectively
c       (b-x,y,z) is inside or not inside the sphere of projectile
        r1=sqrt((b-x)*(b-x)+y*y+z*z)
        if(r1.lt.r0p)goto 100
c        c17(ii,1)=x
c        c17(ii,2)=y
c        c17(ii,3)=z
        endif
        if(jj.eq.1)then   ! ii in projectile (+b/2.)
c       relative to target center, they are x+b, y, and z, respectively
c       (x+b,y,z) is inside or not inside the sphere of target
        r1=sqrt((x+b)*(x+b)+y*y+z*z)
        if(r1.lt.r0t)goto 100
c        c17(ii,1)=x
c        c17(ii,2)=y
c        c17(ii,3)=z
	endif
200	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine samp(xf,i)
c       arrange i-th particle on the surface of sphere with radius xf
c100821 sampling on the surface of a sphere with radius xf
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
	PARAMETER (kszj=80000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	cita=2*pyr(1)-1.
	fi=2.*pio*pyr(1)
	sita=dsqrt(1.-cita**2)
	c17(i,1)=xf*sita*dcos(fi)
	c17(i,2)=xf*sita*dsin(fi)
	c17(i,3)=xf*cita
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sampi(xf,i)   ! 100821
c       sampling in a sphere with radius xf
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=80000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
100     continue        
        x=2.*pyr(1)-1.
        y=2.*pyr(1)-1.
        z=2.*pyr(1)-1.
        rr=x*x+y*y+z*z
        if(rr.gt.1)goto 100
	c17(i,1)=xf*x
	c17(i,2)=xf*y
	c17(i,3)=xf*z
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine filt
c       filter out particles wanted to study and make them in  
c       the order of proton,neutron,pba,nba,pi+,pi-,pi0,k-,k0-,sigma0,
c       sigma-,sigma+,sigma0ba,sigma-ba,sigma+ba,lamda,lamdaba,k0,k+,
c       cascade-,cascade-ba,cascade0,cascade0ba,omega-,omega+,Delta-,
c       Delta0,Delta+,Delta++,rho+,rho-,rho0,J/Psi,Psi',x0c,x1c,x2c,
c	D,Dba,D0,D0ba,lamdac+,sigmac0,sigmac+,sigmac++,omega,k*+,K*0,
c	D*,D*ba,D*0,D*0ba (52 kind of particle altogether)
c060813 120214 in case of lepton+A, one images lepton as a initial 
c	 projectile proton
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
      PARAMETER (kszj=80000)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
	common/sa6/kfmaxi,nwhole
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	COMMON/PYCIDAT2/KFMAXT,NONCI2,PARAM(20),WEIGH(600)
	SAVE/PYCIDAT2/
	iii=0
	jjj=0
	do i=1,kfmax
	kf=kfaco(i)
	do j=iii+1,n
	call ord(jjj,j,kf)
	enddo
	iii=jjj
	numbs(i)=jjj
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ord(ipi,j,kf)
c	order particles according to flavor code
c	j: the particle needed to order
c	ipi: j-th particle should order after ipi
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
	parameter(kszj=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	dimension pp(5),vv(5),kk(5)
	ik=k(j,2)
	if(ik.eq.kf)then
	ipi=ipi+1
	do jj=1,5
	kk(jj)=k(ipi,jj)
	pp(jj)=p(ipi,jj)
	vv(jj)=v(ipi,jj)
	enddo
	do jj=1,5
	k(ipi,jj)=k(j,jj)
	p(ipi,jj)=p(j,jj)
	v(ipi,jj)=v(j,jj)
	enddo
	do jj=1,5
	k(j,jj)=kk(jj)
	p(j,jj)=pp(jj)
	v(j,jj)=vv(jj)
	enddo
	endif
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine prt_pyj(nn,cc)
c	print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/pyjets/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(mstu(11),*)i,ksa(i,2),(psa(i,j),j=1,4)
c       enddo
	call psum(psa,1,nsa,peo)
	ich1=0.
	do i1=1,nn
	kf=ksa(i1,2)
	ich1=ich1+pychge(kf)
	enddo
        cc=ich1/3.
        write(22,*)'pyj nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   ! 
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbh(nn,cc)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(mstu(11),*)i,kbh(i,2),(pbh(i,j),j=1,4)
c	write(9,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sbh nn=',nn
        write(mstu(11),*)'c & p sum=',cc,peo   !
c	write(9,*)peo,ich1/3   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sa2(nn,cc)
c       print particle list and sum of momentum and energy
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sa2 nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_sbe(nn,cc)   ! 220110
c       print particle list and sum of momentum and energy
        parameter (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sbe/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'sbe nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
c       do i=1,nn
c       write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prt_saf(nn,cc)   ! 220110
c       print particle list and sum of momentum and energy
        parameter (kszj=80000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/saf/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension peo(4)
c       do i=1,nn
c       write(22,*)i,kbh(i,2),(pbh(i,j),j=1,4)
c       enddo
        call psum(pbh,1,nbh,peo)
        ich1=0.
        do i1=1,nn
        kf=kbh(i1,2)
        ich1=ich1+pychge(kf)
        enddo
        cc=ich1/3.
        write(22,*)'saf nn=',nn
        write(22,*)'c & p sum=',cc,peo   !
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine psum(pei,il,ih,peo)
c       calculate sum of momentum and energy
c       pei: two dimension array of input momentum and energy
c       il and ih: lower and upper limits of sum
c       peo : one dimension array of output momentum and energy  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
        parameter (kszj=80000)
        dimension pei(kszj,5),peo(4)
        do i=1,4
        peo(i)=0.
        enddo
        do i=il,ih
        do j=1,4
        peo(j)=peo(j)+pei(i,j)
        enddo
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine scat(time,lc,tc,tw,win,parp21,parp22,psno,ijk,
     c	 ipau,irecon,gamt)   ! 021207
c	administrate a nucleus-nucleus collision !060813 120214
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP  
      PARAMETER (kszj=80000)
        parameter(nsize=280000)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa7/ispmax,isdmax,iflmax,ispkf(20),non7,asd(5),
     c   afl(20,5,2)
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23 
	common/sa15/nps,npsi,pps(5000,5),ppsi(5000,5)
 	common/sa16/dtt,dni(10),dpi(10),edi(10),bmin,bmax
     &   ,bar(10),abar(10),barf(10),abarf(10)   
     &   ,emin(10),eminf(10),eplu(10),epluf(10)   
        common/sa18/tdh,itnum,non18,cptl,cptu,cptl2,cptu2,snum(4,20),
     &	 v1(4,20),v2(4,20),v12(4,20),v22(4,20)
	common/sa21/pincl(5),pscal(5),pinch(5),vnu,fq2,w2l,yyl,zl,xb,pph
     c	 ,vnlep   ! 260314
	common/sa23/kpar,knn,kpp,knp,kep   ! 060813   
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa25/mstj1_1,mstj1_2,para1_1,para1_2   
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 220110
        common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3,
     c   parj21,parj4,adiv,gpmax,nnc   !   070417 010518
        common/sa28/nstr,nstra(kszj),nstrv(kszj),nstr0,
     c   nstr1,nstr1a(kszj),nstr1v(kszj)   ! 030620   
        common/sa30/vneump,vneumt,mstptj   ! 230722
        common/sa34/itorw,iikk,cp0,cr0,kkii   ! 060617 010418 010518 040920
        common/sbe/nbe,nonbe,kbe(kszj,5),pbe(kszj,5),vbe(kszj,5)
	common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
	common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5) 
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121
        common/sppb/nppb,non3,kppb(1000,5),pppb(1000,5),vppb(1000,5) ! 281121
        dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension pi(4),pj(4),pii(4),pjj(4),peo(4),pint(4)	
	dimension nni(10),ndi(10),npi(10)
	dimension pkk(kszj,4)   
	dimension cc(5),b(3),bkk(3),pl(100,5)   ! 260314
c	arraies in 'pyjets' are given after calling 'pythia' 
c	arraies in 'sa2' are used in the collision processes 
c       arraies in 'sbh' are used to store hadron after calling 'pythia'
c	numbs(i) is is given in 'filt', updated with transport processes, and 
c        numbs(i)->numb(i) in the initiation of nucleus-nucleus collisin only
c        numb(i) is updated with transport processes
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	0. is the hard distance between two pions
c	0.5 is the hard distance between two nucleons
c	0. is the hard distance between pion and nucleon
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       lc(i,1) and lc(i,2) are the line # of colliding particles 1 and 2 of i-th 
c        collision pair in particle list, respectively
c       lc(i,3) and lc(i,4) are the flavor codes of scattered particles 3 and 4
c        of i-th collision, respectively
c       lc(i,5) identifies the different inelastic processes,
c       lc(i,5)=592 refers to the process calling 'pythia'
c       tc(i) is the collision time of i-th colli.
c       tw(i) is the cross section ratio of (i-th inelas.)/tot
c	array 'sbe' stores cumulatively parton (q,qq,g and their anti-particle) 
c	 configuration before breaking the diquarks 
c	array 'saf' stores cumulatively parton (q,g and their anti-particle) 
c	 configuration after breaking the diquarks
c	idi: counts cumunatively the number of diquark (anti-diquark)
c	idio: value of idi after last nn collision
c       ndiq(j): = 0 if j is quark (antiquark) 
c                = idi if j is diquark (anti-diquark) 
c       note: j is line number in 'sbe' ('saf')
c220110 ifcom(idi): line number of first component of idi-th diquark
c       npt(idi): line number of second component of idi-th diquark 
c        (anti-diquark) in 'sbe' ('saf')
c       nstr: statitics of number of strings in a nucleus-nucleus collision 
c        when fragmentation string-by-string
c       nstra(i): line number of first component of i-th string
c       nstrv(i): line number of last component of i-th string
c220110 nstr0: number of strings after call break
        adj140=adj1(40)   ! 180520
        ijk=0
	do i=1,10
	nni(i)=0
	ndi(i)=0
	npi(i)=0
	dni(i)=0.
	dpi(i)=0.
	edi(i)=0.
c033101
	bar(i)=0.
	abar(i)=0.
	barf(i)=0.
	abarf(i)=0.
	emin(i)=0.
	eplu(i)=0.
	eminf(i)=0.
	epluf(i)=0.
c033101
	enddo

csa**********************************************************
        do i1=1,20
        do i2=1,4
        snum(i2,i1)=0.
        v1(i2,i1)=0.
        v2(i2,i1)=0.
        v12(i2,i1)=0.
        v22(i2,i1)=0.
        enddo
        enddo
c100821
	m1=numb(1)
        m2=numb(2)
	m3=numb(3)
        m4=numb(4)
        m7=numb(7)
c       write(9,*)'m1-m4=',m1,m2,m3,m4
c100821
        nctl0=nctl
        nctlm=nctl0   ! 100821
c060805	mstj(1)=mstj1_1   ! 221203

c       loop over hadron-hadron collisions in a nucleus-nucleus collision
	iii=1   ! 220110, iii-th hadron-hadron collis.
10 	if(iii.eq.1)goto 1000
101	call copl(time)
c	find out the binary colli. with minimum collsion time
1000	call find(icp,tcp,lc,tc,tw,1)
	if(icp.eq.0)goto 100
c	icp=0 means the collision list is empty
	l=lc(icp,1)
	l1=lc(icp,2)
cm	write(9,*)'af find, iii,icp,l,l1,tcp=',iii,icp,l,l1,tcp   ! sa
c131019 writing initialization and differential cross section maximum
        if(iiii.eq.1 .and. iii.eq.1)then
        mstp(122)=1
        else
        mstp(122)=0
        endif
c131019
	time0=time
	kfa=ksa(l,2)
	kfb=ksa(l1,2)
        ikfa=iabs(kfa)   ! 070417
        ikfb=iabs(kfb)   ! 070417
        kfaab=iabs(kfa)   ! 060813 120214
        kfbab=iabs(kfb)   ! 060813 120214
	time=tcp
c	record this collision time

cc	tlco(l,4)=tcp
cc	tlco(l1,4)=tcp
20	continue
	ilo=0
	pi(4)=psa(l,4)
	pj(4)=psa(l1,4)
	if(pi(4).lt.1.e-20)pi(4)=1.e-20   ! 041204
	if(pj(4).lt.1.e-20)pj(4)=1.e-20   ! 041204
	do i=1,3
	pi(i)=psa(l,i)
	pj(i)=psa(l1,i)
c	if(pi(4).lt.1.e-5.or.pj(4).lt.1.e-5)then
c     &	pi(4),pj(4),nsa,l,l1,icp,nctl,tcp,tc(icp)
c	do iop=1,nctl
c	write(*,*)lc(iop,1),lc(iop,2),tc(iop)
c	enddo
c	endif
	b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
	enddo
c200601
        pti=dsqrt(pi(1)**2+pi(2)**2)
        ptj=dsqrt(pj(1)**2+pj(2)**2)
c200601
c080818
        if((ipden.eq.0 .and. itden.eq.0 .and. ifram.eq.1) .or. ipden
     c   .ge.11)then
        ss=win
        goto 1067
        endif
c080818
c       write(9,*)'bf. loren iiii,iii,pi,pj(4)=',iiii,iii,pi(4),pj(4)!100821
c	boost to CMS frame of colliding pair
	call lorntz(ilo,b,pi,pj)
	ss=pi(4)+pj(4)
        if(ss.lt.1.e-18)ss=1.e-18

c	calculate the angular 'seta' of the momenta pi and pj
	ctai=pyangl(pi(3),dsqrt(pi(1)**2+pi(2)**2))
	ctaj=pyangl(pj(3),dsqrt(pj(1)**2+pj(2)**2))
	cctai=dcos(ctai)
	cctaj=dcos(ctaj)
	if(cctai.gt.0.)then
c       calculate the 'orentation' of the vector pi
	call codi(pi,cfi1,sfi1,ccta1,scta1)
	else
	call codi(pj,cfi1,sfi1,ccta1,scta1)
	endif

c	perform classical Newton motion   
1067	continue
c       write(9,*)'af. loren pi,pj(4)=',pi(4),pj(4)   ! 100821
        call his(time,lc,tc,tw,istop)   ! 080818
c       write(9,*)'af. his iiii,iii,istop=',iiii,iii,istop   ! 100821
	if(istop.eq.1)goto 100
c	istop=1 means all particles have get out of considered volume

        if(win.ge.parp21)then   !! if 110123, high energy loop
c060805	if((l.le.m4 .and. l1.le.m4) .and. ss.ge.parp21)then   ! if 1
c241110	if(((l.le.m2 .and. l1.le.m2).or.(kfa.eq.2212.and.kfb.eq.-2212)
c241110     c   .or.(kfb.eq.2212.and.kfa.eq.-2212)) .and. ss.ge.parp21)
c241110     c   then   ! if 1
        if((ipden.lt.2.and.(l.le.m2 .and. l1.le.m2) .and. ss.ge.parp21)
     c	 .or.(ipden.gt.2.and.(((kfaab.ge.11.and.kfaab.le.16).and.l1.le.
     c	 m2).or.(l.le.m2.and.(kfbab.ge.11.and.kfbab.le.16))).and.ss.ge.
     c	 parp21))then   ! if 1 011210 060813 120214
c060813	m7 to m2
c171022 Now long-written statements are replaced as calling the subroutine xevent.
c171022 The statements are re-wrote for better readability of code.
c       Executes collision event.
        call xevent(ifram,kfa,kfb,ss,pti,ptj,cctai,cctaj)   ! 171022 Lei

500     continue
        if((ipden.eq.0 .and. itden.eq.0 .and. ifram.eq.1) .or. 
     c	 ipden.ge.11)goto 1066   ! 080818
	do j=1,n
	do j1=1,4
	pint(j1)=p(j,j1)
	enddo
c     c	 write(9,*)'before pint=',(pint(i1),i1=1,3)   ! sa

c	if(cctai.gt.0.99)goto 1002

c	'cctai.gt.0.99' means pi (or pj) nearly on the z axis, don't need 
c	  rotation
c	perform the rotate for produced particle from calling 'pythia'
	call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
	do j1=1,4
	p(j,j1)=pint(j1)
	enddo
	enddo

1002	continue
c	boost back to Lab.
	ilo=1
	do j=1,n,2
	if(j.eq.n)then
	do j1=1,4
	pi(j1)=p(j,j1)
	enddo
	call lorntz(ilo,b,pi,pi)
	do j1=1,4
	p(j,j1)=pi(j1)
	enddo
	goto 510
	endif
	do j1=1,4
	pi(j1)=p(j,j1)
	pj(j1)=p(j+1,j1)
	enddo
	call lorntz(ilo,b,pi,pj)
	do j1=1,4
	p(j,j1)=pi(j1)
	p(j+1,j1)=pj(j1)
	enddo
510	enddo
1066	continue
c131019
c       remove gamma from 'pyjets' to 'sgam'
        n66=0
        do j=1,n
        kf=k(j,2)
        if(kf.eq.22)then
        k(j,2)=66
        n66=n66+1
        endif
        enddo
c       move "66" from 'pyjets' to 'sgam'
        if(n66.gt.0)call remo_gam(66)

c       checking charge conservation for current hadron-hadron (hh) collision        
        chai=pychge(kfa)+pychge(kfb)
        chaf=0.
        do i3=1,n
        ik=k(i3,2)
        chaf=chaf+pychge(ik)
        enddo                                                          
        if(abs(chaf-chai).gt.0.1d0)then   !! charge was not conserved
        siww=siww+1.

c       re-generate the current hh collision for pp and lepton incidence 
        if((ipden.eq.0 .and. itden.eq.0) .or. ipden.ge.11)then
        iiii=iiii-1
        ijk=1
        return
        endif

c       otherwise, remove current hh collision pair from collision list
        do j1=icp+1,nctl   
        j=j1-1   
        tc(j)=tc(j1)   
        tw(j)=tw(j1)   
        do m=1,5   
        lc(j,m)=lc(j1,m)   
        enddo   
        enddo   
        nctl=nctl-1   
        iii=iii+1   
        goto 10
        endif   !!
c131019      
c-------------------------------------------------------------------
c	give four position to the particles after calling pyevnt
c110517 for particles in pyjets
        call ptcre(l,l1,time)
c       arrange particles (quark,diquark, and gluon mainly) after
c        calling pyevnt into the overlap region randomly
c061207	call ptcre(l,l1,time,gamt)
c--------------------------------------------------------------------
        noinel(592)=noinel(592)+1   ! 280722
c	592-th scattering process is referred to calling 'pythia'
c	write(9,*)'# of calling pythia=',noinel(592)   ! 070802

        if(mstptj.eq.1)goto 997   ! toward the case PYTHIA-like 230722

c260314	statistics of number of leptons studied, identify scattered lepton,  
c	 and fill up pscal(5) 
	if(ipden.ge.11.and.ipden.le.16)then   !
c	identify the studied leptons
	kfl=ipden
	if(nzp.gt.0.)kfl=-ipden
	nlep=0
	do j=1,n
	ikl=k(j,2)
	if(ikl.eq.kfl)then
	nlep=nlep+1
	pl(nlep,1)=p(j,1)
	pl(nlep,2)=p(j,2)
	pl(nlep,3)=p(j,3)
	pl(nlep,4)=p(j,4)
	pl(nlep,5)=p(j,5)
	endif
	enddo
c	find the scattered lepton (with largest energy among studied leptons)
	if(nlep.gt.1)then   !!
	vnlep=vnlep+nlep
	elep=1.d0
	do j1=1,nlep
	plj14=pl(j1,4)
	if(plj14.ge.elep)then
	elep=plj14
	jj=j1
	endif
	enddo
	do j2=1,5
	pscal(j2)=pl(jj,j2)
	enddo
	elseif(nlep.eq.1)then   !!
	vnlep=vnlep+nlep
	do j2=1,5
        pscal(j2)=pl(nlep,j2)
        enddo
	else   !!
	endif   !!
c	write(9,*)'pscal=',(pscal(i1),i1=1,5)
c	calculate kinematic variables relevant to incident and scattered 
c	 lepton only, in cms
	pdotk=pinch(4)*pincl(4)-pinch(1)*pincl(1)-pinch(2)*pincl(2)   
     c   -pinch(3)*pincl(3)   ! P.k  
	q11=pincl(1)-pscal(1)
	q22=pincl(2)-pscal(2)
	q33=pincl(3)-pscal(3)
	q44=pincl(4)-pscal(4)
	q112=q11*q11
	q222=q22*q22
	q332=q33*q33
	q442=q44*q44
	pdotq=pinch(4)*q44-pinch(1)*q11-pinch(2)*q22-pinch(3)*q33   ! P.q
	vnu=pdotq/pinch(5)   ! \nu
	fq2=-(q442-q112-q222-q332)   ! Q^2=-q^2
	w2l=(pinch(4)+q44)**2-(pinch(1)+q11)**2-(pinch(2)+q22)**2-
     c	 (pinch(3)+q33)**2   ! W^2
	pdotk=dmax1(pdotk,1.d-20)
	yyl=pdotq/pdotk   ! y
	pdotq=dmax1(pdotq,1.d-20)
	xb=fq2/2./pdotq   ! x_b
	endif   !
c260314
	igq=0
	do j1=1,n
	kfj1=iabs(k(j1,2))
        if(kfj1.le.8.or.kfj1.eq.2101.or.kfj1.eq.3101.or.kfj1.eq.3201
     c   .or.kfj1.eq.1103.or.kfj1.eq.2103.or.kfj1.eq.2203.or.kfj1.eq.
     c   3103.or.kfj1.eq.3203.or.kfj1.eq.3303.or.kfj1.eq.21)igq=igq+1! 140805
	enddo
	if(igq.eq.0)then   ! no q, diquark, and g at all
c080818
        if((ipden.eq.0 .and. itden.eq.0) .or. ipden.ge.11)then
        iiii=iiii-1
        ijk=1
        return
        endif
c080818
c	remove current nn collision pair from collision list
	do j1=icp+1,nctl   ! active on 241110
	j=j1-1   ! 241110
	tc(j)=tc(j1)   ! 241110
        tw(j)=tw(j1)   ! 241110
        do m=1,5   ! 241110
        lc(j,m)=lc(j1,m)   ! 241110
        enddo   ! 241110   
        enddo   ! 241110
        nctl=nctl-1   ! 241110
        iii=iii+1   ! 060805 241110
        goto 10   ! 241110
	endif

c       removes hadrons from 'pyjets' to 'sbh' and truncate 'pyjets'
c	 correspondingly
	call remo   ! 010418,161021 removed from after 'recons' to before
c100223
        if(adj12.eq.2)then
        call stahad
        goto 997   ! similar to pythia-like
        endif
c100223
        if(ipden.lt.11)call pyedit(2)   
        if(ipden.ge.11)call pyedit(1)  
c260620 write(22,*)'af. remo iii,n=',iii,n
c260620 call pylist(1)

c161021 reconstruct nucleon (anti-nucleon) in order to increase leading 
c        proton effect
c161021 if(nap.ne.1.and.nat.ne.1)then
        irecon=irecon+1
c161021 write(22,*)'be. recons iiii,iii=',iiii,iii
c161021 call pyedit(2)
c161021 call pylist(1)
c161021 write(22,*)'af. recons ='
        nppb=0   ! 281121
        if((ipden.eq.0 .and. itden.eq.1) .or. 
     c   (ipden.eq.1 .and. itden.eq.0))then
c       write(9,*)'bf. recons iiii,iii=',iiii,iii   !!
c       reconstructs leading particle according to probability distribution 
c        resulted from impact parameter density distribution of f(b)=b
c100322 Lei
        r_max = max(suppc,suptc)
        call bp_prob(bp,r_max,probb)  
c140322 nlead = probb*4.+1   ! '4': a model parameter
        if(pyr(1).ge.probb)then   ! 140322
c100322 Lei
        call recons(irecon,l,l1,ss,nlead,time,iii)   ! 150322   
c150322 call recons_g(irecon,l,l1,ss,nlead,time,iii)   
c       call recons_gg (irecon,l,l1,ss,nlead,time,iii)
        endif   ! 140322
        endif  
c161021 call pyedit(2)
c161021 call pylist(1)
c161021 endif

c080104
c	'pyjets' to 'sbe'. etc.
	if(n.ge.1)then   ! 1
	do i1=1,n
	i3=i1+nbe
	kf=k(i1,2)
        kfab=iabs(kf)
c150520 identifies diquarks        
        if(kfab.eq.2101 .or. kfab.eq.3101 .or. kfab.eq.3201 .or. kfab
     c   .eq.1103 .or. kfab.eq.2103 .or. kfab.eq.2203 .or. kfab.eq.3103
     c   .or. kfab.eq.3203 .or. kfab.eq.3303)then   ! 2
c     c   .or. kfab.eq.3203 .or. kfab.eq.3303 .or. kfab.eq.21)then   ! 2
        idi=idi+1
        ndiq(i1+naf)=idi
	endif   ! 2
c150520 'pyjets' to 'sbe'        
        do i2=1,5
        kbe(i3,i2)=k(i1,i2)
        pbe(i3,i2)=p(i1,i2)
        vbe(i3,i2)=v(i1,i2)
        enddo
	enddo
	nbeo=nbe   ! 190204
	nbe=i3
	endif   ! 1
c	write(9,*)'af. fill sbe n,nbe,irecon,iii,event=',n,nbe,irecon,
c     c   iii,iiii   ! sa
c	goto 200
c080104
c       break up diquark and give four momentum and four position
c        to the broken quarks (working in 'pyjets')
	call break
        if(ipden.lt.11)call pyedit(2)   
        if(ipden.ge.11)call pyedit(1)  
c241110
c        if(iiii.eq.5 .and. (iii.ge.12.and.iii.lt.15))then
c        write(22,*)'af break'
c        call pylist(1)
c        call prt_sbh(nbh,cc)
c        endif
c241110
c030620
c       find number of strings and line number of first and last components
c        of each string
        nstr1=0
        jb=0
10000   do i1=jb+1,n
        if(k(i1,1).eq.2)then   ! i1 is 'A'
        do i2=i1+1,n
        if(k(i2,1).eq.1)then   ! i2 is 'V'
        nstr1=nstr1+1
        nstr1a(nstr1)=i1   ! line number of first component of nstr-th string
        nstr1v(nstr1)=i2   ! line number of first component of nstr-th string
        jb=i2
        if(jb.lt.n)goto 10000
        if(jb.eq.n)goto 20000
        endif
        enddo
        endif
        enddo
20000   continue
        nstr0=nstr1   ! 090620
c       nstr1: number of strings after call break
c260620 write(22,*)'af. break nstr1,iii,n=',nstr1,iii,n
c260220 call pylist(1)
c030620
        goto 777   ! without Pauli blocking 230121
c191202
c       Pauli effect (working in 'pyjets' (current), in 'saf' (past))
        tpaul=1.
c       tpaul: product of the unoccupation probabilities
        do i1=1,n   ! current
        kfp=k(i1,2)
        kfp=iabs(kfp)
        ppaul=1.
        if(kfp.eq.1 .or. kfp.eq.2 .or. kfp.eq.3)then
        
	call pauli(i1,ppaul)
c       ppaul: the unoccupation probability of particle i1
        if(ppaul .lt. 0.)then   ! over occupation, should be blocked
        tpaul=0.
        goto 666
        endif
        endif
        tpaul=tpaul*ppaul
        enddo
666     if(pyr(1) .ge. tpaul)then   ! blocked "1"
c080104 
	ipau=1
c190204
c	remove "current part of 'sbe'" from 'sbe' and truncate 'ndiq' and 
c	 'npt' correspondingly
        do i1=nbeo+1,nbe
        do i2=1,5
        kbe(i1,i2)=0
        pbe(i1,i2)=0.
        vbe(i1,i2)=0.
        enddo
	ndiq(i1+naf)=0
        enddo
	do i1=idio,idi
	npt(i1)=0
        ifcom(i1)=0   ! 220110
	enddo
	nbe=nbeo
c190204
	idi=idio    
c080104
c	remove current nn collision pair from collision list
	do j1=icp+1,nctl
	j=j1-1
	tc(j)=tc(j1)
        tw(j)=tw(j1)
        do m=1,5
        lc(j,m)=lc(j1,m)
        enddo
	enddo
	nctl=nctl-1
	goto 10
        endif   ! "1"
777	continue   ! unblocked
c080104
c       'pyjets' to 'saf'. etc.
        if(n.ge.1)then   
        do i1=1,n
        naf=naf+1
	if(naf.gt.kszj)then
	write(9,*)'iiii,naf,kszj=',iiii,naf,kszj   ! sa
	stop 11111
	endif
        do i2=1,5
        kaf(naf,i2)=k(i1,i2)
        paf(naf,i2)=p(i1,i2)
        vaf(naf,i2)=v(i1,i2)
        enddo
        enddo
        endif
c080104
200	continue
	idio=idi   ! 080104
c080104 241110
c       if(iiii.eq.5 .and. iii.eq.13)then
c       write(22,*)'be updpip'
c       write(9,*)'be updpip iiii,iii=',iiii,iii
c       call prt_sa2(nsa,cc)
c       call prt_sbh(nbh,cc)
c       endif
c080104 241110
c140414

997     continue   ! 230722

c       add CME charge separation for u d s,c   ! 042021 She
c        print*,"n=",n
        icme=adj1(23)   ! 221022
        if(icme.eq.0)goto 902
        if((nap.eq.nat).and.(nzp.eq.nzt).and.(icme.eq.1)) 
     c    call chargecme(win)
902     continue

c230722	if(adj140.eq.5)then
        if(mstptj.eq.1)then   ! 230722
c	if(ipden.lt.11)call pyedit(2)   
c	if(ipden.ge.11)call pyedit(1)
c	'pyjets' to 'sbh'
        if(n.eq.0)goto 5001
        do li=1,n
        do m=1,5
        kbh(li,m)=k(li,m)
        pbh(li,m)=p(li,m)
        vbh(li,m)=v(li,m)
        enddo
        enddo
        nbh=n
        do i=n+1,kszj
        do j=1,5
        kbh(i,j)=0
        pbh(i,j)=0.
        vbh(i,j)=0.
        enddo
        enddo
5001    continue
c	n=0
c	do i=1,kszj
c        do j=1,5
c        k(i,j)=0
c        p(i,j)=0.
c        v(i,j)=0.
c        enddo
c        enddo
	endif
c       write(9,*)'be. updpip iiii,iii,nctl=',iiii,iii,nctl0,nctlm,nctl ! 100821
c140414
c281121 update hadron list 'sa2' after calling pythia ('SBH' to 'sa2'), 
c        remove collision pair composed of l and/or l1, remove l (l1)
c        from 'sa2'
	call updpip(l,l1,icp,lc,tc,tw,time,iii)   
c241110
c       if(iiii.eq.5 .and. (iii.ge.12.and.iii.lt.15))then
c       write(9,*)'af. updpip iiii,iii,nctl=',iiii,iii,nctl0,nctlm,nctl ! 100821
c       call prt_sa2(nsa,cc)
c       call prt_saf(naf,cc)
c       endif
c241110
c011204	l=lc(icp,1)
c011204	l1=lc(icp,2)
c       update collision list after calling 'pythia'
        call updtlp(time,lc,tc,tw,iii)
c       write(9,*)'af. updtlp iiii,iii,nctl=',iiii,iii,nctl0,nctlm,nctl   ! sa
cm	call prt_sa2(nsa,cc)   
cm	do i=1,nctl
cm	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
cm	enddo
	if(nctl.eq.0)goto 100   ! 021204
c170121
c240121 noinel(592)=noinel(592)+1  
c       noinel(592): statistics of # of nn collition calling pythia
c170121
	goto 300   ! ss is enough to call pythia	
	endif   ! if 1

c010223 if ss is not enough to call pythia or current hadron-hadron
c        collision pair is not in the plan of calling pythia, then it is 
c       treat by two-body collison kinematic (ela. and inela.) 
        noinel(593)=noinel(593)+1   ! 140820
c140820 noinel(593): statistics of # of hadron-hadron collition which energy 
c        is not enough to call pythia or current hadron-hadron
c        collision pair is not in the plan of calling pythia

        ww=rcsit
c       the cross section ratio of (ela.)/tot =1- rcsit
c       ww=0.   ! active  for ela. scattering only
c       ww=1.   ! active  for inela. scattering only
        rrlu=pyr(1)

        if(rrlu.gt.ww)then    !!1
        call coelas(l,l1,ss,pi,pj)
c       update the particle list for elastic scattering, pi and pj have been
c       boosted back to Lab frame or cms of nucleus-nucleus collision 
        call updple(l,l1,b,pi,pj)
c	statistics of the number of ela. hadron-hadron collisions for 
c        high energy channel        
        noinel(1)=noinel(1)+1
c       update the collision list after ela. scattering
        call updatl(l,l1,time,lc,tc,tw,iii)
c       note: CME is not included for ela. scattering
        if(nctl.eq.0)goto 100
        goto 300

        elseif(rrlu.le.ww)then   !!1 
c       perform hadron-hadron inelastic scattering                
c     p: 2212; n:2112; delta0: 2114; delta+: 2214; delta-: 1114; delta++: 2224
c       consider following 2->2 nucleon-nucleon ielas. channels 
c     1       p + p to delta+ + p
c     2       p + p to delta++ + n
c     3       p + n to delta+ + n
c     4       p + n to delta0 + p
c     5       n + n to delta0 + n
c     6       n + n to delta- + p
c       reverse scattering is not considered

c       give flavor code to inelastically scattered particles
        rpy=pyr(1)
        if(kfa.eq.2212.and.kfb.eq.2212)then   ! pp   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2212
c       give four momentum to inelastically scattered particles 
        call coinel_nn(l,l1,2214,2212,pi,pj)
c       decay of unstable hadron        
        call padecy(l,l1,time)
c       output particles of 'padecy' compose 'pyjets'
        jorn=3   ! delta+ -> p + pi0
c       jorn=0 (3): inelastically scattered particles not join (join) 
c        reconstruction of hadronic collision pair
        goto 302
        else
        ksa(l,2)=2224   ! delta++
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2224,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta++ -> p + pi+
        goto 302
        endif

        elseif(kfa.eq.2212.and.kfb.eq.2112)then   ! pn   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2214,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta+ -> p + pi0
        goto 302
        else
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=2212 
        call coinel_nn(l,l1,2114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 302
        endif

        elseif(kfa.eq.2112.and.kfb.eq.2212)then   ! np   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2214,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta+ -> p + pi0
        goto 302
        else
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=2212
        call coinel_nn(l,l1,2114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 302
        endif

        elseif(kfa.eq.2112.and.kfb.eq.2112)then   ! nn   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2114   !  delta0
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2114,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 302
        else
        ksa(l,2)=1114   ! delta-
        ksa(l1,2)=2212
        call coinel_nn(l,l1,1114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta- -> n + pi-
        goto 302
        endif

        else   !!
c       otherwise treat as ela.
c       perform hadron-hadron ela. scattering and rotation
        call coelas(l,l1,ss,pi,pj)
c       boosted back to Lab frame (cms of nucleus-nucleus collision) &
c        then update the particle list for ela. scattering
c       boosted back to Lab frame or cms of nucleus-nucleus collision
        call updple(l,l1,b,pi,pj)
c       update the collision list after elastic scattering,
        call updatl(l,l1,time,lc,tc,tw,iii)
        if(nctl.eq.0)goto 100
        goto 300

        endif   !!

302     continue
c       prepare rotation        
        do j=1,n
        do j1=1,4
        pint(j1)=p(j,j1)
        enddo
c       perform the rotate for inelas. scattered particles in 'pythia'
        call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
        do j1=1,4
        p(j,j1)=pint(j1)
        enddo
        enddo
c       write(22,*)'af. rosa n=',iii,n
c       call prt_pyj(n,c)

c       boost back to Lab.
        ilo=1
        do j=1,n,2
        if(j.eq.n)then
        do j1=1,4
        pi(j1)=p(j,j1)
        enddo
        call lorntz(ilo,b,pi,pi)
        do j1=1,4
        p(j,j1)=pi(j1)
        enddo
        goto 5200
        endif
        do j1=1,4
        pi(j1)=p(j,j1)
        pj(j1)=p(j+1,j1)
        enddo
        call lorntz(ilo,b,pi,pj)
        do j1=1,4
        p(j,j1)=pi(j1)
        p(j+1,j1)=pj(j1)
        enddo
5200    enddo
c       write(22,*)'af. lor n=',n
c       call prt_pyj(n,c)

c       give four position to the particles after after inelas. scattering
c        in 'pyjets'
        call ptcre(l,l1,time)   
c       write(22,*)'af. ptcre n=',n
c       call prt_pyj(n,c)

c       add CME charge separation for u d s,c   ! 042021 She
c        print*,"n=",n
        icme=adj1(23)   ! 221022
        if(icme.eq.0)goto 9023
        if((nap.eq.nat).and.(nzp.eq.nzt).and.(icme.eq.1))
     c    call chargecme(win)
9023     continue
c       write(22,*)'af. CME n=',n
c       call prt_pyj(n,c)

c       'pyjets' to 'sbh'
        nbh=0
        do i=1,kszj
        do j=1,5
        kbh(i,j)=0
        pbh(i,j)=0.
        vbh(i,j)=0.
        enddo
        enddo
        if(n.eq.0)goto 5003
        nbh=n
        do li=1,nbh
        do m1=1,5
        kbh(li,m1)=k(li,m1)
        pbh(li,m1)=p(li,m1)
        vbh(li,m1)=v(li,m1)
        enddo
        enddo
5003    continue   
c       write(22,*)'pyjets to sbh nbh=',nbh
c       call prt_sbh(nbh,c)
c       write(22,*)'flavor=',(kbh(i1,2),i1=1,nbh)     
c       call prt_sa2(nsa,c)                
c       update hadron list 'sa2' after inela. scattering  ('sbh' to 'sa2')
        nsa0=nsa   ! 250123
        call updpip_nn(l,l1,icp,lc,tc,tw,time)   ! 250123
c       write(22,*)'af. updpip nsa=',nsa
c       call prt_sa2(nsa,c)

c       update collision time list after inela. scattering
        call updatl_nn(l,l1,time,lc,tc,tw,jorn,nsa0)   ! 250123
c       write(22,*)'end of loop iii,nsa=',iii,nsa
c       call prt_sa2(nsa,c)
        if(nctl.eq.0)goto 100
        goto 300
        else   !!1
        endif   !!1

c010223        
        endif   !! if 110123, high energy loop end

c161222
        if(win.lt.parp21)then    !!! if 110123, low energy loop
        ww=rcsit
c       ww=0.   ! active  for ela. scattering only
c       ww=1.   ! active  for inela. scattering only
c       if(iiii.eq.1.and.iii.eq.1)write(9,*)'inela./tot=rcsit=',rcsit
c       the cross section ratio of (ela.)/tot =1- rcsit
        rrlu=pyr(1)

        if(rrlu.gt.ww)then   !1 ela.
c       perform nucleon-nucleon (NN) elastic scattering and rotation            
        call coelas(l,l1,ss,pi,pj)
c       boosted back to Lab frame (cms of nucleus-nucleus collision) & 
c        then update the particle list for elastic scattering
c       boosted back to Lab frame or cms of nucleus-nucleus collision 
        call updple(l,l1,b,pi,pj)
c       update the collision list after elastic scattering, 
        call updatl(l,l1,time,lc,tc,tw,iii)
        if(nctl.eq.0)goto 100
        goto 300
        
        elseif(rrlu.le.ww)then   !1 inela.
c       perform NN inelastic scattering                
c     p: 2212; n:2112; delta0: 2114; delta+: 2214; delta-: 1114; delta++: 2224
c       consider following 2->2 nucleon-nucleon ielas. channels 
c     1       p + p to delta+ + p
c     2       p + p to delta++ + n
c     3       p + n to delta+ + n
c     4       p + n to delta0 + p
c     5       n + n to delta0 + n
c     6       n + n to delta- + p
c       reverse scattering is not considered
c       perform piN inelastic scattering
c     pi+: 211; pi-: -211; pi0: 111 
c       consider following 2->2 pi-nucleon ielas. channels
c     7	      pion- + p to delta- + pion+
c     8	      pion- + p to rho0 + n
c     9	      pion- + p to rho- + p
c     10      pion- + p to delta+ + pion-
c     11      pion- + p to delta0 + pion0
c     12      pion- + n to delta- + pion0
c     13      pion- + n to rho- + n
c     14      pion- + n to delta0 + pion-
c     15      pion+ + p to delta++ + pion0
c     16      pion+ + p to delta+ + pion+
c     17      pion+ + p to rho+ + p
c     18      pion+ + n to delta++ + pion-
c     19      pion+ + n to delta0 + pion+
c     20      pion+ + n to delta+ + pion0
c     21      pion+ + n to rho0 + p
c     22      pion+ + n to rho+ + n
c       assume pion0 decay instantly, it is not able to reconstrut
c        collion pair, following processes is not considered 
c     23      pion0 + p to delta0 + pion+
c     24      pion0 + p to delta++ + pion-
c     25      pion0 + p to rho+ + n
c     26      pion0 + p to rho0 + p
c     27      pion0 + p to delta+ + pion0
c     28      pion0 + n to delta+ + pion-
c     29      pion0 + n to delta- + pion+
c     30      pion0 + n to delta0 + pion0

c       write(22,*)'enter loop iii,nsa,kfa,kfb=',iii,nsa,kfa,kfb
c       give flavor code to inelastically scattered particles
        rpy=pyr(1)
        if(kfa.eq.2212.and.kfb.eq.2212)then   ! pp   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2212
c       give four momentum to inelastically scattered particles 
        call coinel_nn(l,l1,2214,2212,pi,pj)
c       decay of unstable hadron        
        call padecy(l,l1,time)
c       output particles of 'padecy' compose 'pyjets'
        jorn=3   ! delta+ -> p + pi0
c       jorn=0 (3): inelastically scattered particles not join (join) 
c        reconstruction of hadronic collision pair
        goto 301
        else
        ksa(l,2)=2224   ! delta++
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2224,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta++ -> p + pi+
        goto 301
        endif

        elseif(kfa.eq.2212.and.kfb.eq.2112)then   ! pn   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2214,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta+ -> p + pi0
        goto 301
        else
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=2212 
        call coinel_nn(l,l1,2114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 301
        endif

        elseif(kfa.eq.2112.and.kfb.eq.2212)then   ! np   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2214,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta+ -> p + pi0
        goto 301
        else
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=2212
        call coinel_nn(l,l1,2114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 301
        endif

        elseif(kfa.eq.2112.and.kfb.eq.2112)then   ! nn   !!
        if(rpy .gt. 0.5)then
        ksa(l,2)=2114   !  delta0
        ksa(l1,2)=2112
        call coinel_nn(l,l1,2114,2112,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta0 -> n + pi0
        goto 301
        else
        ksa(l,2)=1114   ! delta-
        ksa(l1,2)=2212
        call coinel_nn(l,l1,1114,2212,pi,pj)
        call padecy(l,l1,time)
        jorn=3   ! delta- -> n + pi-
        goto 301
        endif

        elseif(kfa.eq.-211.and.kfb.eq.2212)then   ! (pi-)p   !!
        if(rpy .le. 0.2)then        
        ksa(l,2)=1114   ! delta-
        ksa(l1,2)=211
        call coinel_nn(l,l1,1114,211,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.2 .and. rpy .le. 0.4)then
        ksa(l,2)=113   ! rho0
        ksa(l1,2)=2112
        call coinel_nn(l,l1,113,2112,pi,pj)
        call sa2pyj(l,l1) 
        jorn=0 
        goto 301        
        elseif(rpy .gt. 0.4 .and. rpy .le. 0.6)then
        ksa(l,2)=-213   ! rho-
        ksa(l1,2)=2212
        call coinel_nn(l,l1,-213,2212,pi,pj)
        call sa2pyj(l,l1)
        jorn=0   
        goto 301
        elseif(rpy .gt. 0.6 .and. rpy .le. 0.8)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=-211
        call coinel_nn(l,l1,2214,-211,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.8)then
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=111
        call coinel_nn(l,l1,2114,111,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        endif

        elseif(kfa.eq.-211.and.kfb.eq.2112)then   ! (pi-)n   !!
        if(rpy .le. 0.3333)then
        ksa(l,2)=1114   ! delta-
        ksa(l1,2)=111
        call coinel_nn(l,l1,1114,111,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.3333 .and. rpy .le. 0.6666)then
        ksa(l,2)=-213   ! rho-
        ksa(l1,2)=2112
        call coinel_nn(l,l1,-213,2112,pi,pj)
        call sa2pyj(l,l1)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.6666)then
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=-211
        call coinel_nn(l,l1,2114,-211,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        endif        
                
        elseif(kfa.eq.211.and.kfb.eq.2212)then   ! (pi+)p   !!
        if(rpy .le. 0.3333)then
        ksa(l,2)=2224   ! delta++
        ksa(l1,2)=111
        call coinel_nn(l,l1,2224,111,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.3333 .and. rpy .le. 0.6666)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=211
        call coinel_nn(l,l1,2214,211,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.6666)then
        ksa(l,2)=213   ! rho+
        ksa(l1,2)=2212
        call coinel_nn(l,l1,213,2212,pi,pj)
        call sa2pyj(l,l1)
        jorn=0
        goto 301        
        endif
 
        elseif(kfa.eq.211.and.kfb.eq.2112)then   ! (pi+)n   !!
        if(rpy .le. 0.2)then        
        ksa(l,2)=2224   ! delta++
        ksa(l1,2)=-211
        call coinel_nn(l,l1,2224,-211,pi,pj)
        call padecy(l,l1,time)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.2 .and. rpy .le. 0.4)then
        ksa(l,2)=2114   ! delta0
        ksa(l1,2)=211
        call coinel_nn(l,l1,2114,211,pi,pj)
        call padecy(l,l1,time)
        jorn=0 
        goto 301        
        elseif(rpy .gt. 0.4 .and. rpy .le. 0.6)then
        ksa(l,2)=2214   ! delta+
        ksa(l1,2)=111
        call coinel_nn(l,l1,2214,111,pi,pj)
        call padecy(l,l1,time)
        jorn=0   
        goto 301
        elseif(rpy .gt. 0.6 .and. rpy .le. 0.8)then
        ksa(l,2)=113   ! rho0
        ksa(l1,2)=2212
        call coinel_nn(l,l1,113,2212,pi,pj)
        call sa2pyj(l,l1)
        jorn=0
        goto 301
        elseif(rpy .gt. 0.8)then
        ksa(l,2)=213   ! rho+
        ksa(l1,2)=2112
        call coinel_nn(l,l1,213,2112,pi,pj)
        call sa2pyj(l,l1)
        jorn=0
        goto 301
        endif

        else   !!
c010223                
c       otherwise treat as ela. 
c       perform hadron-hadron ela. scattering and rotation
        call coelas(l,l1,ss,pi,pj)
c       boosted back to Lab frame (cms of nucleus-nucleus collision) &
c        then update the particle list for ela. scattering
c       boosted back to Lab frame or cms of nucleus-nucleus collision
        call updple(l,l1,b,pi,pj)
c       update the collision list after elastic scattering,
        call updatl(l,l1,time,lc,tc,tw,iii)
        if(nctl.eq.0)goto 100
        goto 300
c010223

        endif   !!

301     continue       
c       write(22,*)'af. inela. =',
c       call prt_pyj(n,c)        
c       prepare rotation        
        do j=1,n
        do j1=1,4
        pint(j1)=p(j,j1)
        enddo
c       perform the rotate for inelas. scattered particles in 'pythia'
        call rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,pint)
        do j1=1,4
        p(j,j1)=pint(j1)
        enddo
        enddo
c       write(22,*)'af. rosa n=',iii,n
c       call prt_pyj(n,c)

c       boost back to Lab.
        ilo=1
        do j=1,n,2
        if(j.eq.n)then
        do j1=1,4
        pi(j1)=p(j,j1)
        enddo
        call lorntz(ilo,b,pi,pi)
        do j1=1,4
        p(j,j1)=pi(j1)
        enddo
        goto 5100
        endif
        do j1=1,4
        pi(j1)=p(j,j1)
        pj(j1)=p(j+1,j1)
        enddo
        call lorntz(ilo,b,pi,pj)
        do j1=1,4
        p(j,j1)=pi(j1)
        p(j+1,j1)=pj(j1)
        enddo
5100    enddo
c       write(22,*)'af. lor n=',n
c       call prt_pyj(n,c)

c       give four position to the particles after after inelas. scattering
c        in 'pyjets'
        call ptcre(l,l1,time)   
c       write(22,*)'af. ptcre n=',n
c       call prt_pyj(n,c)

c       add CME charge separation for u d s,c   ! 042021 She
c        print*,"n=",n
        icme=adj1(23)   ! 221022
        if(icme.eq.0)goto 9022
        if((nap.eq.nat).and.(nzp.eq.nzt).and.(icme.eq.1))
     c    call chargecme(win)
9022     continue
c       write(22,*)'af. CME n=',n
c       call prt_pyj(n,c)

c       'pyjets' to 'sbh'
        nbh=0
        do i=1,kszj
        do j=1,5
        kbh(i,j)=0
        pbh(i,j)=0.
        vbh(i,j)=0.
        enddo
        enddo
        if(n.eq.0)goto 5002
        nbh=n
        do li=1,nbh
        do m1=1,5
        kbh(li,m1)=k(li,m1)
        pbh(li,m1)=p(li,m1)
        vbh(li,m1)=v(li,m1)
        enddo
        enddo
5002    continue   
c       write(22,*)'pyjets to sbh nbh=',nbh
c       call prt_sbh(nbh,c)
c       write(22,*)'flavor=',(kbh(i1,2),i1=1,nbh)     
c       call prt_sa2(nsa,c)                
c       update hadron list 'sa2' after inela. scattering  ('sbh' to 'sa2')
        nsa0=nsa   ! 250123
        call updpip_nn(l,l1,icp,lc,tc,tw,time)   ! 250123
c       write(22,*)'af. updpip nsa=',nsa
c       call prt_sa2(nsa,c)

c       update collision time list after inela. scattering
        call updatl_nn(l,l1,time,lc,tc,tw,jorn,nsa0)   ! 250123
c       write(22,*)'end of loop iii,nsa=',iii,nsa
c       call prt_sa2(nsa,c)
        if(nctl.eq.0)goto 100
        goto 300

        else   !1         
        endif  !1

        endif   !!! if 110123   !!! if 110123, low energy loop end

300	continue
        iii=iii+1
	if(iii.gt.100*(nctl0))then
        write(9,*)'infinite loop may have happened in'
        write(9,*)'subroutine scat iiii=',iiii
c10/08/98       stop 'infinite loop occurs'
        iiii=iiii-1   ! 10/08/98
        ijk=1  ! 10/08/98
        return   ! 10/08/98
        endif

        if(nctl.gt.nctlm)nctlm=nctl   ! 180121
c       write(9,*)'iiii,iii,nctl0,nctlm=',iiii,iii,nctl0,nctlm   ! 100821 
	goto 10
100	continue
	call copl(time)
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine xevent(ifram,kf_a,kf_b,ss,pt_a,pt_b,ccta_a,ccta_b)
c171022 A new subroutine to execute the binary NN and lN collision.
c       It replaces the previous long-written statements.   ! 171022 Lei
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c   nap,nat,nzp,nzt,pio
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
c       common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
c     &   iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa23/kpar,knn,kpp,knp,kep   ! 060813
        common/sa30/vneump,vneumt,mstptj   ! 230722
        common/sa34/itorw,iikk,cp0,cr0,kkii   ! 060617 010418 010518 040920

        character name_a*16, name_b*16, name_frame*16, name_x*16

        mstp(111)=mstptj   ! =0 230722
c       Gets name of particles a and b.
        call pyname(kf_a,name_a)
        call pyname(kf_b,name_b)
c       Sets name of frame. The collider frame is the default.
        if( ifram.eq.0 )then
            name_frame = "FIXT"
        elseif( ifram.eq.1 )then
            name_frame = "CMS"
        else
            name_frame = "CMS"
        endif
c       ccta_a: the cosine of the angular 'seta' of 
!        the momentum for particle a.
        if( (kf_a.ne.kf_b) .AND. (ccta_a.lt.0.) )then
            name_x = name_a
            name_a = name_b
            name_b = name_x
        endif
c       Initilizes the colllision.
        call pyinit( TRIM(ADJUSTL(name_frame)),
     &               TRIM(ADJUSTL(name_a)),
     &               TRIM(ADJUSTL(name_b)),
     &               ss )
c       Executes the collision. Calling PYEVNW is the default.
        if( itorw.eq.1 )then
            call pyevnt
        elseif( itorw.eq.2 )then
            call pyevnw
        else
            call pyevnw
        endif
c       Removes unnecessary entries in PYJET.
        if( ipden.lt.11 ) call pyedit(2)
        if( ipden.ge.11 ) call pyedit(1)

c       Records the numbers of the participant nucleons (kpar), 
c        binary pp (kpp), pn/np (knp), and nn collisions (knn), 
c        or lepton-nucleon collisions (kep).
        if( (ipden.lt.2) .AND. (itden.lt.2) )then
            if( (kf_a.eq.kf_b) .AND. (kf_a.eq.2212) ) kpp = kpp + 1
            if( (kf_a.eq.kf_b) .AND. (kf_a.eq.2112) ) knn = knn + 1
            if(  kf_a.ne.kf_b  )                      knp = knp + 1
            if(  pt_a.le.1D-4 ) kpar = kpar + 1
            if(  pt_b.le.1D-4 ) kpar = kpar + 1
        elseif( (ipden.ge.11) .AND. (itden.ge.11) )then
            kep = kep + 1
        endif


        return
        end



c************************************************************
        subroutine chargecme(win)
c       The CME-induced charge initial charge separation by switching the 
c        py values of a fraction of the downward(upward) moving(u,d,s,c)quarks 
c        for symmetrical collision systems,i.e., Ru&Ru Zr&Zr at RHIC and LHC.
c       Here in symmetrical systems, nap=nat,nzp=nzt, and the fraction and
c        magnetic field function is A*bp-B*bp^3 type.  by shezl 2021
c       She and Lei For CME.
        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(KSZJ=80000)
        COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     &   nap,nat,nzp,nzt,pio
        common/schuds/schun,schudn,schudsn,sfra
        common/sa24/adj1(40),nnstop,non24,zstop
        dimension numk(kszj)
        real(kind=8) p2u,erhic,erela,rerzcp,ruzcp

        erhic=200.                     ! RHIC energy 200
        erela=0.45+0.55*(win/erhic)  !RHIC energy as a base
        rnzp=real(nzp)
        rnap=real(nap)
        rerz=rnzp/rnap
        ruzcp=((96./42.)*rerz)**(0.667) !isobar Zr Ru(96,42)as a base

        sfra=3.1*(2448.135*nap**(-1.667)*bp-160.810*nap**(-2.333)
     &           *bp**3.)*erela*ruzcp*0.01

c        print*,"erela,ruzcp,nap,nzp,bp,sfra",erela,ruzcp,nap,nzp,bp,sfra

        do i=1,n
        if(abs(k(i,2)).eq.1.or.abs(k(i,2)).eq.2.or.abs(k(i,2)).eq.3
     &   .or.(k(i,2)).eq.4)then
        schun=schun+1
        if(pyr(1).gt.0..and.pyr(1).le.sfra)then
        numk(i)=0
        schudn=schudn+1
        do ii=1,n
        if(numk(ii).eq.1) cycle
        do jj=ii+1,n
        if(numk(jj).eq.1) cycle
        if((k(ii,2)+k(jj,2)).eq.0.and.(k(ii,2)*p(ii,2).lt.0).and.
     &   (k(jj,2)*p(jj,2).lt.0))then

        p2u=p(ii,2)
        p(ii,2)=p(jj,2)
        p(jj,2)=p2u
        schudsn=schudsn+1
        numk(ii)=1
        numk(jj)=1
        endif
        enddo
        enddo
        endif
        endif
        enddo


        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine lorntz(ilo,b,pi,pj)
c	perform Lorentz (or inverse Lorentz) transformation
c	implicit real*8 (a-h,o-z)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	dimension pi(4),pj(4),b(3),dpi(4),dpj(4)
	bb=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
	DB=DSQRT(bb)
        eps1=1d0-1d-12   ! 121108
c121108	IF(DB.GT.0.99999999D0) THEN
        IF(DB.GT.eps1) THEN   ! 121108 
	do i=1,3
c       rescale boost vector if too close to unity. 
c121108	b(i)=b(i)*(0.99999999D0/DB)
        b(i)=b(i)*(eps1/DB)   ! 121108 	
	enddo
c121108	DB=0.99999999D0
        DB=eps1   ! 121108 
	bb=DB**2
	endif
	bbb=1d0-bb
c	if(bbb.le.1.d-10)bbb=1.d-10
	gam=1d0/dsqrt(bbb)
	ga=gam*gam/(gam+1d0)
	do i=1,4
	dpi(i)=pi(i)
	dpj(i)=pj(i)
	enddo
	if(ilo.eq.1) goto 100
c	Lorentz transformation
	pib=dpi(1)*b(1)+dpi(2)*b(2)+dpi(3)*b(3)
	pjb=dpj(1)*b(1)+dpj(2)*b(2)+dpj(3)*b(3)
	do i=1,3
	pi(i)=dpi(i)+b(i)*(ga*pib-gam*dpi(4))
	pj(i)=dpj(i)+b(i)*(ga*pjb-gam*dpj(4))
	enddo
	pi(4)=gam*(dpi(4)-pib)
	pj(4)=gam*(dpj(4)-pjb)
	return
100	continue
c	inverse Lorentz transformation
	pib=dpi(1)*b(1)+dpi(2)*b(2)+dpi(3)*b(3)
	pjb=dpj(1)*b(1)+dpj(2)*b(2)+dpj(3)*b(3)
	do i=1,3
	pi(i)=dpi(i)+b(i)*(ga*pib+gam*dpi(4))
	pj(i)=dpj(i)+b(i)*(ga*pjb+gam*dpj(4))
	enddo
	pi(4)=gam*(dpi(4)+pib)
	pj(4)=gam*(dpj(4)+pjb)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine his(t1,lc,tc,tw,istop)
c	classical Newton motion in Lab. system
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
	parameter(nsize=280000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa4/tau(kszj),tlco(kszj,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	istop=1
	in=0
	do 200 i=1,nsa
	r0=rao*dmax1(rnt,rnp)   ! 060813
c	if(t1.le.tau(i))goto 100
c	do move particles which have not produced
	if(ishp(i).eq.1) goto 10
c	pp4=psa(i,4)
c	do j=1,3
c	vp=psa(i,j)/pp4
c	vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))	
c	enddo
	in=in+1
	goto 200   ! 100 271004
10	aa=0.
	pp4=psa(i,4)
c	due to the fast speed of bayons, we could not use a limited interaction
c	region
c060813	r0=rao*dmax1(rnt,rnp)
c	if(iabs(k(i,2)).gt.1000)r0=1.E+10*r0
	do j=1,3
	vp=psa(i,j)/pp4
	vsa(i,j)=vsa(i,j)+vp*(t1-vsa(i,4))
	aa=aa+(vsa(i,j)-coor(j))**2
	enddo
c251004	vsa(i,4)=t1
	aa=dsqrt(aa)
	if(aa.lt.r0) goto 100
c	if freeze-out occurs deduct the distance between the last collision 
c	and now
	do j=1,3
	vp=psa(i,j)/pp4
	vsa(i,j)=vsa(i,j)-vp*(t1-vsa(i,4))
	enddo
	ishp(i)=0
	do il=1,nctl
	if(lc(il,1).eq.i.or.lc(il,2).eq.i) tc(il)=0.
	enddo
	goto 200   ! 271004
100	continue
	vsa(i,4)=t1   ! 251004
200	continue
	if(in.eq.nsa) return
	istop=0
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ptcre(l,l1,time)   ! 110517   
c	give four position to the particles after calling pythia  
c	l and l1 are colliding particles 060813 120214
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (kszj=80000)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	do i=1,n
c060813 if(ipden.ne.0 .or. itden.ne.0)then
	rl=pyr(1)
	do m=1,3
c	write(*,*)'v=',v(i,m),k(i,2),time
	v(i,m)=v(i,m)+vsa(l,m)*rl+vsa(l1,m)*(1.-rl)
        enddo
c060813 endif
c210921 generated particles are distributed on the surface with unit radius 
	if((ipden.eq.0 .and. itden.eq.0) .or. (ipden.eq.2 .and. 
     c   itden.eq.2))then   ! 180921 yan
	cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=dsqrt(1.-cita**2)
	v(i,1)=sita*dcos(fi)
	v(i,2)=sita*dsin(fi)
	v(i,3)=cita
	endif
	v(i,4)=time
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ptcre_n(l,l1,time,gamt)   ! 021207
c	arrange particles (quark,diquark, and gluon mainly) after 
c	 calling pythia into the overlap region randomly  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (kszj=80000)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	if(ipden.ne.0 .or. itden.ne.0)then
        b=bp/r0t
	do i=1,n
        iii=0
54      iii=iii+1
        if(iii.eq.10000)then
        write(22,*)'difficult to arrange produced particles in'
        write(22,*)'subroutine ptcre,infinitive loop may occur'
        endif
c       sample a point in the unit sphere of target
        x=1.-2.*pyr(1)
        y=1.-2.*pyr(1)
        z=1.-2.*pyr(1)
        rr=dsqrt(x*x+y*y+z*z)
        if(rr.gt.1) goto 54
c       x and y components of that point in the system of unit sphere of
c        projectile are x and y-b, respectively. Adjudge that does (x,y-b) is 
c        in the sphere of projectile
        r1=r0p*dsqrt(x*x+(b-y)*(b-y))
        if(r1.gt.r0p)goto 54
        xx=x*r0t
        yy=y*r0t
        zz=z*r0t/gamt
        v(i,1)=xx
        v(i,2)=yy
        v(i,3)=zz  
c       write(5,*)xx,yy,zz       
        enddo
	endif
	if(ipden.eq.0 .and. itden.eq.0)then
        do i=1,n
	cita=2*pyr(1)-1.
        fi=2.*pio*pyr(1)
        sita=dsqrt(1.-cita**2)
	v(i,1)=sita*dcos(fi)
	v(i,2)=sita*dsin(fi)
	v(i,3)=cita
        v(i,3)=v(i,3)/gamt
        enddo
	endif
        do i=1,n
	v(i,4)=time   ! 230805
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updpip(l,l1,icp,lc,tc,tw,time,iii)   
c281121 update hadron list 'sa2' after calling pythia ('SBH' to 'sa2') 
c        remove collision pair composed of l and/or l1, remove l (l1)
c        from 'sa2'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=80000)
        parameter(nsize=280000)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        COMMON/SBH/N,NONBH,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)   ! 080104   
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121 
	common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa4/tau(kszj),tlco(kszj,4)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa6/kfmaxi,nwhole
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa14/ipyth(2000),idec(2000),iwide(2000)
        common/sa30/vneump,vneumt,mstptj   ! 290123
        common/sppb/nppb,non3,kppb(1000,5),pppb(1000,5),vppb(1000,5) ! 281121
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio   ! 060813
c060813	ipyth: stord line number of produced hadron in hadron list (sa2),101221
        dimension lc(nsize,5),tc(nsize),tw(nsize)  
	dimension peo(4)
        do m=1,2000
        ipyth(m)=0
        enddo
c281121 
c       'sppb' (reconstructed hadrons) to 'SBH' 
        if(((ipden.eq.0 .and. itden.eq.1) .or. (ipden.eq.1 .and. 
     c   itden.eq.0)) .and. (nppb.ge.1 .and. mstptj.eq.0))then   ! 290123
        do i1=1,nppb
        n=n+1
        do i2=1,5
        k(n,i2)=kppb(i1,i2)
        p(n,i2)=pppb(i1,i2)
        v(n,i2)=vppb(i1,i2)
        enddo
        enddo
        do i1=n+1,kszj
        do i2=1,5
        k(i1,i2)=0
        p(i1,i2)=0.
        v(i1,i2)=0.
        enddo
        enddo
        endif
c281121
c       'SBH' to 'sa2' (i.e. produced hadrons-> hadron list 'sa2'),101221 
c241110
c        if(iiii.eq.5)then
c        write(22,*)'in updpip iiii,iii,nsa,l,l1,nbh=',
c     c   iiii,iii,nsa,l,l1,n
c       write(22,*)'kfaco=',(kfaco(i1),i1=1,52)
c        write(22,*)'numb=',(numb(i1),i1=1,52)
c        call prt_sa2(nsa,cc)
c        call prt_sbh(n,cc)
c        endif
c241110
        ll=l
        ll1=l1
        if(n.eq.0)goto 200   ! 241110
        do 500 i=1,n
	kf=k(i,2)
        do 600 j=1,kfmax
        if(kf.ne.kfaco(j))goto 600
        jj=numb(j)+1
c        if(iiii.eq.5 .and. iii.eq.13)write(22,*)'j,jj=',j,jj
c       update particle list etc.
	do m=nsa,jj,-1
        mm=m+1
c080104	ksa(mm,2)=ksa(m,2)
c080104	ksa(mm,1)=1
c080104	ksa(mm,3)=ksa(m,3)
        do m1=1,5
        ksa(mm,m1)=ksa(m,m1)   ! 080104
        psa(mm,m1)=psa(m,m1)
	vsa(mm,m1)=vsa(m,m1)
	enddo
        ishp(mm)=ishp(m)
        tau(mm)=tau(m)
	enddo
        do m=1,2000
        ipym=ipyth(m)
        if(ipym.ge.jj)ipyth(m)=ipym+1
        enddo
        if(ll.ge.jj)ll=ll+1
        if(ll1.ge.jj)ll1=ll1+1
c       update the values of lc(m,1-2) with value.ge.jj
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
c       give proper values to particle jj.
c221203	ksa(jj,2)=kf
c221203	ksa(jj,1)=1
c221203	ksa(jj,3)=0
        do m=1,5
        ksa(jj,m)=k(i,m)   ! 221203
        psa(jj,m)=p(i,m)
	vsa(jj,m)=v(i,m)
	enddo
        ishp(jj)=1
        tau(jj)=time+t0*p(i,4)/p(i,5)
c	the values of 'ishp' and 'tau' for hadrons from 'pythia' 
c	 are given here, the proper formation time of 'pythia' particle 
c	 is assume to be equal to t0 fm/c, except nucleon and j/psi
	if(kf.eq.2212 .or. kf.eq.2112)then
	tau(jj)=time+t0*p(i,4)/p(i,5)*taup
        elseif(kf.eq.443.or.kf.eq.30443)then
        tau(jj)=time+t0*p(i,4)/p(i,5)*taujp
        endif
	ipyth(i)=jj
	do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
	nsa=nsa+1
        goto 500
600     enddo   ! 040223
c040223 if produced hadron is not in given hadron classification
	nsa=nsa+1
        do m=1,5
        ksa(nsa,m)=k(i,m)
        psa(nsa,m)=p(i,m)
        vsa(nsa,m)=v(i,m)
        enddo
        ishp(nsa)=0
        tau(nsa)=0.
	ipyth(i)=nsa  
500	enddo   ! 040223
200     continue   ! 241110
c241110
c        if(iiii.eq.5)then
c        write(22,*)'updpip af remove iiii,iii,l,l1,nbh=',iiii,iii,l,l1,n
c        write(22,*)'numb=',(numb(i1),i1=1,52)
c        call prt_sa2(nsa,cc)
c        call prt_sbh(n,cc)
c        endif
c241110
        l=ll
        l1=ll1
c       remove colli. pair composed of l or l1  
	jj=0
	do 300 ii=1,nctl
	i1=lc(ii,1)
        j1=lc(ii,2)
        if(i1.eq.l .or. i1.eq.l1)goto 300
        if(j1.eq.l .or. j1.eq.l1)goto 300
	jj=jj+1
        tc(jj)=tc(ii)
        tw(jj)=tw(ii)
        do m=1,5
        lc(jj,m)=lc(ii,m)
        enddo	
300	continue
	do ii=jj+1,nctl+1
        tc(ii)=0.0
        tw(ii)=0.0
        do m=1,5
        lc(ii,m)=0
        enddo
        enddo
	nctl=jj
cm	write(9,*)'updpip af remove collis. pairs nctl=',nctl   ! sa
cm	do i=1,nctl
cm	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
cm	enddo
c	remove hadrons l and l1 from 'sa2'
	kf1=ksa(l,2)   
        kf2=ksa(l1,2)
c	write(9,*)'l,l1,kf1,kf2=',l,l1,kf1,kf2   !!
	kf=kf1
        ll=l
        do 700 i=1,2
        if(ll.eq.nsa)then   ! 
	do i1=1,kfmax
	if(kf.ne.kfaco(i1))goto 400
c241110	numbm=numb(i1)
c241110	do i2=1,i1
c241110	if(numb(i2).eq.numbm)numb(i2)=numb(i2)-1
c241110	enddo
        do m=i1,kfmax
        numb(m)=numb(m)-1
        enddo
        if(i1.gt.1)then
        numba=numb(i1)
        do m=1,i1-1
        if(numb(m).eq.numba)numb(m)=numb(m)-1
        enddo
        endif
c241110
	goto 100
400	enddo
        endif   !
        do j=ll+1,nsa
        jj=j-1
c080504	ksa(jj,2)=ksa(j,2)
c080504	ksa(jj,1)=1
c080504	ksa(jj,3)=ksa(j,3)
        do m=1,5
	ksa(jj,m)=ksa(j,m)   ! 080504
        psa(jj,m)=psa(j,m)
	vsa(jj,m)=vsa(j,m)
	enddo
        ishp(jj)=ishp(j)
	tau(jj)=tau(j)
	enddo
	if(nctl.eq.0)goto 900
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.gt.ll)lc(m,1)=lc1-1
        if(lc2.gt.ll)lc(m,2)=lc2-1
        enddo
900	do 800 j=1,kfmax
        if(kf.ne.kfaco(j))goto 800
        do m=j,kfmax
        numb(m)=numb(m)-1
        enddo
	if(j.gt.1)then
	numba=numb(j)
	do m=1,j-1
	if(numb(m).eq.numba)numb(m)=numb(m)-1
        enddo
	endif
        goto 100
800     continue
100     continue
        nsa=nsa-1
        if(l1.gt.ll)l1=l1-1
        do m=1,2000
        ipym=ipyth(m)
        if(ipym.gt.ll)ipyth(m)=ipym-1
        enddo
        if(i.eq.2)goto 700
        ll=l1
        kf=kf2
700     continue
c       if(((ipden.eq.0 .and. itden.eq.1) .or.
c       c   (ipden.eq.1 .and. itden.eq.0)) .and. nppb.ge.1)then
c       n=n-nppb
c       do i1=n+1,kszj
c       do i2=1,5
c       k(i1,i2)=0
c       p(i1,i2)=0.
c       v(i1,i2)=0.
c       enddo
c       enddo
c       endif
	return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coelas(ic,jc,eij,pi,pj)
c	perform elastic scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arraies in 'sa2' in this subroutine
	dimension pi(4),pj(4)
	iic=k(ic,2)
	jjc=k(jc,2)
	d=3.65*(eij-pmas(pycomp(iic),1)-pmas(pycomp(jjc),1))
	if(d.lt.1.e-10)return
	pt=0.2
	a=dmin1(10.3d0,1.d0/(1.12d0*pt)/(1.12d0*pt))
	d6=d**6
	b=d6*a/(1.+d6)
	if(b.lt.1.e-20)then
	b=1.e-20
	endif
	pm2=pi(1)**2+pi(2)**2+pi(3)**2
	pm=dsqrt(pm2)
	t0=-4.*pm2
	if(dabs(t0).lt.1.d-20)then
	cctas=1.
	goto 100
	endif
	cc=pyr(1)
	if(dabs(b*t0).lt.0.0001d0)then
	abt=1.
c	elseif(b*t0.lt.-50.)then
c	abt=0.
	else
	abt=dexp(dmax1(-7.0D2,dble(b*t0)))
	endif
	tt1=dlog(cc+(1.-cc)*abt)
	if(dabs(tt1).lt.1.d-30 .and. b.le.1.d-20)then
	cctas=1.
	goto 100
	endif
	tt=tt1/b
	if(dabs(tt).lt.1.d-20)then
	cctas=1.
	goto 100
	endif
	cctas=1.-tt*2./t0
	if(dabs(cctas).gt.1.d0)then
	cctas=dsign(1.d0,cctas)
	endif
100	continue
	sctas=dsqrt(1.-cctas**2)
	fis=2.*3.1416*pyr(1)
	cfis=dcos(fis)
	sfis=dsin(fis)
	call rotate(cctas,sctas,cfis,sfis,pm,pi,pj)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rotate(cctas,sctas,cfis,sfis,pp3,pi,pj)
c	perform rotation
c       pi,pj: input, four momentum of colliding pair before scattering
c              output,four momentum of scattered particles after rotation
c       pp3: momentum modulus of pi or pj, both are equal in their cms,
c        after scattering
c       cctas,sctas,cfis,sfis: direction cosines of momentum of one of 
c        scattered particle relative to the momentum
c        of corresponding particle before scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	dimension pi(4),pj(4)
c	fi1=atan2(pi(2),pi(1))
c	cta1=atan2(dsqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1=pyangl(pi(1),pi(2))
	cta1=pyangl(pi(3),dsqrt(pi(1)**2+pi(2)**2))
	cfi1=dcos(fi1)
	sfi1=dsin(fi1)
	ccta1=dcos(cta1)
	scta1=dsin(cta1)
	pi(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
	pi(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
	pi(3)=ccta1*cctas-scta1*sctas*cfis
	pi(1)=pp3*pi(1)
	pi(2)=pp3*pi(2)
	pi(3)=pp3*pi(3)
	do i=1,3
	pj(i)=0.-pi(i)
	enddo
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updple(ic,jc,b,pi,pj)
c	update particle list for elastic scattering 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
c       note the name of the arrays in 'sa2'
	dimension pi(4),pj(4),b(3)
	ilo=1
c       ilo=1 for inverse Lorentz transformation
	call lorntz(ilo,b,pi,pj)
	do i=1,4
	p(ic,i)=pi(i)
	p(jc,i)=pj(i)
	enddo
	return
	end



c********************************************************************
        subroutine tran_saf
c       'saf' to 'pyjets' 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 080104
        do l=1,nsa
        do m=1,5
        k(l,m)=ksa(l,m)
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        enddo
        n=nsa
	do l=n+1,kszj
	do m=1,5
	k(l,m)=0
	p(l,m)=0.
	v(l,m)=0.
	enddo
	enddo
        return
        end

c********************************************************************
        subroutine tran_sbe
c       'sbe' to 'pyjets' 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        PARAMETER (kszj=80000)
        COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbe/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)   ! 080104
        do l=1,nsa
        do m=1,5
        k(l,m)=ksa(l,m)
        p(l,m)=psa(l,m)
        v(l,m)=vsa(l,m)
        enddo
        enddo
        n=nsa
	do l=n+1,kszj
	do m=1,5
	k(l,m)=0
	p(l,m)=0.
	v(l,m)=0.
	enddo
	enddo
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	subroutine conse(np,pp,ps,ii,jj)
c	keep four momentum conservation
c	np : the # of particles
c	ps : four momentum to which the four momenta of particles should 
c            conserve
c	pp : four momenta of particles
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
      	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	dimension pp(250,5),ps(4),ff(250),pxyz(3),arp(3)
	ps4=ps(4)
	do i=1,3
	pxyz(i)=0.
	enddo
	jj=0
100	es=0.
	do i=1,np
	es=es+pp(i,4)
	enddo
	fr=es/ps4
	if(dabs(1.-fr) .le. dep)goto 200
	do i=1,np
	ppm=pp(i,4)/0.938
	ppf=ppm/fr
	ff(i)=dsqrt(dabs(ppf*ppf-1.)/(ppm*ppm-1.))
	do j=1,3
	ppp=ff(i)*pp(i,j)
	pp(i,j)=ppp
	pxyz(j)=pxyz(j)+ppp
	enddo
	enddo
	do i=1,3
	arp(i)=dabs(1.-pxyz(i)/ps(i))
	pxyz(i)=pxyz(i)-ps(i)
	enddo
	if(dabs(1.-fr).le.dep .and.arp(1).le.dep .and. arp(2).le.dep  
     c   .and. arp(3).le.dep) goto 200
	do i=1,3
	pxyz(i)=pxyz(i)/np
	enddo
	do i=1,np
	do j=1,3
	pp(i,j)=pp(i,j)-pxyz(j)
	enddo
	pp(i,4)=dsqrt(0.880+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
c	0.880 = 0.938*0.938
	enddo
	jj=jj+1
	if(jj.eq.4000)then
	write(9,*)'infinitive loop may occur in subroutine conse(),'
        write(9,*)'which means four-momentum conservation' 
        write(9,*)'needed is hard to be achieved,check value' 
        write(9,*)'of PARAM(9)'
	return
	endif
	goto 100
200	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine codi(pis,cfi1,sfi1,ccta1,scta1)
c	calculate the 'orientation' of the vector pis
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	dimension pis(4),pi(4)
c	do i=1,4
c	pi(i)=pis(i)
c	enddo
c	if(pi(1).lt.1.d-15)pi(1)=1.d-15 
c	fi1=atan2(pi(2),pi(1))
c        cta1=atan2(dsqrt(pi(1)**2+pi(2)**2),pi(3))
	fi1s=pyangl(pis(1),pis(2))
        cta1s=pyangl(pis(3),dsqrt(pis(1)**2+pis(2)**2))
	fi1=fi1s
	cta1=cta1s
        cfi1=dcos(fi1)
        sfi1=dsin(fi1)
        ccta1=dcos(cta1)
        scta1=dsin(cta1)
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rosa(cfi1,sfi1,ccta1,scta1,cfis,sfis,cctas,sctas,
     c	 pis)
c       perform rotate for produced particles from 'pythia'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        dimension pis(4),pi(4)
	do i=1,4
	pi(i)=pis(i)
	enddo
	pp=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
	pp=dsqrt(pp)
	call codi(pis,cfis,sfis,cctas,sctas)
        pi(1)=pp*(cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis)
        pi(2)=pp*(sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis)
        pi(3)=pp*(ccta1*cctas-scta1*sctas*cfis)
	do i=1,4
	pis(i)=pi(i)
	enddo
        return
        end



cLei20220223ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bp_prob(bp,r_max,prob)
c       impact parameter density distribution f(b)db=bdb  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      
      b = bp
      
!     integrates f(b) over b, from 0 to b_max. (r_max here)
          sum_b2_max = 0.5 * r_max * r_max
      
!     integrates f(b_i) over b_i, from 0 to b.
      sum_b2 = 0.5 * b * b
      
      prob = (1. - sum_b2 / sum_b2_max)
      
      return
      end 
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine recons(irecon,l,l1,ss,nlead,time,lll)    ! 281121
c161021 a model to reconstruct diquark-quark (quark-diquark) 'A and V' pair 
c        into proton directely to increase leading proton effect
c       l and l1: line # 0f current nucleon-nucleon collision pair in 'sa2',
c        ss: total energy of that pair
c       lll: # of loops over nucleon-nucleon collision in a
c        nucleus-nucleus collision
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sppb/nppb,non3,kppb(1000,5),pppb(1000,5),vppb(1000,5) ! 281121
	dimension ps(4),rs(4),pp(20,5),isuc(1000)! 230407
	dimension pk(5),vk(5),rr(3),kk(5)
	delte=0.
        deltx=0.
        delty=0.
        deltz=0.
        imc=adj1(13)
        ibc=adj1(14)
        do j1=1,1000
        isuc(j1)=0
        enddo

c150322 remove junctions
        jb=0
2010    do i1=jb+1,n  ! i1 loop      
        kf=k(i1,2)
        if(kf.ne.88)then
        jb=jb+1
        goto 2020
        endif
c       move particle list 'pyjets' one step downward since i1+1 to n
        do j=i1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 2010
2020    enddo   ! i1 loop        
c150132        
c       write(9,*)'enter recons iiii,lll,nsa,nlead=',iiii,lll,nsa,nlead !!

600	jjj=0
        jjjj=0   ! 161021, # of reconstracted nucleon (anti-nucleon) 
	iii=0   ! 281121 from 1 -> 0    
c	find out string composed of dd_1,uu_1,ud_0,ud_1,d,dbar,u,or ubar 
        do 500 i2=iii+1,n   ! 281121 from iii -> iii+1   ! 2 
	kf=k(i2,2)
        kfab=iabs(kf)
        if(kfab.ne.1103 .and. kfab.ne.2203 .and. kfab.ne.2101
     c   .and. kfab.ne.2103 .and. kfab.gt.2)then   ! composed of u, d only
	iii=iii+1
	goto 500
	endif
	k1=k(i2,1)
	if(k1.eq.2)then   ! k1=2 means 'A'  if 1
	do 501 i3=i2+1,n   ! 3
	kf4=k(i3,2)
	kf4ab=iabs(kf4)
	k2=k(i3,1)
	if(k2.eq.1.and.((kfab.le.2.and.(kf4ab.eq.1103.or.kf4ab.eq.2203
     c   .or.kf4ab.eq.2101.or.kf4ab.eq.2103)).or.(kf4ab.le.2.and.(kfab
     c   .eq.1103.or.kfab.eq.2203.or.kfab.eq.2101.or.
     c   kfab.eq.2103))))then   ! k2=1 means 'V'  230407 if 2
	jjj=jjj+1   ! # of found diquark-quark (quark-diquark) 'A-V' pair
	p1x=p(i2,1)
        p1y=p(i2,2)
        p1z=p(i2,3)
	p1e=p(i2,4)
	p2x=p(i3,1)
        p2y=p(i3,2)
        p2z=p(i3,3)
	p2e=p(i3,4)
	p12x=p1x+p2x
	p12y=p1y+p2y
	p12z=p1z+p2z
	p12e=p1e+p2e
	cm2=p12e*p12e-p12x*p12x-p12y*p12y-p12z*p12z
	if(cm2.gt.1.d40)cm2=1.d40
	if(cm2.lt.1.d-40)cm2=1.d-40
	cm=dsqrt(cm2) 
c16101  cm: invariant mass of found diquark-quark (quark-diquark) 'A-V' pair

	if(kf.gt.10)then   ! i2 is diquark
	kfbb=kf/1000
	kf1=kfbb
	kf2=(kf-kfbb*1000)/100
	kf3=kf4
	sdir=dsign(1d0,p(i2,3))   ! sign of third momentm of i2
	else   ! i2 is quark
        kfbb=kf4/1000
	kf1=kfbb
        kf2=(kf4-kfbb*1000)/100
        kf3=kf
	sdir=dsign(1d0,p(i3,3))   ! sign of third momentm of i3 
	endif
c161021 cf. pythia manual p.71 for flavor code of diquark

c   found diquark-quark (quark-diquark) 'A-V' pair can be p (pbar,n,or nbar) ?
	call tabhb
	if(kf1.gt.0.and.kf2.gt.0.and.kf3.gt.0)then
	call findb(kf1,kf2,kf3,cm,kfii,amasi,isucc,1)
	elseif(kf1.lt.0.and.kf2.lt.0.and.kf3.lt.0)then
	call findb(-kf1,-kf2,-kf3,cm,kfii,amasi,isucc,-1)   ! 020605 Tan
	else
	endif
	kiab=iabs(kfii)
c281121 if(isucc.eq.1 .and. (kiab.eq.2212 .or. kiab.eq.2112))then ! if 3 161021
        if(isucc.eq.1 .and. kiab.eq.2212)then ! if 3 161021
        jjjj=jjjj+1   ! 161021, # of reconstracted nucleon (anti-nucleon)
c       write(9,*)'# of reconstracted nucleon,jjjj,nlead=',jjjj,nlead !!
c140322 if(jjjj.gt.nlead)return
        if(jjjj.gt.1)return   ! 140322
	isuc(jjj)=1

c	set reconstracted nucleon (anti-nucleon) on line i2 temporarily
c        and give proper variables to that nucleon (anti-nucleon)
	k(i2,1)=1   
	k(i2,2)=kfii
	k(i2,3)=0
	p(i2,5)=amasi
c281121 
        p(i2,1)=0.
        p(i2,2)=0.
        p(i2,4)=0.5*ss
        emp=p(i2,4)*p(i2,4)-amasi*amasi
        if(emp.gt.1.d40)emp=1.d40
        if(emp.lt.1.d-40)emp=1.d-40
        p(i2,3)=dsqrt(emp)   ! pA
        if(ipden.eq.1 .and. itden.eq.0)p(i2,3)=-p(i2,3)   ! Ap
        delte=delte+(p12e-p(i2,4))
        deltx=deltx+p12x
        delty=delty+p12y
        deltz=deltz+(p12z-p(i2,3))
        rl=pyr(1)
        do m=1,3
        v(i2,m)=v(i2,m)+vsa(l,m)*rl+vsa(l1,m)*(1.-rl)
        enddo
        v(i2,4)=time+2.*ddt

c       'pyjets' to 'sppb'        
        nppb=nppb+1
        do ii=1,5
        kppb(nppb,ii)=k(i2,ii)
        pppb(nppb,ii)=p(i2,ii)
        vppb(nppb,ii)=v(i2,ii)
        enddo
c       write(9,*)'jjjj,nppb',jjjj,nppb   !!

c161021 i2+1 ('V')->i3 ('V')
        do jj=1,5
        k(i3,jj)=k(i2+1,jj)
        p(i3,jj)=p(i2+1,jj)
        v(i3,jj)=v(i2+1,jj)
        enddo
c       new (i3-1)-i3 is a 'A-V' pair        

c161021 move 'pyjets',one step downward since i2+2 to n (i.e. throw away i2+1)
        do j=i2+2,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1

c161021 move 'pyjets',one step downward since i2+1 to n (i.e. throw away i2)
        do j=i2+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1

	goto 400   ! success
c281121 goto 888   ! 240805
	endif   ! if 3
c281121 888	iii=i3+1   ! 240805
	endif   ! if 2
c161021 if(k2.eq.1)then
c161021 iii=i3+1
c161021 goto 400
c161021 endif   
c281121 fail, procced
501     enddo   ! 3
	endif   ! if 1
400	continue
500	enddo   ! 2

        goto 889
c	share energy in 'delte' to particles in 'pyjets'
       if(n.gt.0)then

        del=delte   ! 161021
        del=del/dfloat(n)
c       write(9,*)'share n,del=',n,del
        do j3=1,n
        p(j3,4)=p(j3,4)+del
        if(del.lt.0.)then
        if(p(j3,4).lt.0.)p(j3,4)=p(j3,4)-del
        pabs=dabs(p(j3,3))
        if(pabs.ge.p(j3,4))p(j3,4)=p(j3,4)-del
        endif
        enddo
c281121
        del=deltx
        del=del/dfloat(n)
        do j3=1,n
        p(j3,1)=p(j3,1)+del
        enddo
        del=delty
        del=del/dfloat(n)
        do j3=1,n
        p(j3,2)=p(j3,2)+del
        enddo
        del=deltz
        del=del/dfloat(n)
        do j3=1,n
        p(j3,3)=p(j3,3)+del
        enddo
c281121 

        endif        
889     continue
c	write(22,*)'out of recons'
c	call pylist(1)
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine recons_g(irecon,l,l1,ss,nlead,time,lll)   ! 281121
c161021 a model to reconstruct diquark-quark (quark-diquark) 'A and V' pair 
c        into proton to increase leading proton effect
c       break-up gluon randumly first 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sppb/nppb,non3,kppb(1000,5),pppb(1000,5),vppb(1000,5) ! 281121
	dimension ps(4),rs(4),pp(20,5),isuc(1000)!230407
	dimension iiglu(kszj),pk(5),vk(5),rr(3),kk(5)
c	isuc(i): =1 if i-th quark-diquark 'A and V' pair can compose into
c	 nucleon, Delta(0), and Delta(+), otherwise =0
c	iiglu(i): line number of i-th gluon in 'pyjets'
        do j1=1,kszj
	iiglu(j1)=0
	enddo
        amd=pymass(1)   ! 161021 (constituant quark mass)
        amu=pymass(2)
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

c	count number of gluons
	jjj=0
	do j=1,n
	ik=k(j,2)
	if(ik.eq.21)then
	jjj=jjj+1
	iiglu(jjj)=j   ! line number of jjj-th gluon
	endif
	enddo
	jglu=jjj   ! number of gluons
c	write(9,*)'jglu,iiglu=',jglu,(iiglu(j),j=1,jglu)

	if(jglu.eq.0)goto 600

        if(jglu.eq.1)then   ! 1
c	break that gluon forcely
	ii1=iiglu(1)
	ps(1)=p(ii1,1)
	ps(2)=p(ii1,2)
	ps(3)=p(ii1,3)
	ps(4)=p(ii1,4)
	rs(1)=v(ii1,1)
	rs(2)=v(ii1,2)
	rs(3)=v(ii1,3)
	rs(4)=v(ii1,4)
	eg=ps(4)   ! 110322

        if(eg.lt.amdd)then   ! 2 thrown away that gluon
        delte=delte+eg   
c       move particle list,'pyjets',one step downward since ii1+1
        do j=ii1+1,n
        do jj=1,5
        k(j-1,jj)=k(j,jj)
        p(j-1,jj)=p(j,jj)
        v(j-1,jj)=v(j,jj)
        enddo
        enddo
        n=n-1
        goto 600
        endif

        call break_f(eg,kf,amq)   ! which is in coales_23.h 161022
        kf1=kf 
        kf2=-kf
        am1=amq
        am2=amq

c	exchange that gluon with the parton ahead
	j1=ii1+1
        do jj=1,5
	kk(jj)=k(j1,jj)
	pk(jj)=p(j1,jj)
	vk(jj)=v(j1,jj)
        k(j1,jj)=k(ii1,jj)
        p(j1,jj)=p(ii1,jj)
        v(j1,jj)=v(ii1,jj)
	k(ii1,jj)=kk(jj)
	p(ii1,jj)=pk(jj)
	v(ii1,jj)=vk(jj)
	enddo
c161021 j1 is now line # of that gluon in 'pyjets'
c	write(22,*)'af. exchange event=',iiii
c	call pylist(1)

c       move particle list,'pyjets',one step forward since j1+1
        do j=n,j1+1,-1
        j2=j+1
        do jj=1,5
        k(j2,jj)=k(j,jj)
        p(j2,jj)=p(j,jj)
        v(j2,jj)=v(j,jj)
        enddo
        enddo
        n=n+1
c	write(22,*)'af. move forward event=',iiii
c	call pylist(1)

c       assume the breaked q and qbar as a string
	k(j1,1)=2   ! A
	k(j1,2)=kf1
	k(j1,3)=0
        k(j1+1,1)=1   ! V
        k(j1+1,2)=kf2
        k(j1+1,3)=0

c       give four momentum to breaked quarks

c110322 subtract 2*am1 (am1: mass of splited quark) from gluon energy
c        and reduce gluon three momentum correspondingly 
        ps4o=ps(4)
        ps(4)=ps4o-2.*am1
        rati=ps(4)/ps4o   ! times of g energy reduction
        ps(1)=ps(1)*rati
        ps(2)=ps(2)*rati
        ps(3)=ps(3)*rati
c110322
	decsuc=1   ! c1
	call decmom(ps,pp,am1,am2,decsuc)   ! c1
c161021 pp(1,5):four momenta and mass of breaked q (obtained from decmom)
c161021 pp(2,5):four momenta and mass of breaked qbar (obtained from 
c        decmom)        
c	as mass of gluon from 'pyjets' may be negative it may be better
c	 (from energy conservation point of view) not using 'decmom' but
c	 random three momentum method if square root s less than 0.1
	if(decsuc.eq.0)then   ! c1
        do i4=1,3  
        pi=pyr(1)*ps(i4)
        pp(2,i4)=ps(i4)-pi
        pp(1,i4)=pi
        enddo
        pp11=pp(1,1)
        pp12=pp(1,2)
        pp13=pp(1,3)
	ampp=am1*am1+pp11*pp11+pp12*pp12+pp13*pp13
	if(ampp.gt.1.d40)ampp=1.d40
	if(ampp.lt.1.d-40)ampp=1.d-40
        pp(1,4)=dsqrt(ampp)
        pp21=pp(2,1)
        pp22=pp(2,2)
        pp23=pp(2,3)
	ampp=am2*am2+pp21*pp21+pp22*pp22+pp23*pp23
	if(ampp.gt.1.d40)ampp=1.d40
        if(ampp.lt.1.d-40)ampp=1.d-40
        pp(2,4)=dsqrt(ampp)
	endif   ! c1
	p(j1,1)=pp(1,1)
	p(j1,2)=pp(1,2)
	p(j1,3)=pp(1,3)
	p(j1,4)=pp(1,4)
        p(j1,5)=am1   ! 161021
	p(j1+1,1)=pp(2,1)
	p(j1+1,2)=pp(2,2)
	p(j1+1,3)=pp(2,3)
	p(j1+1,4)=pp(2,4)
        p(j1+1,5)=am2   ! 161021
c       give four coordinate to breaked quarks
c        first breaked quark takes the four coordinate of diquark
c        second breaked quark is arranged around first ones within
c        0.5 fm randumly in each of three coordinates and has same
c        fourth coordinate as diquark
	v(j1,1)=rs(1)
	v(j1,2)=rs(2)
	v(j1,3)=rs(3)
	v(j1,4)=rs(4)
	do j4=1,3
	rr(j4)=pyr(1)*0.5
	v(j1+1,j4)=rs(j4)+rr(j4)
	if(pyr(1).gt.0.5)v(j1+1,j4)=rs(j4)-rr(j4)
	enddo
	v(j1+1,4)=rs(4)

	delte=delte+(eg-p(j1,4)-p(j1+1,4))   ! 161021 
c	write(9,*)'ii1,kf1,kf2,eg,delte=',iiglu(1),kf1,kf2,eg,delte
c	write(9,*)'pp(1,)=',(pp(1,ii),ii=1,4)
c	write(9,*)'pp(2,)=',(pp(2,ii),ii=1,4)
c	write(9,*)'sum=',(pp(1,ii)+pp(2,ii),ii=1,4)
c	write(9,*)'ps=',(ps(ii),ii=1,4)
c	write(9,*)'v(j1,)=',(rs(ii),ii=1,4)
c	write(22,*)'af. break g event=',iiii
c	call pylist(1)
	goto 600
	endif   ! 1

c161021 proceed for case of gluon # > 1        
c	move particle list, 'pyjets', jglu steps forward since n to 1
800	do j=n,1,-1
	j1=j+jglu
	do jj=1,5
	k(j1,jj)=k(j,jj)
        p(j1,jj)=p(j,jj)
        v(j1,jj)=v(j,jj)
	enddo
	enddo
	n=n+jglu
	do j1=1,jglu
	iiglu(j1)=iiglu(j1)+jglu
	enddo

c	move g to the beginning of 'pyjets'
	jjj=0
	do j=jglu+1,n   ! do 1
	ik=k(j,2)
	if(ik.eq.21)then
	jjj=jjj+1
	do jj=1,5
	k(jjj,jj)=k(j,jj)
	p(jjj,jj)=p(j,jj)
        v(jjj,jj)=v(j,jj)
	enddo
	endif
	enddo   ! do 1

c161021 delete original gluon lines in 'pyjets'
	do j2=1,jglu   ! do 2
	j11=iiglu(j2)
	do j1=j11+1,n
	do jj=1,5
        k(j1-1,jj)=k(j1,jj)
        p(j1-1,jj)=p(j1,jj)
        v(j1-1,jj)=v(j1,jj)
        enddo
	enddo
	n=n-1
	if(j2.lt.jglu)then
	do j3=j2+1,jglu
	iiglu(j3)=iiglu(j3)-1
	enddo
	endif
	enddo   ! do 2
c090505
c        mglu=jglu-1   ! jglu is odd
c        if(mod(jglu,2).eq.0)mglu=jglu   ! jglu is even
c        do j=1,mglu,2
c        k(j,1)=2   ! A
c        k(j1+1,1)=1   ! V
c        enddo
c        if(mod(jglu,2).ne.0)then
c        k(jglu-2,1)=2   ! A
c        k(jglu-1,1)=2   ! A
c        k(jglu,1)=1   ! V
c        endif
c090505

c161021 treat those gluons as a string
	k(jglu,1)=1   ! V
c161021 note: gluon in 'pyjest' always has k(i,1)=2 (i=1,2,...,jglu), 
c        so k(1,1)=2   ! A
c	write(22,*)'jglu,jjj,k(jglu,1),event=',jglu,jjj,k(jglu,1),iiii
c	call pylist(1)

600     continue
c161021 all are simple string (without gluon) upto here
        call recons(irecon,l,l1,ss,nlead,time,lll)
c	write(22,*)'out of recons'
c	call pylist(1)
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine recons_gg(irecon,l,l1,ss,nlead,time,lll)   ! 281121
c161021 a model to reconstruct diquark-quark (quark-diquark) 'A and V' pair 
c        into proton to increase leading proton effect
c       move gluons to 'gu' first and throw gluons randumly into strings 
c        at last
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYSUBS/MSEL,MSUB(500),KFIN(2,-40:40),NON,CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa24/adj1(40),nnstop,non24,zstop
        common/sa4_c/kqh(80,2),kfh(80,2),proh(80,2),amash(80,2),imc
        common/sa5_c/kqb(80,3),kfb(80,2),prob(80,2),amasb(80,2),ibc
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        common/sppb/nppb,non3,kppb(1000,5),pppb(1000,5),vppb(1000,5) ! 281121
        dimension pgu(kszj,5),vgu(kszj,5),kgu(kszj,5)
        dimension nstra(kszj),nstrv(kszj)
	dimension ps(4),rs(4),pp(20,5),isuc(1000) ! 230407
	dimension pk(5),vk(5),rr(3),kk(5)

c       move gluon from 'pyjets' to 'gu' & count number of gluons
        ngu=0
        do jj=1,n   ! do
        ik=k(jj,2)
        if(ik.eq.21)then   ! if 
        ngu=ngu+1
        do j1=1,5
        kgu(ngu,j1)=k(jj,j1)
        pgu(ngu,j1)=p(jj,j1)
        vgu(ngu,j1)=v(jj,j1)
        enddo
c       move particle list 'pyjets' one step downward since jj+1
        do j=jj+1,n
        do j1=1,5
        k(j-1,j1)=k(j,j1)
        p(j-1,j1)=p(j,j1)
        v(j-1,j1)=v(j,j1)
        enddo
        enddo
        n=n-1
        endif   ! if 
        enddo   ! do        

c161021 all are simple string (without gluon) upto here

        call recons(irecon,l,l1,ss,nlead,time,lll)

c       find string & line number of its first ('A') & last ('V')) components
        nstr=0
        jj=0
503     do j1=jj+1,n
        if(k(j1,1).eq.2)then   ! j1 is 'A'
        do j2=j1+1,n
        if(k(j2,1).eq.1)then   ! j2 is 'V'
        nstr=nstr+1
        nstra(nstr)=j1   ! line number of first component of nstr-th string
        nstrv(nstr)=j2   ! line number of last component of nstr-th string
        if(j2.eq.n)then
        goto 504
        else
        jj=j2
        goto 503
        endif
        endif   ! j2
        enddo
        endif   ! j1
        enddo
504     continue

c       arrange gluons into string randumly
        do 502 j1=1,ngu
        j2=int(nstr*pyr(1)+1)   ! j1-th gluon fall in j2-th string
        jj=nstrv(j2)   ! line # of 'V' in j2-th string
c       move 'pyjets' one step forward from n down to jj
        do j=n,jj,-1
        do j3=1,5
        k(j+1,j3)=k(j,j3)
        p(j+1,j3)=p(j,j3)
        v(j+1,j3)=v(j,j3)
        enddo
        enddo
        n=n+1
c       update string
        do j3=j2,nstr
        if(j3.eq.j2)then
        nstrv(j3)=nstrv(j3)+1
        elseif(j3.gt.j2)then
        nstra(j3)=nstra(j3)+1
        nstrv(j3)=nstrv(j3)+1
        else
        endif
        enddo
c       'gu' to 'pyjets'
        do j3=1,5
        k(jj,j3)=kgu(j1,j3)
        p(jj,j3)=pgu(j1,j3)
        v(jj,j3)=vgu(j1,j3)
        enddo
502     enddo                

c	write(22,*)'out of recons'
c	call pylist(1)
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine remo   ! 110517 010418 040223
c       moves hadron (lepton) and junction from 'pyjets' to 'sbh'   060813
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
        common/sa24/adj1(40),nnstop,non24,zstop   ! 060223
	nbh=0
	do i1=1,kszj
        do j1=1,5
        kbh(i1,j1)=0
        pbh(i1,j1)=0.
        vbh(i1,j1)=0.
        enddo
        enddo
        adj12=adj1(12)   ! 060223
        jb=0
201     do 202 i1=jb+1,n   ! do loop 090122 060223
        kf=k(i1,2)
        kfab=iabs(kf)
c060223        

        if(adj12.eq.0)then   ! keep junction in 'pyjets'
        if(kfab.le.8 .or. kfab.eq.2101 .or. kfab.eq.3101
     c   .or. kfab.eq.3201 .or. kfab.eq.1103 .or. kfab.eq.2103
     c   .or. kfab.eq.2203 .or. kfab.eq.3103 .or. kfab.eq.3203
     c   .or. kfab.eq.3303 .or. kfab.eq.21 .or. kfab.eq.88)then
        jb=jb+1
        goto 202
        endif
        goto 204
        endif   !

        if(adj12.eq.1 .or. adj12.eq.2)then   !! 
c100223  remove junction from 'pyjets' to 'sbh'
        if(kfab.le.8 .or. kfab.eq.2101 .or. kfab.eq.3101
     c   .or. kfab.eq.3201 .or. kfab.eq.1103 .or. kfab.eq.2103
     c   .or. kfab.eq.2203 .or. kfab.eq.3103 .or. kfab.eq.3203                
     c   .or. kfab.eq.3303 .or. kfab.eq.21)then
c060223        
        jb=jb+1
        goto 202
        endif
        endif   !!
c	write(9,*)'n,i1,jb=',n,i1,jb   ! sa
204     continue   ! 060223
        nbh=nbh+1
        do i2=1,5
        kbh(nbh,i2)=k(i1,i2)
        pbh(nbh,i2)=p(i1,i2)
        vbh(nbh,i2)=v(i1,i2)
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
        goto 201
202     enddo   ! do loop 090122
203     continue
	return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine break
c       breaks up diquark (anti-diquark), gives four momenta 
c	 and four positions to the broken objects
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa24/adj1(40),nnstop,non24,zstop   ! 170205
        common/sa26/ndiq(kszj),npt(kszj),ifcom(kszj),idi,idio   ! 080104 220110
        common/saf/naf,nonaf,kaf(kszj,5),paf(kszj,5),vaf(kszj,5)
	jb=0
	ii=idio   ! 080104

100     do i1=jb+1,n
        kf=k(i1,2)
	kfab=iabs(kf)
        if(kfab.ne.2101 .and. kfab.ne.3101
     c   .and. kfab.ne.3201 .and. kfab.ne.1103 .and. kfab.ne.2103
     c   .and. kfab.ne.2203 .and. kfab.ne.3103 .and. kfab.ne.3203
     c   .and. kfab.ne.3303)then
c     c   .and. kfab.ne.3303 .and. kfab.ne.21)then
        jb=jb+1
        goto 300
        endif

	if(kf.eq.2101)then
	kf1=2
	kf2=1
	goto 200
	endif
	if(kf.eq.3101)then
        kf1=3
        kf2=1
        goto 200
	endif
	if(kf.eq.3201)then
        kf1=3
        kf2=2
        goto 200
	endif
	if(kf.eq.1103)then
        kf1=1
        kf2=1
        goto 200
	endif
	if(kf.eq.2103)then
        kf1=2
        kf2=1
        goto 200
	endif
	if(kf.eq.2203)then
        kf1=2
        kf2=2
        goto 200
	endif
	if(kf.eq.3103)then
        kf1=3
        kf2=1
        goto 200
	endif
	if(kf.eq.3203)then
        kf1=3
        kf2=2
        goto 200
	endif
	if(kf.eq.3303)then
        kf1=3
        kf2=3
        goto 200
	endif
c251103
	if(kf.eq.-2101)then
	kf1=-2
	kf2=-1
	goto 200
	endif
	if(kf.eq.-3101)then
        kf1=-3
        kf2=-1
        goto 200
	endif
	if(kf.eq.-3201)then
        kf1=-3
        kf2=-2
        goto 200
	endif
	if(kf.eq.-1103)then
        kf1=-1
        kf2=-1
        goto 200
	endif
	if(kf.eq.-2103)then
        kf1=-2
        kf2=-1
        goto 200
	endif
	if(kf.eq.-2203)then
        kf1=-2
        kf2=-2
        goto 200
	endif
	if(kf.eq.-3103)then
        kf1=-3
        kf2=-1
        goto 200
	endif
	if(kf.eq.-3203)then
        kf1=-3
        kf2=-2
        goto 200
	endif
	if(kf.eq.-3303)then
        kf1=-3
        kf2=-3
        goto 200
	endif
c251103
200	k(i1,2)=kf1
        k(n+1,2)=kf2
c221203	k(i1,1)=1
        k(n+1,1)=1
c221203	k(i1,3)=0
        k(n+1,3)=0
c221203
        k(n+1,4)=0
        k(n+1,5)=0
c221203
c080104
	ii=ii+1
        npt(ii)=n+1+naf
        ifcom(ii)=i1+naf   ! 220110
c080104
c       give four momentum to the breaked quarks
	call bream(i1,kf1,kf2)
c       give four coordinate to the breaked quarks
        call coord(i1)
	if(i1.eq.n)then
	n=n+1
	goto 400
	endif
        n=n+1
        goto 100

300     enddo
400	return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bream(ii,kf1,kf2)
c       give four momentum to the broken quarks
c       ii: line number of diquark in 'pyjets'
c       kf1,kf2: flavor codes of broken quarks
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   ! 260503
        am1=pymass(kf1)
        am2=pymass(kf2)
        pp(1,5)=am1
        pp(2,5)=am2
c       pp : four momenta & mass of broken quarks, local variable 
        do i1=1,4
        ps(i1)=p(ii,i1)
        enddo
c       ps : four momentum of diquark, local variable 
	goto 400   ! activate it for 'decay method'
c       broken quarks share out diquark four momentum randomly,
c        denoted as 'random four momentum method'
c       do i1=1,4   ! activate it for 'random four momentum method'
c       broken quarks share out diquark three momentum randomly,
c        denoted as 'random three momentum method'
401	do i1=1,3   ! activate it for 'random three momentum method'
        pi(i1)=pyr(1)*p(ii,i1)
        pp(2,i1)=ps(i1)-pi(i1)
        pp(1,i1)=pi(i1)
        enddo
c	goto 300   ! activate it for 'random four momentum method'
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
400	continue
	decsuc=1
	call decmom(ps,pp,am1,am2,decsuc)
	if(decsuc.eq.0)goto 401   ! return to random three momentum method
300	continue
c       adjust four momentum conservation by iteration,no more than
c        4000 iterations
c	call conser(2,pp,ps)   
c260503
c        write(9,*)'after	'
c        do i=1,2
c        write(9,*)(pp(i,j),j=1,5)
c        enddo
c260503
        do i1=1,4
        p(ii,i1)=pp(1,i1)
        enddo
        p(ii,5)=am1
        do i1=1,4
        p(n+1,i1)=pp(2,i1)
        enddo
        p(n+1,5)=am2
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decmom(ps,pp,am1,am2,decsuc)
c	calculate four momentum of decayed particles
c	ps: four momentum of decaying particle
c       pp: four momentum of decayed particles   ! 090922
c	am1 (am2): mass of decayed particle   ! 090922
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension pi(4),pj(4),ps(4),pp(20,5),bb(3)   
c       calculate the E and |p| of broken quark in rest frame of diquark
        sm2=ps(4)*ps(4)-ps(1)*ps(1)-ps(2)*ps(2)-ps(3)*ps(3)
c       one problem here is that 'sm2' may not equal to square of diquark 
c	 (gluon) rest mass,'bream' is called for spliting g especially
c030603
c1	if(sm2.lt.1.d-10)then
c1	sm2=1.d-10
c1	endif
c       write(9,*)'in decmom sm2=',sm2   ! sa
c       write(9,*)'ps=',ps   ! sa
c030603
	if(sm2.lt.0.005)then   ! 110211
	decsuc=0   ! go back to random three momentum method
	return
	endif
        sm=dsqrt(sm2)   ! M (should be diquark mass)
c       pp(1,4)=(sm2-am2*am2+am1*am1)/2./sm
c       pp(2,4)=(sm2-am1*am1+am2*am2)/2./sm
        ppp=(sm2-(am1+am2)*(am1+am2))*(sm2-(am1-am2)*(am1-am2))
c161204	ppp=dabs(ppp)   ! 030603 ?
	if(ppp.lt.1.d-28)ppp=1.d-28   !161204
        ppp=dsqrt(ppp)/2./sm
c110211 goto 500   ! activate it for exponential cos(seta) distribution
c       the direction of broken quark is sampled isotropically in '4pi'
        coset=1.-2.*pyr(1)
        if(dabs(coset).gt.1.d0)then
        coset=coset/dabs(coset)
        endif
c021005
        siset=1.-coset*coset
        if(siset.lt.1.d-28)siset=1.d-28
c021005
        siset=dsqrt(siset)   ! 021005
100     cosi1=pyr(1)
        cosi12=cosi1*cosi1
        eta2=2.*pyr(1)-1.
        eta22=eta2*eta2
        coseta=cosi12+eta22
        if(coseta.gt.1.)goto 100
        if(coseta.lt.1.d-28)coseta=1.d-28
        cofi=(cosi12-eta22)/coseta
        sifi=2.*cosi1*eta2/coseta
        goto 600
500     continue
c       cos(seta) is sampled from exponential distribution when
c        0<seta<pi/2 and its absolute value is assumed to be symmetry
c        about seta=pi/2. 'fi' is assumed to be isotropic in 2pi
        coset=dlog(1.d0+1.7183*pyr(1))
        if(pyr(1).lt.0.5d0)coset=-coset
c021005
        siset=1.-coset*coset
        if(siset.lt.1.d-28)siset=1.d-28
c021005
        siset=dsqrt(siset)
        fi=2.*3.1416*pyr(1)
        cofi=dcos(fi)
        sifi=dsin(fi)
600     continue
        pi(1)=ppp*siset*cofi
        pi(2)=ppp*siset*sifi
        pi(3)=ppp*coset
c021005
        pi4=ppp*ppp+am1*am1
        if(pi4.lt.1.d-28)pi4=1.d-28
c021005
        pi(4)=dsqrt(pi4)
        pj(1)=-pi(1)
        pj(2)=-pi(2)
        pj(3)=-pi(3)
c021005
        pj4=ppp*ppp+am2*am2
        if(pj4.lt.1.d-28)pj4=1.d-28
c021005
        pj(4)=dsqrt(pj4)
c       write(9,*)'before rotation'   ! sa
c       write(9,*)ppp,(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c050603
c       rotate to the frame where diquark (gluon), ps, is described
c       calculate the direction cosines of ps
        fi1=pyangl(ps(1),ps(2))
c021005
        ps12=ps(1)*ps(1)+ps(2)*ps(2)
        if(ps12.lt.1.d-28)ps12=1.d-28
        ps12=dsqrt(ps12)
c021005
        cta1=pyangl(ps(3),ps12)
        cfi1=dcos(fi1)
        sfi1=dsin(fi1)
        ccta1=dcos(cta1)
        scta1=dsin(cta1)
        sctas=siset
        cctas=coset
        sfis=sifi
        cfis=cofi
        pi(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
        pi(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
        pi(3)=ccta1*cctas-scta1*sctas*cfis
        pi(1)=ppp*pi(1)
        pi(2)=ppp*pi(2)
        pi(3)=ppp*pi(3)
        do i=1,3
        pj(i)=0.-pi(i)
        enddo
c021005
        pi4=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)+am1*am1
        if(pi4.lt.1.e-28)pi4=1.e-28
        pj4=pj(1)*pj(1)+pj(2)*pj(2)+pj(3)*pj(3)+am2*am2
        if(pj4.lt.1.d-28)pj4=1.d-28
        pi(4)=dsqrt(pi4)
        pj(4)=dsqrt(pj4)
c021005
c       write(9,*)'after rotation'   ! sa
c       write(9,*)(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c       boost to moving frame of diquark
        ee=ps(4)
        if(ee.lt.1.d-14)ee=1.d-14   ! 021005
        do i1=1,3
        bb(i1)=ps(i1)/ee
        enddo
c       write(9,*)'b=',(bb(i),i=1,3)
        call lorntz(1,bb,pi,pj)
c       write(9,*)'after boost back, ps=',(ps(i),i=1,4)   ! sa
c       write(9,*)(pi(i),i=1,4)   ! sa
c       write(9,*)(pj(i),i=1,4)   ! sa
c       write(9,*)(pi(i)+pj(i),i=1,4)   ! sa
c050603
        pp(1,1)=pi(1)
        pp(1,2)=pi(2)
        pp(1,3)=pi(3)
        pp(1,4)=pi(4)
        pp(2,1)=pj(1)
        pp(2,2)=pj(2)
        pp(2,3)=pj(3)
        pp(2,4)=pj(4)
c050603
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine conser(np,pp,ps)
c       adjust four momentum conservation
c       np: the # of particles
c       pp: four momenta of particles have to be conserved
c       ps: the four momentum should be conserved to
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        dimension pp(20,5),ps(4),ff(20),pxyz(3),arp(3)
        ps4=ps(4)
        do i=1,3
        pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=1,np
        es=es+pp(i,4)
        enddo
        fr=es/ps4
        if(dabs(1.-fr) .le. dep)goto 200
        do i=1,np
        amas=pp(i,5)
        ppm=pp(i,4)/amas
        ppf=ppm/fr
        ff(i)=dsqrt(dabs(ppf*ppf-1.d0)/(ppm*ppm-1.))
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=dabs(1.-pxyz(i)/ps(i))
        pxyz(i)=pxyz(i)-ps(i)
        enddo
        if(dabs(1.-fr).le.dep .and. arp(1).le.dep .and. arp(2).le.dep
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)/np
        enddo
        do i=1,np
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
        pp5=pp(i,5)
        pp52=pp5*pp5
        pp(i,4)=dsqrt(pp52+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
        enddo
        jj=jj+1
        if(jj.eq.4000)then
        write(9,*)'infinitive loop may occur in subroutine conser(),'
        write(9,*)'which means four-momentum conservation'
        write(9,*)'needed is hard to be achieved,check value'
        write(9,*)'of PARAM(9)'
        return
        endif
        goto 100
200     return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine coord(ii)
c       give four position to broken quarks
c       first broken quark takes the four position of diquark
c       second broken quark is arranged around first ones within
c        0.5 fm randumly in each of three position coordinates and has same
c        fourth position coordinate as diquark
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSZJ=80000)
      COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        dimension rr(3)
        do i1=1,3
c261002        rr(i1)=pyr(1)*v(ii,i1)
        rr(i1)=pyr(1)*0.5   ! 261002
        v(n+1,i1)=v(ii,i1)+rr(i1)
        if(pyr(1).gt.0.5d0)v(n+1,i1)=v(ii,i1)-rr(i1)
        enddo
        v(n+1,4)=v(ii,4)
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine copl(tt)
c	calculate coordinate of center of mass of non-freeze-out system
c	position of a particle, checking is it freezes out or not, is 
c	 calculated with respect to this origin.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
        COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
     	 COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
	common/sa4/tau(kszj),tlco(kszj,4)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	do ii=1,3
	coor(ii)=0.
	enddo
	samass=0.
	do 110 ii=1,n
c	if(tau(ii).gt.tt)goto 110
	if(ishp(ii).eq.0)goto 110
	kf=k(ii,2)
	amass=pmas(pycomp(kf),1)
	if(iabs(kf).eq.213 .or. kf.eq.113)amass1=p(ii,5)   ! 010600 
	if((iabs(kf).eq.213 .or. kf.eq.113) .and. dabs(amass-amass1)
     &	 .gt.0.001d0)amass=amass1   ! 010600 
	samass=samass+amass
	do 100 jj=1,3
	coor(jj)=coor(jj)+amass*v(ii,jj)
100	continue
110	continue
	do ii=1,3
	coor(ii)=coor(ii)/dmax1(0.14d0,samass)
	enddo
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ctlcre(lc,tc,tw)
c	create initial collision list  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000,nsize=280000)   ! 280722
      	COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non10,
     c   disbe(100,100)
        common/sa35/ncpart,ncpar(kszj)   ! 280722
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121 
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c	nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)
c081010	time=0.
c280722
        ncpart=0
        do i1=1,kszj
        ncpar(i1)=0
        enddo
c280722
	nctl=1
	dminf=100.
	nzpab=iabs(nzp)   ! in order to consider ppbar or pbarp
	nztab=iabs(nzt)
	nzpt=nzpab+nztab
	napt=nap+nat
	do 10 l=1,nzpab   ! projectile proton or lepton   ! 060813
	do l1=nzpab+1,nzpt   ! target proton
	tc(nctl)=0.
	mtc=0
	call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,iif,jf)
	if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
	enddo
	do l1=nap+nztab+1,napt   ! target neutron
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,iif,jf)
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
10      continue
        do 20 l=nzpt+1,nap+nztab   ! projectile neutron
	do l1=nzpab+1,nzpt   ! target proton	
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,iif,jf)
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
	do l1=nap+nztab+1,napt   ! target neutron
        tc(nctl)=0.
	mtc=0
        call coij(l,l1,nctl,lc,tc,tw,mtc,dminf,iif,jf)
        if(mtc.gt.0)then
	nctl=nctl+1
	mtc=0
	endif
        enddo
20	continue
        if(mtc.eq.0)nctl=nctl-1
c	do iop=1,nctl
c	write(*,*)lc(iop,1),lc(iop,2),tc(iop)
c	enddo	
	if(nctl.eq.0)then
c 	at least one collision should occur. this collision has the smallest 
c	 'least approaching distance', that is guaranteed by the variable 
c	 'dminf'
	lc(1,1)=iif
	lc(1,2)=jf
	tc(1)=0.02
	nctl=1
	endif
	do i=nctl+1,nsize
	do m=1,5
	lc(i,m)=0
	enddo
	tc(i)=0.
	tw(i)=0.
	enddo
c280722
        do i1=1,kszj
        if(ncpar(i1).eq.1)ncpart=ncpart+1
        enddo
c       print*,'in ctlcre ncpart=',ncpart
c280722
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coij(i,j,icp,lc,tc,tw,mtc,dminf,iif,jf)
c       calculate collision time & fill up lc(i,1-2) as well as tc(i)
c	 for creating the initial collsion list 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	PARAMETER (kszj=80000)
	parameter(nsize=280000)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sa35/ncpart,ncpar(kszj)   ! 280722
        dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension dr(3),db(3),dv(3),px(4),py(4),vi(4),vj(4)
     c	,pi(4),pj(4),b(3)
	ki=ksa(i,2)
	kj=ksa(j,2)
	pi(4)=psa(i,4)
	if(pi(4).lt.1.e-20)pi(4)=1.e-20
	pj(4)=psa(j,4)
	if(pj(4).lt.1.e-20)pj(4)=1.e-20
	deno6=pi(4)+pj(4)
	do k=1,3
	pi(k)=psa(i,k)
	pj(k)=psa(j,k)
	b(k)=(pi(k)+pj(k))/deno6
	enddo
	ilo=0
	call lorntz(ilo,b,pi,pj)
	do l=1,3
	px(l)=vsa(i,l)
	py(l)=vsa(j,l)
	enddo
	px(4)=0.
	py(4)=0.
	call lorntz(ilo,b,px,py)
	rb=0.
	bb=0.
	rr=0.
	rtai=0.
	do k=1,3
	vi(k)=pi(k)/pi(4)
	vj(k)=pj(k)/pj(4)
	enddo
	do k=1,3
	dr(k)=px(k)-py(k)-(vi(k)*px(4)-vj(k)*py(4))
	db(k)=vi(k)-vj(k)
	dv(k)=db(k)
	rb=rb+dr(k)*db(k)
	bb=db(k)**2+bb
	rr=rr+dr(k)*dr(k)
	enddo
	if(bb.le.1.e-10) return
	tcol=0.-rb/bb
c        if(tcol-px(4) .le. ddt)return
c        if(tcol-py(4) .le. ddt)return
c       for collision to occur,time must one step ahead
csa	if(tcol.lt.1.0e-7)return
	do iik=1,3
	dr(iik)=px(iik)-py(iik)-(vi(iik)*
     &	px(4)-vj(iik)*py(4))+tcol*db(iik)
	rtai=rtai+dr(iik)*dr(iik)
	enddo
	sg=rtai
c	sg=rr+tcol*rb
c	if(sg.lt.0)then
c	write(*,*)'sg=',sg   !
c	return
c	endif
	dmin=dsqrt(sg)
	if(dmin.lt.dminf)then
	dminf=dmin
	iif=i
	jf=j
	endif
	if(ipden.lt.2 .and. dmin.gt.ecsnn)return   ! 060813 120214
        if(ipden.gt.2 .and. dmin.gt.ecsen)return   ! 060813 120214
c	distance between the two particles should be smaller than ecsnn (ecsen)
c        060813
	do ik=1,3
	px(ik)=px(ik)+vi(ik)*(tcol-px(4))
	py(ik)=py(ik)+vj(ik)*(tcol-py(4))
	enddo
c	move along Newton trajectory in CMS
	px(4)=tcol
	py(4)=tcol
c	write(*,*)'CMStcol=',tcol
	ilo=1

	call lorntz(ilo,b,px,py)

c	transform back to Lab.
	if(px(4).gt.py(4)) px(4)=py(4)
	tcol=px(4)
c041204
	drmax=rao*dmax1(rnt,rnp)
	if(tcol.le.drmax)goto 180   ! 041204
        return   ! 041204
c041204
180	tc(icp)=tcol
	mtc=1
        lc(icp,1)=i
        lc(icp,2)=j
c	write(*,*)'LABtcol=',tcol
c280722
        if(ncpar(i).eq.0)ncpar(i)=1
        if(ncpar(j).eq.0)ncpar(j)=1
c280722        
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find(icp,tcp,lc,tc,tw,ico)
c	find out the binary collision with minimum collision time
	parameter(kszj=80000,nsize=280000)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP   
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121 
	dimension lc(nsize,5),tc(nsize),tw(nsize)
	icp=0
	tcp=20000.
	do i=1,nctl
	if(ico.eq.0)goto 100
	if(tc(i).le.1.0e-7) goto 241
100	if(tcp.lt.tc(i))  goto 241
	icp=i
	tcp=tc(i)
241	continue
	enddo
	if(nctl.eq.0)icp=0
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updtlp(time,lc,tc,tw,iii)
c	update collision list after calling 'pythia' successfully
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000)
        parameter(nsize=280000)
        common/sa1/kjp21,non1,bp,iiii,neve,nout,nosc
	common/sa2/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp
	common/sa14/ipyth(2000),idec(2000),iwide(2000)
c010530        common/sa19/kji   ! 16/09/99
        common/sa30/vneump,vneumt,mstptj   ! 010223
        COMMON/SBH/N,NON,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)  
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c  ,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c	ipyth: store line number of produced hardon in hadron list 'sa2',101221 
c	write(9,*)'in updtlp nctl=',nctl
c	do i=1,nctl
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo
c       loop over old colliding pairs
        j=0
	if(nctl.eq.0)goto 200
        do 400 i=1,nctl
        if((tc(i)-time).le.ddt) goto 400
c       through away the pair whih tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
400     continue
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
c010223
        if(mstptj.eq.1)then   ! for (c) channel (pythia-like channel)
        nctl=j
        return
        endif           
c010223        
c	write(9,*)'updtlp rmove nctl,nbh,m4=',j,n,numb(4)   ! sa
c	do i=1,j
c	write(9,*)'i,lci,lcj,t=',i,lc(i,1),lc(i,2),tc(i)
c	enddo

200	nctl=j+1

	m2=numb(2)   ! 060813
        m4=numb(4)
	m7=numb(7)   ! 241110
c	m9=numb(9)
c	m17=numb(17)
c	m19=numb(19)
c	m25=numb(25)
c	m29=numb(29)
c	m32=numb(32)
c	m34=numb(34)
c        m34=numb(kfmax-11)
c       subtract 11, since we do not consider the rescattering of x0c, etc
c101221 note: # of produced hadrons equal to zero (n=0) after call 'pythia'
c        in case of w/o reconstruction leading proton
c101221 proceed for case of with reconstruction leading nucleon
c101221 constract hadron collision pair composed of one from produced hadrons
c        and another one in 'sa2'
c       loop over produced hadrons in 'sbh'   ! 101221
	do j11=1,n
	j1=ipyth(j11)
        kfjab=iabs(ksa(j1,2))   ! 060813 120214   
c	write(9,*)'n,j11,j1,kf,m4=',n,j11,j1,kfj,m4   
	if(kfjab.ne.2212.and.kfjab.ne.2112.and.kfjab.ne.11.and.kfjab
     c	 .ne.12.and.kfjab.ne.13.and.kfjab.ne.14.and.kfjab.ne.15
     c	 .and.kfjab.ne.16)goto 300   ! 241110 m7 to m2 060813 120214
c060813	consider only the reinteraction among nucleons & nucleon with lepton
c060813 loop over particle list ('sa2')
c	mm=m34   
c060813	mm=m7      ! 241110
        mm=m2   ! 130913 m7 to m2
c060813	120214 consider only the reinteraction of j11 with nucleons
	do i=1,mm   

        if(nctl.gt.nsize)then
	write(9,*)'iiii,nsize,n,nctl=',iiii,nsize,n,nctl   ! sa
	stop 22222
	endif
c010600
	do j22=1,n
	j2=ipyth(j22)
	if(i.eq.j2)goto 600   ! avoid particle collide with itself
	enddo
c010600

	i1=i
	kfi=ksa(i1,2)
c	write(9,*)'j1,i1,kfj,kfi=',j1,i1,kfj,kfi   ! sa
        iflag=0
        call rsfilt(j1,i1,iflag)
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
c011204	call tcolij(i1,j1,time,nctl,lc,tc,tw)
        call tcolij(j1,i1,time,nctl,lc,tc,tw)
        if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100     continue
600     enddo   ! loop for i
300	enddo   ! loop for j11
700     if(tc(nctl).le.1.e-7) nctl=nctl-1
        do i=nctl+1,nsize
        do m=1,5
        lc(i,m)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo
        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine rsfilt(l,l1,iflag)
c	 subroutine rsfilt plays the role of first range filter 
c	 subroutine intdis plays the role of second range filter
c       collision pairs not interested can not filter through both of rsfilt 
c        and intdis
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)
	kl1ab=iabs(kl1)
	if(l.eq.l1) goto 10
	if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10

c060813 consider nn collision 
	if(klab.eq.2212.and.(kl1ab.eq.2112.or.kl1ab.eq.2212))goto 11
	if(klab.eq.2112.and.(kl1ab.eq.2112.or.kl1ab.eq.2212))goto 11
c060813 120214 consider interaction between lepton and nucleon
        if((klab.ge.11.and.klab.le.16).and.(kl1.eq.2112.or.kl1.eq.2212))
     c	 goto 11
        if((klab.eq.2112.or.klab.eq.2212).and.(kl1ab.ge.11.and.kl1ab.le.
     c	 16))goto 11
c060813
	goto 10

11	iflag=1
10	continue
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine tcolij(l,l1,time,icp,lc,tc,tw)
c	calculate collision time & fill up lc(i,1-2) as well as tc(i) 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
	parameter(nsize=280000)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
	common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram,
     &	iabsb,iabsm,non10,ajpsi,csspn,csspm,csen   ! 060813
        common/sa4/tau(kszj),tlco(kszj,4)
      	COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/wz/c17(500,3),ishp(kszj),tp(500),coor(3),p17(500,4)
	COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)

	dimension lc(nsize,5),tc(nsize),tw(nsize)
	dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
	dimension ri(4),rj(4),rfi(4),rfj(4),b(3)
       	 pel=p(l,4)
	pel1=p(l1,4)
	if(pel.lt.1.e-20)pel=1.e-20   ! 041204
	if(pel1.lt.1.e-20)pel1=1.e-20   ! 041204
	pi(4)=pel
	pj(4)=pel1
	do i=1,3
c	write(*,*)'p,psa=',p(l,i)
c	write(*,*)'p1,psa1=',p(l1,i)
	pi(i)=p(l,i)
	pj(i)=p(l1,i)
	b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
	enddo
	ilo=0
	call lorntz(ilo,b,pi,pj)
c	perform Lorentz transf. to CMS frame for momentum.
	bta=dsqrt(b(1)**2+b(2)**2+b(3)**2)
c	if boost is too violent,put particles on mass shell by hand.
	if(bta.gt.0.99999d+0)then
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)
	kl1ab=iabs(kl1)
	bmi=pmas(pycomp(kl),1)
	bmj=pmas(pycomp(kl1),1)
	pi(4)=dsqrt(bmi**2+pi(1)**2+pi(2)**2+pi(3)**2)
	pj(4)=dsqrt(bmj**2+pj(1)**2+pj(2)**2+pj(3)**2)
	endif
	ss=pi(4)+pj(4)
c	do not pair into the collision list if the threshold is too small.
c	if(((klab.eq.2211.or.klab.eq.2112).and.
c     &	(kl1ab.eq.2112.or.kl1ab.eq.2212)).and.ss.le.parp(2))goto 10

	do i=1,4
	ri(i)=v(l,i)
	rj(i)=v(l1,i)
c	write(*,*)'v,vsa=',v(l,i)
c	write(*,*)'v1,vsa1=',v(l1,i)
	enddo
cc	ri(4)=time
cc	rj(4)=time
	call lorntz(ilo,b,ri,rj)
c	perform Lorentz transf. to CMS frame for coordinate.
	rb=0.
	bb=0.
	rr=0.
	rtai=0.
	kflag=0
	do ik=1,3
	vi(ik)=pi(ik)/pi(4)
	vj(ik)=pj(ik)/pj(4)
	enddo

	do i=1,3
	rfi(i)=v(l,i)+(tau(l)-time)*(p(l,i)/p(l,4))
	rfj(i)=v(l1,i)+(tau(l1)-time)*(p(l1,i)/p(l1,4))
	enddo
	rfi(4)=tau(l)
	rfj(4)=tau(l1)
	call lorntz(ilo,b,rfi,rfj)
c	gamli=p(l,4)/p(l,5)
c	gamlj=p(l1,4)/p(l1,5)
 	ctaui=rfi(4)
  	ctauj=rfj(4)
 	tcol=ctaui
	if(ctaui.lt.ctauj)tcol=ctauj
 	do ik=1,3
 	db(ik)=(vi(ik)-vj(ik))*tcol
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+db(ik)
	rtai=rtai+dr(ik)*dr(ik)	
	enddo
  	dot=0.
	do ik=1,3
	dot=dr(ik)*pi(ik)+dot
	enddo
c	dot=-1
	if(dot.ge.0.)then
	kflag=1
	if(tcol.le.ri(4) )goto 10
	if(tcol.le.rj(4) )goto 10	
	else
	rtai=0.
	do ik=1,3
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))
	db(ik)=vi(ik)-vj(ik)
	rb=rb+dr(ik)*db(ik)
	bb=bb+db(ik)*db(ik)
	rr=rr+dr(ik)*dr(ik)
	enddo
	if(bb .le. 1.e-10)goto 10
	tcol=0.-rb/bb
	if(tcol.le.ri(4) )goto 10
	if(tcol.le.rj(4) )goto 10
	if(tcol-ctaui .le. 0.)goto 10
	if(tcol-ctauj .le. 0.)goto 10
c	for collision to occur,time must one step ahead
cTai
	do ik=1,3
	dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+tcol*db(ik)
	rtai=rtai+dr(ik)*dr(ik)
	enddo
c	gamai=pi(4)/pmas(pycomp(k(l,2)),1)
c	gamaj=pj(4)/pmas(pycomp(k(l1,2)),1)

C TAIAN

c	when collision happens,particles should already be produced
c	we give a zero formation time for particles produced after
c       calling 'pythia'
	sg1=rr+tcol*rb
c		write(*,*)'sar',sg1
ctai
	endif
	sg=rtai
c		write(*,*)'tair',sg
c	if(sg1.lt.0.and.iii.le.50)then
c	write(*,*)'sar,tair=',sg1,sg
c	dmin=0.
c	tcol=-rr/rb
c	goto 20
c	endif

	dmin=dsqrt(sg)
c	calculate the interaction distance between particles l & l1.
20	call intdis(l,l1,ss,rsig)
c	distance between the two particles should be smaller than rsig
	if(dmin.gt.rsig)goto 10
c	move along Newton trajectory in CMS
	do ik=1,3
	ri(ik)=ri(ik)+vi(ik)*(tcol-ri(4))
	rj(ik)=rj(ik)+vj(ik)*(tcol-rj(4))
	enddo
	ri(4)=tcol
	rj(4)=tcol
	ilo=1

	call lorntz(ilo,b,ri,rj)

c	transform back to Lab.
	tcol1=ri(4)
	tcol2=rj(4)
	if(kflag.eq.0)then
	if(tcol1-tau(l).lt.0.) goto 10
	if(tcol2-tau(l1).lt.0.) goto 10
	else
	if(tcol1-tau(l).lt.-1.E-4) goto 10
	if(tcol2-tau(l1).lt.-1.E-4) goto 10
	endif
	if(ri(4).gt.rj(4)) ri(4)=rj(4)
	tcol=ri(4)
	if(tcol.le.time)goto 10
c	collision happens in the future
c	if(ifram.eq.0)coor(3)=coor(3)+rnt
	do i=1,3
	ri(i)=v(l,i)+p(l,i)*(tcol-time)/pel-coor(i)
	rj(i)=v(l1,i)+p(l1,i)*(tcol-time)/pel1-coor(i)
	enddo
	rri=dsqrt(ri(1)*ri(1)+ri(2)*ri(2)+ri(3)*ri(3))
	rrj=dsqrt(rj(1)*rj(1)+rj(2)*rj(2)+rj(3)*rj(3))
c	the rnt in rao*max(rnt,rnp)+rnt is due to the fact that
c        we could not know the postion of the mass-center in the future
	rrr=rao*dmax1(rnt,rnp)
	if(ifram.eq.0)rrr=rao*dmax1(rnt,rnp)+rnt
c	if(dabs(k(l,2)).gt.1000.or.dabs(k(l1,2)).gt.1000)rrr=1.E+10*rrr
	if(rri.gt.rrr)goto 10
	if(rrj.gt.rrr)goto 10
c	particles under consideration must be still within considered region
c	 when the collision happens
	if(tcol.le.rrr)goto 18   ! 041204
        return   ! 041204
18	tc(icp)=tcol
	lc(icp,1)=l
	lc(icp,2)=l1
10	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine intdis(l,l1,ss,rsig)
c	calculate interaction distance between particles l and l1.
c	It plays also the role of second range filter
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
        COMMON/SA2/N,NON2,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
	common/sa13/kjp20,non13,vjp20,vjp21,vjp22,vjp23
	rsig=0.
	kl=k(l,2)
	kl1=k(l1,2)
	klab=iabs(kl)   ! 060813 120214
	kl1ab=iabs(kl1)   ! 060813 120214

	if(iabs(kl).eq.2212 .or. iabs(kl).eq.2112)idpl=1
        if(iabs(kl1).eq.2212 .or. iabs(kl1).eq.2112)idpl1=1
	if(klab.ge.11.and.klab.le.16)idpl=8   ! lepton 060813 120214
	if(kl1ab.ge.11.and.kl1ab.le.16)idpl1=8   ! lepton 060813 120214
        if(idpl.eq.1 .and. idpl1.eq.1)rsig=ecsnn
	if(idpl.eq.8 .and. idpl1.eq.1)rsig=ecsen   ! 060813
        if(idpl.eq.1 .and. idpl1.eq.8)rsig=ecsen   ! 060813
	return
	end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updatl(ic,jc,time,lc,tc,tw,iii)
c	update collision time list after elastic scattering
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
	parameter(nsize=280000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp  
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121 
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)

c	loop over old colliding pairs
	j=0
	do i=1,nctl
	i1=lc(i,1)
	j1=lc(i,2)
c	ia=(i1-ic)*(j1-jc)*(i1-jc)*(j1-ic)
c	if(ia.eq.0) goto 400
	if(i1.eq.ic .or. i1.eq.jc)goto 400
	if(j1.eq.ic .or. j1.eq.jc)goto 400
	if((tc(i)-time).le.ddt) goto 400
c	through away the pair whih tc<= time
	j=j+1
	tc(j)=tc(i)
	tw(j)=tw(i)
	do m=1,5
	lc(j,m)=lc(i,m)
	enddo
400	continue
	enddo
	do i=j+1,nctl+1
	tc(i)=0.0
	tw(i)=0.0
	do m=1,5
	lc(i,m)=0
	enddo
	enddo

	nctl=j+1
c	loop over particle list

        m2=numb(2)   ! 130913
        m4=numb(4)
        m7=numb(7)   ! 241110
	j1=ic
	do ik=1,2

c010530
	kfab=iabs(ksa(j1,2))   ! 060813 120214
	if(kfab.ne.2212.and.kfab.ne.2112.and.kfab.ne.11.and.kfab
     c   .ne.12.and.kfab.ne.13.and.kfab.ne.14.and.kfab.ne.15
     c   .and.kfab.ne.16)goto 300   ! 241110 m7 to m2 060813 120214
c130913	consider only the reinteraction among nucleons & nucleon with lepton 
c	 060813 120214
	mm=m2   ! 241110 060813 m7 to m2 
	do i=1,mm
	if(j1.eq.ic .and. i.eq.jc)goto 600 
	if(j1.eq.jc .and. i.eq.ic)goto 600
c	forbiden scattered particles colliding with each other
	if(nctl.gt.nsize)then
        write(MSTU(11),*)'size of array "nsize" needs to be extended'
        write(MSTU(11),*)'error is serious,stop running'
        stop 22222
        endif

 	i1=i
	iflag=0
	call rsfilt(j1,i1,iflag)
	if(iflag.eq.0)goto 100
	tc(nctl)=0.0
	call tcolij(i1,j1,time,nctl,lc,tc,tw)
	if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100	continue		
600	enddo
300	if(ik.eq.2)goto 500
	j1=jc
500	enddo
700	if(tc(nctl).le.1.e-7) nctl=nctl-1
	do i=nctl+1,nsize
	do m=1,5
	lc(i,m)=0
	enddo
	tc(i)=0.
	tw(i)=0.
	enddo
	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pauli(ii,ppaul)
c	calculate the unoccupation probability (ppaul) of particle ii
c	 in 'pyjets'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter(kszj=80000)
	COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/saf/nsa,nonsa,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	dimension rkk(kszj,4),pkk(kszj,4),rr(4),pp(4),b(3)
	kf=k(ii,2)   ! new
	xxii=v(ii,1) 
	yyii=v(ii,2) 
	zzii=v(ii,3) 
	ttii=v(ii,4)
	pxii=p(ii,1) 
	pyii=p(ii,2) 
	pzii=p(ii,3)
	eeii=p(ii,4) 
	b(1)=pxii/eeii
	b(2)=pyii/eeii
	b(3)=pzii/eeii
c	pick up the partons with same flavor as ii from 'pyjets' and 
c	 'saf' 
	nkk=0 
c	the new produced partons, except ii, in 'pyjets' are also to be
c	 considered
	do i=1,n   ! loop over new
	if(i.eq.ii)goto 100
	kfi1=k(i,2)
        if(kfi1.eq.kf)then
	nkk=nkk+1
	do j=1,4
	rkk(nkk,j)=v(i,j)
	pkk(nkk,j)=p(i,j)
	enddo
	endif
100	enddo
	do i1=1,nsa   ! loop over old
        kfi1=ksa(i1,2)
        if(kfi1.eq.kf)then
	nkk=nkk+1
        do j=1,4
        rkk(nkk,j)=vsa(i1,j)
        pkk(nkk,j)=psa(i1,j)
        enddo
	endif
        enddo
c	boost to the rest frame of ii
        ilo=0
        do 200 j2=1,nkk
	do j1=1,4
        rr(j1)=rkk(j2,j1)
        pp(j1)=pkk(j2,j1)
        enddo
        call lorntz(ilo,b,rr,pp)
        do j1=1,4
        rkk(j2,j1)=rr(j1)
        pkk(j2,j1)=pp(j1)
        enddo
200	enddo
	rr(1)=xxii
	rr(2)=yyii
	rr(3)=zzii
	rr(4)=ttii
	call lorntz(ilo,b,rr,rr)
	xxii=rr(1)	
	yyii=rr(2)	
	zzii=rr(3)	
	ttii=rr(4)	
c	calculate the number of partons occupied in or on the surface of 
c	 six dimension cub (around ii): (dr*dp)**3=h**3, dr*dp=h, h=1.24 
c	 GeV*fm/c, if dr=1.0 fm then dp=1.24 GeV/c
	anq=0   ! statistics of ocuupation number in (dr*dp)**3=h**3
	do i1=1,nkk   
	dxx=xxii-rkk(i1,1)
	dyy=yyii-rkk(i1,2)
	dzz=zzii-rkk(i1,3)
c	following three staments for without boost
cc	dpx=pxii-pkk(i1,1)
cc	dpy=pyii-pkk(i1,2)
cc	dpz=pzii-pkk(i1,3)
c       following three staments for with boost
	dpx=pkk(i1,1)
	dpy=pkk(i1,2)
	dpz=pkk(i1,3)
	dxx=dabs(dxx)
	dyy=dabs(dyy)
	dzz=dabs(dzz)
	dpx=dabs(dpx)
	dpy=dabs(dpy)
	dpz=dabs(dpz)
	if(dxx.le.0.5.and.dyy.le.0.5.and.dzz.le.0.5.and.
     c	 dpx.le.0.62.and.dpy.le.0.62.and.dpz.le.0.62)anq=anq+1.
	enddo
	proba=anq/6.
c	6=2*3, spin and colour degeneracies of quark (antiquark)
	ppaul=1.-proba
	return
	end


c060813 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine crosep(rots,csen)
!!exponential interpolation for ep total cross section (in mbarn) 
!! calculated with herafitter by Xing-Long Li on 10/Dec./2013
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (rotsMin=2.d0) !!the lower limit of rots (GeV)
      parameter (rotsMax=1003) !!the upper limit of rots (GeV)
      parameter (Ndat=315)
!!sect: ep cross section data for sqrtS range: 2~1003GeV
!!the unit is mbarn
!!when i=1,sqrtS=2.0 GeV; sqrtS=2.0*(1.02)^(i-1) GeV otherwise
      dimension sect(Ndat)
      parameter (sect=[
     & 0.125377D-04,0.136063D-04,0.147199D-04,0.158792D-04,0.170824D-04,
     & 0.183314D-04,0.196228D-04,0.209579D-04,0.223386D-04,0.237574D-04,
     & 0.252260D-04,0.267283D-04,0.282729D-04,0.298683D-04,0.314933D-04,
     & 0.331676D-04,0.348724D-04,0.366252D-04,0.384099D-04,0.402340D-04,
     & 0.421004D-04,0.439997D-04,0.459366D-04,0.479080D-04,0.499176D-04,
     & 0.519633D-04,0.540432D-04,0.561548D-04,0.583026D-04,0.604863D-04,
     & 0.626984D-04,0.649432D-04,0.672178D-04,0.695343D-04,0.718663D-04,
     & 0.742366D-04,0.766461D-04,0.790697D-04,0.815331D-04,0.840216D-04,
     & 0.865379D-04,0.890822D-04,0.916557D-04,0.942586D-04,0.968830D-04,
     & 0.995384D-04,0.102212D-03,0.104924D-03,0.107655D-03,0.110411D-03,
     & 0.113187D-03,0.115993D-03,0.118824D-03,0.121673D-03,0.124552D-03,
     & 0.127461D-03,0.130369D-03,0.133312D-03,0.136294D-03,0.139258D-03,
     & 0.142289D-03,0.145318D-03,0.148358D-03,0.151438D-03,0.154537D-03,
     & 0.157648D-03,0.160778D-03,0.163928D-03,0.167102D-03,0.170296D-03,
     & 0.173515D-03,0.176740D-03,0.179988D-03,0.183258D-03,0.186551D-03,
     & 0.189855D-03,0.193175D-03,0.196515D-03,0.199880D-03,0.203265D-03,
     & 0.206648D-03,0.210060D-03,0.213497D-03,0.216952D-03,0.220408D-03,
     & 0.223888D-03,0.227392D-03,0.230909D-03,0.234442D-03,0.237992D-03,
     & 0.241560D-03,0.245143D-03,0.248743D-03,0.252364D-03,0.256010D-03,
     & 0.259659D-03,0.263324D-03,0.267009D-03,0.270715D-03,0.274444D-03,
     & 0.278182D-03,0.281933D-03,0.285693D-03,0.289480D-03,0.293325D-03,
     & 0.297112D-03,0.300951D-03,0.304824D-03,0.308704D-03,0.312598D-03,
     & 0.316505D-03,0.320433D-03,0.324381D-03,0.328354D-03,0.332347D-03,
     & 0.336349D-03,0.340362D-03,0.344403D-03,0.348459D-03,0.352536D-03,
     & 0.356641D-03,0.360763D-03,0.364893D-03,0.369052D-03,0.373232D-03,
     & 0.377425D-03,0.381629D-03,0.385859D-03,0.390121D-03,0.394418D-03,
     & 0.398725D-03,0.403041D-03,0.407380D-03,0.411739D-03,0.416122D-03,
     & 0.420532D-03,0.424966D-03,0.429429D-03,0.433908D-03,0.438403D-03,
     & 0.442923D-03,0.447466D-03,0.452044D-03,0.456646D-03,0.461274D-03,
     & 0.465905D-03,0.470575D-03,0.475274D-03,0.480000D-03,0.484746D-03,
     & 0.489511D-03,0.494317D-03,0.499153D-03,0.504004D-03,0.508874D-03,
     & 0.513782D-03,0.518727D-03,0.523706D-03,0.528696D-03,0.533711D-03,
     & 0.538757D-03,0.543840D-03,0.548953D-03,0.554107D-03,0.559291D-03,
     & 0.564487D-03,0.569720D-03,0.574975D-03,0.580275D-03,0.585612D-03,
     & 0.590979D-03,0.596380D-03,0.601802D-03,0.607256D-03,0.612757D-03,
     & 0.618283D-03,0.623860D-03,0.629465D-03,0.635094D-03,0.640766D-03,
     & 0.646474D-03,0.652220D-03,0.658001D-03,0.663819D-03,0.669669D-03,
     & 0.675552D-03,0.681488D-03,0.687464D-03,0.693475D-03,0.699515D-03,
     & 0.705602D-03,0.711726D-03,0.717890D-03,0.724097D-03,0.730344D-03,
     & 0.736634D-03,0.742995D-03,0.749397D-03,0.755772D-03,0.762216D-03,
     & 0.768727D-03,0.775280D-03,0.781876D-03,0.788516D-03,0.795202D-03,
     & 0.801928D-03,0.808704D-03,0.815513D-03,0.822366D-03,0.829276D-03,
     & 0.836237D-03,0.843243D-03,0.850306D-03,0.857418D-03,0.864568D-03,
     & 0.871777D-03,0.879027D-03,0.886321D-03,0.893664D-03,0.901059D-03,
     & 0.908518D-03,0.916017D-03,0.923583D-03,0.931201D-03,0.938865D-03,
     & 0.946581D-03,0.954347D-03,0.962169D-03,0.970049D-03,0.977973D-03,
     & 0.985961D-03,0.993995D-03,0.100210D-02,0.101027D-02,0.101850D-02,
     & 0.102677D-02,0.103510D-02,0.104349D-02,0.105192D-02,0.106040D-02,
     & 0.106897D-02,0.107761D-02,0.108635D-02,0.109510D-02,0.110390D-02,
     & 0.111274D-02,0.112163D-02,0.113063D-02,0.113971D-02,0.114886D-02,
     & 0.115806D-02,0.116732D-02,0.117663D-02,0.118600D-02,0.119543D-02,
     & 0.120495D-02,0.121453D-02,0.122419D-02,0.123391D-02,0.124371D-02,
     & 0.125356D-02,0.126348D-02,0.127347D-02,0.128353D-02,0.129367D-02,
     & 0.130387D-02,0.131415D-02,0.132451D-02,0.133493D-02,0.134544D-02,
     & 0.135602D-02,0.136666D-02,0.137737D-02,0.138816D-02,0.139902D-02,
     & 0.140998D-02,0.142099D-02,0.143209D-02,0.144327D-02,0.145453D-02,
     & 0.146587D-02,0.147729D-02,0.148878D-02,0.150035D-02,0.151200D-02,
     & 0.152373D-02,0.153554D-02,0.154745D-02,0.155945D-02,0.157155D-02,
     & 0.158371D-02,0.159591D-02,0.160822D-02,0.162062D-02,0.163313D-02,
     & 0.164572D-02,0.165839D-02,0.167115D-02,0.168400D-02,0.169693D-02,
     & 0.170996D-02,0.172306D-02,0.173628D-02,0.174958D-02,0.176297D-02,
     & 0.177646D-02,0.179003D-02,0.180371D-02,0.181748D-02,0.183135D-02,
     & 0.184530D-02,0.185935D-02,0.187350D-02,0.188775D-02,0.190209D-02
     & ])
!!check if rots is in range of data set,if not in range 2~1003,return -1
      if(rots.lt.rotsMin.or.rots.gt.rotsMax)then
            csen=-1.d0
            return
      endif
!!calculate csen by Interpolation
      x=dlog(rots/2.d0)/dlog(1.02d0)+1.d0
      i=floor(x)
      csen=(sect(i+1)-sect(i))*(x-i)+sect(i)

      return
      end



c060813 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ltof(ii)
c       move ii-th particle (lepton) in pyjets to first position 060813 120214
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=80000)
      COMMON/PYJETS/N,NONJ,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c  disbe(100,100)
        dimension kk(5),pp(5),vv(5)
        do jj=1,5
        kk(jj)=k(ii,jj)
        pp(jj)=p(ii,jj)
        vv(jj)=v(ii,jj)
        enddo
c       move particle list (pyjets) one step forward from ii-1 to 1
        do j1=ii-1,1,-1
        do jj=1,5
        k(j1+1,jj)=k(j1,jj)
        p(j1+1,jj)=p(j1,jj)
        v(j1+1,jj)=v(j1,jj)
        enddo
        enddo
        do jj=1,5
        k(1,jj)=kk(jj)
        p(1,jj)=pp(jj)
        v(1,jj)=vv(jj)
        enddo
        do i1=1,kfmax
        numbs(i1)=numbs(i1)+1   ! now first particle in pyjets is e-
        enddo
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc161222
        subroutine coinel_nn(l,l1,kl,kl1,pi,pj)   ! 110123
c       simulate four momentum of scattered particles in inelastic 
c        nucleon-nucleon collision according formula of (4. 33)
c        in book of Ben-Hao Sa, etc., 'Simulation Physics for High Energy 
c        Nucleus Collisions'
c       l (l1): order number of scattering and/or scattered prticle in 'sa2'
c       kl (kl1): flavor code of scattered particle  
c       pi (pj): four momentum of scattering particle before scatering
c       pi (pj): four momentum of scattered particle after scatering also

        IMPLICIT DOUBLE PRECISION(A-H, O-Z)
        IMPLICIT INTEGER(I-N)
        INTEGER PYK,PYCHGE,PYCOMP
        parameter(kszj=80000)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        dimension pi(4),pj(4)
        am3=pymass(kl)   ! mass of one scattered particle
        am4=pymass(kl1)   ! mass of another one scattered particle
        ss=pi(4)+pj(4)
        pp=(ss*ss-(am3+am4)**2)*(ss*ss-(am3-am4)**2)/(4.*ss*ss)
        if(pp.lt.0.)pp=1.e-10
        pp=sqrt(pp)
        pi(4)=(ss*ss+am3**2-am4**2)/(2.*ss)
c       energy of one particle (between two) after scattering
        fis=2.*3.1415926*pyr(1)
        cfis=cos(fis)
        sfis=sin(fis)
        csita=2*pyr(1)-1.
        ssita=dsqrt(1.-csita**2)
        pi(1)=pp*ssita*cfis
        pi(2)=pp*ssita*sfis
        pi(3)=pp*csita
        do i=1,3
        pj(i)=0.-pi(i)
        enddo
        pj(4)=ss-pi(4)

        do i2=1,4
        psa(l,i2)=pi(i2)   ! four momentum of one collided particle
        psa(l1,i2)=pj(i2)   ! four momentum of another one collided particle
        enddo

        return
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine padecy(l,l1,time)   ! 161222
c	deacy of a hardron
c       l: order # of colliding hadron in 'sa2'
c       l1: order # of another colliding hadron in 'sa2'
c       time: current time
c       output hadrons compose 'pyjets'
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (kszj=80000)
        common/pyjets/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
c       write(22,*)'enter nsa,l,l1,k=',nsa,l,l1,ksa(l,2),ksa(l1,2)
c       write(22,*)'pl=',psa(l,1),psa(l,2),psa(l,3)
c       write(22,*)'pl1=',psa(l1,1),psa(l1,2),psa(l1,3)
c       write(22,*)'e=',psa(l,4),psa(l1,4)
c       write(22,*)'enter loop nsa =',nsa
c       call prt_sa2(nsa,c)
       
        n=0
        do i1=1,kszj
        do i2=1,5
        k(i1,i2)=0
        p(i1,i2)=0.
        v(i1,i2)=0
        enddo
        enddo

        n=1
        k(1,1)=1
        kf=ksa(l,2)
        k(1,2)=kf
        do i2=1,4
        p(1,i2)=psa(l,i2)
        enddo
        p(1,5)=pymass(kf)
        mdcy(pycomp(kf),1)=1
c       call pyedit(2)
c       write(22,*)'af. fill first'
c       call pylist(1)

        call pydecy(1)
        call pyedit(2)
c       write(22,*)'af. call pydecy'
c       call pylist(1)

        n=3
        k(3,1)=1
        k(3,2)=ksa(l1,2)
        do i1=1,5
        p(3,i1)=psa(l1,i1)
        enddo
c       write(22,*)'af. padecy n=',n
c       call prt_pyj(n,c)

        return
	end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sa2pyj(l,l1)   ! 161222
c	part of 'sa2' to 'pyjets'
c       l: order # of colliding hadron in 'sa2'
c       l1: order # of another colliding hadron in 'sa2'
c       output hadrons compose 'pyjets'
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      parameter (kszj=80000)
        common/pyjets/n,npad,k(kszj,5),p(kszj,5),v(kszj,5)
        COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        n=0
        do i1=1,kszj
        do i2=1,5
        k(i1,i2)=0
        p(i1,i2)=0.
        v(i1,i2)=0
        enddo
        enddo

        n=2
        k(1,1)=1
        k(1,2)=ksa(l,2)
        do i1=1,5
        p(1,i1)=psa(l,i1)
        enddo
        k(2,1)=1
        k(2,2)=ksa(l1,2)
        do i1=1,5
        p(2,i1)=psa(l1,i1)
        enddo

        return
        end       


 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updpip_nn(l,l1,icp,lc,tc,tw,time)   ! 250123
c       update particle list 'sa2' ('sbh' to 'sa2') after inela. scattering 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
        parameter (kszj=80000,nsize=280000)   ! 080223
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
        common/sbh/nbh,nonbh,kbh(kszj,5),pbh(kszj,5),vbh(kszj,5)
        dimension lc(nsize,5),tc(nsize),tw(nsize)         


c       update particle list 'sa2' ('sbh' to 'sa2')      
        if(nbh.eq.3)then   ! from call 'padecay'
        do i2=1,5
        ksa(l,i2)=kbh(1,i2)
        psa(l,i2)=pbh(1,i2)
        vsa(l,i2)=vbh(1,i2)
        enddo
        do i2=1,5
        ksa(l1,i2)=kbh(3,i2)
        psa(l1,i2)=pbh(3,i2)
        vsa(l1,i2)=vbh(3,i2)
        enddo
        nsa=nsa+1
        do i2=1,5
        ksa(nsa,i2)=kbh(2,i2)
        psa(nsa,i2)=pbh(2,i2)
        vsa(nsa,i2)=vbh(2,i2)
        enddo
        endif

        if(nbh.eq.2)then   ! from 'call sa2pyj' 
        do i2=1,5
        ksa(l,i2)=kbh(1,i2)
        psa(l,i2)=pbh(1,i2)
        vsa(l,i2)=vbh(1,i2)
        enddo
        do i2=1,5
        ksa(l1,i2)=kbh(2,i2)
        psa(l1,i2)=pbh(2,i2)
        vsa(l1,i2)=vbh(2,i2)
        enddo
        endif

        return 
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updatl_nn(ic,jc,time,lc,tc,tw,ik,nsa0)   ! 250123
c	update collision time list after inela. scattering
c       ic (jc): order number of inela. scattering patticle in 'sa2'
c       ik: number of hadrons joining reconstrution of hadronic colli. pair 
c       time: current time
c       nsa0: last oder number in 'sa2' before call 'updpip_nn'
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
	parameter (kszj=80000)
	parameter(nsize=280000)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa2/nsa,non2,ksa(kszj,5),psa(kszj,5),vsa(kszj,5)
	common/sa5/kfmax,kfaco(100),numb(100),numbs(100),non5,
     c   disbe(100,100)
        common/sa12/ppsa(5),nchan,nsjp,sjp,taup,taujp  
	common/ctllist/nctl,noinel(600),nctl0,nctlm   ! 180121 230121 
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn,ecspsn,ecspsm
     c 	,rnt,rnp,rao,rou0,vneu,vneum,ecsspn,ecsspm,ecsen   ! 060813
        common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t,
     c  nap,nat,nzp,nzt,pio
	dimension lc(nsize,5),tc(nsize),tw(nsize)

c	loop over old colliding pairs
	j=0
	do i=1,nctl
	i1=lc(i,1)
	j1=lc(i,2)
	if(i1.eq.ic .or. i1.eq.jc)goto 400
	if(j1.eq.ic .or. j1.eq.jc)goto 400
	if((tc(i)-time).le.ddt) goto 400
c	through away the pair whih tc<= time
	j=j+1
	tc(j)=tc(i)
	tw(j)=tw(i)
	do m=1,5
	lc(j,m)=lc(i,m)
	enddo
400	continue
	enddo
	do i=j+1,nctl+1
	tc(i)=0.0
	tw(i)=0.0
	do m=1,5
	lc(i,m)=0
	enddo
	enddo
        if(ik.eq.0)then
        nctl=j        
        return
        endif

	nctl=j+1
c       loop over inelastically scatered (new generated) hadrons

        do 500 ii=1,2   
c       loop over new generated hadron which has filled in 'sa2' at ic (jc)
        if(ii.eq.1)j1=ic
        if(ii.eq.2)j1=jc
        kfab=iabs(ksa(j1,2))   ! 060813 120214
c       if(kfab.ne.2212.and.kfab.ne.2112.and.kfab.ne.11.and.kfab
c       c   .ne.12.and.kfab.ne.13.and.kfab.ne.14.and.kfab.ne.15
c       c   .and.kfab.ne.16)goto 300   ! 241110 m7 to m2 060813 120214
        if(kfab.ne.2212.and.kfab.ne.2112.and.kfab.ne.211
     c   .and.kfab.ne.11.and.kfab.ne.12.and.kfab.ne.13
     c   .and.kfab.ne.14.and.kfab.ne.15.and.kfab.ne.16)goto 300 ! not join
!        reconstruction of hadronic collision pair        
c       consider only the reinteraction between NN, Npi+, Npi-, & Nlepton 
c	loop over old particle list
        do i=1,nsa   ! do 600
	if(j1.eq.ic .and. i.eq.jc)goto 600 
	if(j1.eq.jc .and. i.eq.ic)goto 600
c	forbiden scattered particles colliding with each other
	if(nctl.gt.nsize)then
        write(MSTU(11),*)'size of array "nsize" needs to be extended'
        write(MSTU(11),*)'error is serious,stop running'
        stop 22222
        endif

 	i1=i
	iflag=0
	call rsfilt(j1,i1,iflag)
	if(iflag.eq.0)goto 100
	tc(nctl)=0.0
	call tcolij(i1,j1,time,nctl,lc,tc,tw)
	if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100	continue		
600	enddo
300	continue
500	enddo

        do 800 ii=nsa0+1,nsa   ! loop over other generated hadrons
        j1=ii
        kfab=iabs(ksa(j1,2))
        if(kfab.ne.2212.and.kfab.ne.2112.and.kfab.ne.211
     c   .and.kfab.ne.11.and.kfab.ne.12.and.kfab.ne.13
     c   .and.kfab.ne.14.and.kfab.ne.15.and.kfab.ne.16)goto 301 ! not join
!        reconstruction of hadronic collision pair        
c       consider only the reinteraction between NN, Npi+, Npi-, & Nlepton 
c	loop over old particle list
        do i=1,nsa0   ! do 600
	if(j1.eq.ic .and. i.eq.jc)goto 601 
	if(j1.eq.jc .and. i.eq.ic)goto 601
c	forbiden scattered particles colliding with each other
	if(nctl.gt.nsize)then
        write(MSTU(11),*)'size of array "nsize" needs to be extended'
        write(MSTU(11),*)'error is serious,stop running'
        stop 22222
        endif

 	i1=i
	iflag=0
	call rsfilt(j1,i1,iflag)
	if(iflag.eq.0)goto 101
	tc(nctl)=0.0
	call tcolij(i1,j1,time,nctl,lc,tc,tw)
	if(tc(nctl).gt.1.0e-7) nctl=nctl+1
101	continue		
601	enddo
301	continue
800	enddo

700	if(tc(nctl).le.1.e-7) nctl=nctl-1
	do i=nctl+1,nsize
	do m=1,5
	lc(i,m)=0
	enddo
	tc(i)=0.
	tw(i)=0.
	enddo
	return
	end



ccccccccccccccccccccccccccccccccccccc  end  cccccccccccccccccccccccccccc
