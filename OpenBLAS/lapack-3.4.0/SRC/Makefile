include ../make.inc

#######################################################################
#  This is the makefile to create a library for LAPACK.
#  The files are organized as follows:
#       ALLAUX  -- Auxiliary routines called from all precisions
#
#       SCLAUX  -- Auxiliary routines called from both REAL and COMPLEX.
#       DZLAUX  -- Auxiliary routines called from both DOUBLE and COMPLEX*16.
#
#       DSLASRC -- Double-single mixed precision real routines called from
#                  single, single-extra and double precision real LAPACK
#                  routines (i.e. from SLASRC, SXLASRC, DLASRC).
#       ZCLASRC -- Double-single mixed precision complex routines called from
#                  single, single-extra and double precision complex LAPACK
#                  routines (i.e. from CLASRC, CXLASRC, ZLASRC).
#
#       SLASRC  -- Single precision real LAPACK routines
#       SXLASRC -- Single precision real LAPACK routines using extra
#                  precision.
#       CLASRC  -- Single precision complex LAPACK routines
#       CXLASRC -- Single precision complex LAPACK routines using extra
#                  precision.
#       DLASRC  -- Double precision real LAPACK routines
#       DXLASRC -- Double precision real LAPACK routines using extra
#                  precision.
#       ZLASRC  -- Double precision complex LAPACK routines
#       ZXLASRC -- Double precision complex LAPACK routines using extra
#                  precision.   
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       lapack.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#  ***Note***
#  The functions lsame, second, dsecnd, slamch, and dlamch may have
#  to be installed before compiling the library.  Refer to the
#  installation guide, LAPACK Working Note 41, for instructions.
#
#######################################################################

ALLAUX = ilaenv.$(SUFFIX) ieeeck.$(SUFFIX) lsamen.$(SUFFIX) xerbla_array.$(SUFFIX) iparmq.$(SUFFIX)	\
    ilaprec.$(SUFFIX) ilatrans.$(SUFFIX) ilauplo.$(SUFFIX) iladiag.$(SUFFIX) chla_transtype.$(SUFFIX) \
    ../INSTALL/ilaver.$(SUFFIX) 

SCLAUX = \
   sbdsdc.$(SUFFIX) \
   sbdsqr.$(SUFFIX) sdisna.$(SUFFIX) slabad.$(SUFFIX) slacpy.$(SUFFIX) sladiv.$(SUFFIX) slae2.$(SUFFIX)  slaebz.$(SUFFIX) \
   slaed0.$(SUFFIX) slaed1.$(SUFFIX) slaed2.$(SUFFIX) slaed3.$(SUFFIX) slaed4.$(SUFFIX) slaed5.$(SUFFIX) slaed6.$(SUFFIX) \
   slaed7.$(SUFFIX) slaed8.$(SUFFIX) slaed9.$(SUFFIX) slaeda.$(SUFFIX) slaev2.$(SUFFIX) slagtf.$(SUFFIX) \
   slagts.$(SUFFIX) slamrg.$(SUFFIX) slanst.$(SUFFIX) \
   slapy2.$(SUFFIX) slapy3.$(SUFFIX) slarnv.$(SUFFIX) \
   slarra.$(SUFFIX) slarrb.$(SUFFIX) slarrc.$(SUFFIX) slarrd.$(SUFFIX) slarre.$(SUFFIX) slarrf.$(SUFFIX) slarrj.$(SUFFIX) \
   slarrk.$(SUFFIX) slarrr.$(SUFFIX) slaneg.$(SUFFIX) \
   slartg.$(SUFFIX) slaruv.$(SUFFIX) slas2.$(SUFFIX)  slascl.$(SUFFIX) \
   slasd0.$(SUFFIX) slasd1.$(SUFFIX) slasd2.$(SUFFIX) slasd3.$(SUFFIX) slasd4.$(SUFFIX) slasd5.$(SUFFIX) slasd6.$(SUFFIX) \
   slasd7.$(SUFFIX) slasd8.$(SUFFIX) slasda.$(SUFFIX) slasdq.$(SUFFIX) slasdt.$(SUFFIX) \
   slaset.$(SUFFIX) slasq1.$(SUFFIX) slasq2.$(SUFFIX) slasq3.$(SUFFIX) slasq4.$(SUFFIX) slasq5.$(SUFFIX) slasq6.$(SUFFIX) \
   slasr.$(SUFFIX)  slasrt.$(SUFFIX) slassq.$(SUFFIX) slasv2.$(SUFFIX) spttrf.$(SUFFIX) sstebz.$(SUFFIX) sstedc.$(SUFFIX) \
   ssteqr.$(SUFFIX) ssterf.$(SUFFIX) slaisnan.$(SUFFIX) sisnan.$(SUFFIX) \
   slartgp.$(SUFFIX) slartgs.$(SUFFIX) \
    ../INSTALL/second_$(TIMER).$(SUFFIX)

DZLAUX = \
   dbdsdc.$(SUFFIX) \
   dbdsqr.$(SUFFIX) ddisna.$(SUFFIX) dlabad.$(SUFFIX) dlacpy.$(SUFFIX) dladiv.$(SUFFIX) dlae2.$(SUFFIX)  dlaebz.$(SUFFIX) \
   dlaed0.$(SUFFIX) dlaed1.$(SUFFIX) dlaed2.$(SUFFIX) dlaed3.$(SUFFIX) dlaed4.$(SUFFIX) dlaed5.$(SUFFIX) dlaed6.$(SUFFIX) \
   dlaed7.$(SUFFIX) dlaed8.$(SUFFIX) dlaed9.$(SUFFIX) dlaeda.$(SUFFIX) dlaev2.$(SUFFIX) dlagtf.$(SUFFIX) \
   dlagts.$(SUFFIX) dlamrg.$(SUFFIX) dlanst.$(SUFFIX) \
   dlapy2.$(SUFFIX) dlapy3.$(SUFFIX) dlarnv.$(SUFFIX) \
   dlarra.$(SUFFIX) dlarrb.$(SUFFIX) dlarrc.$(SUFFIX) dlarrd.$(SUFFIX) dlarre.$(SUFFIX) dlarrf.$(SUFFIX) dlarrj.$(SUFFIX) \
   dlarrk.$(SUFFIX) dlarrr.$(SUFFIX) dlaneg.$(SUFFIX) \
   dlartg.$(SUFFIX) dlaruv.$(SUFFIX) dlas2.$(SUFFIX)  dlascl.$(SUFFIX) \
   dlasd0.$(SUFFIX) dlasd1.$(SUFFIX) dlasd2.$(SUFFIX) dlasd3.$(SUFFIX) dlasd4.$(SUFFIX) dlasd5.$(SUFFIX) dlasd6.$(SUFFIX) \
   dlasd7.$(SUFFIX) dlasd8.$(SUFFIX) dlasda.$(SUFFIX) dlasdq.$(SUFFIX) dlasdt.$(SUFFIX) \
   dlaset.$(SUFFIX) dlasq1.$(SUFFIX) dlasq2.$(SUFFIX) dlasq3.$(SUFFIX) dlasq4.$(SUFFIX) dlasq5.$(SUFFIX) dlasq6.$(SUFFIX) \
   dlasr.$(SUFFIX)  dlasrt.$(SUFFIX) dlassq.$(SUFFIX) dlasv2.$(SUFFIX) dpttrf.$(SUFFIX) dstebz.$(SUFFIX) dstedc.$(SUFFIX) \
   dsteqr.$(SUFFIX) dsterf.$(SUFFIX) dlaisnan.$(SUFFIX) disnan.$(SUFFIX) \
   dlartgp.$(SUFFIX) dlartgs.$(SUFFIX) \
   ../INSTALL/dsecnd_$(TIMER).$(SUFFIX)

SLASRC = \
   sgbbrd.$(SUFFIX) sgbcon.$(SUFFIX) sgbequ.$(SUFFIX) sgbrfs.$(SUFFIX) sgbsv.$(SUFFIX)  \
   sgbsvx.$(SUFFIX) sgbtf2.$(SUFFIX) sgbtrf.$(SUFFIX) sgbtrs.$(SUFFIX) sgebak.$(SUFFIX) sgebal.$(SUFFIX) sgebd2.$(SUFFIX) \
   sgebrd.$(SUFFIX) sgecon.$(SUFFIX) sgeequ.$(SUFFIX) sgees.$(SUFFIX)  sgeesx.$(SUFFIX) sgeev.$(SUFFIX)  sgeevx.$(SUFFIX) \
   sgegs.$(SUFFIX)  sgegv.$(SUFFIX)  sgehd2.$(SUFFIX) sgehrd.$(SUFFIX) sgelq2.$(SUFFIX) sgelqf.$(SUFFIX) \
   sgels.$(SUFFIX)  sgelsd.$(SUFFIX) sgelss.$(SUFFIX) sgelsx.$(SUFFIX) sgelsy.$(SUFFIX) sgeql2.$(SUFFIX) sgeqlf.$(SUFFIX) \
   sgeqp3.$(SUFFIX) sgeqpf.$(SUFFIX) sgeqr2.$(SUFFIX) sgeqr2p.$(SUFFIX) sgeqrf.$(SUFFIX) sgeqrfp.$(SUFFIX) sgerfs.$(SUFFIX) \
   sgerq2.$(SUFFIX) sgerqf.$(SUFFIX) sgesc2.$(SUFFIX) sgesdd.$(SUFFIX) sgesv.$(SUFFIX)  sgesvd.$(SUFFIX) sgesvx.$(SUFFIX) \
   sgetc2.$(SUFFIX) sgetri.$(SUFFIX) \
   sggbak.$(SUFFIX) sggbal.$(SUFFIX) sgges.$(SUFFIX)  sggesx.$(SUFFIX) sggev.$(SUFFIX)  sggevx.$(SUFFIX) \
   sggglm.$(SUFFIX) sgghrd.$(SUFFIX) sgglse.$(SUFFIX) sggqrf.$(SUFFIX) \
   sggrqf.$(SUFFIX) sggsvd.$(SUFFIX) sggsvp.$(SUFFIX) sgtcon.$(SUFFIX) sgtrfs.$(SUFFIX) sgtsv.$(SUFFIX)  \
   sgtsvx.$(SUFFIX) sgttrf.$(SUFFIX) sgttrs.$(SUFFIX) sgtts2.$(SUFFIX) shgeqz.$(SUFFIX) \
   shsein.$(SUFFIX) shseqr.$(SUFFIX) slabrd.$(SUFFIX) slacon.$(SUFFIX) slacn2.$(SUFFIX) \
   slaein.$(SUFFIX) slaexc.$(SUFFIX) slag2.$(SUFFIX)  slags2.$(SUFFIX) slagtm.$(SUFFIX) slagv2.$(SUFFIX) slahqr.$(SUFFIX) \
   slahrd.$(SUFFIX) slahr2.$(SUFFIX) slaic1.$(SUFFIX) slaln2.$(SUFFIX) slals0.$(SUFFIX) slalsa.$(SUFFIX) slalsd.$(SUFFIX) \
   slangb.$(SUFFIX) slange.$(SUFFIX) slangt.$(SUFFIX) slanhs.$(SUFFIX) slansb.$(SUFFIX) slansp.$(SUFFIX) \
   slansy.$(SUFFIX) slantb.$(SUFFIX) slantp.$(SUFFIX) slantr.$(SUFFIX) slanv2.$(SUFFIX) \
   slapll.$(SUFFIX) slapmt.$(SUFFIX) \
   slaqgb.$(SUFFIX) slaqge.$(SUFFIX) slaqp2.$(SUFFIX) slaqps.$(SUFFIX) slaqsb.$(SUFFIX) slaqsp.$(SUFFIX) slaqsy.$(SUFFIX) \
   slaqr0.$(SUFFIX) slaqr1.$(SUFFIX) slaqr2.$(SUFFIX) slaqr3.$(SUFFIX) slaqr4.$(SUFFIX) slaqr5.$(SUFFIX) \
   slaqtr.$(SUFFIX) slar1v.$(SUFFIX) slar2v.$(SUFFIX) ilaslr.$(SUFFIX) ilaslc.$(SUFFIX) \
   slarf.$(SUFFIX)  slarfb.$(SUFFIX) slarfg.$(SUFFIX) slarfgp.$(SUFFIX) slarft.$(SUFFIX) slarfx.$(SUFFIX) slargv.$(SUFFIX) \
   slarrv.$(SUFFIX) slartv.$(SUFFIX)  \
   slarz.$(SUFFIX)  slarzb.$(SUFFIX) slarzt.$(SUFFIX) slasy2.$(SUFFIX) slasyf.$(SUFFIX) \
   slatbs.$(SUFFIX) slatdf.$(SUFFIX) slatps.$(SUFFIX) slatrd.$(SUFFIX) slatrs.$(SUFFIX) slatrz.$(SUFFIX) slatzm.$(SUFFIX) \
   sopgtr.$(SUFFIX) sopmtr.$(SUFFIX) sorg2l.$(SUFFIX) sorg2r.$(SUFFIX) \
   sorgbr.$(SUFFIX) sorghr.$(SUFFIX) sorgl2.$(SUFFIX) sorglq.$(SUFFIX) sorgql.$(SUFFIX) sorgqr.$(SUFFIX) sorgr2.$(SUFFIX) \
   sorgrq.$(SUFFIX) sorgtr.$(SUFFIX) sorm2l.$(SUFFIX) sorm2r.$(SUFFIX) \
   sormbr.$(SUFFIX) sormhr.$(SUFFIX) sorml2.$(SUFFIX) sormlq.$(SUFFIX) sormql.$(SUFFIX) sormqr.$(SUFFIX) sormr2.$(SUFFIX) \
   sormr3.$(SUFFIX) sormrq.$(SUFFIX) sormrz.$(SUFFIX) sormtr.$(SUFFIX) spbcon.$(SUFFIX) spbequ.$(SUFFIX) spbrfs.$(SUFFIX) \
   spbstf.$(SUFFIX) spbsv.$(SUFFIX)  spbsvx.$(SUFFIX) \
   spbtf2.$(SUFFIX) spbtrf.$(SUFFIX) spbtrs.$(SUFFIX) spocon.$(SUFFIX) spoequ.$(SUFFIX) sporfs.$(SUFFIX) sposv.$(SUFFIX)  \
   sposvx.$(SUFFIX) spotri.$(SUFFIX)  spstrf.$(SUFFIX) spstf2.$(SUFFIX) \
   sppcon.$(SUFFIX) sppequ.$(SUFFIX) \
   spprfs.$(SUFFIX) sppsv.$(SUFFIX)  sppsvx.$(SUFFIX) spptrf.$(SUFFIX) spptri.$(SUFFIX) spptrs.$(SUFFIX) sptcon.$(SUFFIX) \
   spteqr.$(SUFFIX) sptrfs.$(SUFFIX) sptsv.$(SUFFIX)  sptsvx.$(SUFFIX) spttrs.$(SUFFIX) sptts2.$(SUFFIX) srscl.$(SUFFIX)  \
   ssbev.$(SUFFIX)  ssbevd.$(SUFFIX) ssbevx.$(SUFFIX) ssbgst.$(SUFFIX) ssbgv.$(SUFFIX)  ssbgvd.$(SUFFIX) ssbgvx.$(SUFFIX) \
   ssbtrd.$(SUFFIX) sspcon.$(SUFFIX) sspev.$(SUFFIX)  sspevd.$(SUFFIX) sspevx.$(SUFFIX) sspgst.$(SUFFIX) \
   sspgv.$(SUFFIX)  sspgvd.$(SUFFIX) sspgvx.$(SUFFIX) ssprfs.$(SUFFIX) sspsv.$(SUFFIX)  sspsvx.$(SUFFIX) ssptrd.$(SUFFIX) \
   ssptrf.$(SUFFIX) ssptri.$(SUFFIX) ssptrs.$(SUFFIX) sstegr.$(SUFFIX) sstein.$(SUFFIX) sstev.$(SUFFIX)  sstevd.$(SUFFIX) sstevr.$(SUFFIX) \
   sstevx.$(SUFFIX) ssycon.$(SUFFIX) ssyev.$(SUFFIX)  ssyevd.$(SUFFIX) ssyevr.$(SUFFIX) ssyevx.$(SUFFIX) ssygs2.$(SUFFIX) \
   ssygst.$(SUFFIX) ssygv.$(SUFFIX)  ssygvd.$(SUFFIX) ssygvx.$(SUFFIX) ssyrfs.$(SUFFIX) ssysv.$(SUFFIX)  ssysvx.$(SUFFIX) \
   ssytd2.$(SUFFIX) ssytf2.$(SUFFIX) ssytrd.$(SUFFIX) ssytrf.$(SUFFIX) ssytri.$(SUFFIX) ssytri2.$(SUFFIX) ssytri2x.$(SUFFIX) \
   ssyswapr.$(SUFFIX) ssytrs.$(SUFFIX) ssytrs2.$(SUFFIX) ssyconv.$(SUFFIX) stbcon.$(SUFFIX) \
   stbrfs.$(SUFFIX) stbtrs.$(SUFFIX) stgevc.$(SUFFIX) stgex2.$(SUFFIX) stgexc.$(SUFFIX) stgsen.$(SUFFIX) \
   stgsja.$(SUFFIX) stgsna.$(SUFFIX) stgsy2.$(SUFFIX) stgsyl.$(SUFFIX) stpcon.$(SUFFIX) stprfs.$(SUFFIX) stptri.$(SUFFIX) \
   stptrs.$(SUFFIX) \
   strcon.$(SUFFIX) strevc.$(SUFFIX) strexc.$(SUFFIX) strrfs.$(SUFFIX) strsen.$(SUFFIX) strsna.$(SUFFIX) strsyl.$(SUFFIX) \
   strtrs.$(SUFFIX) stzrqf.$(SUFFIX) stzrzf.$(SUFFIX) sstemr.$(SUFFIX) \
   slansf.$(SUFFIX) spftrf.$(SUFFIX) spftri.$(SUFFIX) spftrs.$(SUFFIX) ssfrk.$(SUFFIX) stfsm.$(SUFFIX) stftri.$(SUFFIX) stfttp.$(SUFFIX) \
   stfttr.$(SUFFIX) stpttf.$(SUFFIX) stpttr.$(SUFFIX) strttf.$(SUFFIX) strttp.$(SUFFIX) \
   sgejsv.$(SUFFIX)  sgesvj.$(SUFFIX)  sgsvj0.$(SUFFIX)  sgsvj1.$(SUFFIX) \
   sgeequb.$(SUFFIX) ssyequb.$(SUFFIX) spoequb.$(SUFFIX) sgbequb.$(SUFFIX) \
   sbbcsd.$(SUFFIX) slapmr.$(SUFFIX) sorbdb.$(SUFFIX) sorcsd.$(SUFFIX) \
   sgeqrt.$(SUFFIX) sgeqrt2.$(SUFFIX) sgeqrt3.$(SUFFIX) sgemqrt.$(SUFFIX) \
   stpqrt.$(SUFFIX) stpqrt2.$(SUFFIX) stpmqrt.$(SUFFIX) stprfb.$(SUFFIX)

DSLASRC = spotrs.$(SUFFIX) 

ifdef USEXBLAS
SXLASRC = sgesvxx.$(SUFFIX) sgerfsx.$(SUFFIX) sla_gerfsx_extended.$(SUFFIX) sla_geamv.$(SUFFIX)		\
   sla_gercond.$(SUFFIX) sla_gerpvgrw.$(SUFFIX) ssysvxx.$(SUFFIX) ssyrfsx.$(SUFFIX)			\
   sla_syrfsx_extended.$(SUFFIX) sla_syamv.$(SUFFIX) sla_syrcond.$(SUFFIX) sla_syrpvgrw.$(SUFFIX)	\
   sposvxx.$(SUFFIX) sporfsx.$(SUFFIX) sla_porfsx_extended.$(SUFFIX) sla_porcond.$(SUFFIX)		\
   sla_porpvgrw.$(SUFFIX) sgbsvxx.$(SUFFIX) sgbrfsx.$(SUFFIX) sla_gbrfsx_extended.$(SUFFIX)		\
   sla_gbamv.$(SUFFIX) sla_gbrcond.$(SUFFIX) sla_gbrpvgrw.$(SUFFIX) sla_lin_berr.$(SUFFIX) slarscl2.$(SUFFIX)	\
   slascl2.$(SUFFIX) sla_wwaddw.$(SUFFIX)
endif

CLASRC = \
   cbdsqr.$(SUFFIX) cgbbrd.$(SUFFIX) cgbcon.$(SUFFIX) cgbequ.$(SUFFIX) cgbrfs.$(SUFFIX) cgbsv.$(SUFFIX)  cgbsvx.$(SUFFIX) \
   cgbtf2.$(SUFFIX) cgbtrf.$(SUFFIX) cgbtrs.$(SUFFIX) cgebak.$(SUFFIX) cgebal.$(SUFFIX) cgebd2.$(SUFFIX) cgebrd.$(SUFFIX) \
   cgecon.$(SUFFIX) cgeequ.$(SUFFIX) cgees.$(SUFFIX)  cgeesx.$(SUFFIX) cgeev.$(SUFFIX)  cgeevx.$(SUFFIX) \
   cgegs.$(SUFFIX)  cgegv.$(SUFFIX)  cgehd2.$(SUFFIX) cgehrd.$(SUFFIX) cgelq2.$(SUFFIX) cgelqf.$(SUFFIX) \
   cgels.$(SUFFIX)  cgelsd.$(SUFFIX) cgelss.$(SUFFIX) cgelsx.$(SUFFIX) cgelsy.$(SUFFIX) cgeql2.$(SUFFIX) cgeqlf.$(SUFFIX) cgeqp3.$(SUFFIX) \
   cgeqpf.$(SUFFIX) cgeqr2.$(SUFFIX) cgeqr2p.$(SUFFIX) cgeqrf.$(SUFFIX) cgeqrfp.$(SUFFIX) cgerfs.$(SUFFIX) \
   cgerq2.$(SUFFIX) cgerqf.$(SUFFIX) cgesc2.$(SUFFIX) cgesdd.$(SUFFIX) cgesv.$(SUFFIX)  cgesvd.$(SUFFIX)  \
   cgesvx.$(SUFFIX) cgetc2.$(SUFFIX) cgetri.$(SUFFIX) \
   cggbak.$(SUFFIX) cggbal.$(SUFFIX) cgges.$(SUFFIX)  cggesx.$(SUFFIX) cggev.$(SUFFIX)  cggevx.$(SUFFIX) cggglm.$(SUFFIX) \
   cgghrd.$(SUFFIX) cgglse.$(SUFFIX) cggqrf.$(SUFFIX) cggrqf.$(SUFFIX) \
   cggsvd.$(SUFFIX) cggsvp.$(SUFFIX) \
   cgtcon.$(SUFFIX) cgtrfs.$(SUFFIX) cgtsv.$(SUFFIX)  cgtsvx.$(SUFFIX) cgttrf.$(SUFFIX) cgttrs.$(SUFFIX) cgtts2.$(SUFFIX) chbev.$(SUFFIX)  \
   chbevd.$(SUFFIX) chbevx.$(SUFFIX) chbgst.$(SUFFIX) chbgv.$(SUFFIX)  chbgvd.$(SUFFIX) chbgvx.$(SUFFIX) chbtrd.$(SUFFIX) \
   checon.$(SUFFIX) cheev.$(SUFFIX)  cheevd.$(SUFFIX) cheevr.$(SUFFIX) cheevx.$(SUFFIX) chegs2.$(SUFFIX) chegst.$(SUFFIX) \
   chegv.$(SUFFIX)  chegvd.$(SUFFIX) chegvx.$(SUFFIX) cherfs.$(SUFFIX) chesv.$(SUFFIX)  chesvx.$(SUFFIX) chetd2.$(SUFFIX) \
   chetf2.$(SUFFIX) chetrd.$(SUFFIX) \
   chetrf.$(SUFFIX) chetri.$(SUFFIX) chetri2.$(SUFFIX) chetri2x.$(SUFFIX) cheswapr.$(SUFFIX) \
   chetrs.$(SUFFIX) chetrs2.$(SUFFIX) chgeqz.$(SUFFIX) chpcon.$(SUFFIX) chpev.$(SUFFIX)  chpevd.$(SUFFIX) \
   chpevx.$(SUFFIX) chpgst.$(SUFFIX) chpgv.$(SUFFIX)  chpgvd.$(SUFFIX) chpgvx.$(SUFFIX) chprfs.$(SUFFIX) chpsv.$(SUFFIX)  \
   chpsvx.$(SUFFIX) \
   chptrd.$(SUFFIX) chptrf.$(SUFFIX) chptri.$(SUFFIX) chptrs.$(SUFFIX) chsein.$(SUFFIX) chseqr.$(SUFFIX) clabrd.$(SUFFIX) \
   clacgv.$(SUFFIX) clacon.$(SUFFIX) clacn2.$(SUFFIX) clacp2.$(SUFFIX) clacpy.$(SUFFIX) clacrm.$(SUFFIX) clacrt.$(SUFFIX) cladiv.$(SUFFIX) \
   claed0.$(SUFFIX) claed7.$(SUFFIX) claed8.$(SUFFIX) \
   claein.$(SUFFIX) claesy.$(SUFFIX) claev2.$(SUFFIX) clags2.$(SUFFIX) clagtm.$(SUFFIX) \
   clahef.$(SUFFIX) clahqr.$(SUFFIX) \
   clahrd.$(SUFFIX) clahr2.$(SUFFIX) claic1.$(SUFFIX) clals0.$(SUFFIX) clalsa.$(SUFFIX) clalsd.$(SUFFIX) clangb.$(SUFFIX) clange.$(SUFFIX) clangt.$(SUFFIX) \
   clanhb.$(SUFFIX) clanhe.$(SUFFIX) \
   clanhp.$(SUFFIX) clanhs.$(SUFFIX) clanht.$(SUFFIX) clansb.$(SUFFIX) clansp.$(SUFFIX) clansy.$(SUFFIX) clantb.$(SUFFIX) \
   clantp.$(SUFFIX) clantr.$(SUFFIX) clapll.$(SUFFIX) clapmt.$(SUFFIX) clarcm.$(SUFFIX) claqgb.$(SUFFIX) claqge.$(SUFFIX) \
   claqhb.$(SUFFIX) claqhe.$(SUFFIX) claqhp.$(SUFFIX) claqp2.$(SUFFIX) claqps.$(SUFFIX) claqsb.$(SUFFIX) \
   claqr0.$(SUFFIX) claqr1.$(SUFFIX) claqr2.$(SUFFIX) claqr3.$(SUFFIX) claqr4.$(SUFFIX) claqr5.$(SUFFIX) \
   claqsp.$(SUFFIX) claqsy.$(SUFFIX) clar1v.$(SUFFIX) clar2v.$(SUFFIX) ilaclr.$(SUFFIX) ilaclc.$(SUFFIX) \
   clarf.$(SUFFIX)  clarfb.$(SUFFIX) clarfg.$(SUFFIX) clarft.$(SUFFIX) clarfgp.$(SUFFIX) \
   clarfx.$(SUFFIX) clargv.$(SUFFIX) clarnv.$(SUFFIX) clarrv.$(SUFFIX) clartg.$(SUFFIX) clartv.$(SUFFIX) \
   clarz.$(SUFFIX)  clarzb.$(SUFFIX) clarzt.$(SUFFIX) clascl.$(SUFFIX) claset.$(SUFFIX) clasr.$(SUFFIX)  classq.$(SUFFIX) \
   clasyf.$(SUFFIX) clatbs.$(SUFFIX) clatdf.$(SUFFIX) clatps.$(SUFFIX) clatrd.$(SUFFIX) clatrs.$(SUFFIX) clatrz.$(SUFFIX) \
   clatzm.$(SUFFIX) cpbcon.$(SUFFIX) cpbequ.$(SUFFIX) cpbrfs.$(SUFFIX) cpbstf.$(SUFFIX) cpbsv.$(SUFFIX)  \
   cpbsvx.$(SUFFIX) cpbtf2.$(SUFFIX) cpbtrf.$(SUFFIX) cpbtrs.$(SUFFIX) cpocon.$(SUFFIX) cpoequ.$(SUFFIX) cporfs.$(SUFFIX) \
   cposv.$(SUFFIX)  cposvx.$(SUFFIX) cpotri.$(SUFFIX) cpstrf.$(SUFFIX) cpstf2.$(SUFFIX) \
   cppcon.$(SUFFIX) cppequ.$(SUFFIX) cpprfs.$(SUFFIX) cppsv.$(SUFFIX)  cppsvx.$(SUFFIX) cpptrf.$(SUFFIX) cpptri.$(SUFFIX) cpptrs.$(SUFFIX) \
   cptcon.$(SUFFIX) cpteqr.$(SUFFIX) cptrfs.$(SUFFIX) cptsv.$(SUFFIX)  cptsvx.$(SUFFIX) cpttrf.$(SUFFIX) cpttrs.$(SUFFIX) cptts2.$(SUFFIX) \
   crot.$(SUFFIX)   cspcon.$(SUFFIX) cspmv.$(SUFFIX)  cspr.$(SUFFIX)   csprfs.$(SUFFIX) cspsv.$(SUFFIX)  \
   cspsvx.$(SUFFIX) csptrf.$(SUFFIX) csptri.$(SUFFIX) csptrs.$(SUFFIX) csrscl.$(SUFFIX) cstedc.$(SUFFIX) \
   cstegr.$(SUFFIX) cstein.$(SUFFIX) csteqr.$(SUFFIX) csycon.$(SUFFIX) csymv.$(SUFFIX) \
   csyr.$(SUFFIX)   csyrfs.$(SUFFIX) csysv.$(SUFFIX)  csysvx.$(SUFFIX) csytf2.$(SUFFIX) csytrf.$(SUFFIX) csytri.$(SUFFIX) csytri2.$(SUFFIX) csytri2x.$(SUFFIX) \
   csyswapr.$(SUFFIX) csytrs.$(SUFFIX) csytrs2.$(SUFFIX) csyconv.$(SUFFIX) ctbcon.$(SUFFIX) ctbrfs.$(SUFFIX) ctbtrs.$(SUFFIX) ctgevc.$(SUFFIX) ctgex2.$(SUFFIX) \
   ctgexc.$(SUFFIX) ctgsen.$(SUFFIX) ctgsja.$(SUFFIX) ctgsna.$(SUFFIX) ctgsy2.$(SUFFIX) ctgsyl.$(SUFFIX) ctpcon.$(SUFFIX) \
   ctprfs.$(SUFFIX) ctptri.$(SUFFIX) \
   ctptrs.$(SUFFIX) ctrcon.$(SUFFIX) ctrevc.$(SUFFIX) ctrexc.$(SUFFIX) ctrrfs.$(SUFFIX) ctrsen.$(SUFFIX) ctrsna.$(SUFFIX) \
   ctrsyl.$(SUFFIX) ctrtrs.$(SUFFIX) ctzrqf.$(SUFFIX) ctzrzf.$(SUFFIX) cung2l.$(SUFFIX) cung2r.$(SUFFIX) \
   cungbr.$(SUFFIX) cunghr.$(SUFFIX) cungl2.$(SUFFIX) cunglq.$(SUFFIX) cungql.$(SUFFIX) cungqr.$(SUFFIX) cungr2.$(SUFFIX) \
   cungrq.$(SUFFIX) cungtr.$(SUFFIX) cunm2l.$(SUFFIX) cunm2r.$(SUFFIX) cunmbr.$(SUFFIX) cunmhr.$(SUFFIX) cunml2.$(SUFFIX) \
   cunmlq.$(SUFFIX) cunmql.$(SUFFIX) cunmqr.$(SUFFIX) cunmr2.$(SUFFIX) cunmr3.$(SUFFIX) cunmrq.$(SUFFIX) cunmrz.$(SUFFIX) \
   cunmtr.$(SUFFIX) cupgtr.$(SUFFIX) cupmtr.$(SUFFIX) icmax1.$(SUFFIX) scsum1.$(SUFFIX) cstemr.$(SUFFIX) \
   chfrk.$(SUFFIX) ctfttp.$(SUFFIX) clanhf.$(SUFFIX) cpftrf.$(SUFFIX) cpftri.$(SUFFIX) cpftrs.$(SUFFIX) ctfsm.$(SUFFIX) ctftri.$(SUFFIX) \
   ctfttr.$(SUFFIX) ctpttf.$(SUFFIX) ctpttr.$(SUFFIX) ctrttf.$(SUFFIX) ctrttp.$(SUFFIX) \
   cgeequb.$(SUFFIX) cgbequb.$(SUFFIX) csyequb.$(SUFFIX) cpoequb.$(SUFFIX) cheequb.$(SUFFIX) \
   cbbcsd.$(SUFFIX) clapmr.$(SUFFIX) cunbdb.$(SUFFIX) cuncsd.$(SUFFIX) \
   cgeqrt.$(SUFFIX) cgeqrt2.$(SUFFIX) cgeqrt3.$(SUFFIX) cgemqrt.$(SUFFIX) \
   ctpqrt.$(SUFFIX) ctpqrt2.$(SUFFIX) ctpmqrt.$(SUFFIX) ctprfb.$(SUFFIX)

ifdef USEXBLAS
CXLASRC = cgesvxx.$(SUFFIX) cgerfsx.$(SUFFIX) cla_gerfsx_extended.$(SUFFIX) cla_geamv.$(SUFFIX) \
   cla_gercond_c.$(SUFFIX) cla_gercond_x.$(SUFFIX) cla_gerpvgrw.$(SUFFIX) \
   csysvxx.$(SUFFIX) csyrfsx.$(SUFFIX) cla_syrfsx_extended.$(SUFFIX) cla_syamv.$(SUFFIX) \
   cla_syrcond_c.$(SUFFIX) cla_syrcond_x.$(SUFFIX) cla_syrpvgrw.$(SUFFIX) \
   cposvxx.$(SUFFIX) cporfsx.$(SUFFIX) cla_porfsx_extended.$(SUFFIX) \
   cla_porcond_c.$(SUFFIX) cla_porcond_x.$(SUFFIX) cla_porpvgrw.$(SUFFIX) \
   cgbsvxx.$(SUFFIX) cgbrfsx.$(SUFFIX) cla_gbrfsx_extended.$(SUFFIX) cla_gbamv.$(SUFFIX) \
   cla_gbrcond_c.$(SUFFIX) cla_gbrcond_x.$(SUFFIX) cla_gbrpvgrw.$(SUFFIX) \
   chesvxx.$(SUFFIX) cherfsx.$(SUFFIX) cla_herfsx_extended.$(SUFFIX) cla_heamv.$(SUFFIX) \
   cla_hercond_c.$(SUFFIX) cla_hercond_x.$(SUFFIX) cla_herpvgrw.$(SUFFIX) \
   cla_lin_berr.$(SUFFIX) clarscl2.$(SUFFIX) clascl2.$(SUFFIX) cla_wwaddw.$(SUFFIX)
endif

ZCLASRC = cpotrs.$(SUFFIX) 

DLASRC = \
   dgbbrd.$(SUFFIX) dgbcon.$(SUFFIX) dgbequ.$(SUFFIX) dgbrfs.$(SUFFIX) dgbsv.$(SUFFIX)  \
   dgbsvx.$(SUFFIX) dgbtf2.$(SUFFIX) dgbtrf.$(SUFFIX) dgbtrs.$(SUFFIX) dgebak.$(SUFFIX) dgebal.$(SUFFIX) dgebd2.$(SUFFIX) \
   dgebrd.$(SUFFIX) dgecon.$(SUFFIX) dgeequ.$(SUFFIX) dgees.$(SUFFIX)  dgeesx.$(SUFFIX) dgeev.$(SUFFIX)  dgeevx.$(SUFFIX) \
   dgegs.$(SUFFIX)  dgegv.$(SUFFIX)  dgehd2.$(SUFFIX) dgehrd.$(SUFFIX) dgelq2.$(SUFFIX) dgelqf.$(SUFFIX) \
   dgels.$(SUFFIX)  dgelsd.$(SUFFIX) dgelss.$(SUFFIX) dgelsx.$(SUFFIX) dgelsy.$(SUFFIX) dgeql2.$(SUFFIX) dgeqlf.$(SUFFIX) \
   dgeqp3.$(SUFFIX) dgeqpf.$(SUFFIX) dgeqr2.$(SUFFIX) dgeqr2p.$(SUFFIX) dgeqrf.$(SUFFIX) dgeqrfp.$(SUFFIX) dgerfs.$(SUFFIX) \
   dgerq2.$(SUFFIX) dgerqf.$(SUFFIX) dgesc2.$(SUFFIX) dgesdd.$(SUFFIX) dgesv.$(SUFFIX)  dgesvd.$(SUFFIX) dgesvx.$(SUFFIX)  \
   dgetc2.$(SUFFIX) dgetri.$(SUFFIX) \
   dggbak.$(SUFFIX) dggbal.$(SUFFIX) dgges.$(SUFFIX)  dggesx.$(SUFFIX) dggev.$(SUFFIX)  dggevx.$(SUFFIX) \
   dggglm.$(SUFFIX) dgghrd.$(SUFFIX) dgglse.$(SUFFIX) dggqrf.$(SUFFIX) \
   dggrqf.$(SUFFIX) dggsvd.$(SUFFIX) dggsvp.$(SUFFIX) dgtcon.$(SUFFIX) dgtrfs.$(SUFFIX) dgtsv.$(SUFFIX)  \
   dgtsvx.$(SUFFIX) dgttrf.$(SUFFIX) dgttrs.$(SUFFIX) dgtts2.$(SUFFIX) dhgeqz.$(SUFFIX) \
   dhsein.$(SUFFIX) dhseqr.$(SUFFIX) dlabrd.$(SUFFIX) dlacon.$(SUFFIX) dlacn2.$(SUFFIX) \
   dlaein.$(SUFFIX) dlaexc.$(SUFFIX) dlag2.$(SUFFIX)  dlags2.$(SUFFIX) dlagtm.$(SUFFIX) dlagv2.$(SUFFIX) dlahqr.$(SUFFIX) \
   dlahrd.$(SUFFIX) dlahr2.$(SUFFIX) dlaic1.$(SUFFIX) dlaln2.$(SUFFIX) dlals0.$(SUFFIX) dlalsa.$(SUFFIX) dlalsd.$(SUFFIX) \
   dlangb.$(SUFFIX) dlange.$(SUFFIX) dlangt.$(SUFFIX) dlanhs.$(SUFFIX) dlansb.$(SUFFIX) dlansp.$(SUFFIX) \
   dlansy.$(SUFFIX) dlantb.$(SUFFIX) dlantp.$(SUFFIX) dlantr.$(SUFFIX) dlanv2.$(SUFFIX) \
   dlapll.$(SUFFIX) dlapmt.$(SUFFIX) \
   dlaqgb.$(SUFFIX) dlaqge.$(SUFFIX) dlaqp2.$(SUFFIX) dlaqps.$(SUFFIX) dlaqsb.$(SUFFIX) dlaqsp.$(SUFFIX) dlaqsy.$(SUFFIX) \
   dlaqr0.$(SUFFIX) dlaqr1.$(SUFFIX) dlaqr2.$(SUFFIX) dlaqr3.$(SUFFIX) dlaqr4.$(SUFFIX) dlaqr5.$(SUFFIX) \
   dlaqtr.$(SUFFIX) dlar1v.$(SUFFIX) dlar2v.$(SUFFIX) iladlr.$(SUFFIX) iladlc.$(SUFFIX) \
   dlarf.$(SUFFIX)  dlarfb.$(SUFFIX) dlarfg.$(SUFFIX) dlarfgp.$(SUFFIX) dlarft.$(SUFFIX) dlarfx.$(SUFFIX) \
   dlargv.$(SUFFIX) dlarrv.$(SUFFIX) dlartv.$(SUFFIX)  \
   dlarz.$(SUFFIX)  dlarzb.$(SUFFIX) dlarzt.$(SUFFIX) dlasy2.$(SUFFIX) dlasyf.$(SUFFIX) \
   dlatbs.$(SUFFIX) dlatdf.$(SUFFIX) dlatps.$(SUFFIX) dlatrd.$(SUFFIX) dlatrs.$(SUFFIX) dlatrz.$(SUFFIX) dlatzm.$(SUFFIX) \
   dopgtr.$(SUFFIX) dopmtr.$(SUFFIX) dorg2l.$(SUFFIX) dorg2r.$(SUFFIX) \
   dorgbr.$(SUFFIX) dorghr.$(SUFFIX) dorgl2.$(SUFFIX) dorglq.$(SUFFIX) dorgql.$(SUFFIX) dorgqr.$(SUFFIX) dorgr2.$(SUFFIX) \
   dorgrq.$(SUFFIX) dorgtr.$(SUFFIX) dorm2l.$(SUFFIX) dorm2r.$(SUFFIX) \
   dormbr.$(SUFFIX) dormhr.$(SUFFIX) dorml2.$(SUFFIX) dormlq.$(SUFFIX) dormql.$(SUFFIX) dormqr.$(SUFFIX) dormr2.$(SUFFIX) \
   dormr3.$(SUFFIX) dormrq.$(SUFFIX) dormrz.$(SUFFIX) dormtr.$(SUFFIX) dpbcon.$(SUFFIX) dpbequ.$(SUFFIX) dpbrfs.$(SUFFIX) \
   dpbstf.$(SUFFIX) dpbsv.$(SUFFIX)  dpbsvx.$(SUFFIX) \
   dpbtf2.$(SUFFIX) dpbtrf.$(SUFFIX) dpbtrs.$(SUFFIX) dpocon.$(SUFFIX) dpoequ.$(SUFFIX) dporfs.$(SUFFIX) dposv.$(SUFFIX)  \
   dposvx.$(SUFFIX) dpotri.$(SUFFIX) dpotrs.$(SUFFIX) dpstrf.$(SUFFIX) dpstf2.$(SUFFIX) \
   dppcon.$(SUFFIX) dppequ.$(SUFFIX) \
   dpprfs.$(SUFFIX) dppsv.$(SUFFIX)  dppsvx.$(SUFFIX) dpptrf.$(SUFFIX) dpptri.$(SUFFIX) dpptrs.$(SUFFIX) dptcon.$(SUFFIX) \
   dpteqr.$(SUFFIX) dptrfs.$(SUFFIX) dptsv.$(SUFFIX)  dptsvx.$(SUFFIX) dpttrs.$(SUFFIX) dptts2.$(SUFFIX) drscl.$(SUFFIX)  \
   dsbev.$(SUFFIX)  dsbevd.$(SUFFIX) dsbevx.$(SUFFIX) dsbgst.$(SUFFIX) dsbgv.$(SUFFIX)  dsbgvd.$(SUFFIX) dsbgvx.$(SUFFIX) \
   dsbtrd.$(SUFFIX)  dspcon.$(SUFFIX) dspev.$(SUFFIX)  dspevd.$(SUFFIX) dspevx.$(SUFFIX) dspgst.$(SUFFIX) \
   dspgv.$(SUFFIX)  dspgvd.$(SUFFIX) dspgvx.$(SUFFIX) dsprfs.$(SUFFIX) dspsv.$(SUFFIX)  dspsvx.$(SUFFIX) dsptrd.$(SUFFIX) \
   dsptrf.$(SUFFIX) dsptri.$(SUFFIX) dsptrs.$(SUFFIX) dstegr.$(SUFFIX) dstein.$(SUFFIX) dstev.$(SUFFIX)  dstevd.$(SUFFIX) dstevr.$(SUFFIX) \
   dstevx.$(SUFFIX) dsycon.$(SUFFIX) dsyev.$(SUFFIX)  dsyevd.$(SUFFIX) dsyevr.$(SUFFIX) \
   dsyevx.$(SUFFIX) dsygs2.$(SUFFIX) dsygst.$(SUFFIX) dsygv.$(SUFFIX)  dsygvd.$(SUFFIX) dsygvx.$(SUFFIX) dsyrfs.$(SUFFIX) \
   dsysv.$(SUFFIX)  dsysvx.$(SUFFIX) \
   dsytd2.$(SUFFIX) dsytf2.$(SUFFIX) dsytrd.$(SUFFIX) dsytrf.$(SUFFIX) dsytri.$(SUFFIX) dsytri2.$(SUFFIX) dsytri2x.$(SUFFIX) \
   dsyswapr.$(SUFFIX) dsytrs.$(SUFFIX) dsytrs2.$(SUFFIX) dsyconv.$(SUFFIX) dtbcon.$(SUFFIX) \
   dtbrfs.$(SUFFIX) dtbtrs.$(SUFFIX) dtgevc.$(SUFFIX) dtgex2.$(SUFFIX) dtgexc.$(SUFFIX) dtgsen.$(SUFFIX) \
   dtgsja.$(SUFFIX) dtgsna.$(SUFFIX) dtgsy2.$(SUFFIX) dtgsyl.$(SUFFIX) dtpcon.$(SUFFIX) dtprfs.$(SUFFIX) dtptri.$(SUFFIX) \
   dtptrs.$(SUFFIX) \
   dtrcon.$(SUFFIX) dtrevc.$(SUFFIX) dtrexc.$(SUFFIX) dtrrfs.$(SUFFIX) dtrsen.$(SUFFIX) dtrsna.$(SUFFIX) dtrsyl.$(SUFFIX) \
   dtrtrs.$(SUFFIX) dtzrqf.$(SUFFIX) dtzrzf.$(SUFFIX) dstemr.$(SUFFIX) \
   dsgesv.$(SUFFIX) dsposv.$(SUFFIX) dlag2s.$(SUFFIX) slag2d.$(SUFFIX) dlat2s.$(SUFFIX) \
   dlansf.$(SUFFIX) dpftrf.$(SUFFIX) dpftri.$(SUFFIX) dpftrs.$(SUFFIX) dsfrk.$(SUFFIX) dtfsm.$(SUFFIX) dtftri.$(SUFFIX) dtfttp.$(SUFFIX) \
   dtfttr.$(SUFFIX) dtpttf.$(SUFFIX) dtpttr.$(SUFFIX) dtrttf.$(SUFFIX) dtrttp.$(SUFFIX) \
   dgejsv.$(SUFFIX)  dgesvj.$(SUFFIX)  dgsvj0.$(SUFFIX)  dgsvj1.$(SUFFIX) \
   dgeequb.$(SUFFIX) dsyequb.$(SUFFIX) dpoequb.$(SUFFIX) dgbequb.$(SUFFIX) \
   dbbcsd.$(SUFFIX) dlapmr.$(SUFFIX) dorbdb.$(SUFFIX) dorcsd.$(SUFFIX) \
   dgeqrt.$(SUFFIX) dgeqrt2.$(SUFFIX) dgeqrt3.$(SUFFIX) dgemqrt.$(SUFFIX) \
   dtpqrt.$(SUFFIX) dtpqrt2.$(SUFFIX) dtpmqrt.$(SUFFIX) dtprfb.$(SUFFIX)

ifdef USEXBLAS
DXLASRC = dgesvxx.$(SUFFIX) dgerfsx.$(SUFFIX) dla_gerfsx_extended.$(SUFFIX) dla_geamv.$(SUFFIX)		\
   dla_gercond.$(SUFFIX) dla_gerpvgrw.$(SUFFIX) dsysvxx.$(SUFFIX) dsyrfsx.$(SUFFIX)			\
   dla_syrfsx_extended.$(SUFFIX) dla_syamv.$(SUFFIX) dla_syrcond.$(SUFFIX) dla_syrpvgrw.$(SUFFIX)	\
   dposvxx.$(SUFFIX) dporfsx.$(SUFFIX) dla_porfsx_extended.$(SUFFIX) dla_porcond.$(SUFFIX)		\
   dla_porpvgrw.$(SUFFIX) dgbsvxx.$(SUFFIX) dgbrfsx.$(SUFFIX) dla_gbrfsx_extended.$(SUFFIX)		\
   dla_gbamv.$(SUFFIX) dla_gbrcond.$(SUFFIX) dla_gbrpvgrw.$(SUFFIX) dla_lin_berr.$(SUFFIX) dlarscl2.$(SUFFIX)	\
   dlascl2.$(SUFFIX) dla_wwaddw.$(SUFFIX)
endif

ZLASRC = \
   zbdsqr.$(SUFFIX) zgbbrd.$(SUFFIX) zgbcon.$(SUFFIX) zgbequ.$(SUFFIX) zgbrfs.$(SUFFIX) zgbsv.$(SUFFIX)  zgbsvx.$(SUFFIX) \
   zgbtf2.$(SUFFIX) zgbtrf.$(SUFFIX) zgbtrs.$(SUFFIX) zgebak.$(SUFFIX) zgebal.$(SUFFIX) zgebd2.$(SUFFIX) zgebrd.$(SUFFIX) \
   zgecon.$(SUFFIX) zgeequ.$(SUFFIX) zgees.$(SUFFIX)  zgeesx.$(SUFFIX) zgeev.$(SUFFIX)  zgeevx.$(SUFFIX) \
   zgegs.$(SUFFIX)  zgegv.$(SUFFIX)  zgehd2.$(SUFFIX) zgehrd.$(SUFFIX) zgelq2.$(SUFFIX) zgelqf.$(SUFFIX) \
   zgels.$(SUFFIX)  zgelsd.$(SUFFIX) zgelss.$(SUFFIX) zgelsx.$(SUFFIX) zgelsy.$(SUFFIX) zgeql2.$(SUFFIX) zgeqlf.$(SUFFIX) zgeqp3.$(SUFFIX) \
   zgeqpf.$(SUFFIX) zgeqr2.$(SUFFIX) zgeqr2p.$(SUFFIX) zgeqrf.$(SUFFIX) zgeqrfp.$(SUFFIX) zgerfs.$(SUFFIX) zgerq2.$(SUFFIX) zgerqf.$(SUFFIX) \
   zgesc2.$(SUFFIX) zgesdd.$(SUFFIX) zgesv.$(SUFFIX)  zgesvd.$(SUFFIX) zgesvx.$(SUFFIX) zgetc2.$(SUFFIX) \
   zgetri.$(SUFFIX) \
   zggbak.$(SUFFIX) zggbal.$(SUFFIX) zgges.$(SUFFIX)  zggesx.$(SUFFIX) zggev.$(SUFFIX)  zggevx.$(SUFFIX) zggglm.$(SUFFIX) \
   zgghrd.$(SUFFIX) zgglse.$(SUFFIX) zggqrf.$(SUFFIX) zggrqf.$(SUFFIX) \
   zggsvd.$(SUFFIX) zggsvp.$(SUFFIX) \
   zgtcon.$(SUFFIX) zgtrfs.$(SUFFIX) zgtsv.$(SUFFIX)  zgtsvx.$(SUFFIX) zgttrf.$(SUFFIX) zgttrs.$(SUFFIX) zgtts2.$(SUFFIX) zhbev.$(SUFFIX)  \
   zhbevd.$(SUFFIX) zhbevx.$(SUFFIX) zhbgst.$(SUFFIX) zhbgv.$(SUFFIX)  zhbgvd.$(SUFFIX) zhbgvx.$(SUFFIX) zhbtrd.$(SUFFIX) \
   zhecon.$(SUFFIX) zheev.$(SUFFIX)  zheevd.$(SUFFIX) zheevr.$(SUFFIX) zheevx.$(SUFFIX) zhegs2.$(SUFFIX) zhegst.$(SUFFIX) \
   zhegv.$(SUFFIX)  zhegvd.$(SUFFIX) zhegvx.$(SUFFIX) zherfs.$(SUFFIX) zhesv.$(SUFFIX)  zhesvx.$(SUFFIX) zhetd2.$(SUFFIX) \
   zhetf2.$(SUFFIX) zhetrd.$(SUFFIX) \
   zhetrf.$(SUFFIX) zhetri.$(SUFFIX) zhetri2.$(SUFFIX) zhetri2x.$(SUFFIX) zheswapr.$(SUFFIX) \
   zhetrs.$(SUFFIX) zhetrs2.$(SUFFIX) zhgeqz.$(SUFFIX) zhpcon.$(SUFFIX) zhpev.$(SUFFIX)  zhpevd.$(SUFFIX) \
   zhpevx.$(SUFFIX) zhpgst.$(SUFFIX) zhpgv.$(SUFFIX)  zhpgvd.$(SUFFIX) zhpgvx.$(SUFFIX) zhprfs.$(SUFFIX) zhpsv.$(SUFFIX)  \
   zhpsvx.$(SUFFIX) \
   zhptrd.$(SUFFIX) zhptrf.$(SUFFIX) zhptri.$(SUFFIX) zhptrs.$(SUFFIX) zhsein.$(SUFFIX) zhseqr.$(SUFFIX) zlabrd.$(SUFFIX) \
   zlacgv.$(SUFFIX) zlacon.$(SUFFIX) zlacn2.$(SUFFIX) zlacp2.$(SUFFIX) zlacpy.$(SUFFIX) zlacrm.$(SUFFIX) zlacrt.$(SUFFIX) zladiv.$(SUFFIX) \
   zlaed0.$(SUFFIX) zlaed7.$(SUFFIX) zlaed8.$(SUFFIX) \
   zlaein.$(SUFFIX) zlaesy.$(SUFFIX) zlaev2.$(SUFFIX) zlags2.$(SUFFIX) zlagtm.$(SUFFIX) \
   zlahef.$(SUFFIX) zlahqr.$(SUFFIX) \
   zlahrd.$(SUFFIX) zlahr2.$(SUFFIX) zlaic1.$(SUFFIX) zlals0.$(SUFFIX) zlalsa.$(SUFFIX) zlalsd.$(SUFFIX) zlangb.$(SUFFIX) zlange.$(SUFFIX) \
   zlangt.$(SUFFIX) zlanhb.$(SUFFIX) \
   zlanhe.$(SUFFIX) \
   zlanhp.$(SUFFIX) zlanhs.$(SUFFIX) zlanht.$(SUFFIX) zlansb.$(SUFFIX) zlansp.$(SUFFIX) zlansy.$(SUFFIX) zlantb.$(SUFFIX) \
   zlantp.$(SUFFIX) zlantr.$(SUFFIX) zlapll.$(SUFFIX) zlapmt.$(SUFFIX) zlaqgb.$(SUFFIX) zlaqge.$(SUFFIX) \
   zlaqhb.$(SUFFIX) zlaqhe.$(SUFFIX) zlaqhp.$(SUFFIX) zlaqp2.$(SUFFIX) zlaqps.$(SUFFIX) zlaqsb.$(SUFFIX) \
   zlaqr0.$(SUFFIX) zlaqr1.$(SUFFIX) zlaqr2.$(SUFFIX) zlaqr3.$(SUFFIX) zlaqr4.$(SUFFIX) zlaqr5.$(SUFFIX) \
   zlaqsp.$(SUFFIX) zlaqsy.$(SUFFIX) zlar1v.$(SUFFIX) zlar2v.$(SUFFIX) ilazlr.$(SUFFIX) ilazlc.$(SUFFIX) \
   zlarcm.$(SUFFIX) zlarf.$(SUFFIX)  zlarfb.$(SUFFIX) \
   zlarfg.$(SUFFIX) zlarft.$(SUFFIX) zlarfgp.$(SUFFIX) \
   zlarfx.$(SUFFIX) zlargv.$(SUFFIX) zlarnv.$(SUFFIX) zlarrv.$(SUFFIX) zlartg.$(SUFFIX) zlartv.$(SUFFIX) \
   zlarz.$(SUFFIX)  zlarzb.$(SUFFIX) zlarzt.$(SUFFIX) zlascl.$(SUFFIX) zlaset.$(SUFFIX) zlasr.$(SUFFIX)  \
   zlassq.$(SUFFIX) zlasyf.$(SUFFIX) \
   zlatbs.$(SUFFIX) zlatdf.$(SUFFIX) zlatps.$(SUFFIX) zlatrd.$(SUFFIX) zlatrs.$(SUFFIX) zlatrz.$(SUFFIX) zlatzm.$(SUFFIX) \
   zpbcon.$(SUFFIX) zpbequ.$(SUFFIX) zpbrfs.$(SUFFIX) zpbstf.$(SUFFIX) zpbsv.$(SUFFIX)  \
   zpbsvx.$(SUFFIX) zpbtf2.$(SUFFIX) zpbtrf.$(SUFFIX) zpbtrs.$(SUFFIX) zpocon.$(SUFFIX) zpoequ.$(SUFFIX) zporfs.$(SUFFIX) \
   zposv.$(SUFFIX)  zposvx.$(SUFFIX) zpotri.$(SUFFIX) zpotrs.$(SUFFIX) zpstrf.$(SUFFIX) zpstf2.$(SUFFIX) \
   zppcon.$(SUFFIX) zppequ.$(SUFFIX) zpprfs.$(SUFFIX) zppsv.$(SUFFIX)  zppsvx.$(SUFFIX) zpptrf.$(SUFFIX) zpptri.$(SUFFIX) zpptrs.$(SUFFIX) \
   zptcon.$(SUFFIX) zpteqr.$(SUFFIX) zptrfs.$(SUFFIX) zptsv.$(SUFFIX)  zptsvx.$(SUFFIX) zpttrf.$(SUFFIX) zpttrs.$(SUFFIX) zptts2.$(SUFFIX) \
   zrot.$(SUFFIX)   zspcon.$(SUFFIX) zspmv.$(SUFFIX)  zspr.$(SUFFIX)   zsprfs.$(SUFFIX) zspsv.$(SUFFIX)  \
   zspsvx.$(SUFFIX) zsptrf.$(SUFFIX) zsptri.$(SUFFIX) zsptrs.$(SUFFIX) zdrscl.$(SUFFIX) zstedc.$(SUFFIX) \
   zstegr.$(SUFFIX) zstein.$(SUFFIX) zsteqr.$(SUFFIX) zsycon.$(SUFFIX) zsymv.$(SUFFIX) \
   zsyr.$(SUFFIX)   zsyrfs.$(SUFFIX) zsysv.$(SUFFIX)  zsysvx.$(SUFFIX) zsytf2.$(SUFFIX) zsytrf.$(SUFFIX) zsytri.$(SUFFIX) zsytri2.$(SUFFIX) zsytri2x.$(SUFFIX) \
   zsyswapr.$(SUFFIX) zsytrs.$(SUFFIX) zsytrs2.$(SUFFIX) zsyconv.$(SUFFIX) ztbcon.$(SUFFIX) ztbrfs.$(SUFFIX) ztbtrs.$(SUFFIX) ztgevc.$(SUFFIX) ztgex2.$(SUFFIX) \
   ztgexc.$(SUFFIX) ztgsen.$(SUFFIX) ztgsja.$(SUFFIX) ztgsna.$(SUFFIX) ztgsy2.$(SUFFIX) ztgsyl.$(SUFFIX) ztpcon.$(SUFFIX) \
   ztprfs.$(SUFFIX) ztptri.$(SUFFIX) \
   ztptrs.$(SUFFIX) ztrcon.$(SUFFIX) ztrevc.$(SUFFIX) ztrexc.$(SUFFIX) ztrrfs.$(SUFFIX) ztrsen.$(SUFFIX) ztrsna.$(SUFFIX) \
   ztrsyl.$(SUFFIX) ztrtrs.$(SUFFIX) ztzrqf.$(SUFFIX) ztzrzf.$(SUFFIX) zung2l.$(SUFFIX) \
   zung2r.$(SUFFIX) zungbr.$(SUFFIX) zunghr.$(SUFFIX) zungl2.$(SUFFIX) zunglq.$(SUFFIX) zungql.$(SUFFIX) zungqr.$(SUFFIX) zungr2.$(SUFFIX) \
   zungrq.$(SUFFIX) zungtr.$(SUFFIX) zunm2l.$(SUFFIX) zunm2r.$(SUFFIX) zunmbr.$(SUFFIX) zunmhr.$(SUFFIX) zunml2.$(SUFFIX) \
   zunmlq.$(SUFFIX) zunmql.$(SUFFIX) zunmqr.$(SUFFIX) zunmr2.$(SUFFIX) zunmr3.$(SUFFIX) zunmrq.$(SUFFIX) zunmrz.$(SUFFIX) \
   zunmtr.$(SUFFIX) zupgtr.$(SUFFIX) \
   zupmtr.$(SUFFIX) izmax1.$(SUFFIX) dzsum1.$(SUFFIX) zstemr.$(SUFFIX) \
   zcgesv.$(SUFFIX) zcposv.$(SUFFIX) zlag2c.$(SUFFIX) clag2z.$(SUFFIX) zlat2c.$(SUFFIX) \
   zhfrk.$(SUFFIX) ztfttp.$(SUFFIX) zlanhf.$(SUFFIX) zpftrf.$(SUFFIX) zpftri.$(SUFFIX) zpftrs.$(SUFFIX) ztfsm.$(SUFFIX) ztftri.$(SUFFIX) \
   ztfttr.$(SUFFIX) ztpttf.$(SUFFIX) ztpttr.$(SUFFIX) ztrttf.$(SUFFIX) ztrttp.$(SUFFIX) \
   zgeequb.$(SUFFIX) zgbequb.$(SUFFIX) zsyequb.$(SUFFIX) zpoequb.$(SUFFIX) zheequb.$(SUFFIX) \
   zbbcsd.$(SUFFIX) zlapmr.$(SUFFIX) zunbdb.$(SUFFIX) zuncsd.$(SUFFIX) \
   zgeqrt.$(SUFFIX) zgeqrt2.$(SUFFIX) zgeqrt3.$(SUFFIX) zgemqrt.$(SUFFIX) \
   ztpqrt.$(SUFFIX) ztpqrt2.$(SUFFIX) ztpmqrt.$(SUFFIX) ztprfb.$(SUFFIX)

ifdef USEXBLAS
ZXLASRC = zgesvxx.$(SUFFIX) zgerfsx.$(SUFFIX) zla_gerfsx_extended.$(SUFFIX) zla_geamv.$(SUFFIX)		\
   zla_gercond_c.$(SUFFIX) zla_gercond_x.$(SUFFIX) zla_gerpvgrw.$(SUFFIX) zsysvxx.$(SUFFIX) zsyrfsx.$(SUFFIX)	\
   zla_syrfsx_extended.$(SUFFIX) zla_syamv.$(SUFFIX) zla_syrcond_c.$(SUFFIX) zla_syrcond_x.$(SUFFIX)	\
   zla_syrpvgrw.$(SUFFIX) zposvxx.$(SUFFIX) zporfsx.$(SUFFIX) zla_porfsx_extended.$(SUFFIX)		\
   zla_porcond_c.$(SUFFIX) zla_porcond_x.$(SUFFIX) zla_porpvgrw.$(SUFFIX) zgbsvxx.$(SUFFIX) zgbrfsx.$(SUFFIX)	\
   zla_gbrfsx_extended.$(SUFFIX) zla_gbamv.$(SUFFIX) zla_gbrcond_c.$(SUFFIX) zla_gbrcond_x.$(SUFFIX)	\
   zla_gbrpvgrw.$(SUFFIX) zhesvxx.$(SUFFIX) zherfsx.$(SUFFIX) zla_herfsx_extended.$(SUFFIX)		\
   zla_heamv.$(SUFFIX) zla_hercond_c.$(SUFFIX) zla_hercond_x.$(SUFFIX) zla_herpvgrw.$(SUFFIX)		\
   zla_lin_berr.$(SUFFIX) zlarscl2.$(SUFFIX) zlascl2.$(SUFFIX) zla_wwaddw.$(SUFFIX)
endif

ALLOBJ = $(SLASRC) $(DLASRC) $(DSLASRC) $(CLASRC) $(ZLASRC) $(ZCLASRC) \
   $(SCLAUX) $(DZLAUX) $(ALLAUX)

ALLOBJ_P = $(ALLOBJ:.$(SUFFIX)=.$(PSUFFIX))

ifdef USEXBLAS
ALLXOBJ = $(SXLASRC) $(DXLASRC) $(CXLASRC) $(ZXLASRC)
endif

all: ../$(LAPACKLIB)

lapack_prof: ../$(LAPACKLIB_P)

../$(LAPACKLIB): $(ALLOBJ) $(ALLXOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ) $(ALLXOBJ)
	$(RANLIB) $@

../$(LAPACKLIB_P): $(ALLOBJ_P)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ_P)
	$(RANLIB) $@

single: $(SLASRC) $(DSLASRC) $(SXLASRC) $(SCLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(SLASRC) $(DSLASRC) \
	$(SXLASRC) $(SCLAUX) $(ALLAUX) $(ALLXAUX)
	$(RANLIB) ../$(LAPACKLIB)

complex: $(CLASRC) $(ZCLASRC) $(CXLASRC) $(SCLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(CLASRC) $(ZCLASRC) \
	$(CXLASRC) $(SCLAUX) $(ALLAUX) $(ALLXAUX)
	$(RANLIB) ../$(LAPACKLIB)

double: $(DLASRC) $(DSLASRC) $(DXLASRC) $(DZLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(DLASRC) $(DSLASRC) \
	$(DXLASRC) $(DZLAUX) $(ALLAUX) $(ALLXAUX)
	$(RANLIB) ../$(LAPACKLIB)

complex16: $(ZLASRC) $(ZCLASRC) $(ZXLASRC) $(DZLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(LAPACKLIB) $(ZLASRC) $(ZCLASRC) \
	$(ZXLASRC) $(DZLAUX) $(ALLAUX) $(ALLXAUX)
	$(RANLIB) ../$(LAPACKLIB)

$(ALLAUX): $(FRC)
$(SCLAUX): $(FRC)
$(DZLAUX): $(FRC)
$(SLASRC): $(FRC)
$(CLASRC): $(FRC)
$(DLASRC): $(FRC)
$(ZLASRC): $(FRC)
$(ZCLASRC): $(FRC)
$(DSLASRC): $(FRC)
ifdef USEXBLAS
$(SXLASRC): $(FRC)
$(CXLASRC): $(FRC)
$(DXLASRC): $(FRC)
$(ZXLASRC): $(FRC)
endif

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.$(SUFFIX) *.$(PSUFFIX)

%.$(SUFFIX): %.f
	$(FORTRAN) $(OPTS) -c $< -o $@

%.$(PSUFFIX): %.f
	$(FORTRAN) $(POPTS) -c $< -o $@


slaruv.$(SUFFIX): slaruv.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@
dlaruv.$(SUFFIX): dlaruv.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@
sla_wwaddw.$(SUFFIX): sla_wwaddw.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@
dla_wwaddw.$(SUFFIX): dla_wwaddw.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@
cla_wwaddw.$(SUFFIX): cla_wwaddw.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@
zla_wwaddw.$(SUFFIX): zla_wwaddw.f ; $(FORTRAN) $(NOOPT) -O0 -c $< -o $@

slaruv.$(PSUFFIX): slaruv.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
dlaruv.$(PSUFFIX): dlaruv.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
sla_wwaddw.$(PSUFFIX): sla_wwaddw.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
dla_wwaddw.$(PSUFFIX): dla_wwaddw.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
cla_wwaddw.$(PSUFFIX): cla_wwaddw.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
zla_wwaddw.$(PSUFFIX): zla_wwaddw.f ; $(FORTRAN) $(PNOOPT) -O0 -c $< -o $@
