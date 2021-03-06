# ===========================================================================
#  Makefile urQMD                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile		[uqmd.f]
##			TYPE	=	operating system	['uname']
##                      LINK    =       linking additional libraries
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  


# ----------------------------------------- 

FC := gfortran
FFLAGS = -O3 -cpp

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(FFLAGS)
SYSTEMFILES     =       $(SRCGNU)


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	urqmd.e
endif

SRC 		= 	addpart.f angdis.f anndec.f blockres.f cascinit.f \
	   		coload.f dectim.f delpart.f detbal.f dwidth.f \
	   		getmass.f getspin.f init.f input.f ityp2pdg.f error.f\
	   		make22.f numrec.f output.f paulibl.f jdecay2.f iso.f\
	   		proppot.f saveinfo.f scatter.f siglookup.f string.f \
	   		tabinit.f whichres.f boxprg.f shine.f urqmd.f pdummy.f
INC		= 	colltab.f comnorm.f complot.f coms.f comres.f comwid.f \
			freezeout.f inputs.f newpart.f options.f boxinc.f \
			outcom.f 

SRCGNU		=	genranf.f


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) $(SYSTEMFILES) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC) $(SYSTEMFILES))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.f:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(FC) $(OBJECTS) $(LIBS) $(LDFLAGS) -o $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj*

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)/$(basename $(TARGET))
		strip $(INSTPATH)/$(basename $(TARGET))

# --------------- Dependencies -------------------

./addpart.f:	coms.f newpart.f
./angdis.f:	options.f
./anndec.f: 	comres.f comwid.f options.f
./blockres.f: 	comres.f
./boxprg.f:    boxinc.f coms.f comres.f options.f
./cascinit.f: 	coms.f options.f 
./coload.f: 	colltab.f coms.f comres.f newpart.f options.f
./dectim.f: 	coms.f colltab.f options.f 
./delpart.f:	coms.f comres.f newpart.f
./detbal.f: 	coms.f comres.f comwid.f newpart.f options.f
./dwidth.f:	coms.f comres.f comwid.f options.f
./getmass.f:    comnorm.f comres.f
./getspin.f:	comres.f 
./init.f: 	boxinc.f colltab.f coms.f comres.f freezeout.f inputs.f options.f
./input.f: 	boxinc.f coms.f comres.f inputs.f options.f
./ityp2pdg.f:	comres.f
./iso.f:	comres.f newpart.f
./jdecay2.f:	newpart.f
./make22.f:	comres.f comwid.f newpart.f options.f
./output.f: 	boxinc.f coms.f comres.f freezeout.f inputs.f newpart.f options.f
./paulibl.f: 	coms.f options.f
./proppot.f: 	boxinc.f coms.f options.f
./saveinfo.f: 	coms.f
./scatter.f: 	coms.f comres.f freezeout.f newpart.f options.f
./shine.f:	coms.f comres.f options.f
./siglookup.f:	comres.f
./string.f:	comres.f options.f
./tabinit.f:	comres.f comwid.f
./urqmd.f: 	boxinc.f colltab.f coms.f comres.f inputs.f newpart.f options.f
./whichres.f: 	comres.f options.f
