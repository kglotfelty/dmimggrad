##############################################################################


#MK_TOP = ../../../..

MK_TOP = /export/ciao_from_source/ciao-4.8/src
KJG = /export/ciao


include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmimggrad
LIB_FILES         =
PAR_FILES         = dmimggrad.par
INC_FILES         =
XML_FILES         = dmimggrad.xml

SRCS	= dmimggrad.c t_dmimggrad.c

LOCAL_LIBS = -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I$(MK_TOP)/da/analysis/dmtools/dmimgio/

OBJS	= $(SRCS:.c=.o)

MAKETEST_SCRIPT = dmimggrad.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LINK)
	@echo

announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |             Building dmimggrad DM host tool           | "
	@echo "   \----------------------------------------------------------/ "

kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)
