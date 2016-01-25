include $(DEPTH)/src/templates/cleanup_macros.mk

DIR_PATH := src/apps/Database_Trial_v1
APPNAME := Database_Trial_v1

C_SOURCES := 

CXX_SOURCES := \
	Database_Trial_v1.cpp


EXPORTED_HEADERS   := 

REQUIRED_LIBS := $(SEDRIS_CORE_LIBS)

LOCAL_INCLUDES := -I$(INC_DIR)

SYS_LIBS =

include $(DEPTH)/src/templates/localtargets.mk
