TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on release

HEADERS	+= inc/filehandler.h \
	inc/gui_constants.h

SOURCES	+= tkgeomgui.cpp \
	src/FileHandler.cc

FORMS	= maindialog.ui

unix {
  UI_DIR = .ui
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}

